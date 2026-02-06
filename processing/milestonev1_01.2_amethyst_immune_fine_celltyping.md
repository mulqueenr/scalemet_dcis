Running subclustering on just immune cells for finer cell typing

1. using RNA based marker genes for hypomethylation over genes
2. overlapping dmrs and rna marker genes (fishers exact test)
3. ranking RNA marker genes and DMRs for better markers
4. integration with rliger

Output percent methylation as bedgraph to view in IGV over marker genes as well

```R
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing

#source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
library(amethyst)
library(data.table)
library(dplyr)
library(msigdbr)
library(fgsea)
library(ggplot2)
library(GenomicRanges)
library(Matrix)
library(parallel)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)

system(paste0("mkdir -p ",project_data_directory,"/fine_celltyping"))

dat<-readRDS(file="05_scaledcis.coarse_clusters.amethyst.rds")
rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")

```

# Immune Cells

## Subcluster on just immune cells
```R
broad_celltype="immune"

output_directory=paste0(project_data_directory,"/fine_celltyping/",broad_celltype)
system(paste0("mkdir -p ",output_directory))
print(paste("Subsetting object based on broad celltype:",broad_celltype))
dat<-subsetObject(dat,cells=row.names(dat@metadata[dat@metadata$broad_celltype %in% broad_celltype,]))


celltype_umap<-function(obj=dat,prefix="allcells",dims=12,regressCov=TRUE,k_pheno=50,neigh=25,dist=1e-5,method="cosine"){
  print("Running IRLBA reduction...")
  obj@reductions[[paste(window_name,"irlba",sep="_")]] <- runIrlba(obj, genomeMatrices = c(window_name), dims = dims, replaceNA = c(0))
  
  if(regressCov){
    print("Running regression...")
  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- regressCovBias(obj, reduction = paste(window_name,"irlba",sep="_")) 
  } else {
      print("Skipping regression...")
  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- obj@reductions[[paste(window_name,"irlba",sep="_")]]
  }

  obj <- amethyst::runCluster(obj, k_phenograph = k_pheno, reduction = paste(window_name,"irlba_regressed",sep="_")) 

  print(paste("Running UMAP...",as.character(neigh),as.character(dist),as.character(method)))
  obj <- amethyst::runUmap(obj, neighbors = neigh, dist = dist, method = method, reduction = paste(window_name,"irlba_regressed",sep="_")) 

  outname=paste(prefix,"integrated_celltype",dims,as.character(regressCov),k_pheno,neigh,as.character(dist),method,sep="_")
  print(paste("Plotting...",outname))

  p1 <- dimFeature(obj, colorBy = sample, reduction = "umap") + ggtitle(paste(window_name,"Samples"))
  p2 <- dimFeature(obj, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red"),guide="colourbar") + ggtitle("Coverage distribution")
  p3 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
  p4 <- dimFeature(obj, colorBy = cluster_id, reduction = "umap") + ggtitle(paste(window_name,"Clusters"))
  p5 <- dimFeature(obj, colorBy = Group, reduction = "umap") + ggtitle(paste(window_name," Group"))
  p6<-ggplot()
  ggsave((p1|p2)/(p3|p4)/(p5|p6),file=paste0(outname,"_umap.pdf"),width=20,height=30)  
  return(obj)
  }

window_name="coarse_cluster_dmr_sites"
dat@genomeMatrices[[window_name]] <- dat@genomeMatrices[[window_name]][rowSums(!is.na(dat@genomeMatrices[[window_name]])) >= 45, ]
est_dim<-dimEstimate(dat, genomeMatrices = c(window_name), dims = c(20), threshold = 0.95)
print(est_dim)
#7

dat<-celltype_umap(obj=dat,prefix="06_immune_finecelltyping",dims=8,regressCov=FALSE,k_pheno=50,neigh=5,dist=0.01,method="cosine") 
dat@metadata$fine_cluster_id<-dat@metadata$cluster_id
#happy with this immune cell clustering! looks like good structure!

saveRDS(dat,file="06_scaledcis.immune_finecelltyping.amethyst.rds")
```

## Summarize new immune clusters for marker plotting


```R
celltype500bpwindows <- calcSmoothedWindows(dat, 
                                        type = "CG", 
                                        threads = 100,
                                        step = 500, 
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "fine_cluster_id",
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)

saveRDS(celltype500bpwindows,file=paste0("06_immune_finecelltyping",".500bp_windows.rds"))

dat@genomeMatrices[["cg_immunecells_perc"]] <- celltype500bpwindows[["pct_matrix"]]
saveRDS(dat,file="06_scaledcis.immune_finecelltyping.amethyst.rds")


```

```R


#modified from histograM
histograModified <- function(obj,
    genes = NULL,
    matrix,colors = NULL,
    ncol = length(genes),
    trackOverhang = 5000,arrowOverhang = 3000,trackScale = 1.5,arrowScale = 0.025,colorMax = 100,
    legend = TRUE,removeNA = TRUE,
    order = NULL,baseline = "mean",
    promoter_focus=FALSE,
    cgisland=NULL) {

    if (!is.null(colors)) {pal <- colors} else {pal <- rev(c("#FF0082", "#dbdbdb", "#cccccc", "#999999"))}
    
    genes<-genes[genes %in% obj@ref$gene_name] #ensure genes are in ref
    p <- vector("list", length(genes)) # empty plot list

    for (i in 1:length(genes)) {
            group_count<-ncol(obj@genomeMatrices[[matrix]])-3 #-3 for chr start end
            ngroup<-vector("list",group_count+1)  #+1 for gene track
            if (!is.null(order)) {
                names(ngroup)<-c("gene_track",order)
            }else{
                names(ngroup)<-c("gene_track",colnames(obj@genomeMatrices[[matrix]])[4:ncol(obj@genomeMatrices[[matrix]])])
            }

            ref <- obj@ref |> dplyr::filter(gene_name == genes[i])
            aggregated <- obj@genomeMatrices[[matrix]]

            toplot <- aggregated[c((aggregated$chr == ref$seqid[ref$type == "gene"] &
                                aggregated$start > (ref$start[ref$type == "gene"] - trackOverhang) &
                                aggregated$end < (ref$end[ref$type == "gene"] + trackOverhang))), ]

            promoter <- ref |> dplyr::filter(type == "gene") |>
                  dplyr::mutate(promoter_start = ifelse(strand == "+", (start - 2000), (end + 2000)),
                  promoter_end = ifelse(strand == "+", (promoter_start+2000), (promoter_start-2000)))
            
            if(promoter_focus){
              if(ref[ref$type == "gene",]$strand[1] == "+"){
                plot_left<-promoter$promoter_start-trackOverhang
                plot_right<-promoter$promoter_end+trackOverhang
              }else{
                plot_right<-promoter$promoter_start+trackOverhang
                plot_left<-promoter$promoter_end-trackOverhang
              }
            }
            exon <- ref |> dplyr::filter(type == "exon")
            
            trackHeight <- group_count * trackScale
            toplot <- tidyr::pivot_longer(toplot, cols = c(4:ncol(toplot)), names_to = "group", values_to = "pct_m") |> dplyr::rowwise() |> dplyr::mutate(middle = mean(c(start, end), na.rm = TRUE))
        
        if (removeNA) {toplot <- toplot |> dplyr::filter(!is.na(pct_m))}

        if (baseline == "mean") {
        glob_m <- data.frame(group = colnames(aggregated[, 4:ncol(aggregated)]),glob_m = colMeans(aggregated[, 4:ncol(aggregated)], na.rm = T))
        toplot <- dplyr::left_join(toplot, glob_m, by = "group")}

        if (!is.null(order)) {toplot$group <- factor(toplot$group, levels = order)}
        
        if(!is.null(cgisland)) {cgi <- cgisland[cgisland$chr == as.character(toplot$chr[1]) &
                            cgisland$start > min(toplot$start) &
                            cgisland$end < max(toplot$end),]}

        for (j in names(ngroup)[2:length(ngroup)]){
            toplot_sub<-toplot[toplot$group==j,]
            # plotting histogram
            p[[i]][[j]] <- ggplot2::ggplot() +
            {if (baseline == 0) ggplot2::geom_col(data = toplot_sub, ggplot2::aes(x = middle, y = pct_m, fill = pct_m), width = mean(toplot_sub$end - toplot_sub$start))} +
            {if (baseline == "mean")ggplot2::geom_rect(data = toplot_sub, ggplot2::aes(xmin = start, xmax = end, ymin = glob_m, ymax = pct_m, fill = pct_m))} + 
            {if(promoter_focus) ggplot2::coord_cartesian(xlim=c(plot_left,plot_right),expand=FALSE)}+
            ggplot2::scale_fill_gradientn(colors = pal, limits = c(0,colorMax), oob = scales::squish) + theme_minimal()+ ylab(j) + ggplot2::theme(legend.position="none") +
            ylim(c(0,100))+
            ggplot2::theme(panel.background = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank(), panel.grid.major.y = element_line(color = "#dbdbdb", linetype = "dashed"))
            }
        #last track is annot
        p[[i]][["gene_track"]] <-ggplot()+
        ggplot2::geom_rect(fill = "pink",alpha=0.8, data = promoter, ggplot2::aes(xmin = promoter_start, xmax = promoter_end, ymin = 0, ymax = trackHeight/4)) +
        ggplot2::geom_rect(alpha=0.8,fill = "gray", data = exon,ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = trackHeight/8)) +
        ggplot2::geom_rect(alpha=0.8,fill = "green", data = cgi, ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = trackHeight/4)) +
        ggplot2::geom_segment(data = ref, color=ifelse(ref$strand == "+","gray","black"),
                            aes(x = ifelse(strand == "+", (min(start) - arrowOverhang), (max(end)) + arrowOverhang),
                            y = trackHeight/8,
                            xend = ifelse(strand == "+", (max(end) + arrowOverhang), (min(start)) - arrowOverhang),
                            yend = trackHeight/8), arrow = arrow(length = unit(trackHeight/40, "cm"))) + 
                            labs(title=paste(genes[i]),subtitle=paste(toplot$chr[1],toplot$start[1],"-",toplot$end[nrow(toplot)]),size=3) +
                            theme(axis.title.x = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank(),
                            axis.text.y = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.border = element_blank(),
                            panel.background = element_blank(),
                            title=element_text(size=8))+
                            theme(legend.position="none") #+xlim(c(xmin_pos,xmax_pos))

    if(promoter_focus){ #focus to plot region to promoter + track overhang for either direction
      p[[i]][["gene_track"]]<-p[[i]][["gene_track"]]+ggplot2::coord_cartesian(xlim=c(plot_left,plot_right),expand=FALSE)
    }
    p[[i]]<-wrap_plots(p[[i]],guides='collect')+plot_layout(ncol=1)
    }
    all_genes_plot<-wrap_plots(p)+plot_layout(nrow=1)
    return(all_genes_plot)
}

#more broad cell markers
cell_markers<-list()
cell_markers[["myeloid"]]<-c("HLA-DRA","HLA-DPA1","CD74")
cell_markers[["tcell"]]<-c("PTPRC","IKZF1","IL7R","GNLY")
cell_markers[["bcell"]]=c("CD37","TCL1A","LTB","HLA-DPB1")
cell_markers[["plasma"]]=c("IGHA2","IGHA1","JCHAIN","IGHM","IGHG1","IGHG4","IGHG3","IGHG2")

#https://www.nature.com/articles/s41588-024-01688-9/figures/2
cell_markers[["immune"]]<-c("PTPRC","CXCR4") #both good
cell_markers[["tcell"]]<-c("CD3D","CD3G") #both good
cell_markers[["tcell_CD4"]]<-c("CD4","CCR4") #not great
cell_markers[["tcell_CD8"]]<-c("CD8A","CD8B") #cd8b good
cell_markers[["bcell"]]<-c("MS4A1","CD79A") #both good
cell_markers[["macro"]]<-c("FCER1G","C1QB","APOE","TREM2","LYVE1","FOLR2","NRP1") #FCER1G good, TREM2 FOLR2 is TAM




immune_markers<-list()
immune_markers[["myeloid_mono_classical"]]<-c('SERPINB2', 'AQP9', 'FCN1', 'VCAN', 'AC245128.3', 'EREG', 'TNIP3', 'CD300E', 'S100A8', 'S100A9', 'THBS1', 'CYP1B1', 'SLC39A8', 'MAP3K20', 'EHD1', 'LUCAT1', 'FCAR', 'S100A12', 'PID1', 'ACOD1')
immune_markers[["myeloid_mono_nonclassical"]]<-c('CCDC68', 'CLEC12A', 'CD48', 'OLIG1', 'CDKN1C', 'S100A4', 'KRTAP2-3', 'CD52', 'WARS', 'LINC02432', 'CFP', 'UNC45B', 'APOBEC3A', 'TKT', 'RIPOR2', 'CDA', 'PLAC8', 'FFAR2', 'PAPSS2', 'GK5')
immune_markers[["myeloid_macro_C3"]]<-c('C3', 'GPR158', 'C1QC', 'C1QB', 'AC079760.2', 'FCGR3A', 'SDS', 'A2M', 'AL590705.1', 'MT3', 'TREM2', 'HLA-DQB2', 'OLR1', 'ZBED9', 'IL4I1', 'DOK5', 'IGLC7', 'PMEPA1', 'IBSP', 'CD9')
immune_markers[["myeloid_macro_LYVE1"]]<-c('PDK4', 'SELENOP', 'FOLR2', 'LYVE1', 'LILRB5', 'MAF', 'KLF4', 'RNASE1', 'SLC40A1', 'RHOB', 'F13A1', 'MRC1', 'TSC22D3', 'PLTP', 'STAB1', 'NR4A2', 'MS4A4A', 'DAB2', 'CD163', 'CTSC')
immune_markers[["myeloid_macro_CCL4"]]<-c('CCL4', 'SELENOP', 'FOLR2', 'RNASE1', 'MRC1', 'F13A1', 'LYVE1', 'CD163', 'STAB1', 'CXCL3', 'HMOX1', 'CCL3', 'CXCL2', 'AC010980.2', 'GCLM', 'CD93', 'CTSL', 'PLTP', 'MAFB', 'THBD')
immune_markers[["myeloid_macro_APOC1"]]<-c('ACP5', 'CYP27A1', 'LIPA', 'PLD3', 'TREM2', 'PLA2G7', 'CD36', 'FABP4', 'GPNMB', 'CCL18', 'NCEH1', 'RARRES1', 'FUCA1', 'APOC1', 'LPL', 'OTOA', 'BLVRB', 'HS3ST2', 'NR1H3', 'CHIT1')
immune_markers[["myeloid_macro_interferon"]]<-c( 'EPSTI1', 'CXCL10', 'GBP1', 'STAT1', 'CXCL11', 'CXCL9', 'GBP4', 'IFI44L', 'FCGR1A', 'SERPING1', 'IFI6', 'IGSF6', 'IFIH1', 'IDO1', 'RASSF4', 'XAF1', 'APOL6', 'ISG15', 'FYB1', 'PSME2')
immune_markers[["myeloid_macro_stress"]]<-c( 'HSPA1B', 'BAG3', 'DNAJB1', 'HSPA6', 'MRPL18', 'ZFAND2A', 'CHORDC1', 'SDS', 'CACYBP', 'FKBP4', 'DNAJA4', 'HSPE1', 'RGS1', 'HSPH1', 'HSPA1A', 'OLR1', 'HSPB1', 'CLEC2B', 'CKS2', 'PLEK')
immune_markers[["myeloid_macro_CTSK"]]<-c( 'AP000904.1', 'CTSK', 'ACP5', 'SLC9B2', 'MMP9', 'DGKI', 'SPP1', 'TCIRG1', 'ATP6V0D2', 'CKB', 'CST3', 'CD109', 'TIMP2', 'JDP2', 'SPARC', 'C9orf16', 'GRN', 'SIGLEC15', 'AK5', 'SNX10')
immune_markers[["myeloid_cycling"]]<-c( 'MKI67', 'PCLAF', 'CKS1B', 'DIAPH3', 'TYMS', 'TOP2A', 'HIST1H1B', 'ASPM', 'TK1', 'DHFR', 'NUSAP1', 'NCAPG', 'DLGAP5', 'CENPM', 'TPX2', 'RRM2', 'CEP55', 'MYBL2', 'BIRC5', 'CENPF')
immune_markers[["myeloid_neutrophil"]]<-c('CRISP3', 'ANKRD34B', 'ABCA13', 'CD177', 'ORM1', 'SERPINB4', 'DEFA3', 'DPYS', 'C7', 'AZU1', 'CLC', 'FOXQ1', 'RAB6B', 'SFRP1', 'TWIST1', 'MAP1B', 'SLPI', 'ADH1B', 'FOLR3', 'BTNL9')
immune_markers[["myeloid_mast"]]<-c('CPA3', 'TPSAB1', 'MS4A2', 'TPSB2', 'HDC', 'IL1RL1', 'GATA2', 'HPGDS', 'GCSAML', 'HPGD', 'ADCYAP1', 'KIT', 'KRT1', 'PRG2', 'CTSG', 'CLU', 'TPSD1', 'CALB2', 'RAB27B', 'SLC18A2')
immune_markers[["myeloid_cDC1"]]<-c('CLEC9A', 'XCR1', 'DNASE1L3', 'IDO1', 'LGALS2', 'CPNE3', 'WDFY4', 'C1orf54', 'DAPP1', 'RAB11FIP1', 'PPA1', 'CPVL', 'C9orf135', 'FOXB1', 'AC016717.2', 'CST3', 'PLEKHA5', 'GPR157', 'FOXN2', 'NCOA7')
immune_markers[["myeloid_cDC2"]]<-c('FCER1A', 'IL1R2', 'CLEC10A', 'CD1C', 'CST7', 'GPAT3', 'DAPP1', 'CFP', 'EREG', 'CCL22', 'LGALS2', 'JAML', 'PID1', 'AREG', 'IL7R', 'AC020656.1', 'IL1R1', 'MARCKSL1', 'ADAM19', 'SLC7A11')
immune_markers[["myeloid_mDC"]]<-c('LAMP3', 'BIRC3', 'LACRT', 'NUB1', 'CCR7', 'MARCKSL1', 'IDO1', 'DAPP1', 'POGLUT1', 'LINC01539', 'GPR157', 'IL12B', 'LAD1', 'KIF2A', 'FSCN1', 'IL7R', 'TXN', 'DUSP5', 'FOXD4L1', 'CD200', 'RAB9A')
immune_markers[["myeloid_pDC"]]<-c('CD5', 'LYPD8', 'PRL', 'AC136475.3', 'LINC01087', 'AC006058.1', 'BDKRB2', 'POTEI', 'DSP', 'AC026369.3', 'GZMB', 'TNFSF4', 'CD2', 'TIGIT', 'LTB', 'TMEM45A', 'PGR', 'CD3G', 'AC015936.1', 'ACOT7')

cell_colors=c(
"bcell"="#65cbe4",
"plasma"="#7ecdc2",
"immune"="#8088c2",
"tcell"="#1d87c8",
"tcell_CD8"="#8088c2",
"tcell_CD4"="#8088c2",
"macro"="#8088c2"
)

#prepare cgi
cgisland="/data/rmulqueen/projects/scalebio_dcis/ref/cpgIslandExt.bed"
cgi<-rtracklayer::import(cgisland)
cgi<-as.data.frame(cgi)
colnames(cgi)<-c("chr","start","end","strand")

#rough ordering by what clustered together
order<-c(
"9",
"1","7","2","10",
"6",
"8",
"5",
"4",
"3")
    
plt_list<-lapply(names(cell_markers),
function(celltype){
    genes<-unlist(cell_markers[celltype])
    genes<-genes[genes %in% dat@ref$gene_name]
    print(paste("Plotting:",celltype))
    print(paste("Genes to plot:",genes))
    plt<-histograModified(obj=dat, 
        baseline="mean",
        genes = unlist(cell_markers[celltype]),
        colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
        matrix = "cg_immunecells_perc", arrowScale = .03, trackScale = .5,
        legend = F, cgisland=cgi,order=order) + ggtitle(celltype)
    return(plt)
})

plt_out<-wrap_plots(plt_list,ncol=1) 
ggsave(plt_out ,file="06_immune.finecelltyping.marker.pdf",width=30,height=length(cell_colors)*30,limitsize=F)



celltype_assignment<-c(
"9"="tcell",
"1"="tcell","7"="tcell","2"="tcell","10"="tcell",
"6",
"8",
"5",
"4"="bcell",
"3")
```

hard to tell from these markers, trying liger integration on just immune cells
and trying dmr marker sites

```R

celltype500bpwindows<-readRDS(file=paste0("06_immune_finecelltyping",".500bp_windows.rds"))

pct_mat<-celltype500bpwindows[["pct_matrix"]] 
sum_mat<-celltype500bpwindows[["sum_matrix"]] 

dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
        eachVsAll = TRUE, # If TRUE, each group found in the sumMatrix will be tested against all others
        nminTotal = 3, # Min number observations across all groups to include the region in calculations
        nminGroup = 3) # Min number observations across either members or nonmembers to include the region

saveRDS(dmrs,file=paste0("06_immune_finecelltyping",".500bp_dmrs.rds"))

dmrs <- filterDMR(dmrs, 
            method = "bonferroni", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
            filter = FALSE, # If TRUE, removes insignificant results
            pThreshold = 0.05, # Maxmimum adjusted p value to allow if filter = TRUE
            logThreshold = 1.5) # Minimum absolute value of the log2FC to allow if filter = TRUE

collapsed_dmrs <- collapseDMR(dat, 
                       dmrs, 
                        maxDist = 2000, # Max allowable overlap between DMRs to be considered adjacent
                        minLength = 500, # Min length of collapsed DMR window to include in the output
                        reduce = T, # Reduce results to unique observations (recommended)
                        annotate = T) # Add column with overlapping gene names

saveRDS(collapsed_dmrs,file=paste0("06_immune_finecelltyping",".500bp_dmrs.filt_collapsed.rds"))

rename_dmr_output<-setNames(nm=1:length(unique(collapsed_dmrs$test)),gsub(colnames(sum_mat)[grepl(colnames(sum_mat),pattern="_t$")],pattern="_t",replacement=""))
collapsed_dmrs$type <- rename_dmr_output[collapsed_dmrs$test]
saveRDS(collapsed_dmrs,file=paste0("06_immune_finecelltyping",".500bp_dmrs.filt_collapsed.rds"))

#plot dmr counts per cluster
pal <- colorRampPalette(colors = c("#F9AB60", "#E7576E", "#630661", "#B5DCA5"))
COLS <- pal(length(unique(collapsed_dmrs$type)))

plt<-ggplot(collapsed_dmrs |> dplyr::group_by(type, direction) |> dplyr::summarise(n = n()), 
    aes(y = type, x = n, fill = type)) + 
    geom_col() + 
    facet_grid(vars(direction), scales = "free_y") + 
    scale_fill_manual(values = COLS) + 
    theme_classic()
ggsave(plt,file=paste0("06_immune_finecelltyping",".500bp_dmrs.filt_collapsed.barplot.pdf"))

```

## DMR site overlap with RNA marker genes

```R
#compare dmr sites with rna marker genes overlap
library(Seurat)
library(GenomicRanges)

rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")
rna<-subset(rna,coarse_celltype %in% c("bcell","myeloid","plasma","tcell"))
Idents(rna)<-rna$fine_celltype
rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)
rna <- ScaleData(rna)
rna <- JoinLayers(rna)
rna_markers<-FindAllMarkers(rna,assay="RNA",only.pos=TRUE)

#recluster RNA on just immune too
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna <- RunPCA(rna, features = VariableFeatures(object = rna))

#function for umap clustering
cluster_object<-function(obj=obj,dims=1:15,res=0.2){
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    obj <- RunPCA(obj, features = VariableFeatures(object = obj))
    obj <- FindNeighbors(obj, dims = dims)
    obj <- FindClusters(obj, resolution = res)
    obj <- RunUMAP(obj, dims = dims)
    return(obj)
}
rna<-cluster_object(obj=rna,dims=1:15,res=0.2)
plt<-DimPlot(obj,group.by="sample",raster=T,label=TRUE)
ggsave(plt,file=paste("tenx_dcis","umap",prefix,"sample.pdf",sep="."),width=30,height=20,limitsize=FALSE)
plt<-DimPlot(obj,group.by="fine_celltype",raster=T,label=TRUE)
ggsave(plt,file=paste("tenx_dcis","umap",prefix,"celltype.pdf",sep="."),width=30,height=20,limitsize=FALSE)


#overexpressed genes by cell type
rna_markers_filt<-rna_markers %>% 
                    filter(avg_log2FC > 3) %>% 
                    filter(p_val_adj <0.01) %>% 
                    group_by(cluster)


#filter collapsed dmrs by padj and extract gene names per type
dmr_out<-collapsed_dmrs %>% 
                filter(direction=="hypo") %>% 
                filter(dmr_padj<0.05) %>% 
                filter(abs(dmr_logFC) > 1.5) %>% 
                filter(dmr_length<50000) %>% 
                filter(!is.na(gene_names)) %>% 
                as.data.frame()

dmr_out$met_feature<-paste(dmr_out$chr,dmr_out$dmr_start,dmr_out$dmr_end,sep="_")

#overlap with rna markers
#expand dmr_out so each row is its own gene name
#this will be used to summarize multiple DMRs to their genes
dmrs_to_gene_names<-lapply(1:nrow(dmr_out), function(x) {
  temp<-dmr_out[x,]
  genes<-unlist(temp$gene_names %>% stringr::str_replace_all(" ","") %>% strsplit(split=","))
  return(cbind(genes=genes,type=temp$type,met_feature=temp$met_feature))
})

dmr_genes<-as.data.frame(do.call("rbind",dmrs_to_gene_names))
colnames(dmr_genes)<-c("gene","cluster","met_feature")

dmr_genes<-dmr_genes %>% filter(!duplicated(gene,cluster)) %>% group_by(cluster) 


#count intersect between type in DMRs and cell type in RNA markers
rna_markers_filt
dmr_genes
celltype_dmr_overlaps <- dmr_genes %>%
  inner_join(rna_markers_filt, by = "gene", relationship="many-to-many") %>%
  group_by(cluster.x,) %>%
  count(cluster.y, name = "overlap_count") %>%
  as.data.frame()

colnames(celltype_dmr_overlaps)<-c("cluster","celltype","dmr_rnamarker_overlap")

#add count of dmr features passing filter per cluster
celltype_dmr_overlaps$dmr_cluster_total<-NA
dmr_count<-dmr_genes %>% filter(!duplicated(met_feature)) %>% count() 
dmr_count<-setNames(nm=dmr_count$cluster,dmr_count$n)
celltype_dmr_overlaps$dmr_cluster_total<-NA
celltype_dmr_overlaps$dmr_cluster_total<-dmr_count[celltype_dmr_overlaps$cluster]

#add count of rna markers passing filter per cell type
rna_marker_count<-rna_markers_filt %>% count() 
rna_marker_count<-setNames(nm=rna_marker_count$cluster,rna_marker_count$n)
celltype_dmr_overlaps$rna_celltype_total<-NA
celltype_dmr_overlaps$rna_celltype_total<-rna_marker_count[celltype_dmr_overlaps$celltype]

#perform fisher test on set
row_fisher <- function(i){
    df<-celltype_dmr_overlaps[i,]
    p11<-df["dmr_rnamarker_overlap"]
    p10<-df["dmr_cluster_total"]-df["dmr_rnamarker_overlap"]
    p01<-df["rna_celltype_total"]
    p00<-length(Features(rna))
    mat <- matrix(as.numeric(c(p11,p10,p01,p00)), ncol=2)
    f <- fisher.test(as.table(mat), alt="greater")
    return(f$p.value)
}

fishers <- unlist(lapply(1:nrow(celltype_dmr_overlaps),row_fisher))
celltype_dmr_overlaps$fishers.pval<-NA
celltype_dmr_overlaps$fishers.pval<-fishers
celltype_dmr_overlaps<-as.data.frame(celltype_dmr_overlaps)

dat <- reshape2::dcast(celltype_dmr_overlaps, cluster~celltype,value.var="fishers.pval",fill=1) 

row.names(dat)<-dat$cluster

#negative log10 on dat
dat<-dat[2:ncol(dat)]
dat<- -(log10(dat))
col_fun=colorRamp2(breaks=quantile(unlist(dat),probs=c(0,0.5,0.9)),colors=c("white","grey","#FF00FF"))
plt<-Heatmap(dat,col=col_fun)
pdf(paste0("06_immune_finecelltyping.dmr_rnamarkergene.fisher.pdf"))
print(plt)
dev.off()

```

# Use Liger to integrate DMRs and RNA markers for immune cells

```R
library(rliger)

dat<-readRDS(file="06_scaledcis.immune_finecelltyping.amethyst.rds")
collapsed_dmrs<-readRDS(file=paste0("06_immune_finecelltyping",".500bp_dmrs.filt_collapsed.rds"))

#select dmrs for integration (same way we selected for dmr and gene overlap)
dmrs_for_integration <- collapsed_dmrs %>% 
                          filter(direction=="hypo") %>%
                          filter(dmr_padj < 0.05) %>%
                          filter(abs(dmr_logFC) > 1.5) %>% 
                          filter(dmr_length<50000) %>% 
                          filter(gene_names != "NA") %>%
                          group_by(type) 
table(dmrs_for_integration$type)

genes_for_integration <- dmrs_for_integration %>% 
                          as.data.frame() %>%
                          select(gene_names) %>% 
                          unlist %>%
                          strsplit(split=",") %>% 
                          unlist() %>% 
                          stringr::str_replace_all(" ","") %>% 
                          unique() %>%
                          sort()

#reduce overlapping dmrs regardless of cluster
dmr_out<-dmrs_for_integration %>% 
        as.data.frame() %>% 
        select(chr,dmr_start,dmr_end) %>% 
        distinct(chr,dmr_start,dmr_end) %>% 
        makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>%
        reduce()

dmr_out<-dmr_out[width(dmr_out)<50000,]
summary(width(dmr_out))

#assign genes to merged dmr_out features
reduced_overlap<-findOverlaps(dmr_out,
            makeGRangesFromDataFrame(dmrs_for_integration,keep.extra.columns=TRUE))

#add list of gene names and store original dmr sites in bed
dmr_out$gene_names<-NA
dmr_out$gene_names<-mclapply(unique(queryHits(reduced_overlap)), function(x){
    genes=dmrs_for_integration[subjectHits(reduced_overlap[queryHits(reduced_overlap)==x,]),]$gene_names
    genes <- genes %>% strsplit(split=",") %>% unlist() %>% stringr::str_replace_all(" ","") %>% unique() %>% paste(collapse=",")
    genes<-paste(unlist(genes),sep=",")
    dmr_out[x,]$gene_names<-genes
},mc.cores=100)
dmr_out$met_feature<-paste(seqnames(dmr_out),start(dmr_out),end(dmr_out),sep="_")

#output bed format for making windows
dmr_bed<-data.frame(chr=seqnames(dmr_out),start=start(dmr_out),end=end(dmr_out))

#create new matrix from DMR sites for refined clustering
dat@genomeMatrices[["immune_cluster_dmr_sites"]] <- makeWindows(dat, 
                                                     bed = dmr_bed,
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 50, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 
saveRDS(dat,file="06_scaledcis.immune_finecelltyping.amethyst.rds")

#add filter to DMRs with atleast 10% coverage per cell
met_features<-which(rowSums(!is.na(dat@genomeMatrices[["immune_cluster_dmr_sites"]]))>(ncol(dat@genomeMatrices[["immune_cluster_dmr_sites"]])*0.1))

#expand dmr_out so each row is its own gene name
#this will be used to summarize multiple DMRs to their genes
dmrs_to_gene_names<-lapply(1:length(dmr_out), function(x) {
  temp<-dmr_out[x,]
  genes<-unlist(temp$gene_names %>% unlist() %>% strsplit(split=",")) 
  return(cbind(genes=genes,met_feature=temp$met_feature))
})

dmrs_to_gene_key<-as.data.frame(do.call("rbind",dmrs_to_gene_names))

#set new matrix with average DMR score per gene using the key above
gene_methylation<-mclapply(unique(dmrs_to_gene_key$genes), function(i){
  gene_name=i
  met_features<-dmrs_to_gene_key[dmrs_to_gene_key$genes==gene_name,2]
  return(colMeans(dat@genomeMatrices[["immune_cluster_dmr_sites"]][row.names(dat@genomeMatrices[["immune_cluster_dmr_sites"]]) %in% met_features,],na.rm=TRUE))
},mc.cores=100)

gene_met<-as.data.frame(do.call("rbind",gene_methylation))
row.names(gene_met)<-unique(dmrs_to_gene_key$genes)

# Create from raw score matrices MET
met_dat <- Matrix(as.matrix(gene_met), sparse = TRUE) # Convert the dense matrix to a sparse dgCMatrix

#set NA values to 0 (following amethyst example)
met_dat[which(is.na(met_dat),arr.ind=TRUE)]<-0

system(paste0("mkdir -p ",project_data_directory,"/integration_immune/"))
saveRDS(met_dat,file=paste0(project_data_directory,"/integration_immune/","scaledcis.gene_dmr_methylation.immune.mat.rds"))

table(rna_markers_filt$cluster)

# RNA markers continued from previous section
#get raw RNA values on filtered methylation genes
rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")
rna<-subset(rna,coarse_celltype %in% c("bcell","myeloid","plasma","tcell"))
Idents(rna)<-rna$fine_celltype
rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)
rna <- ScaleData(rna)
rna <- JoinLayers(rna)

#overexpressed genes by fine cell type
rna_markers<-FindAllMarkers(rna,assay="RNA",only.pos=TRUE)
rna_markers_filt<-rna_markers %>% 
                    filter(avg_log2FC > 3) %>% 
                    filter(p_val_adj <0.01) %>% 
                    group_by(cluster)


#filter met_dat to only genes in rna markers
met_dat<-met_dat[row.names(met_dat) %in% unique(rna_markers_filt$gene),]


rna_dat<-LayerData(rna,features=row.names(met_dat),assay="RNA",layer="counts")
saveRDS(rna_dat,file=paste0(project_data_directory,"/integration_immune/","scaledcis.gene_rna.mat.immune.rds"))

rna$broad_celltype<-rna$fine_celltype
rna$fine_cluster_id <- "NA"
rna$mcg_pct <- "NA"

#create shared metadata, add coarse_cluster_id to meta
met_dat<-met_dat[,order(colnames(met_dat))]
meta<-rbind(dat@metadata[colnames(met_dat),c("sample","broad_celltype","fine_cluster_id","mcg_pct")],rna@meta.data[c("sample","broad_celltype","fine_cluster_id","mcg_pct")])

summary(unlist(as.data.frame(met_dat)))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-1.00000  0.00000  0.00000  0.06306  0.00000  1.00000 

summary(unlist(as.data.frame(rna_dat)))
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#    0.000     0.000     0.000     0.374     0.000 15099.000 

met_dat<-met_dat[row.names(met_dat)%in% row.names(rna_dat),]

#make liger object
met.liger <- createLiger(rawData=list(met=met_dat,rna=rna_dat), 
                        modal=c("meth","rna"),cellMeta=meta,
                        removeMissing=TRUE)

#now follow liger for integration
rna.met <- met.liger %>% 
            rliger::normalize() %>%
            rliger::selectGenes(useDatasets = "rna",nGenes=nrow(met_dat)) %>%
            rliger::scaleNotCenter()



#try https://welch-lab.github.io/liger/reference/runUINMF.html on unshared features as well

summary(unlist(as.data.frame(rna.met@datasets$met@scaleData)))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  1.0000  1.0000  0.9002  1.0000  2.0000 

summary(unlist(as.data.frame(rna.met@datasets$rna@scaleData)))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.000   0.000   0.000   0.116   0.000 258.780 

rna.met <- rliger::runIntegration(rna.met, k = 20)
rna.met <- rliger::quantileNorm(rna.met)
rna.met <- runCluster(rna.met)
rna.met <- runUMAP(rna.met)

plt1<-plotDatasetDimRed(rna.met,splitBy="dataset")
ggsave(wrap_plots(plt1),file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","dataset",".integrated.umap.pdf"),width=20,height=10)

plt2<-plotClusterDimRed(rna.met,splitBy="dataset","broad_celltype")
ggsave(wrap_plots(plt2),file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","celltype",".integrated.umap.pdf"),width=20,height=10)

plt3<-plotClusterDimRed(rna.met,splitBy="dataset","sample")
ggsave(wrap_plots(plt3),file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","sample",".integrated.umap.pdf"),width=20,height=10)

plt4<-plotClusterDimRed(rna.met,splitBy="dataset","fine_cluster_id")
ggsave(wrap_plots(plt4),file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","fine_cluster_id",".integrated.umap.pdf"),width=20,height=10)

#plt5<-plotClusterDimRed(rna.met,slot="cellMeta","mcg_pct")
#ggsave(wrap_plots(plt5),file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","mcg",".integrated.umap.pdf"),width=20,height=10)

umap_dims<-rna.met@dimReds$UMAP
colnames(umap_dims)<-c("X","Y")
umap_dims<-umap_dims[paste0("met_",colnames(met_dat)),]
plt<-ggplot(umap_dims,aes(x=X,y=Y,color=rna.met@cellMeta[row.names(umap_dims),]$fine_cluster_id))+geom_point()+theme_void()
ggsave(plt,file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","fine_cluster_id",".integrated.umap.pdf"),width=20,height=10)

plots <- plotGeneDimRed(rna.met, c("PDGFB", "IL3RA"), splitBy = "dataset",
                        titles = c(names(rna.met), names(rna.met)))
plt<-cowplot::plot_grid(plotlist = plots, nrow = 2)
ggsave(plt,file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","fine_cluster_id",".integrated.umap.pdf"),width=20,height=10)

saveRDS(rna.met,file=paste0(project_data_directory,"/integration_immune/","immune.met_rna.integrated.liger.rds"))

rna.met@cellMeta[rna.met@cellMeta$dataset=="met",]$fine_cluster_id == dat@metadata[colnames(met_dat),c("fine_cluster_id")]
gsub("met_","",row.names(rna.met@cellMeta[rna.met@cellMeta$dataset=="met",])) == row.names(dat@metadata[colnames(met_dat),])
colnames(met_dat)==gsub("met_","",row.names(rna.met@cellMeta[rna.met@cellMeta$dataset=="met",])) 
```


Promoter regions on marker genes
Also doesnt work particularly well.

```R


rna_dat<-LayerData(rna,features=unique(rna_markers_filt$gene),assay="RNA",layer="counts")
saveRDS(rna_dat,file=paste0(project_data_directory,"/integration_immune/","scaledcis.gene_rna.mat.immune.rds"))

rna$broad_celltype<-rna$fine_celltype
rna$fine_cluster_id <- "NA"
rna$mcg_pct <- "NA"

#create new matrix from DMR sites for refined clustering
promoter_met <- makeWindows(dat, 
                    genes = unique(rna_markers_filt$gene),
                    promoter=TRUE,
                    type = "CG", 
                    metric = "score", 
                    threads = 50, 
                    index = "chr_cg", 
                    nmin = 2) 

#get the perc matrix before progressing

promoter_met2<-promoter_met[which(rowSums(!is.na(promoter_met))>(ncol(promoter_met)*0.1)),]

rna_dat<-rna_dat[row.names(promoter_met2),]

promoter_met <- Matrix(as.matrix(promoter_met2), sparse = TRUE) # Convert the dense matrix to a sparse dgCMatrix
promoter_met[which(is.na(promoter_met),arr.ind=TRUE)]<-0

#make liger object
met.liger <- createLiger(rawData=list(met=promoter_met,rna=rna_dat), 
                        modal=c("meth","rna"),cellMeta=meta,
                        removeMissing=TRUE)

#now follow liger for integration
rna.met <- met.liger %>% 
            rliger::normalize() %>%
            rliger::selectGenes(useDatasets = "rna") %>%
            rliger::scaleNotCenter()

rna.met <- rliger::runIntegration(rna.met, k = 20)
rna.met <- rliger::quantileNorm(rna.met)
rna.met <- runCluster(rna.met)
rna.met <- runUMAP(rna.met)

plt1<-plotDatasetDimRed(rna.met,splitBy="dataset")
ggsave(wrap_plots(plt1),file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","dataset",".integrated.umap.pdf"),width=20,height=10)

plt2<-plotClusterDimRed(rna.met,splitBy="dataset","broad_celltype")
ggsave(wrap_plots(plt2),file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","celltype",".integrated.umap.pdf"),width=20,height=10)

plt3<-plotClusterDimRed(rna.met,splitBy="dataset","sample")
ggsave(wrap_plots(plt3),file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","sample",".integrated.umap.pdf"),width=20,height=10)

plt4<-plotClusterDimRed(rna.met,splitBy="dataset","fine_cluster_id")
ggsave(wrap_plots(plt4),file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","fine_cluster_id",".integrated.umap.pdf"),width=20,height=10)


```
Trying with bindSC as well
https://htmlpreview.github.io/?https://github.com/KChen-lab/bindSC/blob/master/vignettes/mouse_retina/retina.html
https://htmlpreview.github.io/?https://github.com/KChen-lab/bindSC/blob/master/vignettes/SC_ST/SC_ST.html 
```R
#devtools::install_github('KChen-lab/bindSC')
library(bindSC)




```
old stuff

```R

#set order by umap groupings
order<-c("2","4","8","6","5","1","7","3","9")
cell_markers<-cell_markers
plt_list<-lapply(names(cell_markers),
    function(celltype){
            genes<-unlist(cell_markers[celltype])
            print(paste("Plotting:",celltype))
            print(paste("Genes to plot:",genes))
            plt<-histograModified(obj=dat_fine, 
                baseline="mean",
                genes = unlist(cell_markers[celltype]),
                order=order, promoter_focus=TRUE,trackOverhang=10000,
                matrix = paste0("cg_",broad_celltype,".clusters_tracks"), arrowScale = .03, trackScale = .5,
                legend = F, cgisland=cgi) + ggtitle(celltype)
            ggsave(plt,
            file=paste(project_data_directory,"fine_celltyping",broad_celltype,
            paste("04_scaledcis",broad_celltype,celltype,"broad_celltype.markers.pdf",sep="."),
            sep="/"),
            width=length(genes)*2,
            height=length(order)*2,
            limitsize=F)
        })


immune_markers<-list()
immune_markers[["myeloid_mono_classical"]]<-c('SERPINB2', 'AQP9', 'FCN1', 'VCAN', 'AC245128.3', 'EREG', 'TNIP3', 'CD300E', 'S100A8', 'S100A9', 'THBS1', 'CYP1B1', 'SLC39A8', 'MAP3K20', 'EHD1', 'LUCAT1', 'FCAR', 'S100A12', 'PID1', 'ACOD1')
immune_markers[["myeloid_mono_nonclassical"]]<-c('CCDC68', 'CLEC12A', 'CD48', 'OLIG1', 'CDKN1C', 'S100A4', 'KRTAP2-3', 'CD52', 'WARS', 'LINC02432', 'CFP', 'UNC45B', 'APOBEC3A', 'TKT', 'RIPOR2', 'CDA', 'PLAC8', 'FFAR2', 'PAPSS2', 'GK5')
immune_markers[["myeloid_macro_C3"]]<-c('C3', 'GPR158', 'C1QC', 'C1QB', 'AC079760.2', 'FCGR3A', 'SDS', 'A2M', 'AL590705.1', 'MT3', 'TREM2', 'HLA-DQB2', 'OLR1', 'ZBED9', 'IL4I1', 'DOK5', 'IGLC7', 'PMEPA1', 'IBSP', 'CD9')
immune_markers[["myeloid_macro_LYVE1"]]<-c('PDK4', 'SELENOP', 'FOLR2', 'LYVE1', 'LILRB5', 'MAF', 'KLF4', 'RNASE1', 'SLC40A1', 'RHOB', 'F13A1', 'MRC1', 'TSC22D3', 'PLTP', 'STAB1', 'NR4A2', 'MS4A4A', 'DAB2', 'CD163', 'CTSC')
immune_markers[["myeloid_macro_CCL4"]]<-c('CCL4', 'SELENOP', 'FOLR2', 'RNASE1', 'MRC1', 'F13A1', 'LYVE1', 'CD163', 'STAB1', 'CXCL3', 'HMOX1', 'CCL3', 'CXCL2', 'AC010980.2', 'GCLM', 'CD93', 'CTSL', 'PLTP', 'MAFB', 'THBD')
immune_markers[["myeloid_macro_APOC1"]]<-c('ACP5', 'CYP27A1', 'LIPA', 'PLD3', 'TREM2', 'PLA2G7', 'CD36', 'FABP4', 'GPNMB', 'CCL18', 'NCEH1', 'RARRES1', 'FUCA1', 'APOC1', 'LPL', 'OTOA', 'BLVRB', 'HS3ST2', 'NR1H3', 'CHIT1')
immune_markers[["myeloid_macro_interferon"]]<-c( 'EPSTI1', 'CXCL10', 'GBP1', 'STAT1', 'CXCL11', 'CXCL9', 'GBP4', 'IFI44L', 'FCGR1A', 'SERPING1', 'IFI6', 'IGSF6', 'IFIH1', 'IDO1', 'RASSF4', 'XAF1', 'APOL6', 'ISG15', 'FYB1', 'PSME2')
immune_markers[["myeloid_macro_stress"]]<-c( 'HSPA1B', 'BAG3', 'DNAJB1', 'HSPA6', 'MRPL18', 'ZFAND2A', 'CHORDC1', 'SDS', 'CACYBP', 'FKBP4', 'DNAJA4', 'HSPE1', 'RGS1', 'HSPH1', 'HSPA1A', 'OLR1', 'HSPB1', 'CLEC2B', 'CKS2', 'PLEK')
immune_markers[["myeloid_macro_CTSK"]]<-c( 'AP000904.1', 'CTSK', 'ACP5', 'SLC9B2', 'MMP9', 'DGKI', 'SPP1', 'TCIRG1', 'ATP6V0D2', 'CKB', 'CST3', 'CD109', 'TIMP2', 'JDP2', 'SPARC', 'C9orf16', 'GRN', 'SIGLEC15', 'AK5', 'SNX10')
immune_markers[["myeloid_cycling"]]<-c( 'MKI67', 'PCLAF', 'CKS1B', 'DIAPH3', 'TYMS', 'TOP2A', 'HIST1H1B', 'ASPM', 'TK1', 'DHFR', 'NUSAP1', 'NCAPG', 'DLGAP5', 'CENPM', 'TPX2', 'RRM2', 'CEP55', 'MYBL2', 'BIRC5', 'CENPF')
immune_markers[["myeloid_neutrophil"]]<-c('CRISP3', 'ANKRD34B', 'ABCA13', 'CD177', 'ORM1', 'SERPINB4', 'DEFA3', 'DPYS', 'C7', 'AZU1', 'CLC', 'FOXQ1', 'RAB6B', 'SFRP1', 'TWIST1', 'MAP1B', 'SLPI', 'ADH1B', 'FOLR3', 'BTNL9')
immune_markers[["myeloid_mast"]]<-c('CPA3', 'TPSAB1', 'MS4A2', 'TPSB2', 'HDC', 'IL1RL1', 'GATA2', 'HPGDS', 'GCSAML', 'HPGD', 'ADCYAP1', 'KIT', 'KRT1', 'PRG2', 'CTSG', 'CLU', 'TPSD1', 'CALB2', 'RAB27B', 'SLC18A2')
immune_markers[["myeloid_cDC1"]]<-c('CLEC9A', 'XCR1', 'DNASE1L3', 'IDO1', 'LGALS2', 'CPNE3', 'WDFY4', 'C1orf54', 'DAPP1', 'RAB11FIP1', 'PPA1', 'CPVL', 'C9orf135', 'FOXB1', 'AC016717.2', 'CST3', 'PLEKHA5', 'GPR157', 'FOXN2', 'NCOA7')
immune_markers[["myeloid_cDC2"]]<-c('FCER1A', 'IL1R2', 'CLEC10A', 'CD1C', 'CST7', 'GPAT3', 'DAPP1', 'CFP', 'EREG', 'CCL22', 'LGALS2', 'JAML', 'PID1', 'AREG', 'IL7R', 'AC020656.1', 'IL1R1', 'MARCKSL1', 'ADAM19', 'SLC7A11')
immune_markers[["myeloid_mDC"]]<-c('LAMP3', 'BIRC3', 'LACRT', 'NUB1', 'CCR7', 'MARCKSL1', 'IDO1', 'DAPP1', 'POGLUT1', 'LINC01539', 'GPR157', 'IL12B', 'LAD1', 'KIF2A', 'FSCN1', 'IL7R', 'TXN', 'DUSP5', 'FOXD4L1', 'CD200', 'RAB9A')
immune_markers[["myeloid_pDC"]]<-c('CD5', 'LYPD8', 'PRL', 'AC136475.3', 'LINC01087', 'AC006058.1', 'BDKRB2', 'POTEI', 'DSP', 'AC026369.3', 'GZMB', 'TNFSF4', 'CD2', 'TIGIT', 'LTB', 'TMEM45A', 'PGR', 'CD3G', 'AC015936.1', 'ACOT7')

immune_markers[["tcell_cd4_memory"]]<-c('ADAM23', 'NEFL', 'LINC02273', 'ANTXR2', 'MFHAS1', 'AP3M2', 'GPR183', 'PASK', 'S1PR1', 'FTH1', 'GLIPR1', 'SESN3', 'IL7R', 'CCR7', 'SLC2A3', 'DDIT4', 'CD28', 'ANXA1', 'TRAT1', 'KLF2')
immune_markers[["tcell_cd4_stress"]]<-c('MYADM', 'LMNA', 'ANXA1', 'KLF6', 'HSP90AA1', 'RGCC', 'FAM107B', 'AHNAK', 'HSPH1', 'VIM', 'HSPA8', 'GPR183', 'HSP90AB1', 'S100A10', 'S100A11', 'EZR', 'HSPD1', 'IL7R', 'ANKRD12', 'TUBB4B')
immune_markers[["tcell_cd4_naive"]]<-c('CA6', 'ADTRP', 'GTSCR1', 'LINC00402', 'LINC01481', 'LEF1', 'MAL', 'SELL', 'AL391097.1', 'ANKRD55', 'NEXMIF', 'CCR7', 'LINC00861', 'TESPA1', 'ANK1', 'LTB', 'SESN3', 'TSHZ2', 'PASK', 'LDHB')
immune_markers[["tcell_treg"]]<-c('CD177', 'FANK1', 'LINC02099', 'IL1R2', 'CCR8', 'FOXP3', 'PTGIR', 'LAYN', 'IKZF4', 'F5', 'RTKN2', 'CARD16', 'MAGEH1', 'IL2RA', 'LINC01943', 'BATF', 'IKZF2', 'HTATIP2', 'GK', 'CTLA4')
immune_markers[["tcell_cd4_Tfh"]]<-c('IGFL2', 'CXCL13', 'CPM', 'CD200', 'KSR2', 'NMB', 'HMSD', 'TUBA3D', 'PTPN13', 'GNG4', 'SERPINE2', 'AC008011.2', 'HMHB1', 'LINC01229', 'TSHZ2', 'GK', 'PGM2L1', 'SESN3', 'MFHAS1', 'FKBP5')
immune_markers[["tcell_cd8_TEM"]]<-c('DKK3', 'ENC1', 'GZMK', 'EOMES', 'YBX3', 'GGA2', 'CMC1', 'CD8A', 'CST7', 'CD8B', 'SH2D1A', 'HLA-DPB1', 'CRTAM', 'DUSP2', 'LYST', 'COTL1', 'ZEB2', 'PIK3R1', 'RNF19A', 'TRAT1')
immune_markers[["tcell_cd8_TRMZNF683"]]<-c('CSN2', 'AL590434.1', 'CCL4', 'OR5AU1', 'FGF10-AS1', 'CCL4L2', 'TNFSF9', 'AL133163.2', 'DUSP6', 'FOSB', 'AC020916.1', 'IER5L', 'MZF1-AS1', 'NR4A2', 'IFNG', 'RASGEF1B', 'PLCH1', 'POLQ', 'GZMA', 'GZMK')
immune_markers[["tcell_cd8_TRMAUTS"]]<-c('AUTS2', 'PRSS3', 'GFPT2', 'DAPK2', 'PAGE5', 'ITGA1', 'KLRC1', 'IGFBP3', 'AC020571.1', 'HPSE2', 'CTH', 'COL25A1', 'LINC01871', 'FGF9', 'PARD6G', 'SPRY1', 'SNTG2', 'CD8A', 'ADAMTS6', 'TAC1')
immune_markers[["tcell_cd8_stress"]]<-c( 'HSPA1A', 'DNAJB1', 'HSPH1', 'HSP90AA1', 'HSPB1', 'HSPE1', 'DNAJA1', 'TUBA4A', 'HSP90AB1', 'HSPD1', 'CDKN1C', 'B3GNT5', 'CSF1', 'HSPA8', 'KLF6', 'EGR1', 'POU3F1', 'ANKRD37', 'ZNF80', 'MIR222HG')
immune_markers[["tcell_gdT"]]<-c('KIR2DL3', 'KLRC2', 'KIR2DL1', 'KIR3DL1', 'CHPT1', 'SMC4', 'TRGC2', 'KLRD1', 'XCL1', 'CEMIP2', 'CD7', 'METRNL', 'IKZF2', 'AHI1', 'CTSW', 'KIR2DL4', 'CAST', 'LINC01871', 'ID3', 'PRKX')
immune_markers[["tcell_cycling"]]<-c('DLGAP5', 'E2F8', 'HJURP', 'SPC25', 'RRM2', 'GTSE1', 'AC007240.1', 'BIRC5', 'MYBL2', 'MKI67', 'MCM10', 'KIF4A', 'CENPA', 'SPC24', 'CDC45', 'PCLAF', 'ESCO2', 'KIF20A', 'UBE2C', 'DIAPH3')
immune_markers[["tcell_cd8_interferon"]]<-c( 'IFIT3', 'IFIT1', 'IFI44L', 'HERC5', 'MX2', 'MX1', 'IFI6', 'XAF1', 'ISG15', 'STAT1', 'EIF2AK2', 'LY6E', 'EPSTI1', 'IFI44', 'PARP14', 'OASL', 'OAS2', 'RSAD2', 'SAMD9L', 'MT2A')
immune_markers[["tcell_cd8_exhausted"]]<-c( 'REG4', 'RYR2', 'HAVCR2', 'NPW', 'VCAM1', 'ETV1', 'PNOC', 'KRT86', 'TNS3', 'MYO7A', 'PHEX', 'LRRN3', 'AC243829.4', 'FNDC9', 'BUB1B', 'BEND4', 'AKAP5', 'ASB2', 'NHS', 'RGS13', 'ENTPD1' )
immune_markers[["tcell_nk_FCGR3A"]]<-c('MGAM', 'MYOM2', 'FCGR3A', 'KLRF1', 'PRSS57', 'MLC1', 'AKR1C3', 'FGFBP2', 'SPON2', 'RBPMS2', 'FCER1G', 'KIR3DL1', 'CX3CR1', 'KIR2DL1', 'SH2D1B', 'S1PR5', 'CHST2', 'CCL3', 'CFAP97D2', 'KIR2DL3')
immune_markers[["tcell_nk_cytotoxic"]]<-c('LILRB1', 'CES1', 'PROK2', 'FGFBP2', 'ASCL2', 'CX3CR1', 'FCRL6', 'LINC02384', 'GZMH', 'CFAP97D2', 'LINC01936', 'FCGR3A', 'S1PR5', 'AL590434.1', 'PLEK', 'NKG7', 'ADGRG1', 'AC011933.3', 'PRSS23', 'KIR2DL1')
immune_markers[["tcell_ilc1"]]<-c('ADGRG3', 'CCNJL', 'LINC00996', 'KIR2DL4', 'CSF2', 'B3GNT7', 'FAM166B', 'FES', 'SCT', 'SH2D1B', 'CFAP46', 'AL157895.1', 'KRT86', 'ATP8B4', 'FCER1G', 'GNLY', 'ZNF660', 'TYROBP', 'GOLIM4', 'HOXA6')
immune_markers[["tcell_ilc23"]]<-c('CD300LF', 'KIT', 'SOST', 'LINC00299', 'HOXA6', 'PLA2G4A', 'NPTX2', 'VSTM2L', 'FAM30A', 'SCT', 'XCL1', 'SOX4', 'XCL2', 'AREG', 'TNFRSF18', 'SCN1B', 'TNFSF4', 'MAP3K8', 'SSBP2', 'PLCH1')
immune_markers[["tcell_MAIT"]]<-c('KLRB1', 'RORA', 'KISS1', 'CEBPD', 'TRDV2', 'AC022392.1', 'AC026461.3', 'FEZ1', 'REL', 'GPR65', 'MYBL1', 'LTB', 'MT3', 'SLC4A10', 'NFKB1', 'IGHV3-23', 'CCL20', 'RORC', 'NCR3', 'DPP4')
immune_markers[["tcell_GZMK"]]<-c('ZMAT4', 'KRT13', 'COL24A1', 'CD69', 'COX6A2', 'AL157895.1', 'CCNJL', 'SOST', 'MLC1', 'ADGRG3', 'SMIM24', 'HOXA6', 'SH2D1B', 'TYROBP', 'XCL2', 'GNLY', 'TAL1', 'XCL1', 'AREG', 'FCER1G')

immune_markers[["bcell_naive"]]<-c('IGHD', 'YBX3', 'TCL1A', 'SPRY1', 'CD83', 'BCL7A', 'TXNIP', 'FCER2', 'CLEC2D', 'HLA-DQA2', 'MYB', 'KLF2', 'ESCO2', 'STAG3', 'HLA-DRB5', 'GPR18', 'CMSS1', 'ACSM3', 'GBP4')
immune_markers[["bcell_memory"]]<-c('REL', 'CD83', 'LTB', 'CD69', 'GPR183', 'CCR7', 'CD48', 'ID3', 'BCL2A1', 'KYNU', 'HLA-DRB5', 'PIKFYVE', 'SAMSN1', 'IER5', 'CD80', 'CHORDC1', 'SCIMP', 'DUSP2', 'KLF6', 'NFKB1')
immune_markers[["bcell_interferon"]]<-c('CCL17', 'STAT1', 'IFI44L', 'GBP4', 'ISG15', 'SERPINA9', 'MX1', 'AICDA', 'DUSP4', 'XAF1', 'GBP1', 'ZBED2', 'GBP5', 'GBP2', 'AIM2', 'SAMD9L', 'IFNG', 'MX2', 'ITGAX', 'PLIN2')
immune_markers[["bcell_stress"]]<-c('HSPA1B', 'RRBP1', 'DNAJB1', 'HSPA1A', 'FOS', 'ANKRD28', 'SAMD7', 'H1FX', 'IGKV1-37', 'RGS1', 'HSPB1', 'AC026369.3', 'AMPD1', 'IGKV2-26', 'FAM92B', 'AC021074.3', 'IGHV3-41', 'LINC01405', 'CPEB4', 'IGKV2-4')

immune_markers[["plasma_IGA"]]<-c('JCHAIN', 'IGHA2', 'IGHA1', 'LGALSL', 'NCAM1', 'CPVL', 'CCR10', 'AC106028.4', 'ZNF215', 'IGHV7-81', 'IGKV2-28', 'CADPS2', 'LHX8', 'IGLV2-18', 'IGHV3OR16-8', 'OTP', 'CTSW', 'AC136428.1', 'IGHV3-72')
immune_markers[["plasma_IGG"]]<-c('MZB1', 'IGHG3', 'DERL3', 'FKBP11', 'IGHG4', 'ITM2C', 'IGHG1', 'PRDX4', 'SEC11C', 'XBP1', 'SSR4', 'SSR3', 'JSRP1', 'ERLEC1', 'CYTOR', 'IGKV3D-11', 'SELENOS', 'DPEP1', 'CD38', 'MYDGF')


#set order by umap groupings
cell_markers<-immune_markers
plt_list<-lapply(names(cell_markers),
    function(celltype){
            genes<-unlist(cell_markers[celltype])
            print(paste("Plotting:",celltype))
            print(paste("Genes to plot:",genes))
            plt<-histograModified(obj=dat_fine, 
                baseline="mean",
                genes = unlist(cell_markers[celltype]),
                order=order, promoter_focus=TRUE,trackOverhang=10000,
                matrix = paste0("cg_",broad_celltype,".clusters_tracks"), arrowScale = .03, trackScale = .5,
                legend = F, cgisland=cgi) + ggtitle(celltype)
            ggsave(plt,
            file=paste(project_data_directory,"fine_celltyping",broad_celltype,
            paste("04_scaledcis",broad_celltype,celltype,"fine_celltype.markers.pdf",sep="."),
            sep="/"),
            width=length(genes)*2,
            height=length(order)*2,
            limitsize=F)
        })

#dmr set
dmr<-readRDS(
        paste(project_data_directory,"fine_celltyping",broad_celltype,
        paste(broad_celltype,"dmr",broad_celltype,"clusters.collapsed.rds",sep="."),
        sep="/"))

#rna markers defined by paired 10x RNA analysis
rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.rds")
table(rna$broad_celltype)
rna<-subset(rna,broad_celltype %in% c("tcell","bcell","plasma","myeloid")) #using all stromal types
rna<-subset(rna,fine_celltype %in% c("suspected_doublet"),invert=TRUE)
table(rna$fine_celltype)

rna<-JoinLayers(rna)
markers<-FindAllMarkers(rna,group.by="fine_celltype")

#plot overlap of overexpressed RNA markers for fine cell types and our cluster hypo DMRs
plot_dmr_RNAmarker_overlap(dmr=dmr,
                            rna_markers=markers,
                            celltype=broad_celltype,
                            clus_order=order,
                            plot_value="pval")



#grab cluster dmr hypo markers
clus_markers<-dmr %>% filter(celltype=="7") %>% filter(direction=="hypo") %>% slice_min(dmr_logFC,n=100)
clus_markers<-unlist(
        lapply(
            unlist(
                strsplit(
                    unlist(clus_markers$gene_names),",")),
                function(x) gsub("[[:space:]]", "", x)
            )
        )
#overlap with RNA markers
markers %>% filter(gene %in% clus_markers) %>% filter(avg_log2FC>2) %>% select(cluster) %>% table()



######FINAL IMMUNE CELL TYPING###########
dat_fine@metadata$fine_celltypes<-"tcell"
dat_fine@metadata[dat_fine@metadata[[paste(broad_celltype,"clusters",sep=".")]] %in% c("9"),]$fine_celltypes<-"bcell"
dat_fine@metadata[dat_fine@metadata[[paste(broad_celltype,"clusters",sep=".")]] %in% c("3","7"),]$fine_celltypes<-"myeloid1"
dat_fine@metadata[dat_fine@metadata[[paste(broad_celltype,"clusters",sep=".")]] %in% c("1","8"),]$fine_celltypes<-"myeloid2"


#summarize celltype-dmr clusters over windows
dat_fine<-dmr_and_1kb_window_gen(dat_fine,
    prefix=paste(output_directory,broad_celltype,sep="/"),
    groupBy="fine_celltypes",
    threads=10,step=500)

print(paste("Saving final amethyst object of",broad_celltype))
saveRDS(dat_fine,paste0(output_directory,"/",paste("04_scaledcis",broad_celltype,"fine_celltype.amethyst.rds",sep=".")))



#fix rna marker overlap (i think heatmap doesnt work, but pulling with dplyr does?)
#fix atac overlap of markers





# #plot with atac peak overlaps
# atac_markers<-read.table("/data/rmulqueen/projects/scalebio_dcis/ref/CEDAR_multiome.celltypes.markers.peaks.df.tsv",header=T)
# atac_markers$chr<-unlist(lapply(strsplit(atac_markers$feature,"-"),"[",1))
# atac_markers$start<-unlist(lapply(strsplit(atac_markers$feature,"-"),"[",2))
# atac_markers$end<-unlist(lapply(strsplit(atac_markers$feature,"-"),"[",3))
# clus_order<-order

# plot_dmr_ATACmarker_overlap<-function(dmr,atac_markers,celltype,clus_order,plot_value="Jaccard"){
#     #plot_value can be one of pval, odds.ratio, intersection, union, Jaccard
#     #flow through:
#     # - filter DMRs to HYPO
#     # - filter RNA markers to overexpressed
#     # - table dmr to marker count, and plot odds ratio
#     # - plot gene body methylation hypoM plots
#     output_directory=paste0(project_data_directory,"/fine_celltyping/",celltype)

#     #PLOT HYPO
#     #filter to hypo DMRs (more likely to increase expression)
#     dmr_hypo<- dmr %>% group_by(test) %>% filter(direction=="hypo") 
#     dmr_granges<-makeGRangesFromDataFrame(dmr_hypo,keep.extra.columns=TRUE)
#     dmr_granges<-split(dmr_granges,dmr_granges$test)

#     atac_markers <- atac_markers %>% filter(logFC>0) #look only for increased expression markers
#     atac_markers_granges<-makeGRangesFromDataFrame(atac_markers,keep.extra.columns=TRUE)
#     atac_markers_granges<-split(atac_markers_granges,atac_markers_granges$group)

#     overlap<-as.data.frame(do.call("rbind",lapply(1:length(dmr_granges),calc_overlap_stats)))
#     colnames(overlap)<-c("cluster","celltype","odds.ratio","pval","jaccard")

#     dmr_count<-dmr_hypo %>% group_by(test,.drop=FALSE) %>% tally() 
#     dmr_count<-setNames(dmr_count$n,nm=dmr_count$test)

#     atac_count<-atac_markers %>% group_by(group,.drop=FALSE) %>% tally() 
#     atac_count<-setNames(atac_count$n,nm=atac_count$group)

#     or_mat<-as.data.frame(data.table::dcast(data=as.data.table(overlap),formula=cluster~celltype,value.var=plot_value))
#     row.names(or_mat)<-or_mat$cluster
#     or_mat<-or_mat[,2:ncol(or_mat)]
#     or_mat <- or_mat %>%
#                 mutate_all(~as.numeric(as.character(.)))    
    
#     print("Plotting hypo overlap of DMRs...")
#     pdf(paste0(output_directory,"/","04_scaledcis.",celltype,".dmr_ATAC_overlap.heatmap.pdf"))
#     print(
#         Heatmap(or_mat,
#             bottom_annotation=HeatmapAnnotation(bar=anno_barplot(atac_count)),
#             row_order=clus_order,cluster_rows=FALSE,
#             show_row_names = TRUE, show_column_names = TRUE, row_names_side="right",column_names_side="bottom",
#             column_title = "Hypomethylated and ATAC Peaks",
#             left_annotation=rowAnnotation(bar=anno_barplot(dmr_count)),
#             name=plot_value,
#             cell_fun = function(j, i, x, y, width, height, fill) {
#                                 grid.text(format(round(or_mat,digits=2),nsmall=2)[i, j], x, y, gp = gpar(fontsize = 10))
#                                 }
#             ))
#     dev.off()
# }


# plot_dmr_ATACmarker_overlap(dmr=dmr,
#                             rna_markers=markers,
#                             celltype=broad_celltype,
#                             plot_value="pval")

# #plot percentage cells in fine celltype per patient(rna), and percentage cells in cluster per patient (met)



# # #plot RNA defined markers that overlap with large DMRs
# # plot_dmr_RNAMarker_hypoM(dat=dat_fine,
# #                             dmr=dmr,
# #                             rna_markers=markers,
# #                             broad_celltype=celltype,
# #                             order=c("1","4","2","7","5","9","3","8","6")) #loose ordering by umap grouping

# # #ASSIGN CELL TYPES BY PLOT
# # dat_fine@metadata$fine_celltypes<-"unknown"
# # dat_fine@metadata[dat_fine@metadata[[paste(celltype,"clusters",sep=".")]] %in% c("0"),]$fine_celltypes<-"unknown"

# # #SAVE FINE CELL TYPES
# # saveRDS(dat_fine,paste("04_scaledcis",celltype,"fine_celltype.amethyst.rds",sep="."))

# ```

# #run for lymphoid, myeloid

# ```R

# ######RUNNING###########

# #read in subset amethyst object and dmrs
# dat_fine<-readRDS(
#     paste(project_data_directory,"fine_celltyping",broad_celltype,
#     paste("04_scaledcis",broad_celltype,"fine_celltype.amethyst.rds",sep="."),
#     sep="/"))


# ```


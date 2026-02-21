Running intial clustering to identify cell type compartments (to help with CNV calling)

After object is initiated, cell lines filtered out, and QC (by coverage and percent methylation done), perform clustering.

1. Cells passing QC have methylation summarized over 5000bp windows > 02_scaledcis.5kbpwin.passfilt.amethyst.rds
2. Cells are clustered on 5000bp windows  > 03_scaledcis.5kbpclus_dmrs.amethyst.rds
3. summarize methylation on 5kbp based clusters over 500bp windows
4. define DMRs between clusters 
5. collapse overlapping dmrs
6. recluster on collapsed dmrs and repeat steps 3-5 for final clusters > 04_scaledcis.5kbpclus_dmrs.amethyst.rds

Define cell types per cluster by:
1. Final clustering on coarse_cluster DMRs (500bp windows from 04_scaledcis.5kbpclus_dmrs.amethyst.rds clusters)
2. 500bp windows methylation plotted over marker genes > 05_scaledcis.coarse_clusters.amethyst.rds
3. Hypomethylation DMRs overlapping genes used for integration with RNA (? TBD)

After copynumber calling, also redefine cells with CNVs as cancer cell type.

```R
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

set.seed(111)
options(future.globals.maxSize= 2500000*1024^2) #80gb limit for parallelizing

task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)

system(paste("mkdir -p",wd))
setwd(wd)

#make merged object
plate_obj=paste(sep="/",project_data_directory,"scalemethyl_pipeline_out/amethyst_plate_obj")
amethyst_files=list.files(path=plate_obj,pattern="*.amethyst.rds",recursive=TRUE,full.names=TRUE)

```

## Read in all amethyst files

```R


window_name="cg_100k_score"
dat_list<-mclapply(amethyst_files, function(x) {
    obj<-readRDS(x)
    return(obj)},mc.cores=20)

dat <- combineObject(objList = dat_list, genomeMatrices=window_name)
```


## Run 5kb windows on all samples to cluster

Goal here is to identify cell type compartments (epithelial, stromal, immune).

This will be used to:
1. generate clusters of cells for refined DMR calling
2. help with CNV calling 
3. validate the RNA integration

Using large windows (5kb) so coverage has less of an effect.

```R
dat@genomeMatrices[["initial_cluster_5kb_win"]] <- makeWindows(dat, 
                                                     stepsize = 5000,
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 100, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

saveRDS(dat,file="00_merged.amethyst.rds")
```

### Subset out cell lines

```R
dat_celllines<-subsetObject(dat,cells=row.names(dat@metadata)[dat@metadata$sample %in% c("MCF10A","MCF7","MDA-MB-231")])
saveRDS(dat_celllines,file="01_celllines.amethyst.rds")
```

### Subset out patient samples
Adding patient metadata as well.

```R
#subset to just patient samples
dat<-subsetObject(dat,cells=row.names(dat@metadata)[!(dat@metadata$sample %in% c("MCF10A","MCF7","MDA-MB-231"))])
patient_metadata<-read.csv("/data/rmulqueen/projects/scalebio_dcis/sample_selection/simplified_patient_metadata.csv")
merge_meta<-merge(dat@metadata,patient_metadata,by.x="sample",by.y="Processing_Name",all.x=TRUE)
table(merge_meta$Sample,useNA="ifany")

dat@metadata<-merge_meta
row.names(dat@metadata)<-paste(dat@metadata$cell_id,dat@metadata$plate_info,sep="+")
saveRDS(dat,file="01_scaledcis.amethyst.rds")

```
Proceeding with only patient samples.

## QC AND FILTER
Filtering by coverage and percent methylation

```R
dat@metadata %>% 
    group_by(Sample) %>% 
    dplyr::summarize(cell_count=n(),mean_reads=mean(unique_reads),mean_mcg=mean(cg_cov)) 

#simple filters by cov and CG percentage
dat<-subsetObject(dat,cells=row.names(dat@metadata[dat@metadata$cg_cov>=25000,]))
dat<-subsetObject(dat,cells=row.names(dat@metadata[dat@metadata$mcg_pct>=50,]))

#replace cov with cg specific cov (amethyst requires these metadata columns)
dat@metadata$ch_cov<-dat@metadata$cov
dat@metadata$cov<-dat@metadata$cg_cov

saveRDS(dat,file="02_scaledcis.5kbpwin.passfilt.amethyst.rds")
```

## Cluster on 5kbp windows

### Functions used for iterative clustering
Writing as a function so I can tweak some parameters more easily.

```R
#note that calculate_DMR returns an RDS of all files, but returns a DMR matrix on just the hypomethylated regions (which has improved clustering heuristically)
calculate_DMR<-function(dat=dat,project_data_directory=project_data_directory,groupBy="cluster_id",prefix="5kbp_initial_clusters"){
  #calculate 500bp windows by groupBy
  celltype500bpwindows <- calcSmoothedWindows(dat, 
                                          type = "CG", 
                                          threads = 100,
                                          step = 500, # change to 500 for real data unless you have really low coverage
                                          smooth = 3,
                                          genome = "hg38",
                                          index = "chr_cg",
                                          groupBy = groupBy,
                                          returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                          returnPctMatrix = TRUE)
  #initiate directory for output
  dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")
  dmr_celltype_outdir=paste(sep="/",dmr_outdir,prefix)
  system(paste("mkdir -p", dmr_celltype_outdir))


  #save 500bp window matrix                            
  saveRDS(celltype500bpwindows,file=paste0(dmr_celltype_outdir,"/","dmr_analysis.",prefix,".500bp_windows.rds"))
  celltype500bpwindows<-readRDS(file=paste0(dmr_celltype_outdir,"/","dmr_analysis.",prefix,".500bp_windows.rds"))


  pct_mat<-celltype500bpwindows[["pct_matrix"]] 
  sum_mat<-celltype500bpwindows[["sum_matrix"]] 

  #save clone object for future genome track plotting
  dat@genomeMatrices[[paste0("cg_",prefix,"_tracks")]] <- pct_mat #load it into amethyst object for plotting

  #calculate DMRs
  dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
          eachVsAll = TRUE, # If TRUE, each group found in the sumMatrix will be tested against all others
          nminTotal = 3, # Min number observations across all groups to include the region in calculations
          nminGroup = 3) # Min number observations across either members or nonmembers to include the region
  #save DMRs
  saveRDS(dmrs,file=paste0(dmr_celltype_outdir,"/",prefix,".dmr.rds"))


  #filter and collapse DMRs
  dmrs <- filterDMR(dmrs, 
              method = "bonferroni", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
              filter = FALSE, # If TRUE, removes insignificant results
              pThreshold = 0.05, # Maxmimum adjusted p value to allow if filter = TRUE
              logThreshold = 1.5) # Minimum absolute value of the log2FC to allow if filter = TRUE

  collapsed_dmrs <- collapseDMR(dat, 
                        dmrs, 
                          maxDist = 2000, # Max allowable overlap between DMRs to be considered adjacent
                          minLength = 2000, # Min length of collapsed DMR window to include in the output
                          reduce = T, # Reduce results to unique observations (recommended)
                          annotate = T) # Add column with overlapping gene names



  saveRDS(collapsed_dmrs,file=paste0(dmr_celltype_outdir,"/",prefix,".dmr_filt_collapse.rds"))
  rename_dmr_output<-setNames(nm=1:length(unique(collapsed_dmrs$test)),gsub(colnames(sum_mat)[grepl(colnames(sum_mat),pattern="_t$")],pattern="_t",replacement=""))

  collapsed_dmrs$type <- rename_dmr_output[collapsed_dmrs$test]
  saveRDS(collapsed_dmrs,file=paste0(dmr_celltype_outdir,"/",prefix,".dmr_filt_collapse.rds"))
  collapsed_dmrs<-readRDS(file=paste0(dmr_celltype_outdir,"/",prefix,".dmr_filt_collapse.rds"))

  #plot dmr counts per group
  pal <- colorRampPalette(colors = c("#F9AB60", "#E7576E", "#630661", "#B5DCA5"))
  COLS <- pal(length(unique(collapsed_dmrs$type)))

  plt<-ggplot(collapsed_dmrs |> dplyr::group_by(type, direction) |> dplyr::summarise(n = n()), 
      aes(y = type, x = n, fill = type)) + 
      geom_col() + 
      facet_grid(vars(direction), scales = "free_y") + 
      scale_fill_manual(values = COLS) + 
      theme_classic()
  ggsave(plt,file=paste0(dmr_celltype_outdir,"/",prefix,".dmr_counts.pdf"))

  #filtering dmrs
  #use only hypo dmrs, filter by length and merge any that overlap (across cell types)
  dmrs <- collapsed_dmrs %>% 
    filter(direction=="hypo") %>% 
    filter(dmr_padj<0.05) %>% 
    filter(dmr_length < 20000) %>% 
    group_by(type) %>% 
    #slice_min(n=20000, dmr_logFC) %>%
    filter(abs(dmr_logFC)>1.5) %>%
    as.data.frame() %>%
    GenomicRanges::makeGRangesFromDataFrame() %>%   GenomicRanges::reduce()
  #22,301 dmrs

  summary(width(dmrs))
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #   2001    6001   10501   11231   16001   81001 

  length(dmrs)
  dmr_bed<-data.frame(chr=seqnames(dmrs),start=start(dmrs),end=end(dmrs))

  #create new matrix from DMR sites for refined clustering
  window_name=paste0(prefix,"_dmr_sites")

  dat@genomeMatrices[[window_name]] <- makeWindows(dat, 
                                                      bed = dmr_bed,
                                                      type = "CG", 
                                                      metric = "score", 
                                                      threads = 50, 
                                                      index = "chr_cg", 
                                                      nmin = 2) 
  return(dat)
}

celltype_umap<-function(obj=dat,prefix="allcells",dims=12,regressCov=TRUE,k_pheno=50,neigh=25,dist=1e-5,method="cosine",output_directory,window_name){
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
  ggsave((p1|p2)/(p3|p4)/(p5|p6),file=paste0(output_directory,"/",outname,"_umap.pdf"),width=20,height=30)  
  return(obj)
}

```


## Iterative clustering
1. Initialize clusters on 5kbp windows, and then:
2. summarize methylation over 500bp windows per cluster
3. define DMRs between clusters
4. collapse overlapping dmrs
5. recluster on dmrs
6. repeat on new DMR cluster to define final clusters and DMRS
7. label cell types by marker


*1. Initialize clusters on 5kbp windows, and then:*

```R
dat<-readRDS(file="02_scaledcis.5kbpwin.passfilt.amethyst.rds")

#preexisting window used for first round of clustering (5kb windows)
window_name="initial_cluster_5kb_win"
dat@genomeMatrices[[window_name]] <- dat@genomeMatrices[[window_name]][rowSums(!is.na(dat@genomeMatrices[[window_name]])) >= 45, ]
est_dim<-dimEstimate(dat, genomeMatrices = c(window_name), dims = c(20), threshold = 0.95)
print(est_dim)
#9

dat<-celltype_umap(obj=dat,prefix="02_5kbp_initialclustering",
                  dims=9,
                  window_name=window_name,
                  regressCov=FALSE,
                  k_pheno=200,
                  neigh=10,
                  dist=0.001,
                  method="cosine",
                  output_directory=getwd())

```

2. summarize methylation over 500bp windows per cluster
3. define DMRs between clusters
4. collapse overlapping dmrs
5. recluster on dmrs

```R
#define DMRs on 5kb window initial clusters
prefix="5kbp_initial_clusters"
dat<-calculate_DMR(dat=dat,
                    project_data_directory=project_data_directory,
                    groupBy="cluster_id",
                    prefix=prefix)
saveRDS(dat,file="03_scaledcis.5kbpclus_dmrs.amethyst.rds") #step 4 in iterative clustering

#use DMRs from initial clusters to refine clustering
window_name=paste0(prefix,"_dmr_sites")
dat@genomeMatrices[[window_name]] <- dat@genomeMatrices[[window_name]][rowSums(!is.na(dat@genomeMatrices[[window_name]])) >= 45, ]
dim(dat@genomeMatrices[[window_name]])
est_dim<-dimEstimate(dat, genomeMatrices = c(window_name), dims = c(20), threshold = 0.95)
print(est_dim)
#10

dat<-celltype_umap(obj=dat,
                  prefix="04_dmrclustering",
                  dims=11,
                  regressCov=FALSE,
                  k_pheno=60,
                  neigh=7,
                  dist=0.001,
                  method="cosine",
                  window_name=window_name,
                  output_directory=getwd())

dat@metadata$coarse_cluster_id<-dat@metadata$cluster_id
saveRDS(dat,file="04_scaledcis.coarse_cluster.amethyst.rds") #step 5 in iterative clustering
```

Calculate final coarse cluster DMRs 

```R
prefix="coarse_cluster"
dat<-calculate_DMR(dat=dat,
                    project_data_directory=project_data_directory,
                    groupBy="coarse_cluster_id",
                    prefix=prefix)
saveRDS(dat,file="05_scaledcis.coarse_clusters.amethyst.rds") #step 7 in iterative clustering

```

Plot markers per cluster for initial cell typing

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

    if (!is.null(colors)) {pal <- colors} else {pal <- c("#FF0082", "#dbdbdb", "#cccccc", "#999999")}
    
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

cell_markers<-list()
cell_markers[["tcell"]]<-c("PTPRC","IKZF1","IL7R","CD3D","CD3E","CD3G") 
cell_markers[["bcell"]]<-c("CD79A","CD79B","MS4A1","HLA-DRA","HLA-DPA1","HLA−DPB1") 
cell_markers[["myeloid"]]<-c("LYZ","CD68","CSF1R","FOLR2","FUT4","FCGR3A","TREM2","FCGR1A") #neutrophils monocytes/macrophages DC mast TAMs
cell_markers[["fibro"]]<-c("DCN","APOD","LUM","COL1A2","COL1A1","FAP","PDGFRA") #fibro and caf
cell_markers[["endo"]]=c("CLDN5","BTNL9","PTPRB","VWF")
cell_markers[["lumsec"]]<-c("GABRP","ELF5","CCL28","KRT15","KIT","MMP7","LTF","SLPI")
cell_markers[["basal"]]<-c("CARMN","ACTA2","KRT17","KRT14","KRT5")
cell_markers[["lumhr"]]<-c("EPCAMP","AREG","AZGP1","KRT18","AGR2","PIP","ANKRD30A","FOXA1","ESR1","PGR")

cell_colors=c(
"basal"="#844c9d",
"lumhr"="#e23e96",
"lumsec"="#ff6498",
"fibro"="#f58e90",
"endo"="#5bbb5a",
"perivasc"="#eaba67",
"myeloid"="#8088c2",
"tcell"="#1d87c8",
"lymphocytes"="#1d87c8",
"mast"="#dcd0ff",
"bcell"="#65cbe4",
"plasma"="#7ecdc2",
"adipo"="#b48454")


#prepare cgi
cgisland="/data/rmulqueen/projects/scalebio_dcis/ref/cpgIslandExt.bed"
cgi<-rtracklayer::import(cgisland)
cgi<-as.data.frame(cgi)
colnames(cgi)<-c("chr","start","end","strand")


#rough ordering by what clustered together
order<-c(
"1","4", #tcell
"16", #bcell
"6", #myeloid
"8","3","15", #fibro
"12",#endo
"11","9", #basal
"14","13", #lumsec
"5", #lumhr
"10","2","17","7" #cancer
)

    
mclapply(names(cell_markers),
function(celltype){
    genes<-unlist(cell_markers[celltype])
    genes<-genes[genes %in% dat@ref$gene_name]
    print(paste("Plotting:",celltype))
    print(paste("Genes to plot:",genes))
    plt<-histograModified(obj=dat, 
        baseline="mean",
        genes = unlist(cell_markers[celltype]),
        colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
        matrix = paste0("cg_",prefix,"_tracks"), arrowScale = .03, trackScale = .5,
        legend = F, cgisland=cgi,order=order) + ggtitle(celltype)
    ggsave(plt,
          file=paste0(wd,"/","05_",prefix,".",celltype,".marker.pdf"),
          width=3*length(genes),
          height=ncol(dat@genomeMatrices[[paste0("cg_",prefix,"_tracks")]])*1,
          limitsize=F)
},mc.cores=20)

#assigning broad cell types
dat@metadata$broad_celltype<-NA
dat@metadata[dat@metadata$coarse_cluster_id %in% c("1","4"),]$broad_celltype<-"tcell"
dat@metadata[dat@metadata$coarse_cluster_id %in% c("16"),]$broad_celltype<-"bcell"
dat@metadata[dat@metadata$coarse_cluster_id %in% c("6"),]$broad_celltype<-"myeloid"
dat@metadata[dat@metadata$coarse_cluster_id %in% c("8","3","15"),]$broad_celltype<-"fibroblast"
dat@metadata[dat@metadata$coarse_cluster_id %in% c("12"),]$broad_celltype<-"endothelial"
dat@metadata[dat@metadata$coarse_cluster_id %in% c("11","9"),]$broad_celltype<-"basal"
dat@metadata[dat@metadata$coarse_cluster_id %in% c("14","13"),]$broad_celltype<-"lumsec"
dat@metadata[dat@metadata$coarse_cluster_id %in% c("5"),]$broad_celltype<-"lumhr"
dat@metadata[dat@metadata$coarse_cluster_id %in% c("10","2","17","7"),]$broad_celltype<-"cancer"
saveRDS(dat,file="05_scaledcis.coarse_clusters.amethyst.rds")

#umap plot of celltypes to confirm
p1 <- dimFeature(dat, colorBy = broad_celltype, reduction = "umap")
ggsave(p1,file="05_scaledcis.coarse_clusters.celltype.umap.pdf",width=10,height=10)  


```


```R

#run 3d plotting just for fun
dim3d<-uwot::umap(
  X=dat@reductions[["5kbp_initial_clusters_dmr_sites_irlba"]],
  n_neighbors = 7,
  n_components = 3,
  metric = "cosine",
  seed = 123,
  n_threads=50,
)


dim3d<-as.data.frame(dim3d)
colnames(dim3d)<-c("X","Y","Z")
dim3d$cellname<-row.names(dim3d)
dim3d$celltype<-as.data.frame(dat@metadata[row.names(dim3d),])$broad_celltype
dim3d$hex_color<-cell_colors[dim3d$celltype]
dim3d$r_col<-unlist(col2rgb(dim3d$hex_color)["red",])
dim3d$g_col<-unlist(col2rgb(dim3d$hex_color)["green",])
dim3d$b_col<-unlist(col2rgb(dim3d$hex_color)["blue",])

write.table(dim3d,col.names=T,file="05_scaledcis.broad_celltype.3dumap.csv",sep=",",row.names=F)

```

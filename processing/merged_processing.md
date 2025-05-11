```bash
#run per line 
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif

cd /data/rmulqueen/projects/scalebio_dcis/data
```

```R
library(amethyst)
library(rhdf5)
library(data.table)
library(ggplot2)
library(patchwork)
library(tibble)
library(tidyr)
library(plyr); library(dplyr)
library(future)
library(furrr)
library(purrr)
library(cowplot)
library(pheatmap)
library(parallel)
library(rtracklayer)
library(gridExtra)

library(optparse,lib.loc="/home/users/rmulqueen/R/x86_64-conda-linux-gnu-library/4.4") #add this
set.seed(111)

options(future.globals.maxSize= 50*1024^3) #50gb limit

task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data"
setwd(project_data_directory)
amethyst_files=list.files(path=project_data_directory,pattern="*.amethyst.rds",recursive=TRUE,full.names=TRUE)

system("mkdir -p merged_data")
window_name="cg_100k_score"

########################################
## Read in all amethyst files
########################################

dat_list<-mclapply(amethyst_files, function(x) {
    obj<-readRDS(x)
    return(obj)},mc.cores=20)



#note BCMDCIS81T is TNBC and has no cells. I should add a filter in pipeline to remove it earlier.
#in this case i just deleted the file
########################################
## Run nakshatri atac peaks on all samples to cluster
########################################
nakshatri_peaks<-readRDS("/data/rmulqueen/projects/scalebio_dcis/ref/celltype_peaks.500bp.rds")
nakshatri_bed<-data.frame(
        seqnames=seqnames(nakshatri_peaks),
        starts=start(nakshatri_peaks),
        ends=end(nakshatri_peaks))
dat <- combineObject(objList = dat_list, genomeMatrices=window_name)

nakshatri_bed<-nakshatri_bed[nakshatri_bed$seqnames %in% unique(dat@ref$seqid),]
dat@genomeMatrices[["nakshatri_peaks"]] <- makeWindows(dat, 
                                                     bed = nakshatri_bed,
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 300, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

saveRDS(dat,file="./merged_data/all_samples.amethyst.rds")

#subset to just patient samples
dat<-subsetObject(dat,cells=row.names(dat@metadata)[!(dat@metadata$sample %in% c("MCF10A","MCF7","MDA-MB-231"))])
saveRDS(dat,file="./merged_data/scaledcis.amethyst.rds")

dat<-subsetObject(dat,cells=row.names(dat@metadata[dat@metadata$cg_cov>=50000,]))
#replace cov with cg specific cov 
dat@metadata$ch_cov<-dat@metadata$cov 
dat@metadata$cov<-dat@metadata$cg_cov

cluster_by_windows<-function(obj,
    window_name,
    metric.="score",
    threads.=100,
    neighbors.=50,
    est_dim=10,
    k_pheno=15,
    dist.=0.1,
    outname="clustering"){
  print(paste("Estimating dimensions..."))                                           
  #filter windows by cell coverage
  obj@genomeMatrices[[window_name]] <- obj@genomeMatrices[[window_name]][rowSums(!is.na(obj@genomeMatrices[[window_name]])) >= 45, ]
  #est_dim<-dimEstimate(obj, genomeMatrices = c(window_name), dims = c(10), threshold = 0.95)
  #print(est_dim)

  print("Running IRLBA reduction...")
  obj@reductions[[paste(window_name,"irlba",sep="_")]] <- runIrlba(obj, genomeMatrices = c(window_name), dims = est_dim, replaceNA = c(0))

  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- regressCovBias(obj, reduction = paste(window_name,"irlba",sep="_")) # Optional; helps reduce coverage bias
  print("Clustering on coverage regressed reduction...")

  obj <- runCluster(obj, k_phenograph = 25, reduction = paste(window_name,"irlba_regressed",sep="_")) # consider increasing k_phenograph to 50 for larger datasets
  obj <- runUmap(obj, neighbors = neighbors., dist = dist., method = "euclidean", reduction = paste(window_name,"irlba_regressed",sep="_")) 
  
  print("Plotting...")
  p1 <- dimFeature(obj, colorBy = cluster_id, reduction = "umap") + ggtitle(paste(window_name,"Clusters"))
  p2 <- dimFeature(obj, colorBy = sample, reduction = "umap") + ggtitle(paste(window_name,"Samples"))
  p3 <- dimFeature(obj, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Coverage distribution")
  p4 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
  plt<-plot_grid(p1, p2,p3, p4,ncol=2)
  ggsave(plt,file=paste0(outname,"_umap.pdf"),width=20,height=20)     
  return(obj)                                             
}

dat<-cluster_by_windows(dat,window_name="cg_100k_score",threads.=task_cpus)
dat<-cluster_by_windows(dat,window_name="nakshatri_peaks",threads.=task_cpus,est_dim=7)

#final hyperparamters
dat<-cluster_by_windows(dat,window_name="nakshatri_peaks",threads.=task_cpus,est_dim=12,neighbors.=25,dist.=0.01,k_pheno=15)

saveRDS(dat,file="./merged_data/scaledcis.amethyst.pf.rds")
########################################
## From Nakshatri ATAC peaks, generate list of DMRs and recluster on those.
### Nakshatri is only normal breast tissue, so clustering is hampered for DCIS/IDC epithelia
########################################


###### 1kb window generation, and DMR Calculation ####
dmr_and_1kb_window_gen<-function(
    obj=dat,
    prefix="nakshatri_peaks",
    groupBy="cluster_id",
    threads.=20){
    cluster1kbwindows <- calcSmoothedWindows(obj, 
                                            type = "CG", 
                                            threads = threads.,
                                            step = 1000,
                                            smooth = 3,
                                            index = "chr_cg",
                                            groupBy = groupBy, 
                                            returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                            returnPctMatrix = TRUE)
    obj@genomeMatrices[[paste0("cg_",groupBy,"_tracks")]] <- cluster1kbwindows[["pct_matrix"]]

    pal=c("#E5E6E4","#CFD2CD","#A6A2A2","#847577","#6E44FF")
    dmrs<-testDMR(cluster1kbwindows[["sum_matrix"]], eachVsAll = TRUE, nminTotal = 0, nminGroup = 0) 
    dmrs2<-filterDMR(dmrs, method = "bonferroni", filter = TRUE,pThreshold=0.05,logThreshold=0.5) #add additional columns direction column
    celltype_test<-setNames(nm=unique(dmrs2$test), gsub("_c$","",colnames(cluster1kbwindows[["sum_matrix"]])[grepl(pattern="_c$",colnames(cluster1kbwindows[["sum_matrix"]]))]))
    dmrs2$celltype<-celltype_test[dmrs2$test]
    saveRDS(dmrs2,file=paste0(prefix,".dmr.",groupBy,".rds"))

    collapsed_dmrs <- collapseDMR(obj, dmrs2, maxDist = 4000, minLength = 1000, reduce = T, annotate = T) 
    collapsed_dmrs$celltype<-celltype_test[collapsed_dmrs$test]
    saveRDS(collapsed_dmrs,file=paste0(prefix,".dmr.",groupBy,".collapsed.rds"))

    plt<-ggplot(collapsed_dmrs |> dplyr::group_by(celltype, direction) |> dplyr::summarise(n = dplyr::n()), 
       aes(y = celltype, x = n, fill = celltype)) + geom_col() + 
        facet_grid(vars(direction), scales = "free_y") + 
        scale_fill_manual(values = makePalette(option = 7, n = length(unique(collapsed_dmrs$celltype)) ) ) + 
        theme_classic()
    
    ggsave(plt,file=paste0(prefix,".met_per_dmr.",groupBy,".pdf"))
    median(collapsed_dmrs$dmr_end-collapsed_dmrs$dmr_start)

    top_dmrs <- collapsed_dmrs |> 
    dplyr::group_by(celltype, direction) |> 
    dplyr::arrange(dmr_padj, .by_group = TRUE) |> dplyr::mutate(rank_padj = 1:dplyr::n()) |>
    dplyr::arrange(desc(abs(dmr_logFC)), .by_group = TRUE) |> dplyr::mutate(rank_logFC = 1:dplyr::n()) |>
    rowwise() |> dplyr::mutate(total_rank = sum(rank_padj, rank_logFC)) |> 
    group_by(celltype, direction) |> slice_min(n = 20, order_by = total_rank) |>
    dplyr::mutate(location = paste0(chr, "_", (dmr_start - 2000), "_", (dmr_end + 2000))) |> dplyr::arrange(direction)
    write.table(top_dmrs,file=paste0(prefix,".cluster_dmr.",groupBy,".tsv"))
    return(list(object=obj,dmr=collapsed_dmrs))
}

out<-dmr_and_1kb_window_gen(dat,prefix="nakshatri_peaks",threads.=20)
dmr_out<-out[["dmr"]]
obj<-out[["object"]]

saveRDS(obj,file="./merged_data/scaledcis.amethyst.pf.nakshatri_clus.rds")
saveRDS(dmr_out,file="./merged_data/scaledcis.nakshatri_cluster.dmr.rds")

dmr_out<-readRDS("scaledcis.nakshatri_cluster.dmr.rds")

#filter by significance <0.01 and logFC > 1, merge regions that overlap
dmr_out<-dmr_out |> dplyr::filter(dmr_padj<0.01) |> dplyr::filter(dmr_logFC < -1 | dmr_logFC > 1) |> select(chr,dmr_start,dmr_end) |> distinct(chr,dmr_start,dmr_end) |> makeGRangesFromDataFrame() |> reduce()
dmr_bed<-data.frame(seqnames=seqnames(dmr_out),start=start(dmr_out),end=end(dmr_out))
dat@genomeMatrices[["dmr_sites"]] <- makeWindows(dat, 
                                                     bed = dmr_bed,
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 20, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 


for(d in c(8,10,12,14,16,18,20)){
    for(neigh in c(20,50,75,100)){
        for(q in c(0.1,0.05,0.01)){
            dat<-cluster_by_windows(dat,window_name="dmr_sites",
            outname=paste("clustest",d,neigh,q,sep="_"),
            threads.=task_cpus,est_dim=d,neighbors.=neigh,dist.=q,
            k_pheno=15)
        }
    }
}


#run DMR and 1kb windows on reclustered data
out<-dmr_and_1kb_window_gen(dat,prefix="dmr_sites",threads.=10)
dmr_out<-out[["dmr"]]
obj<-out[["object"]]
saveRDS(dmr_out,file="dmr_cluster.dmr.rds")
saveRDS(obj,file="scaledcis.amethyst.pf.dmr_clus.rds")




########################################
## Make sense of new clusters by cell type assignment and rGREAT for DMR enrichment
## User promoter methylation + into gene body for heatmap plot
########################################

cell_markers<-list()
cell_markers[["lumhr"]]<-c("AREG","AZGP1","KRT18","AGR2")
cell_markers[["lumsec"]]<-c("GABRP","ELF5","CCL28","KRT15","KIT","MMP7")
cell_markers[["basal"]]<-c("CARMN","ACTA2","KRT17")
cell_markers[["adipo"]]=c("PDE3B","ACACB","WDPCP","PCDH9","CLSTN2","ADIPOQ","TRHDE")
cell_markers[["vascular"]]=c("MECOM","BTNL9","MCTP1","PTPRB","VWF","ADGRL4","LDB2")
cell_markers[["lymphatic"]]=c("AL357507.1","PKHD1L1","KLHL4","LINC02147","RHOJ","ST6GALNAC3","MMRN1")
cell_markers[["perivasc"]]=c("RGS6","KCNAB1","COL25A1","ADGRL3","PRKG1","NR2F2-AS1","AC012409.2")
cell_markers[["fibro"]]<-c("DCN","APOD","LUM","COL1A2","COL1A1")
cell_markers[["myeloid"]]<-c("HLA-DRA","HLA-DPA1")
cell_markers[["bcell"]]=c("CD37","TCL1A","LTB","HLA-DPB1","HLA-DRA","HLA-DPA1")
cell_markers[["plasma"]]=c("IGHA2","IGHA1","JCHAIN","IGHM","IGHG1","IGHG4","IGHG3","IGHG2")
cell_markers[["tcell"]]<-c("PTPRC","IKZF1","IL7R","GNLY")

cell_colors=c("lumhr"="#e23e96",
"lumsec"="#ff6498",
"basal"="#844c9d",
"adipo"="#b48454",
"vascular"="#5bbb5a",
"lymphatic"="#b5d564",
"perivasc"="#eaba67",
"fibro"="#f58e90",
"myeloid"="#8088c2",
"bcell"="#65cbe4",
"plasma"="#7ecdc2",
"tcell"="#1d87c8")

cgisland="/data/rmulqueen/projects/scalebio_dcis/ref/cpgIslandExt.bed"

#prepare cgi
cgi<-rtracklayer::import(cgisland)
cgi<-as.data.frame(cgi)
colnames(cgi)<-c("chr","start","end","strand")

#modified from histograM
histograModified <- function(obj,
    genes = NULL,
    matrix,colors = NULL,
    ncol = length(genes),
    trackOverhang = 5000,arrowOverhang = 3000,trackScale = 1.5,arrowScale = 0.025,colorMax = 100,
    legend = TRUE,removeNA = TRUE,
    order = NULL,baseline = "mean",
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
            ggplot2::scale_fill_gradientn(colors = pal, limits = c(0,colorMax), oob = scales::squish) + theme_minimal()+ ylab(j) + ggplot2::theme(legend.position="none") +
            ylim(c(0,100))+
            ggplot2::theme(panel.background = element_blank(), axis.ticks = element_blank(), panel.grid.major.y = element_line(color = "#dbdbdb", linetype = "dashed"))
            }

        #first track is annot
        p[[i]][["gene_track"]] <-ggplot()+
        ggplot2::geom_rect(fill = "pink",alpha=0.8, data = ref |> dplyr::filter(type == "gene") |>
                            dplyr::mutate(promoter_start = ifelse(strand == "+", (start - 1500), (end + 1500)),
                            promoter_end = ifelse(strand == "+", (promoter_start+3000), (promoter_start-3000))),
                            ggplot2::aes(xmin = promoter_start, xmax = promoter_end, ymin = 0, ymax = trackHeight/4)) +
        ggplot2::geom_rect(alpha=0.8,fill = "black", data = ref |> dplyr::filter(type == "exon"),ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = trackHeight/4)) +
        ggplot2::geom_rect(alpha=0.8,fill = "green", data = cgi, ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = trackHeight/4)) +
        ggplot2::geom_segment(data = ref, color="gray",
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
    p[[i]]<-wrap_plots(p[[i]],guides='collect')+plot_layout(ncol=1)
    }
    all_genes_plot<-wrap_plots(p)+plot_layout(nrow=1)
    return(all_genes_plot)
}

plot_histogram_page<-function(celltype){
    plt<-histograModified(obj, 
        baseline="mean",
        genes = unlist(cell_markers[celltype]),
        colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
        matrix = "cg_cluster_id_tracks", arrowScale = .03, trackScale = .5,
        legend = F, cgisland=cgi)
    return(plt)
}

#rough ordering by what clustered together

order<-c('14','9','17','2','5','3','1','11','19','18','12','15','4','13','16','7','10','6','8','20','21')
plt_list<-lapply(names(cell_colors),
function(celltype){
    genes<-unlist(cell_markers[celltype])
    genes<-genes[genes %in% obj@ref$gene_name]
    print(paste("Plotting:",celltype))
    print(paste("Genes to plot:",genes))
    plt<-histograModified(obj, 
        baseline="mean",
        genes = unlist(cell_markers[celltype]),
        colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
        matrix = "cg_cluster_id_tracks", arrowScale = .03, trackScale = .5,
        legend = F, cgisland=cgi,order=order) + ggtitle(celltype)
    return(plt)
})

plt_out<-wrap_plots(plt_list,ncol=1) 
ggsave(plt_out ,file="test.marker.pdf",width=30,height=length(cell_colors)*30,limitsize=F)




#dmr around promoter regions for cell type calling too (add to plot list)
#methyltree output
#scage output with celltypes




methyltree_output<-function(obj=obj,prefix="DCIS-41T",sample="DCIS-41T",filt_min_pct=10,filt_max_pct=80,threads=1){
        dcis<-subsetObject(obj, cells = row.names(obj@metadata[obj@metadata$sample %in% sample,]))
        dcis@metadata$methyltree_group<-"all"
        #make 500bp windows with methylation percentages
        methyltreewindows <- calcSmoothedWindows(dcis, 
                                            type = "CG", 
                                            threads = threads,
                                            step = 500,
                                            smooth = 1,
                                            index = "chr_cg",
                                            groupBy = "methyltree_group", 
                                            returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                            returnPctMatrix = TRUE)
        print(paste("Starting number of windows:",as.character(nrow(methyltreewindows[["pct_matrix"]]))))
        methyltreewindows[["pct_matrix"]]<-methyltreewindows[["pct_matrix"]][methyltreewindows[["pct_matrix"]]$all>=filt_min_pct & methyltreewindows[["pct_matrix"]]$all<=filt_max_pct,]
        #filter to windows to middling methylation values
        print(paste("Filtering by m0 >=", as.character(filt_min_pct), "m1 <=", as.character(filt_max_pct),as.character(nrow(methyltreewindows[["pct_matrix"]]))))
        #merge windows that are touching
        methyltreewindows<-reduce(GenomicRanges::makeGRangesFromDataFrame(methyltreewindows[["pct_matrix"]]))
        print(paste("Filtered window count:",as.character(nrow(as.data.frame((methyltreewindows))))))
        print(paste("Filtered window average width:",as.character(mean(width(methyltreewindows)))))
        print(paste("Total genome covered:",as.character(sum(width(methyltreewindows))/1000000),"Mbp"))
        #make a merged windows percentile matrix per cell for methyltree

        methyltreeoutput<-makeWindows(dcis,
                    type = "CG", 
                    metric = "percent", 
                    bed = as.data.frame(methyltreewindows,col.names=NULL)[,1:3],
                    threads = threads, 
                    index = "chr_cg", 
                    nmin = 1) 
      print(paste("Mean percent cells covered per window:",
            mean((rowSums(!is.na(methyltreeoutput))/ncol(methyltreeoutput))*100)))
      print("Filtering to windows with >10% of cells with coverage")
      methyltreeoutput<-methyltreeoutput[(rowSums(!is.na(methyltreeoutput))/ncol(methyltreeoutput)*100)>=10,]
      methyltreewindows<-data.frame(do.call("rbind",strsplit(row.names(methyltreeoutput),"_")))
      colnames(methyltreewindows)<-c("chr","start","end")
      methyltreewindows<-GenomicRanges::makeGRangesFromDataFrame(methyltreewindows)
      print(paste("Final Filtered window count:",as.character(nrow(as.data.frame((methyltreewindows))))))
      print(paste("Final Filtered window average width:",as.character(mean(width(methyltreewindows)))))
      print(paste("Final Total genome covered:",as.character(sum(width(methyltreewindows))/1000000),"Mbp"))
      methyltreeoutput<-makeWindows(dcis,
                  type = "CG", 
                  metric = "percent", 
                 bed = as.data.frame(methyltreewindows,col.names=NULL)[,1:3],
                  threads = threads, 
                  index = "chr_cg", 
                  nmin = 1) 

      methyltreeoutput$genomic_region_id<-row.names(methyltreeoutput)
      methyltreeoutput <- methyltreeoutput |> 
          pivot_longer(
          cols = !genomic_region_id, 
          names_to = "cell_id",
          values_to = "value",
          values_drop_na = TRUE)
    #make a metadata sheet with cluster info
    out_metadata<-dcis@metadata[,c("pass","cluster_id","cg_cov","mcg_pct","subclones")]
    colnames(out_metadata)<-c("HQ","celltype","nCG","met_rate","large_clone_id") #match names
    out_metadata$sample<-row.names(out_metadata) #sample (cell) names
    out_metadata$met_rate<-out_metadata$met_rate/100 #percentage to rate

    if(file.exists(paste0(prefix,"_methyltree_input.h5"))){
        system(paste0("rm -rf ",prefix,"_methyltree_input.h5"))
    }
      h5createFile(file=paste0(prefix,"_methyltree_input.h5"))
      h5write(methyltreeoutput,file=paste0(prefix,"_methyltree_input.h5"),name="data")
      h5write(out_metadata,file=paste0(prefix,"_methyltree_input.h5"),name="metadata")
    }











#promoter plus into gene body
    promoter_into_genebody_bed <- dat@ref |> dplyr::filter( type == "gene") |>
            dplyr::select(seqid, start, end, gene_name, strand) |>
            dplyr::distinct(gene_name, .keep_all = TRUE) |>
            dplyr::rename("chr" = "seqid", "gene" = "gene_name") |>
            data.table::data.table()
    
    #adding another 5kb into gene body from promoter regions
    promoter_into_genebody_bed[, `:=`(p_start = ifelse(strand == "+", start - 1500, end - 1500),
                    p_end = ifelse(strand == "+", start + 1500 + 5000, end + 1500 + 5000))] # make bed file for promoter cg
        
    promoter_into_genebody_bed[, `:=`(start = NULL, end = NULL, strand = NULL)]
    setnames(promoter_into_genebody_bed, c("chr", "gene", "start", "end"))
    row.names(promoter_into_genebody_bed)<-promoter_into_genebody_bed$gene
    options(future.globals.maxSize= 50*1000*1024^2) #50gb
    dat@genomeMatrices[["promoter_5kb_genebody"]] <- makeWindows(dat, 
                                                        bed = promoter_into_genebody_bed[,c(1,3,4)],
                                                        type = "CG", 
                                                        metric = "score", 
                                                        threads = 5, 
                                                        index = "chr_cg", 
                                                        nmin = 2) 
histograM(dat, 
          genes = c("ELANE", "MPEG1", "SPI1", "CD2", "CD3D"), 
          matrix = "cg_cluster_tracks",
          legend = F)
#bigwig output

bigwig_output<-function(
    obj,
    tracks="cg_cluster_tracks"){
    for(i in 4:ncol(obj@genomeMatrices[[tracks]])){
        cluster=names(obj@genomeMatrices[[tracks]])[i]
        out_bw<-as.data.frame(obj@genomeMatrices[[tracks]])
        out_bw<-out_bw[c("chr","start","end",cluster)]
        out_bw<-GRanges(out_bw[complete.cases(out_bw),])
        names(out_bw@elementMetadata)<-"score"
        out_bw<-out_bw[unique(findOverlaps(out_bw, type = "any", select = "first"))]
        out_bw <- resize(out_bw, width=1000, fix='start') #resize to avoid 1base overlap
        genome(out_bw)<-"hg38"
        hg38_seq_info<-Seqinfo(genome="hg38")
        seqlengths(out_bw)<-as.data.frame(hg38_seq_info)[hg38_seq_info@seqnames %in% out_bw@seqnames,]$seqlengths
        print(paste("Saving bigwig for...",cluster))
        export(out_bw,con=paste(tracks,cluster,"bw",sep="."),format='bigWig')}
}



#add additional methylation windows

dat@genomeMatrices[["cg_promoter"]] <- makeWindows(dat, 
                                                     genes = protein_coding,
                                                     promoter = TRUE, 
                                                     type = "CG", 
                                                     metric = "percent", 
                                                     threads = 300, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 



dat@genomeMatrices[["cg_genebody"]] <- makeWindows(dat, 
                                                     genes = protein_coding,
                                                     promoter = FALSE, 
                                                     type = "CG", 
                                                     metric = "percent", 
                                                     threads = 300, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

Running intial clustering to identify cell type compartments (to help with CNV calling)

After object is initiated, cell lines filtered out, and QC (by coverage and percent methylation done), perform clustering.

1. Cells passing QC have methylation summarized over 5000bp windows > 02_scaledcis.5kbpwin.passfilt.amethyst.rds
2. Cells are clustered on 5000bp windows  > 03_scaledcis.5kbpclus_dmrs.amethyst.rds
3. summarize methylation on clusters over 500bp windows
4. define DMRs between clusters 
5. collapse overlapping dmrs
6. recluster on collapsed dmrs and repeat steps 3-5 for final clusters > 04_scaledcis.5kbpclus_dmrs.amethyst.rds

Define cell types per cluster by:
1. 500bp windows methylation plotted over marker genes > 05_scaledcis.coarse_clusters.amethyst.rds
2. Hypomethylation DMRs overlapping genes used for integration with RNA

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

Writing as a function so I can tweak some parameters more easily.

```R

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
  ggsave((p1|p2)/(p3|p4)/(p5|p6),file=paste0(prefix,"_",outname,"_umap.pdf"),width=20,height=30)  
  return(obj)
  }

window_name="initial_cluster_5kb_win"
dat@genomeMatrices[[window_name]] <- dat@genomeMatrices[[window_name]][rowSums(!is.na(dat@genomeMatrices[[window_name]])) >= 45, ]
est_dim<-dimEstimate(dat, genomeMatrices = c(window_name), dims = c(20), threshold = 0.95)
print(est_dim)
#9

dat<-celltype_umap(obj=dat,prefix="02_5kbp_initialclustering",dims=12,regressCov=FALSE,k_pheno=200,neigh=8,dist=0.001,method="cosine")

```

Take initial clusters on 5kbp windows, and then:
1. summarize methylation over 500bp windows per cluster
2. define DMRs between clusters
3. collapse overlapping dmrs
4. recluster on dmrs
5. repeat on new DMR cluster to define final clusters and DMRS
6. label cell types by marker

```R

celltype500bpwindows <- calcSmoothedWindows(dat, 
                                        type = "CG", 
                                        threads = 100,
                                        step = 500, # change to 500 for real data unless you have really low coverage
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "cluster_id",
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)
#Get DMR per 5kbp clusters
dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")
dmr_celltype_outdir=paste(sep="/",dmr_outdir,"5kbp_initial_clusters")
system(paste("mkdir -p", dmr_celltype_outdir))
                                   
saveRDS(celltype500bpwindows,file=paste0(dmr_celltype_outdir,"/","dmr_analysis.5kbp_clusters.500bp_windows.rds"))
celltype500bpwindows<-readRDS(file=paste0(dmr_celltype_outdir,"/","dmr_analysis.5kbp_clusters.500bp_windows.rds"))

pct_mat<-celltype500bpwindows[["pct_matrix"]] 
sum_mat<-celltype500bpwindows[["sum_matrix"]] 

#save clone object for future genome track plotting
dat@genomeMatrices[["cg_5kbpcluster_tracks"]] <- pct_mat #load it into amethyst object for plotting

dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
        eachVsAll = TRUE, # If TRUE, each group found in the sumMatrix will be tested against all others
        nminTotal = 3, # Min number observations across all groups to include the region in calculations
        nminGroup = 3) # Min number observations across either members or nonmembers to include the region

saveRDS(dmrs,file=paste0(dmr_celltype_outdir,"/","5kbp_clusters",".dmr.rds"))

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

saveRDS(collapsed_dmrs,file=paste0(dmr_celltype_outdir,"/","5kbp_clusters",".dmr_filt_collapse.rds"))

rename_dmr_output<-setNames(nm=1:length(unique(collapsed_dmrs$test)),gsub(colnames(sum_mat)[grepl(colnames(sum_mat),pattern="_t$")],pattern="_t",replacement=""))
collapsed_dmrs$type <- rename_dmr_output[collapsed_dmrs$test]
saveRDS(collapsed_dmrs,file=paste0(dmr_celltype_outdir,"/","5kbp_clusters",".dmr_filt_collapse.rds"))
collapsed_dmrs<-readRDS(file=paste0(dmr_celltype_outdir,"/","5kbp_clusters",".dmr_filt_collapse.rds"))

#plot dmr counts per clone
pal <- colorRampPalette(colors = c("#F9AB60", "#E7576E", "#630661", "#B5DCA5"))
COLS <- pal(length(unique(collapsed_dmrs$type)))

plt<-ggplot(collapsed_dmrs |> dplyr::group_by(type, direction) |> dplyr::summarise(n = n()), 
    aes(y = type, x = n, fill = type)) + 
    geom_col() + 
    facet_grid(vars(direction), scales = "free_y") + 
    scale_fill_manual(values = COLS) + 
    theme_classic()
ggsave(plt,file=paste0(dmr_celltype_outdir,"/","5kbp_clusters",".dmr_counts.pdf"))

#filtering dmrs
#use only hypo dmrs, take top 5000 per celltype, filter by length and merge any that overlap (across cell types)
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
window_name="5kbp_cluster_dmr_sites"


dat@genomeMatrices[[window_name]] <- makeWindows(dat, 
                                                     bed = dmr_bed,
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 50, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

saveRDS(dat,file="03_scaledcis.5kbpclus_dmrs.amethyst.rds") #step 5 in iterative clustering


```

Repeat clustering now with DMRs identified on 5kbp clusters

```R
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

window_name="5kbp_cluster_dmr_sites"
dat@genomeMatrices[[window_name]] <- dat@genomeMatrices[[window_name]][rowSums(!is.na(dat@genomeMatrices[[window_name]])) >= 45, ]
est_dim<-dimEstimate(dat, genomeMatrices = c(window_name), dims = c(20), threshold = 0.95)
print(est_dim)
#10

dat<-celltype_umap(obj=dat,prefix="04_dmrclustering",dims=11,regressCov=FALSE,k_pheno=50,neigh=5,dist=0.01,method="cosine") #current fav
dat@metadata$coarse_cluster_id<-dat@metadata$cluster_id
saveRDS(dat,file="04_scaledcis.5kbpclus_dmrs.amethyst.rds") #step 5 in iterative clustering

```

Calculate final coarse cluster DMRs 

```R
celltype500bpwindows <- calcSmoothedWindows(dat, 
                                        type = "CG", 
                                        threads = 100,
                                        step = 500, # change to 500 for real data unless you have really low coverage
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "coarse_cluster_id",
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)

dat@genomeMatrices[["cg_coarse_cluster_id_perc"]] <- celltype500bpwindows[["pct_matrix"]]

#Get DMR per dmr clusters
dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")
dmr_celltype_outdir=paste(sep="/",dmr_outdir,"dmr_coarse_cluster")

system(paste("mkdir -p", dmr_celltype_outdir))
saveRDS(celltype500bpwindows,file=paste0(dmr_celltype_outdir,"/","dmr_analysis.dmr_coarse_cluster.500bp_windows.rds"))

pct_mat<-celltype500bpwindows[["pct_matrix"]] 
sum_mat<-celltype500bpwindows[["sum_matrix"]] 

#save clone object for future genome track plotting
dat@genomeMatrices[["cg_coarse_cluster_tracks"]] <- pct_mat #load it into amethyst object for plotting

dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
        eachVsAll = TRUE, # If TRUE, each group found in the sumMatrix will be tested against all others
        nminTotal = 3, # Min number observations across all groups to include the region in calculations
        nminGroup = 3) # Min number observations across either members or nonmembers to include the region

saveRDS(dmrs,file=paste0(dmr_celltype_outdir,"/","coarse_cluster",".dmr.rds"))

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

saveRDS(collapsed_dmrs,file=paste0(dmr_celltype_outdir,"/","coarse_cluster",".dmr_filt_collapse.rds"))

rename_dmr_output<-setNames(nm=1:length(unique(collapsed_dmrs$test)),gsub(colnames(sum_mat)[grepl(colnames(sum_mat),pattern="_t$")],pattern="_t",replacement=""))
collapsed_dmrs$type <- rename_dmr_output[collapsed_dmrs$test]
saveRDS(collapsed_dmrs,file=paste0(dmr_celltype_outdir,"/","coarse_cluster",".dmr_filt_collapse.rds"))
collapsed_dmrs<-readRDS(file=paste0(dmr_celltype_outdir,"/","coarse_cluster",".dmr_filt_collapse.rds"))

#plot dmr counts per clone
pal <- colorRampPalette(colors = c("#F9AB60", "#E7576E", "#630661", "#B5DCA5"))
COLS <- pal(length(unique(collapsed_dmrs$type)))

plt<-ggplot(collapsed_dmrs |> dplyr::group_by(type, direction) |> dplyr::summarise(n = n()), 
    aes(y = type, x = n, fill = type)) + 
    geom_col() + 
    facet_grid(vars(direction), scales = "free_y") + 
    scale_fill_manual(values = COLS) + 
    theme_classic()
ggsave(plt,file=paste0(dmr_celltype_outdir,"/","coarse_cluster",".dmr_counts.pdf"))

#filtering dmrs
#use only hypo dmrs, take top 5000 per celltype, filter by length and merge any that overlap (across cell types)
dmrs <- collapsed_dmrs %>% 
  filter(direction=="hypo") %>% 
  filter(dmr_padj<0.05) %>% 
  filter(dmr_length < 20000) %>% 
  group_by(type) %>% 
  #slice_min(n=20000, dmr_logFC) %>%
  filter(abs(dmr_logFC)>1.5) %>%
  as.data.frame() %>%
  GenomicRanges::makeGRangesFromDataFrame() %>%   GenomicRanges::reduce()

#25,808 dmrs

summary(width(dmrs))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    501    6501   11501   12040   16501   91001

length(dmrs)
dmr_bed<-data.frame(chr=seqnames(dmrs),start=start(dmrs),end=end(dmrs))


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
cell_markers[["basal"]]<-c("CARMN","ACTA2","KRT17","KRT14","DST","KRT5")
cell_markers[["lumhr"]]<-c("AREG","AZGP1","KRT18","AGR2","PIP","ANKRD30A")
cell_markers[["lumsec"]]<-c("GABRP","ELF5","CCL28","KRT15","KIT","MMP7","LTF","SLPI")
cell_markers[["fibro"]]<-c("DCN","APOD","LUM","COL1A2","COL1A1","FAP")
cell_markers[["endo"]]=c("CCL21","TFF3","MMRN1","CLDN5","AL357507.1","PKHD1L1","KLHL4","LINC02147","RHOJ","ST6GALNAC3","MMRN1","MECOM","BTNL9","MCTP1","PTPRB","VWF","ADGRL4","LDB2")
cell_markers[["perivasc"]]=c("RGS6","KCNAB1","COL25A1","ADGRL3","PRKG1","MYL9","ADIRF","NR2F2-AS1","AC012409.2")
cell_markers[["myeloid"]]<-c("HLA-DRA","HLA-DPA1","CD74")
cell_markers[["tcell"]]<-c("PTPRC","IKZF1","IL7R","GNLY")
cell_markers[["mast"]]<-c("NTM","SYTL3","SLC24A3","TPSB2","HDC")
cell_markers[["bcell"]]=c("CD37","TCL1A","LTB","HLA-DPB1","HLA-DRA","HLA-DPA1")
cell_markers[["plasma"]]=c("IGHA2","IGHA1","JCHAIN","IGHM","IGHG1","IGHG4","IGHG3","IGHG2")
cell_markers[["adipo"]]=c("PDE3B","ACACB","WDPCP","PCDH9","CLSTN2","ADIPOQ","TRHDE")

cell_colors=c(
"basal"="#844c9d",
"lumhr"="#e23e96",
"lumsec"="#ff6498",
"fibro"="#f58e90",
"endo"="#5bbb5a",
"perivasc"="#eaba67",
"myeloid"="#8088c2",
"tcell"="#1d87c8",
"mast"="#dcd0ff",
"bcell"="#65cbe4",
"plasma"="#7ecdc2",
"adipo"="#b48454")


#prepare cgi
cgisland="/data/rmulqueen/projects/scalebio_dcis/ref/cpgIslandExt.bed"
cgi<-rtracklayer::import(cgisland)
cgi<-as.data.frame(cgi)
colnames(cgi)<-c("chr","start","end","strand")

plot_histogram_page<-function(celltype){
    plt<-histograModified(dat, 
        baseline="mean",
        genes = unlist(cell_markers[celltype]),
        colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
        matrix = "cg_cluster_id_tracks", arrowScale = .03, trackScale = .5,
        legend = F, cgisland=cgi)
    return(plt)
}

#rough ordering by what clustered together
order<-c(
"5","1","7",
"17","18",
"2",
"10","9",
"8",
"19","6","20",
"14",
"12","11","13",
"16","15",
"4","3")
    
plt_list<-lapply(names(cell_colors),
function(celltype){
    genes<-unlist(cell_markers[celltype])
    genes<-genes[genes %in% dat@ref$gene_name]
    print(paste("Plotting:",celltype))
    print(paste("Genes to plot:",genes))
    plt<-histograModified(dat, 
        baseline="mean",
        genes = unlist(cell_markers[celltype]),
        colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
        matrix = "cg_coarse_cluster_tracks", arrowScale = .03, trackScale = .5,
        legend = F, cgisland=cgi,order=order) + ggtitle(celltype)
    return(plt)
})

plt_out<-wrap_plots(plt_list,ncol=1) 
ggsave(plt_out ,file="05_scaledcis.coarse_clusters.marker.pdf",width=30,height=length(cell_colors)*30,limitsize=F)


#assigning coarse cell types
dat@metadata$broad_celltype<-"lumhr" #2, 10, 9, 8, 19, 6, 20, i suspect 10 and 9 are normal, rest are cancer
dat@metadata[dat@metadata$coarse_cluster_id %in% c("12","13","11"),]$broad_celltype<-"basal"
dat@metadata[dat@metadata$coarse_cluster_id %in% c("15","16"),]$broad_celltype<-"lumsec"
dat@metadata[dat@metadata$coarse_cluster_id %in% c("4","3"),]$broad_celltype<-"fibroblast"
dat@metadata[dat@metadata$coarse_cluster_id %in% c("14"),]$broad_celltype<-"endothelial"
dat@metadata[dat@metadata$coarse_cluster_id %in% c("18","17","7","1","5"),]$broad_celltype<-"immune"

saveRDS(dat,file="05_scaledcis.coarse_clusters.amethyst.rds")

```

## Add DMR site information per coarse cluster
This will be used 

```R
#################################
#Get DMR per Celltypes
#################################

dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")
dmr_celltype_outdir=paste(sep="/",dmr_outdir,"dmr_coarse_cluster")
collapsed_dmrs<-readRDS(file=paste0(dmr_celltype_outdir,"/","coarse_cluster",".dmr_filt_collapse.rds"))

#reading in pre-computed windows per cell type, to determine DMRs
#from scalemet_dcis/processing/milestonev1_08_DiffMet_analysis.md

#filter DMRS that are different by cell type to
# significant 
# gene overlap
# hypomethylated
# logFC atleast 1.5

table(collapsed_dmrs$type)
#   1     10     11     12     13     14     15     16     17     18     19 
#184447 185076 182009 177556 175150 173075 184027 167856 178465 157240 173844 
#     2     20      3      4      5      6      7      8      9 
#158315 161403 181188 176171 182814 171735 173040 137075 182240 

dmrs_for_subtyping <- collapsed_dmrs %>% 
                          filter(direction=="hypo") %>%
                          filter(dmr_padj < 0.05) %>%
                          filter(abs(dmr_logFC) > 1.5) %>% 
                          filter(dmr_length<50000) %>% 
                          group_by(type) 

table(dmrs_for_subtyping$type)
#  1   10   11   12   13   14   15   16   17   18   19    2   20    3    4    5 
#1385 3330 3714 4934 4394 3913 5510 6089 2707 2945  445 3208 4929 3931 3056 1426 
#   6    7    8    9 
#4879 1119 2147 3131 

#merge dmrs that overlap regardless of cluster
dmr_out<-dmrs_for_subtyping %>% 
        as.data.frame() %>% 
        select(chr,dmr_start,dmr_end) %>% 
        distinct(chr,dmr_start,dmr_end) %>% 
        makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>%
        reduce()

dmr_out<-dmr_out[width(dmr_out)<50000,]
summary(width(dmr_out))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    501    8501   18501   21103   32501   49501 
#22713 ranges

#convert to bed file for calculation of windows
dmr_bed<-data.frame(chr=seqnames(dmr_out),start=start(dmr_out),end=end(dmr_out))


#create new matrix from DMR sites for refined clustering
dat@genomeMatrices[["coarse_cluster_dmr_sites"]] <- makeWindows(dat, 
                                                     bed = dmr_bed,
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 50, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

saveRDS(dat,file="05_scaledcis.coarse_clusters.amethyst.rds")

# Use these coarse_cluster_dmr_sites methylation sites for subclustering for immune and stromal cells.

```

<!--
```R

#run 3d plotting just for fun
dim3d<-uwot::umap(
  X=dat@reductions$nakshatri_dmr_sites_irlba_regressed,
  n_neighbors = 15,
  n_components = 3,
  metric = "euclidean",
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

write.table(dim3d,col.names=T,file="03_scaledcis.broad_celltypes.3dumap.csv",sep=",",row.names=F)

```
-->

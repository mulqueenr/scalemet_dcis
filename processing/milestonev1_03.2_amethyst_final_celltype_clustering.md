# Final cell typing

Cell type assignment by order of importance:
1. CNVs identified by copy number clones assigned cancer in scalemet_dcis/processing/milestonev1_03.2_copykit_cnv_calling.md
2. Cluster based on DMRs in scalemet_dcis/processing/milestonev1_04.1_liger_scrna_met_integration.md assigned cell types after integration
3. Assign cell types for clusters with clear (>75% of cell type consistently assigned)
4. Identify clusters which group together on umap at lower cluster resolution (assign celltype by shared cluster)
5. Confirm cell types by RNA based marker genes and hypomethylation over gene body

Use integration to finalize cell typing. Note some cell types were mislabelled based on small set of marker genes originally used. 
Using integration data for final classifications.

```R
BiocManager::install("HDF5Array")
install.packages('leidenAlg')
install.packages('rliger')
```


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
options(future.globals.maxSize= 200000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)

obj<-readRDS(file="07_scaledcis.integrated_celltyping.amethyst.rds")

#################################
#Get DMR per Integrated Celltypes
#################################
dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")
system(paste("mkdir -p", dmr_outdir))

#celltype
dmr_celltype_outdir=paste(sep="/",dmr_outdir,"integrated_celltype")
system(paste("mkdir -p", dmr_celltype_outdir))


celltype500bpwindows <- calcSmoothedWindows(obj, 
                                        type = "CG", 
                                        threads = 5,
                                        step = 500, # change to 500 for real data unless you have really low coverage
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "integrated_celltype",
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)
saveRDS(celltype500bpwindows,file=paste0(dmr_celltype_outdir,"/","dmr_analysis.integrated_celltype.500bp_windows.rds"))
celltype500bpwindows<-readRDS(file=paste0(dmr_celltype_outdir,"/","dmr_analysis.integrated_celltype.500bp_windows.rds"))

pct_mat<-celltype500bpwindows[["pct_matrix"]] 
sum_mat<-celltype500bpwindows[["sum_matrix"]] 

#save clone object for future genome track plotting
obj@genomeMatrices[["cg_celltype_tracks"]] <- pct_mat #load it into amethyst object for plotting

dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
        eachVsAll = TRUE, # If TRUE, each group found in the sumMatrix will be tested against all others
        nminTotal = 3, # Min number observations across all groups to include the region in calculations
        nminGroup = 3) # Min number observations across either members or nonmembers to include the region

saveRDS(dmrs,file=paste0(dmr_celltype_outdir,"/","integrated_celltype_allcells",".dmr.rds"))

dmrs <- filterDMR(dmrs, 
            method = "bonferroni", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
            filter = FALSE, # If TRUE, removes insignificant results
            pThreshold = 0.05, # Maxmimum adjusted p value to allow if filter = TRUE
            logThreshold = 1.5) # Minimum absolute value of the log2FC to allow if filter = TRUE

collapsed_dmrs <- collapseDMR(obj, 
                       dmrs, 
                        maxDist = 2000, # Max allowable overlap between DMRs to be considered adjacent
                        minLength = 2000, # Min length of collapsed DMR window to include in the output
                        reduce = T, # Reduce results to unique observations (recommended)
                        annotate = T) # Add column with overlapping gene names

saveRDS(collapsed_dmrs,file=paste0(dmr_celltype_outdir,"/","integrated_celltype_allcells",".dmr_filt_collapse.rds"))

rename_dmr_output<-setNames(nm=1:length(unique(collapsed_dmrs$test)),gsub(colnames(sum_mat)[grepl(colnames(sum_mat),pattern="_t$")],pattern="_t",replacement=""))
collapsed_dmrs$type <- rename_dmr_output[collapsed_dmrs$test]
saveRDS(collapsed_dmrs,file=paste0(dmr_celltype_outdir,"/","integrated_celltype_allcells",".dmr_filt_collapse.rds"))

collapsed_dmrs<-readRDS(file=paste0(dmr_celltype_outdir,"/","integrated_celltype_allcells",".dmr_filt_collapse.rds"))

#plot dmr counts per clone
pal <- colorRampPalette(colors = c("#F9AB60", "#E7576E", "#630661", "#B5DCA5"))
COLS <- pal(length(unique(collapsed_dmrs$type)))

plt<-ggplot(collapsed_dmrs |> dplyr::group_by(type, direction) |> dplyr::summarise(n = n()), 
    aes(y = type, x = n, fill = type)) + 
    geom_col() + 
    facet_grid(vars(direction), scales = "free_y") + 
    scale_fill_manual(values = COLS) + 
    theme_classic()
ggsave(plt,file=paste0(dmr_celltype_outdir,"/","celltype_allcellsm",".dmr_counts.pdf"))

#clear that plasma and epithelial is too few cell count
# rename them when happy with umap projections

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


summary(width(dmrs))
length(dmrs)
dmr_bed<-data.frame(chr=seqnames(dmrs),start=start(dmrs),end=end(dmrs))


#create new matrix from DMR sites for refined clustering
window_name="integrated_celltype_dmr_sites"


obj@genomeMatrices[[window_name]] <- makeWindows(obj, 
                                                     bed = dmr_bed,
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 50, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

obj@genomeMatrices[[window_name]] <- obj@genomeMatrices[[window_name]][rowSums(!is.na(obj@genomeMatrices[[window_name]])) >= 45, ]
est_dim<-dimEstimate(obj, genomeMatrices = c(window_name), dims = c(20), threshold = 0.95)
print(est_dim)

saveRDS(obj,file="08_scaledcis.final_celltyping.amethyst.rds")



celltype_umap<-function(obj=obj,prefix="allcells",dims=12,regressCov=TRUE,k_pheno=50,neigh=25,dist=1e-5,method="cosine"){

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

  p1 <- dimFeature(obj, colorBy = integrated_celltype, reduction = "umap") + ggtitle(paste(window_name,"Cell Types"))
  p2 <- dimFeature(obj, colorBy = sample, reduction = "umap") + ggtitle(paste(window_name,"Samples"))
  p3 <- dimFeature(obj, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red"),guide="colourbar") + ggtitle("Coverage distribution")
  p4 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
  p5 <- dimFeature(obj, colorBy = cluster_id, reduction = "umap") + ggtitle(paste(window_name,"Clusters"))
  p6 <- dimFeature(obj, colorBy = fine_celltype, reduction = "umap") + ggtitle(paste(window_name,"Old Cell Types"))
  p7 <- dimFeature(obj, colorBy = Group, reduction = "umap") + ggtitle(paste(window_name," Group"))
  p8<-ggplot()
  ggsave((p1|p2)/(p3|p4)/(p5|p6)/(p7|p8),file=paste0(prefix,"_",outname,"_umap.pdf"),width=20,height=40)  
  
  #plot barplot of clusterid to integrated_celltype and to fine_celltype
  plt_dat<-obj@metadata %>% select(cluster_id,fine_celltype,integrated_celltype,Group) 

  plt1<-ggplot(plt_dat, aes(fill=fine_celltype, x=cluster_id)) + 
      geom_bar(position="stack", stat="count")
  plt2<-ggplot(plt_dat, aes(fill=integrated_celltype, x=cluster_id)) + 
      geom_bar(position="stack", stat="count")
  plt3<-ggplot(plt_dat, aes(fill=Group, x=cluster_id)) + 
      geom_bar(position="stack", stat="count")

  ggsave(plt1/plt2,file=paste0(prefix,"_cluster_celltype.pdf"),width=20,height=10) 
  return(obj)
  }




#run on all cells
obj_allcells<-obj
#obj<-obj_allcells
obj<-celltype_umap(dims=13,regressCov=FALSE,k_pheno=15,neigh=8,dist=0.001,method="cosine") #13 dims, no regress cov, k pheno to overcluster, modify neighbors to clean?
#neigh 5 looks good
#14 dims starts splitting out cancer clones

celltype_dat<-obj@metadata %>% select(cluster_id,fine_celltype,integrated_celltype,Group) 
celltype_dat<-table(celltype_dat$cluster_id,celltype_dat$integrated_celltype)
celltype_dat<-round(celltype_dat/rowSums(celltype_dat),2)

Group_dat<-obj@metadata %>% select(cluster_id,fine_celltype,integrated_celltype,Group) 
Group_dat<-table(Group_dat$cluster_id,Group_dat$Group)
Group_dat<-round(Group_dat/rowSums(Group_dat),2)

fine_celltype_dat<-obj@metadata %>% select(cluster_id,fine_celltype,integrated_celltype,Group) 
fine_celltype_dat<-table(fine_celltype_dat$cluster_id,fine_celltype_dat$fine_celltype)
fine_celltype_dat<-round(fine_celltype_dat/rowSums(fine_celltype_dat),2)

sample_dat<-obj@metadata %>% select(cluster_id,Sample) 
sample_dat<-table(sample_dat$cluster_id,sample_dat$Sample)
sample_dat<-round(sample_dat/rowSums(sample_dat),2)

plt_dat<-cbind(celltype_dat,Group_dat,fine_celltype_dat,sample_dat)
write.table(plt_dat,file=paste0(prefix,"_cluster_celltype.csv"),row.names=T,col.names=T,sep=",")

#based on these data, we can define epithelial cells well:
#going to subcluster nonepithelial cells for final celltypes
obj@metadata$epithelial<-"nonepithelial"
obj@metadata[obj@metadata$cluster_id %in% c("12","13"),]$epithelial<-"basal"
obj@metadata[obj@metadata$cluster_id %in% c("15","18","16"),]$epithelial<-"lumsec"
obj@metadata[obj@metadata$cluster_id %in% c("2"),]$epithelial<-"lumhr"
obj@metadata[obj@metadata$cluster_id %in% c("10","7","21","8","6","9","22","23","24"),]$epithelial<-"cancer"

#anything with label "cancer" is from cnv calling, so labelling that as well
obj@metadata[obj@metadata$cnv_ploidy_500kb=="aneuploid",]$epithelial<-"cancer"

#from this plot it is clear that epithelial markers are dominating signal
#splitting to cancer and non cancer (or better yet, epithelial and non epithelial)

#non epithelial cells, filter DMRS to non-epithelial cells
obj_nonepi<-subsetObject(obj,cells=row.names(obj@metadata[obj@metadata$epithelial=="nonepithelial",]))
est_dim<-dimEstimate(obj_nonepi, genomeMatrices = c(window_name), dims = c(20), threshold = 0.95)
est_dim
obj_nonepi<-celltype_umap(obj=obj_nonepi,prefix="nonepi",dims=13,regressCov=FALSE,k_pheno=15,neigh=8,dist=1E-5,method="cosine") 


# epithelial cells
obj_epi<-subsetObject(obj,cells=row.names(obj@metadata[obj@metadata$epithelial!="nonepithelial",]))
est_dim<-dimEstimate(obj_epi, genomeMatrices = c(window_name), dims = c(20), threshold = 0.95)
est_dim
obj_epi<-celltype_umap(obj=obj_epi,prefix="epi",dims=10,regressCov=FALSE,k_pheno=15,neigh=8,dist=1E-5,method="cosine") 

#summarize in final object and plot markers


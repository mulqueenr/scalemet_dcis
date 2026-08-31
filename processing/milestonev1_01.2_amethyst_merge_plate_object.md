# Merge plate objects and cluster cells
After plate objects are initiated, cell lines filtered out, perform initial clustering.

1. Merge amethyst objects by plate.
2. Split out cell line from patient objects.
3. Add apriori metadata
4. Perform initial clustering via iterative PCA into UMAP
5. Identify doublets by initial cluster, via artificial doublet generation and random forest model
6. Remove identified doublets.
7. Recluster to final coarse cluster for cell typing.
8. Identify cell types by coarse marker genes via IGV bigwig output.
9. Save final objects.

```R
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
library(tidyverse)
library(irlba)
library(uwot)
library(readr)
library(Matrix)
library(igraph)  


set.seed(111)
options(future.globals.maxSize= 2500000*1024^2) #80gb limit for parallelizing

task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
processing_folder="01_amethyst_initial_object"
wd=paste(sep="/",project_data_directory,processing_folder)

system(paste("mkdir -p",wd))
setwd(wd)

#make merged object
plate_obj=paste(sep="/",project_data_directory,"scalemethyl_pipeline_out/amethyst_plate_obj")
amethyst_files=list.files(path=plate_obj,pattern="*.amethyst.rds",recursive=TRUE,full.names=TRUE)

##########################################################################################
######################1. Merge amethyst objects by plate.#################################
##########################################################################################

# Read in all amethyst files
window_name="cg_5mb_score"
dat_list<-mclapply(amethyst_files, function(x) {
    obj<-readRDS(x)
    return(obj)},mc.cores=50)

dat <- combineObject(objList = dat_list, genomeMatrices=window_name)

# Cluster via methscan processing on VMRs

#matrix columns
#  "row_i": row_i + 1,  # 1-indexed indices, like cellranger output
#"col_i": col_i + 1 + region_n_offset,
#"residuals": residuals,
#"mfracs": mfracs,
#"coverage": coverage,

#https://github.com/satijalab/seurat/issues/9711

dt <- readr::read_delim(
        file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/methscan/VMR_matrix/matrix.mtx.gz",
        delim=' ',
        trim_ws=TRUE,
        num_threads=200,
        col_names=FALSE,
        col_types="iiddd")

dims = c(max(dt$X1), max(dt[,2]))

meth_mtx <- new('dgTMatrix', 
            i = dt$X1 - 1L, 
            j = dt$X2 - 1L, 
            x = dt$X3,
            Dim=dims)


rownames(meth_mtx) <-readLines("/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/methscan/VMR_matrix/barcodes.tsv.gz")
colnames(meth_mtx) <- readLines("/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/methscan/VMR_matrix/features.tsv.gz")

#correct meth_mtx row.names and save
row.names(meth_mtx)<-paste(
  unlist(lapply(strsplit(row.names(meth_mtx),"[.]"),"[",4)),
  unlist(lapply(strsplit(row.names(meth_mtx),"[.]"),"[",1)),
  sep="+")

dim(meth_mtx)
dim(dat@metadata)
row.names(meth_mtx) %in% row.names(dat@metadata) %>% sum()

#those that are not overlapping did not pass filters
meth_mtx<-meth_mtx[row.names(meth_mtx) %in% row.names(dat@metadata),]
saveRDS(meth_mtx,file="01.0.VMR.residuals.meth_mat.rds")
saveRDS(dat,file="01.0.amethyst.initialobj.rds")

#add matrix to amethyst object
dat@genomeMatrices[["vmr_matrix_cg_residuals"]]<-t(meth_mtx[row.names(meth_mtx) %in% row.names(dat@metadata),])


##########################################################################################
######################2. Split out cell line from patient objects.########################
##########################################################################################

#split out cell lines and tissue as separate objects
#cell line
dat_cellline<-subsetObject(dat,cells=row.names(dat@metadata[dat@metadata$sample %in% c("MDA-MB-231","MCF7","MCF10A"),]))
dat_cellline@metadata$Sample<-dat_cellline@metadata$sample
dat_cellline@metadata$cellid<-row.names(dat_cellline@metadata)
saveRDS(dat_cellline,file="01.0.cellline.amethyst.rds")
dat_cellline<-readRDS(file="01.0.cellline.amethyst.rds")

plt<-ggplot(dat_cellline@metadata,aes(x=Sample,y=mcg_pct,color=Sample),)+geom_jitter()+geom_boxplot(outlier.shape=NA,color="black")

ggsave(plt,file="./plot/cellline_methylation.pdf")
dat_cellline@metadata %>% group_by(Sample) %>% summarize(cg_cov=mean(cg_cov),cg_perc=mean(mcg_pct))

##########################################################################################
######################3. Add apriori metadata#############################################
##########################################################################################

#patient data
dat_patient<-subsetObject(dat,cells=row.names(dat@metadata[!(dat@metadata$sample %in% c("MDA-MB-231","MCF7","MCF10A")),]))

#correct sample names for consistency across experiments and RNA
#and add patient specific metadata
sample_names<-c(
  'BCMDCIS05T'='BCMDCIS05T',
  'BCMDCIS07T'='BCMDCIS07T',
  'BCMDCIS102T-4h'='BCMDCIS102T_24hTis',
  'BCMDCIS124T'='BCMDCIS124T',
  'BCMDCIS22T'='BCMDCIS22T',
  'BCMDCIS28T'='BCMDCIS28T',
  'BCMDCIS32T'='BCMDCIS32T',
  'BCMDCIS35T-3h'='BCMDCIS35T',
  'BCMDCIS41T'='BCMDCIS41T',
  'BCMDCIS49T-24hTis'='BCMDCIS49T',
  'BCMDCIS52T'='BCMDCIS52T',
  'BCMDCIS65T'='BCMDCIS65T',
  'BCMDCIS66T'='BCMDCIS66T',
  'BCMDCIS70T'='BCMDCIS70T',
  'BCMDCIS74T'='BCMDCIS74T',
  'BCMDCIS80T'='BCMDCIS80T_24hTis',
  'BCMDCIS82T-24hTis'='BCMDCIS82T_24hTis',
  'BCMDCIS94T-24hTis'='BCMDCIS94T_24hTis',
  'BCMDCIS97T'='BCMDCIS97T',
  'BCMDCIS99T'='BCMDCIS99T',
  'BCMHBCA03R'='BCMHBCA03R',
  'BCMHBCA04R'='BCMHBCA04R',
  'BCMHBCA09R-3h'='BCMHBCA09R-3h',
  'BCMHBCA12R-3h'='BCMHBCA12R-3h',
  'BCMHBCA22R-4h'='BCMHBCA22R-4h',
  'BCMHBCA26L-24hTis-4h'='BCMHBCA26L-24hTis-4h',
  'BCMHBCA29L-2h'='BCMHBCA29L-2h',
  'BCMHBCA38L-3h'='BCMHBCA38L-3h',
  'BCMHBCA85L-3h'='BCMHBCA85L-3h',
  'DCIS-66T'='BCMDCIS66T',
  'DCIS-79T'='BCMDCIS79T_24hTis_DCIS',
  'DCIS-92T'='BCMDCIS92T_24hTis',
  'ECIS25T'='ECIS25T',
  'ECIS26T'='ECIS26T',
  'ECIS36T'='ECIS36T',
  'ECIS48T'='ECIS48T',
  'ECIS57T'='ECIS57T',
  'HBCA-16R'='BCMHBCA16R-3h',
  'HBCA-17T'='BCMHBCA17R-3h',
  'HBCA-19T'='BCMHBCA19R-4h',
  'HBCA-83L'='BCMHBCA83L-3h',
  'IDC-79T'='BCMDCIS79T_24hTis_IDC')

dat_patient@metadata$Sample<-sample_names[dat_patient@metadata$sample]
dat_patient@metadata$cellid<-row.names(dat_patient@metadata)
patient_metadata<-read.table("/data/rmulqueen/projects/scalebio_dcis/sample_selection/simplified_patient_metadata.csv",header=T,sep=",")
patient_metadata<-patient_metadata[!duplicated(patient_metadata$Sample),]

dat_patient@metadata<-base::merge(dat_patient@metadata,patient_metadata,by="Sample",all.x=T)
row.names(dat_patient@metadata)<-dat_patient@metadata$cellid

saveRDS(dat_patient,file="01.0.patient.amethyst.rds")
dat_patient<-readRDS(file="01.0.patient.amethyst.rds")




##########################################################################################
######################4-9###########
##########################################################################################

# PCA that iteratively imputes missing values

```R
prcomp_iterative <- function(x, n=10, n_iter=50, min_gain=0.001, ...) {
  mse <- rep(NA, n_iter)
  na_loc <- is.na(x)
  x[na_loc] = 0  # zero is our first guess

  for (i in 1:n_iter) {
    prev_imp <- x[na_loc]  # what we imputed in the previous round
    # PCA on the imputed matrix
    pr <- prcomp_irlba(x, center = F, scale. = F, n = n, ...)
    # impute missing values with PCA
    new_imp <- (pr$x %*% t(pr$rotation))[na_loc]
    x[na_loc] <- new_imp
    # compare our new imputed values to the ones from the previous round
    mse[i] = mean((prev_imp - new_imp) ^ 2)
    # if the values didn't change a lot, terminate the iteration
    gain <- mse[i] / max(mse, na.rm = T)
    if (gain < min_gain) {
      message(paste(c("\n\nTerminated after ", i, " iterations.")))
      break
    }
  }
  pr$mse_iter <- mse[1:i]
  pr
}

vmr_cluster<-function(obj=dat_cellline,suffix="cellline",
                      outdir="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/01_amethyst_initial_object/",
                      npc=15,
                      reduction_name="irlba",
                      leiden_cluster_resolution=3e-6,
                      pheno_k=200){

  system(paste0("mkdir -p ",outdir,"/plot"))
  #run clustering on doublet object, code same as initial clustering code
  print("Calculating PCA...")
  pca <- t(obj@genomeMatrices[["vmr_matrix_cg_residuals"]]) %>%
      scale(center = T, scale = F) %>%
      prcomp_iterative(n = npc)  # this is number of principle components (npc in args)

  pca_dims <- pca$x %>% 
    magrittr::set_rownames(colnames(obj@genomeMatrices[["vmr_matrix_cg_residuals"]]))
  
  print("Running UMAP...")
  umap_obj <- uwot::umap(X=as.matrix(pca_dims), 
                          metric="cosine", 
                          min_dist=0.05, 
                          n_neighbors=5, 
                          seed=2, 
                          ret_nn=T)

  umap_tbl <- umap_obj$embedding %>% as.data.frame() %>% 
    mutate(cell=row.names(umap_obj$embedding)) %>%
    magrittr::set_colnames(c("UMAP1", "UMAP2","cell")) 

  # get the edges of the neighbor graph from the UMAP object
  neighbor_graph_edges <- 
    tibble(from = rep(1:nrow(umap_obj$nn$cosine$idx), times=ncol(umap_obj$nn$cosine$idx)),
          to = as.vector(umap_obj$nn$cosine$idx),
          weight = as.vector(umap_obj$nn$cosine$dist)) %>%
    filter(from != to) %>%
    mutate(from = colnames(obj@genomeMatrices[["vmr_matrix_cg_residuals"]])[from],
          to = colnames(obj@genomeMatrices[["vmr_matrix_cg_residuals"]])[to])

  # run Leiden clustering
  print("Leiden clustering...")
  clust_obj <- neighbor_graph_edges %>%
    igraph::graph_from_data_frame(directed=F) %>% 
    igraph::cluster_leiden(resolution_parameter = leiden_cluster_resolution)  # adjust the resolution parameter to your needs, overclustering to identify more doublet clusters

  # put the clustering results into a data frame (tibble) for plotting
  clust_tbl <- tibble(
    leiden_cluster = as.character(clust_obj$membership),
    cell= clust_obj$names) %>% 
    full_join(umap_tbl, by="cell")

  print("Adding PCA to reduction slot...")
  obj@reductions[[reduction_name]]<-pca_dims

  print("Adding clusters to metadata...")
  coarse_cluster_leidenclus<-setNames(nm=clust_tbl$cell,clust_tbl$leiden_cluster)
  coarse_cluster_UMAP1<-setNames(nm=clust_tbl$cell,clust_tbl$UMAP1)
  coarse_cluster_UMAP2<-setNames(nm=clust_tbl$cell,clust_tbl$UMAP2)

  obj@metadata$coarse_cluster_UMAP_X<-coarse_cluster_UMAP1[row.names(obj@metadata)]
  obj@metadata$coarse_cluster_UMAP_Y<-coarse_cluster_UMAP2[row.names(obj@metadata)]
  obj@metadata$coarse_cluster_leidenclus<-coarse_cluster_leidenclus[row.names(obj@metadata)]

  #cluster on umap as well
  umap_clus<-obj@metadata %>% 
    select(c("coarse_cluster_UMAP_X","coarse_cluster_UMAP_Y")) %>% 
    magrittr::set_rownames(row.names(obj@metadata))

  umap_clus<-Rphenograph::Rphenograph(umap_clus,k=pheno_k)

  obj@metadata$coarse_cluster_phenograph<-as.character(unlist(as.list(igraph::membership(umap_clus[[2]]))))
  print("Plotting...")
  if("coarse_cluster_leidenclus" %in% colnames(obj@metadata)){
  plt_clus<-obj@metadata %>% 
    ggplot(aes(x = coarse_cluster_UMAP_X, y = coarse_cluster_UMAP_Y, color = coarse_cluster_leidenclus)) +
    geom_point() +
    coord_fixed()
      ggsave(plt_clus,file=paste0(outdir,"/plot/","01.0.VMR_umap.",suffix,".clus.leiden.pdf"),width=20,height=20)
  }

  if("coarse_cluster_phenograph" %in% colnames(obj@metadata)){
  plt_clus<-obj@metadata %>% 
    ggplot(aes(x = coarse_cluster_UMAP_X, y = coarse_cluster_UMAP_Y, color = coarse_cluster_phenograph)) +
    geom_point() +
    coord_fixed()
      ggsave(plt_clus,file=paste0(outdir,"/plot/","01.0.VMR_umap.",suffix,".clus.pheno.pdf"),width=20,height=20)
  }

  if("Sample" %in% colnames(obj@metadata)){
  plt_sample<-obj@metadata %>% 
    ggplot(aes(x = coarse_cluster_UMAP_X, y = coarse_cluster_UMAP_Y, color = Sample)) +
    geom_point() +
    coord_fixed()
      ggsave(plt_sample,file=paste0(outdir,"/plot/","01.0.VMR_umap.",suffix,".sample.pdf"),width=20,height=20)
  }

  if("mcg_pct" %in% colnames(obj@metadata)){
  plt_mcg<-obj@metadata %>% 
    ggplot(aes(x = coarse_cluster_UMAP_X, y = coarse_cluster_UMAP_Y, color = mcg_pct)) +
    geom_point() +
    coord_fixed()
  ggsave(plt_mcg,file=paste0(outdir,"/plot/","01.0.VMR_umap.",suffix,".mcg.pdf"),width=20,height=20)
  }

  if("unique_reads" %in% colnames(obj@metadata)){
  plt_reads<-obj@metadata %>% 
    ggplot(aes(x = coarse_cluster_UMAP_X, y = coarse_cluster_UMAP_Y, color = log10(unique_reads))) +
    geom_point() +
    coord_fixed()
  ggsave(plt_reads,file=paste0(outdir,"/plot/","01.0.VMR_umap.",suffix,".reads.pdf"),width=20,height=20)
  }

  if("doublet_score" %in% colnames(obj@metadata)){
  plt_reads<-obj@metadata %>% 
    ggplot(aes(x = coarse_cluster_UMAP_X, y = coarse_cluster_UMAP_Y, color = doublet_score)) +
    geom_point() +
    coord_fixed()
  ggsave(plt_reads,file=paste0(outdir,"/plot/","01.0.VMR_umap.",suffix,".doublet_score.pdf"),width=20,height=20)
  }
  return(obj)
}

generate_doublet_score<-function(obj=dat_patient,
  suffix="patient",
  reduction_name="irlba",
  outdir="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/01_amethyst_initial_object/",
  cpus=100){
  ##make simulated doublets for detection
  print("Making artificial doublets...")
  dbobj <- makeDoubletObject(obj, simFraction=0.25, threads = cpus, genomeMatrices=c("vmr_matrix_cg_residuals"))
  dbobj@metadata$cellid<-row.names(dbobj@metadata)

  print("Clustering with artificial doublets...")
  dbobj<-vmr_cluster(obj=dbobj) #generate reduction on doublet object

  print("Predicting true doublets...")
  result <- buildDoubletModel(dbobj = dbobj, method="rf", reduction = reduction_name)
  dbobj <- predictDoubletScores(dbobj = dbobj, model = result$model, reduction = reduction_name)

  print("Plotting doublets...")
  plt<-ggplot(dbobj@metadata,aes(x=doublet_info,y=doublet_score,color=doublet_info))+geom_violin(fill=NA)+geom_jitter(size=0.2)+theme_minimal()
  ggsave(plt,file=paste0(outdir,"/plot/","01.0.VMR_predicted.doublet_score.",suffix,".violin.pdf"))

  doubletscore<-setNames(nm=row.names(dbobj@metadata),dbobj@metadata$doublet_score)
  obj@metadata$doublet_score<-unname(doubletscore[row.names(obj@metadata)])

  plt<-ggplot(obj@metadata,aes(x=coarse_cluster_leidenclus,y=doublet_score,color=coarse_cluster_leidenclus))+geom_violin(fill=NA)+geom_jitter(size=0.2)+theme_minimal()
  ggsave(plt,file=paste0(outdir,"/plot/","01.0.VMR_doublet_score_per_cluster.",suffix,".violin.pdf"))
  return(obj)
}

write_binned_bigwigs<-function(celltype_tracks=celltype_tracks,
                                outdir=celltype_outdir,
                                i){
    #split bw into 4 files per cell type by methylation/average methylation
    out_dat<-celltype_tracks %>% select(chr,start,end,i) 

    hg38_seq_info<-Seqinfo(genome="hg38")
    out_dat<-GRanges(out_dat[complete.cases(out_dat),]) #filter NA
    out_dat<-out_dat[out_dat@seqnames %in% hg38_seq_info@seqnames,] #filter chr
    out_dat<-resize(out_dat,width=500)
    names(out_dat@elementMetadata)<-"score"
    mean_score<-mean(mcols(out_dat)$score)

    #bin to 100-mean, mean-50, 50-20, 20-0  
    #color black, grey, lightgrey, celltypecol
    #333333, #444444, #bcbcbc, celltypecol
    #subtract score-meanscorevalue

    out_dat_hypermet <- out_dat %>% 
                        as.data.frame() %>% 
                        filter(mcols(out_dat)$score > mean_score) %>% 
                        mutate(score=score-mean_score) %>% 
                        GRanges() %>% unique()
    names(out_dat_hypermet@elementMetadata)<-"score"
    genome(out_dat_hypermet)<-"hg38"
    seqlengths(out_dat_hypermet)<-as.data.frame(hg38_seq_info)[hg38_seq_info@seqnames %in% out_dat_hypermet@seqnames,]$seqlengths #filter by seqlengths

    out_dat_met_mid <- out_dat %>% 
                      as.data.frame() %>% 
                      filter(mcols(out_dat)$score <= mean_score & mcols(out_dat)$score > 50) %>% 
                      mutate(score=score-mean_score) %>% 
                      GRanges() %>% unique()
    names(out_dat_met_mid@elementMetadata)<-"score"
    genome(out_dat_met_mid)<-"hg38"
    seqlengths(out_dat_met_mid)<-as.data.frame(hg38_seq_info)[hg38_seq_info@seqnames %in% out_dat_met_mid@seqnames,]$seqlengths #filter by seqlengths

    out_dat_met_low <- out_dat %>% 
                        as.data.frame() %>% 
                        filter(mcols(out_dat)$score <= 50 & mcols(out_dat)$score > 20) %>% 
                        mutate(score=score-mean_score) %>% 
                        GRanges() %>% unique()
    names(out_dat_met_low@elementMetadata)<-"score"
    genome(out_dat_met_low)<-"hg38"
    seqlengths(out_dat_met_low)<-as.data.frame(hg38_seq_info)[hg38_seq_info@seqnames %in% out_dat_met_low@seqnames,]$seqlengths #filter by seqlengths

    out_dat_met_hypomet <- out_dat %>% 
                            as.data.frame() %>% filter(mcols(out_dat)$score <= 20) %>%
                            mutate(score=score-mean_score) %>% 
                            GRanges() %>% unique()
    names(out_dat_met_hypomet@elementMetadata)<-"score"
    genome(out_dat_met_hypomet)<-"hg38"
    seqlengths(out_dat_met_hypomet)<-as.data.frame(hg38_seq_info)[hg38_seq_info@seqnames %in% out_dat_met_hypomet@seqnames,]$seqlengths #filter by seqlengths

    print(paste("Saving bedgraphs for...",i))
    rtracklayer::export(out_dat_hypermet,con=paste0(outdir,"/",i,".hypermet.bw"))
    rtracklayer::export(out_dat_met_mid,con=paste0(outdir,"/",i,".midmet.bw"))
    rtracklayer::export(out_dat_met_low,con=paste0(outdir,"/",i,".lowmet.bw"))
    rtracklayer::export(out_dat_met_hypomet,con=paste0(outdir,"/",i,".hypomet.bw"))
}

generate_bigwig<-function(obj=dat_cellline,
                      suffix="cellline",
                      groupBy="coarse_cluster_leidenclus",
                      outdir="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/01_amethyst_initial_object/"){

  bigwig_output_dir=paste0(outdir,"/bigwig_output_",suffix)
  system(paste("mkdir -p ", bigwig_output_dir))

  obj@h5paths$path<-obj@h5paths$paths
  obj@h5paths$barcode<-row.names(obj@h5paths)

  #update to new smoothed windows
  celltype500bpwindows <- calcSmoothedWindows(obj, 
                                          type = "CG", 
                                          threads = 200,
                                          step = 500, 
                                          smooth = 3,
                                          genome = "hg38",
                                          index = "chr_cg",
                                          groupBy = groupBy,
                                          returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                          returnPctMatrix = TRUE)
  saveRDS(celltype500bpwindows,file=paste0(bigwig_output_dir,"/","01.0.VMR_umap.",suffix,".coarse_cluster.500bp_windows.rds"))
  obj@genomeMatrices[[paste0("cg_",suffix,"_cells_perc")]] <- celltype500bpwindows[["pct_matrix"]]
  #output tracks as bigwig
  groups<-colnames(celltype500bpwindows[["pct_matrix"]])[4:ncol(celltype500bpwindows[["pct_matrix"]])]
  groups<-groups[groups!="NA"]
  lapply(groups,function(i) {
    write_binned_bigwigs(celltype_tracks=celltype500bpwindows[["pct_matrix"]],
                          outdir=bigwig_output_dir,
                          i=i)})
  return(obj)
}

```

Cell line
```R

#4. Perform initial clustering via iterative PCA into UMAP
dat_cellline<-vmr_cluster(obj=dat_cellline,suffix="cellline.predoubletfilter") 

#5. Identify doublets by initial cluster, via artificial doublet generation and random forest model
dat_cellline<-generate_doublet_score(obj=dat_cellline,suffix="cellline.predoubletfilter",reduction_name="irlba") 

#6. Remove identified doublets.
#filter by leiden clus (vmr filter)
dat_cellline<-subsetObject(dat_cellline,cells=row.names(dat_cellline@metadata[!is.na(dat_cellline@metadata$coarse_cluster_leidenclus),]))
#filter by doubletscore
singlets <- rownames(dat_cellline@metadata[dat_cellline@metadata$doublet_score < 0.5, ])
dat_cellline <- subsetObject(dat_cellline, cells = singlets)

#7. Recluster to final coarse cluster for cell typing.
dat_cellline<-vmr_cluster(obj=dat_cellline,suffix="cellline")

#8. Identify cell types by coarse marker genes via IGV bigwig output.
dat_cellline<-generate_bigwig(obj=dat_cellline,suffix="cellline")

#9. Save final objects.
saveRDS(dat_cellline,file="01.0.cellline.filt.amethyst.rds")

```




Patient
````R
#4. Perform initial clustering via iterative PCA into UMAP

dat_patient<-vmr_cluster(obj=dat_patient,suffix="patient.predoubletfilter") 

#5. Identify doublets by initial cluster, via artificial doublet generation and random forest model
dat_patient<-generate_doublet_score(obj=dat_patient,suffix="patient.predoubletfilter",reduction_name="irlba") #mamma mia this uses a lotta memory

#6. Remove identified doublets.
#filter by leiden clus (vmr filter)
dat_patient<-subsetObject(dat_patient,cells=row.names(dat_patient@metadata[!is.na(dat_patient@metadata$coarse_cluster_leidenclus),]))
#filter by doubletscore
singlets <- rownames(dat_patient@metadata[dat_patient@metadata$doublet_score < 0.5, ])
dat_patient <- subsetObject(dat_patient, cells = singlets)

#7. Recluster to final coarse cluster for cell typing.
dat_patient<-vmr_cluster(obj=dat_patient,suffix="patient",npc=20,pheno_k=180)

#8. Identify cell types by coarse marker genes via IGV bigwig output.
dat_patient<-generate_bigwig(obj=dat_patient,suffix="patient",groupBy="coarse_cluster_phenograph")

#assigning cell lineage and initial coarse cell types per cluster

clus_lineage<-c(
'8'='immune',
'1'='immune',
'12'='immune',
'5'='immune',
'3'='immune',
'27'='stromal',
'2'='stromal',
'15'='stromal',
'4'='stromal',
'11'='stromal',
'22'='stromal',
'26'='stromal',
'16'='epithelial',
'18'='epithelial',
'21'='epithelial',
'17'='epithelial',
'9'='epithelial',
'19'='epithelial',
'20'='epithelial',
'24'='epithelial',
'10'='epithelial',
'6'='epithelial',
'14'='epithelial',
'35'='epithelial',
'34'='epithelial',
'33'='epithelial',
'32'='epithelial',
'31'='epithelial',
'30'='epithelial',
'29'='epithelial',
'28'='epithelial',
'25'='epithelial',
'23'='epithelial',
'13'='epithelial',
'7'='epithelial')

clus_celltype<-c(
'8'='tcell',
'1'='tcell',
'12'='tcell',
'5'='bcell',
'3'='myeloid',
'27'='endothelial',
'2'='endothelial',
'15'='fibroblast',
'4'='fibroblast',
'11'='peri_VSMC',
'22'='peri_VSMC',
'26'='adipocyte',
'16'='basal',
'18'='basal',
'21'='basal',
'17'='lumhr', #changed after finding cancer cells in cluster
'9'='lumsec',
'19'='lumsec',
'20'='lumsec',
'24'='lumsec',
'10'='lumhr',
'6'='lumhr',
'14'='lumhr',
'35'='lumhr',
'34'='lumhr',
'33'='lumhr',
'32'='lumhr',
'31'='lumhr',
'30'='lumhr',
'29'='lumhr',
'28'='lumhr',
'25'='lumhr',
'23'='lumhr',
'13'='lumhr',
'7'='lumhr')

dat_patient@metadata$celltype_lineage<-unname(clus_lineage[dat_patient@metadata$coarse_cluster_phenograph])
dat_patient@metadata$coarse_celltype<-unname(clus_celltype[dat_patient@metadata$coarse_cluster_phenograph])
#9. Save final objects.
saveRDS(dat_patient,file="01.0.patient.filt.amethyst.rds")
```
# next up! call CNVs and define cancer cells.
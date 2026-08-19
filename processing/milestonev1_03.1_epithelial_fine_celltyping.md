# Running final cell typing
1. check major lineages against reclustered data
2. subset out major cell lineages
3. perform reclustering on vmrs for just those cells
4. use multiple resolutions of clustering to combine ambiguous clusters
5. output bigwigs in IGV and use RNA markers for pairwise discrimination of cell states
6. use lineage to sublineage to cell type markers to define cells
# Prepare RNA marker genes per cell type


# Read in methylation data and additional libraries
From processing/milestonev1_01_amethyst_fine_celltyping.md
```R

set.seed(111)
options(future.globals.maxSize= 500000*1024^2) #80gb limit for parallelizing
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
library(GeneOverlap)
library(RPhenograph)
library(matrixStats)
library(sparseMatrixStats)
library(irlba)
library(rtracklayer)
library(cowplot)
library(patchwork)

project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"

#read in object from directory
task_cpus=300
processing_folder="03_fine_celltyping"
wd=paste(sep="/",project_data_directory,processing_folder)
system(paste0("mkdir -p ",wd))
setwd(wd)
```


From processing/milestonev1_00.2_seurat_scrna_copykat.md seurat object.
```R

#compare dmr sites with rna marker genes overlap
library(Seurat)
library(GenomicRanges)

rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")
#just run on all cell types
Idents(rna)<-rna$coarse_celltype
table(Idents(rna))
rna<-subset(rna,coarse_celltype %in% c("lumsec","lumhr","cancer","basal"))

rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)
rna <- ScaleData(rna)
rna <- JoinLayers(rna)

res=0.2
dims=1:50 #originally 1:50
#umap here is just for visualization, not used for redefining cell types
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 10000)
rna <- RunPCA(rna, features = VariableFeatures(object = rna),npcs=max(dims))
plt<-ElbowPlot(rna, reduction = "pca",ndims=max(dims))
ggsave(plt,file="test.elbow.rna.pdf",width=30)

rna <- FindNeighbors(rna, dims = 1:50)
rna <- FindClusters(rna, resolution = 0.2)
rna <- RunUMAP(rna, dims = 1:50)

rna@meta.data$fine_cluster_UMAP_X<-rna@reductions$umap@cell.embeddings[,1]
rna@meta.data$fine_cluster_UMAP_Y<-rna@reductions$umap@cell.embeddings[,2]
rna@meta.data$fine_clusters<-rna@meta.data$seurat_clusters

plt<-DimPlot(rna,group.by=c("fine_clusters","sample","coarse_celltype","rna_ploidy"),raster=T,label=TRUE)
ggsave(plt,file="test.rna.pdf",width=20,height=20)

#assigning cancer to clusters those with majority aneuploid signatures
aneuploid_clust_table<-table(rna$rna_ploidy,rna$fine_clusters)
pct_table <- sweep(aneuploid_clust_table, 2, colSums(aneuploid_clust_table), FUN = "/") %>% 
  t() %>% 
  as.data.frame() %>% 
  tidyr::pivot_wider(id_cols=Var2,names_from=Var1,values_from=Freq) %>% 
  tibble::column_to_rownames(var="Var2") %>% 
  mutate(Var2=NULL)

#majority aneuploid clusters defined as cancer
aneuploid_clusters<-colnames(pct_table)[pct_table["aneuploid",]>0.5]
rna@meta.data[rna@meta.data$fine_clusters %in% aneuploid_clusters,]$coarse_celltype<-"cancer"


#also assigning to those with majority single sample signatures, consistent with cancer clones
#this just grabs cluster 10 and 19
sample_clust_table<-table(rna$sample,rna$fine_clusters)
pct_table <- sweep(sample_clust_table, 2, colSums(aneuploid_clust_table), FUN = "/") %>% 
  t() %>% 
  as.data.frame() %>% 
  tidyr::pivot_wider(id_cols=Var2,names_from=Var1,values_from=Freq) %>% 
  tibble::column_to_rownames(var="Var2") %>% 
  mutate(Var2=NULL)

#majority aneuploid clusters defined as cancer
sample_based_cancer<-which(pct_table>0.9,arr.ind=T)
sample_cancer_col<-names(pct_table)[sample_based_cancer[,2]]
rna@meta.data[rna@meta.data$fine_clusters %in% sample_cancer_col,]$coarse_celltype<-"cancer"


saveRDS(rna,file="03_rna.fine_celltyping.epithelial.rds")
Idents(rna)<-rna$fine_celltype
rna_fine_markers<-FindAllMarkers(rna,assay="RNA",only.pos=TRUE)
rna_markers<-rna$coarse_celltype
rna_fine_markers<-FindAllMarkers(rna,assay="RNA",only.pos=TRUE)
rna<-readRDS(file="03_rna.fine_celltyping.epithelial.rds")


```



Functions
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

write_binned_bigwigs<-function(celltype_tracks=celltype_tracks,
                                outdir=celltype_outdir,
                                i){
    #split bw into 4 files per cell type by methylation/average methylation
    out_dat<-celltype_tracks %>% distinct(chr,start,end, .keep_all=TRUE) %>% select(chr,start,end,i) 

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
                        GRanges() 
    names(out_dat_hypermet@elementMetadata)<-"score"
    genome(out_dat_hypermet)<-"hg38"
    seqlengths(out_dat_hypermet)<-as.data.frame(hg38_seq_info)[hg38_seq_info@seqnames %in% out_dat_hypermet@seqnames,]$seqlengths #filter by seqlengths

    out_dat_met_mid <- out_dat %>% 
                      as.data.frame() %>% 
                      filter(mcols(out_dat)$score <= mean_score & mcols(out_dat)$score > 50) %>% 
                      mutate(score=score-mean_score) %>% 
                      GRanges() 
    names(out_dat_met_mid@elementMetadata)<-"score"
    genome(out_dat_met_mid)<-"hg38"
    seqlengths(out_dat_met_mid)<-as.data.frame(hg38_seq_info)[hg38_seq_info@seqnames %in% out_dat_met_mid@seqnames,]$seqlengths #filter by seqlengths

    out_dat_met_low <- out_dat %>% 
                        as.data.frame() %>% 
                        filter(mcols(out_dat)$score <= 50 & mcols(out_dat)$score > 20) %>% 
                        mutate(score=score-mean_score) %>% 
                        GRanges() 
    names(out_dat_met_low@elementMetadata)<-"score"
    genome(out_dat_met_low)<-"hg38"
    seqlengths(out_dat_met_low)<-as.data.frame(hg38_seq_info)[hg38_seq_info@seqnames %in% out_dat_met_low@seqnames,]$seqlengths #filter by seqlengths

    out_dat_met_hypomet <- out_dat %>% 
                            as.data.frame() %>% filter(mcols(out_dat)$score <= 20) %>%
                            mutate(score=score-mean_score) %>% 
                            GRanges() 
    names(out_dat_met_hypomet@elementMetadata)<-"score"
    genome(out_dat_met_hypomet)<-"hg38"
    seqlengths(out_dat_met_hypomet)<-as.data.frame(hg38_seq_info)[hg38_seq_info@seqnames %in% out_dat_met_hypomet@seqnames,]$seqlengths #filter by seqlengths

    print(paste("Saving bedgraphs for...",i))
    export.bw(out_dat_hypermet,con=paste0(outdir,"/",i,".hypermet.bw"))
    export.bw(out_dat_met_mid,con=paste0(outdir,"/",i,".midmet.bw"))
    export.bw(out_dat_met_low,con=paste0(outdir,"/",i,".lowmet.bw"))
    export.bw(out_dat_met_hypomet,con=paste0(outdir,"/",i,".hypomet.bw"))
  #make a ucsc/igv trackhub
  col_out=c("hypermet"="#333333","midmet"="#444444","lowmet"="#bcbcbc","hypomet"="#FF13F0")
  
  #write out multiwig trackhub
  multiwig<-c(
      paste(c("track", paste0(i)),collapse=" "),
      paste(c("container", "multiWig"),collapse=" "),
      paste(c("shortLabel",  paste0(i)),collapse=" "),
      paste(c("longLabel",  paste0(i)),collapse=" "),
      paste(c("visibility", "full"),collapse=" "),
      paste(c("aggregate", "transparentOverlay"),collapse=" "),
      paste(c("showSubtrackColorOnUi", "on"),collapse=" "))
  track_lines<-lapply(names(col_out),function(met){
    return(c(
      paste(c("track", paste0(i,".",met,".bw")),collapse=" "),
      paste(c("shortLabel", paste0(i,".",met)),collapse=" "),
      paste(c("longLabel",  paste0(i,".",met)),collapse=" "),
      paste(c("parent",  paste0(i)),collapse=" "),
      paste(c("type", "bigWig"),collapse=" "),
      paste(c("group", i),collapse=" "),
      paste(c("autoscale", "group"),collapse=" "),
      paste(c("graphTypeDefault", "bar"),collapse=" "),
      paste(c("color", paste(col2rgb(col_out[met])[,1],collapse=","),collapse=" "))))
  })
  return(c(list(multiwig),track_lines))
  
}

generate_bigwig<-function(obj=immune,
                      suffix="immune",
                      groupBy="fine_cluster_phenograph",
                      outdir="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/03_fine_celltyping/"){
  bigwig_output_dir=paste0(outdir,"/bigwig_output_",suffix)
  system(paste("mkdir -p ", bigwig_output_dir))
  system(paste("mkdir -p ", paste0(bigwig_output_dir,"/","hg38")))

  obj@h5paths$path<-obj@h5paths$paths
  obj@h5paths$barcode<-row.names(obj@h5paths)

  #update to new smoothed windows
  celltype500bpwindows <- calcSmoothedWindows(obj, 
                                          type = "CG", 
                                          threads = 300,
                                          step = 500, 
                                          smooth = 1,
                                          genome = "hg38",
                                          index = "chr_cg",
                                          groupBy = groupBy,
                                          returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                          returnPctMatrix = TRUE)
  saveRDS(celltype500bpwindows,file=paste0(bigwig_output_dir,"/","03.1.VMR_umap.",suffix,".fine_cluster.500bp_windows.rds"))
  obj@genomeMatrices[[paste0("cg_",suffix,"_cells_perc")]] <- celltype500bpwindows[["pct_matrix"]]
  #output tracks as bigwig
  groups<-colnames(celltype500bpwindows[["pct_matrix"]])[4:ncol(celltype500bpwindows[["pct_matrix"]])]
  groups<-groups[groups!="NA"]
  
  #make track hub.txt
  writeLines(c(
    paste(c("hub", groupBy),collapse=" "),
    paste(c("shortLabel", groupBy),collapse=" "),
    paste(c("longLabel", groupBy),collapse=" "),
    paste(c("genomesFile", "genomesFile.txt"),collapse=" "),
    paste(c("email ryan"),collapse=" ")),
  paste0(bigwig_output_dir,"/","hub.txt"))

  #make track genomesFile.txt
  writeLines(c(
    paste(c("genome", "hg38"),collapse=" "),
    paste(c("trackDb", "hg38/trackDb.txt"),collapse=" ")),
  paste0(bigwig_output_dir,"/","genomesFile.txt"))
  
  #write out bigwigs and trackhub data
  trackhub_dat<-lapply(groups,function(i) {
    write_binned_bigwigs(celltype_tracks=celltype500bpwindows[["pct_matrix"]],
                          outdir=paste0(bigwig_output_dir,"/","hg38"),
                          i=i)})
  file_conn <- file(paste0(bigwig_output_dir,"/hg38/","trackDb.txt"), open = "a")
  for(i in trackhub_dat){writeLines(unlist(i), con = file_conn)}
  close(file_conn)

  return(obj)
}

cell_subtyping_clustering<-function(suffix="stromal",
                                    outdir="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/03_fine_celltyping/",
                                    init_npc=50,
                                    pc_var_explained=50,
                                    reduction_name="irlba_immune",
                                    leiden_cluster_resolution=5e-5,
                                    pheno_k=200,
                                    min_dist=1e-6,
                                    var_features_count=20000,
                                    n_neighbors=5){
  
  plot_dir=paste0(outdir,"/","plot_",suffix)
  system(paste0("mkdir -p ",plot_dir))

  print(paste0("Saving RNA umap to ",plot_dir))

  for(i in c("coarse_celltype","fine_celltype","Group")){
      plt_dim<-DimPlot(rna,group.by=i,label=TRUE,raster=FALSE)
      ggsave(plt_dim,file=paste0(plot_dir,"/","tenx_dcis.pf.final.",suffix,".",i,".umap.pdf"),width=10)
      }

  print(paste0("Saving RNA Markers to ",plot_dir))
  Idents(rna)<-rna$coarse_celltype
  rna_markers<-FindAllMarkers(rna,assay="RNA",only.pos=TRUE)
  saveRDS(rna_markers,file=paste0(plot_dir,"/tenx_dcis.rna_markers.",suffix,".rds"))

  print(paste0("Saving plots to ",plot_dir))

  print("Calculating PCA...")

  #filter features to top variance
  feature_variance<-apply(t(obj@genomeMatrices[["vmr_matrix_cg_residuals"]]), 2, var, na.rm=T)

  top_var_features<-feature_variance %>% as.data.frame() %>% setNames("var") %>% filter(is.finite(var))%>% slice_max(var, n=var_features_count)

  var_plt<-ggplot()+geom_violin(aes(x=1,y=feature_variance))+theme_minimal()+geom_hline(yintercept=min(top_var_features$var),color="red")+ggtitle("VMR Feature Variance")
  ggsave(var_plt,file=paste0(plot_dir,"/vmr.feature_variance.",suffix,".pdf"))

  pca<- t(obj@genomeMatrices[["vmr_matrix_cg_residuals"]][row.names(top_var_features),]) %>% 
      scale(center = T, scale = F) %>% 
      prcomp_iterative(n = init_npc)  # this is number of initial run principle components (init_npc in args) to get a sense of variance explained

  # Compute variance explained
  variance_explained <- pca$sdev^2
  total_variance <- sum(variance_explained)
  percent_variance <- (variance_explained / total_variance) * 100
  
  #cut off pca at 90% variance explained
  pca_to_use<-min(which(cumsum(percent_variance)>pc_var_explained))

  pca_plt<-ggplot()+geom_point(aes(x=1:length(percent_variance),y=percent_variance)) + geom_vline(xintercept=pca_to_use,color="red")+theme_minimal() + ggtitle("PCA elbow with 75% variance explained")
  ggsave(pca_plt,file=paste0(plot_dir,"/pca_variance_explained.",suffix,".pdf"))
  print(paste("PCs to explain",as.character(pc_var_explained),"%", "of captured variance:",as.character(pca_to_use)))

  pca_dims <- pca$x[,1:pca_to_use] %>% 
    magrittr::set_rownames(colnames(obj@genomeMatrices[["vmr_matrix_cg_residuals"]]))


  print("Running UMAP...")
  umap_obj <- uwot::umap(X=as.matrix(pca_dims), 
                          metric="cosine", 
                          min_dist=min_dist, 
                          n_neighbors=n_neighbors, 
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
    igraph::cluster_leiden(resolution_parameter = leiden_cluster_resolution) 
  table(clust_obj$membership)

  # put the clustering results into a data frame (tibble) for plotting
  clust_tbl <- tibble(
    leiden_cluster = as.character(clust_obj$membership),
    cell= clust_obj$names) %>% 
    full_join(umap_tbl, by="cell")

    print("Adding PCA to reduction slot...")
    obj@reductions[[reduction_name]]<-pca_dims

    print("Adding clusters to metadata...")
    fine_cluster_leidenclus<-setNames(nm=clust_tbl$cell,clust_tbl$leiden_cluster)
    fine_cluster_UMAP1<-setNames(nm=clust_tbl$cell,clust_tbl$UMAP1)
    fine_cluster_UMAP2<-setNames(nm=clust_tbl$cell,clust_tbl$UMAP2)

    obj@metadata$fine_cluster_UMAP_X<-NA
    obj@metadata$fine_cluster_UMAP_Y<-NA
    obj@metadata$fine_cluster_leidenclus<-NA

    obj@metadata$fine_cluster_UMAP_X<-fine_cluster_UMAP1[row.names(obj@metadata)]
    obj@metadata$fine_cluster_UMAP_Y<-fine_cluster_UMAP2[row.names(obj@metadata)]
    obj@metadata$fine_cluster_leidenclus<-paste0(suffix,"_",fine_cluster_leidenclus[row.names(obj@metadata)])

    #cluster on umap as well
    umap_clus<-obj@metadata %>% 
      select(c("fine_cluster_UMAP_X","fine_cluster_UMAP_Y")) %>% 
      magrittr::set_rownames(row.names(obj@metadata))

    umap_clus<-Rphenograph::Rphenograph(umap_clus,k=pheno_k)

    obj@metadata$fine_cluster_phenograph<-paste0(suffix,"_",as.character(unlist(as.list(igraph::membership(umap_clus[[2]])))))
    print("Plotting...")
    

    print("Plotting...")
    lapply(c("fine_cluster_leidenclus","fine_cluster_phenograph","Sample","mcg_pct","unique_reads","Group"), 
      function(groupby){
        if(groupby %in% colnames(obj@metadata)){
          plt_clus<-obj@metadata %>% 
            ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = obj@metadata[[groupby]])) +
              geom_point() +
              coord_fixed()
            ggsave(plt_clus,file=paste0(plot_dir,"/","03.1.VMR_umap.",suffix,".",as.character(pca_to_use),".",as.character(min_dist),".","clus.",as.character(groupby),".pdf"),width=20,height=20)
          }
      })
    return(obj)
}

plot_umap_panels_met<-function(x=epithelial,outname="epithelial"){
  #plot celltype
  plt_celltype<-x@metadata %>% 
      ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = celltype)) +
      geom_point() +
      coord_fixed() + 
      theme_minimal() + xlab("UMAP X") + ylab("UMAP Y")+ ggtitle("cell type")+
      theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))+
      scale_color_manual(values=celltype_col)
  plt_celltype_legend <- get_legend(plt_celltype)
  plt_celltype_umap<-plt_celltype+ theme(legend.position='none')

  #plot group
  plt_group<-x@metadata %>% 
      ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = Group)) +
      geom_point() +
      coord_fixed() + 
      theme_minimal() + xlab("UMAP X") + ylab("UMAP Y")+ ggtitle("read counts")+
      theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))+
      scale_color_manual(values=group_col)
  plt_group_legend <- get_legend(plt_group)
  plt_group_umap<-plt_group+ theme(legend.position='none')


  #plot Sample
  plt_sample<-x@metadata %>% 
      ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = Sample)) +
      geom_point() +
      coord_fixed() + 
      theme_minimal() + xlab("UMAP X") + ylab("UMAP Y")+ ggtitle("percent methylation")+
      theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black")) +
      scale_color_manual(values=sample_col)
  plt_sample_legend <- get_legend(plt_sample)
  plt_sample_umap<-plt_sample+ theme(legend.position='none')

  #plot %met
  plt_cgperc<-x@metadata %>% 
      ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = mcg_pct)) +
      geom_point() +
      coord_fixed() + 
      theme_minimal() + xlab("UMAP X") + ylab("UMAP Y")+ ggtitle("diagnostic group")+
      theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black")) +
      scale_color_gradient2(low="#ff70ff",mid="#CCCCCC",high="#000000",midpoint=median(x@metadata$mcg_pct))
  plt_cgperc_legend <- get_legend(plt_cgperc)
  plt_cgperc_umap<-plt_cgperc+ theme(legend.position='none')


  #plot %met
  plt_read<-x@metadata %>% 
      ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = log10(unique_reads))) +
      geom_point() +
      coord_fixed() + 
      theme_minimal() + xlab("UMAP X") + ylab("UMAP Y")+ ggtitle("sample")+
      theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black")) +
      scale_color_gradient2(low="#cccccc",mid="#CCCCFF",high="#000066",midpoint=median(log10(x@metadata$unique_reads)))
  plt_read_legend <- get_legend(plt_read)
  plt_read_umap<-plt_read+ theme(legend.position='none')

  #plot ploidy
  plt_ploidy<-x@metadata %>% 
      ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = cnv_ploidy_500kb)) +
      geom_point() +
      coord_fixed() + 
      theme_minimal() + xlab("UMAP X") + ylab("UMAP Y")+ ggtitle("cnv based ploidy")+
      theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black")) +
      scale_color_manual(values=ploidy_col)
  plt_ploidy_legend <- get_legend(plt_ploidy)
  plt_ploidy_umap<-plt_ploidy+ theme(legend.position='none')

  layout<-"
  AABDF
  AACEG"

  umap_plts<-plt_celltype_umap+plt_read_umap+plt_cgperc_umap+plt_group_umap+plt_sample_umap+plt_ploidy_umap+ggplot() + plot_layout(design = layout)+theme(plot.margin = margin(l = 3, r = 3))
  legends_plts<-plot_grid(plt_celltype_legend, plt_read_legend, plt_cgperc_legend, plt_group_legend,plt_sample_legend,plt_ploidy_legend, nrow = 1, align = "h")

  ggsave(umap_plts,
        file=paste0("03.1.",outname,".celltype.met.umap.pdf"),
        width=60,
        height=40,
        limitsize=FALSE)

  ggsave(legends_plts,
          file=paste0("03.1.",outname,".celltype.umap.met.legends.pdf"),
          width=60,
          height=10,
          limitsize=FALSE)
}

plot_umap_panels_rna<-function(x=rna,outname="epithelial"){
  #plot celltype
  plt_celltype<-x@meta.data %>% 
      ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = coarse_celltype)) +
      geom_point() +
      coord_fixed() + 
      theme_minimal() + xlab("UMAP X") + ylab("UMAP Y")+ ggtitle("cell type")+
      theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))+
      scale_color_manual(values=celltype_col)
  plt_celltype_legend <- get_legend(plt_celltype)
  plt_celltype_umap<-plt_celltype+ theme(legend.position='none')

  #plot group
  plt_group<-x@meta.data %>% 
      ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = Group)) +
      geom_point() +
      coord_fixed() + 
      theme_minimal() + xlab("UMAP X") + ylab("UMAP Y")+ ggtitle("read counts")+
      theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))+
      scale_color_manual(values=group_col)
  plt_group_legend <- get_legend(plt_group)
  plt_group_umap<-plt_group+ theme(legend.position='none')


  #plot Sample
  plt_sample<-x@meta.data %>% 
      ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = sample)) +
      geom_point() +
      coord_fixed() + 
      theme_minimal() + xlab("UMAP X") + ylab("UMAP Y")+ ggtitle("percent methylation")+
      theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black")) +
      scale_color_manual(values=sample_col)
  plt_sample_legend <- get_legend(plt_sample)
  plt_sample_umap<-plt_sample+ theme(legend.position='none')

  #plot ploidy
  plt_ploidy<-x@meta.data %>% 
      ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = rna_ploidy)) +
      geom_point() +
      coord_fixed() + 
      theme_minimal() + xlab("UMAP X") + ylab("UMAP Y")+ ggtitle("cnv based ploidy")+
      theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black")) +
      scale_color_manual(values=ploidy_col)
  plt_ploidy_legend <- get_legend(plt_ploidy)
  plt_ploidy_umap<-plt_ploidy+ theme(legend.position='none')

  layout<-"
  AABD
  AACE"

  umap_plts<-plt_celltype_umap+plt_group_umap+plt_sample_umap+plt_ploidy_umap+ggplot() + plot_layout(design = layout)+theme(plot.margin = margin(l = 3, r = 3))
  legends_plts<-plot_grid(plt_celltype_legend, plt_group_legend,plt_sample_legend,plt_ploidy_legend, nrow = 1, align = "h")

  ggsave(umap_plts,
        file=paste0("03.1.",outname,".celltype.rna.umap.pdf"),
        width=60,
        height=40,
        limitsize=FALSE)

  ggsave(legends_plts,
          file=paste0("03.1.",outname,".celltype.rna.umap.legends.pdf"),
          width=60,
          height=10,
          limitsize=FALSE)
}


```


#subcluster on epithelial cells

```R
obj<-readRDS(file=paste(project_data_directory,"02_copykit_cnv_calling","02_scaledcis.cnv_clones.amethyst.rds",sep="/"))

#filter to immune cells by fine_celltype defined lineages
obj<-subsetObject(obj,cells=row.names(obj@metadata[obj@metadata$celltype_lineage %in% c("epithelial"),]))
dim(obj@metadata)

#run once with 50 npcs
#cut back based on pca variance explained in second run
epithelial<-cell_subtyping_clustering(suffix="epithelial",
                                    outdir="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/03_fine_celltyping/",
                                    init_npc=50,
                                    pc_var_explained=90,
                                    reduction_name="irlba_epithelial",
                                    leiden_cluster_resolution=1e-5,
                                    pheno_k=200,
                                    min_dist=1e-6,
                                    var_features_count=25000,
                                    n_neighbors=10)

epithelial<-generate_bigwig(obj=epithelial,
                        suffix="epithelial",
                        groupBy="fine_cluster_phenograph",
                        outdir=getwd())
#assign cell types via IGV tracks

celltype_assignment=c(
  'epithelial_11'='lumsec',
  'epithelial_8'='lumsec',
  'epithelial_9'='lumsec',
  'epithelial_6'='lumsec',
  'epithelial_18'='lumsec',
  'epithelial_14'='lumsec',

  'epithelial_13'='basal',
  'epithelial_2'='basal',
  'epithelial_7'='basal',
  
    'epithelial_3'='lumhr',
    'epithelial_4'='lumhr',
    'epithelial_17'='lumhr',
    'epithelial_10'='lumhr',
    'epithelial_5'='lumhr',
    'epithelial_15'='lumhr',
    'epithelial_1'='lumhr',

    'epithelial_12'='lumhr',
    'epithelial_16'='lumhr',
    'epithelial_19'='lumhr',
    'epithelial_20'='lumhr',
    'epithelial_21'='lumhr',
    'epithelial_22'='lumhr',
    'epithelial_23'='lumhr',
    'epithelial_24'='lumhr',
    'epithelial_25'='lumhr',
    'epithelial_26'='lumhr',
    'epithelial_27'='lumhr',
    'epithelial_28'='lumhr',
    'epithelial_29'='lumhr',
    'epithelial_30'='lumhr')



#get clusters with cancer clones assigned, assign the whole cluster as cancer
epithelial@metadata$celltype<-"NA"
epithelial@metadata$celltype<-celltype_assignment[epithelial@metadata$fine_cluster_phenograph]

 plt_celltype<-epithelial@metadata %>% 
    ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = celltype)) +
    geom_point() +
    coord_fixed()
  ggsave(plt_celltype,file="03.1.epithelial.celltype.umap.pdf",width=20,height=20)


 plt_celltype<-epithelial@metadata %>% 
    ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = fine_cluster_phenograph)) +
    geom_point() +
    coord_fixed()
  ggsave(plt_celltype,file="03.1.epithelial.cluster.umap.pdf",width=20,height=20)

saveRDS(epithelial,file="03_scaledcis.fine_celltyping.epithelial.rds")

#now assign cancer to cells in the same cluster and CNV called cancer cells
hm<-table(epithelial@metadata$cnv_ploidy_500kb,epithelial@metadata$fine_cluster_phenograph)
hm_count<-as.data.frame(hm) #convert it back into a data frame, now with the counts
hm_scale<-as.data.frame(scale(hm,center=F))
colnames(hm_count)<-c("group","cluster","count")
hm_count$scaled <- hm_scale$Freq #combine count and scaled data to plot both
plt<-ggplot(hm_count, aes(x=group, y=cluster, fill=scaled,label=count)) + 
  geom_tile() + 
  geom_text(color="white") + 
  theme(axis.text.x = element_text(angle = 90))
ggsave(plt,file="03_scaledcis.fine_celltyping.epithelial.cancer_cluster.heatmap.pdf")

#assign celltype as cancer to those with >80% aneuploid cells within cluster and those from single sample origin (which tracks with cancer clustering)
#this is to pick up any cells with NA calls for aneuploidy
cancer_assignment<-c("epithelial_30","epithelial_29","epithelial_27","epithelial_26","epithelial_25","epithelial_24","epithelial_23","epithelial_22","epithelial_21","epithelial_19","epithelial_16","epithelial_12","epithelial_1","epithelial_28","epithelial_20")


#get clusters with cancer clones assigned, assign the whole cluster as cancer
epithelial@metadata$celltype<-"NA"
epithelial@metadata$celltype<-celltype_assignment[epithelial@metadata$fine_cluster_phenograph]
epithelial@metadata[epithelial@metadata$fine_cluster_phenograph %in% cancer_assignment,]$celltype<-"cancer"
#and assign any aneuploidy cells as cancer
epithelial@metadata[epithelial@metadata$cnv_ploidy_500kb=="aneuploid",]$celltype<-"cancer"

saveRDS(epithelial,file="03_scaledcis.fine_celltyping.epithelial.rds")
epithelial<-readRDS(file="03_scaledcis.fine_celltyping.epithelial.rds")
```

Assign cancer by RNA clusters as well
```R

```

```R
###Generate UMAP plots for supp figure###
#colors taken from manuscript_colors.md in ref folder
#source(manuscript_colors.r)

plot_umap_panels_met(x=epithelial,outname="epithelial")
plot_umap_panels_rna(x=rna,outname="epithelial")

table(rna$rna_clonename,rna$rna_ploidy,useNA="ifany")
rna@meta.data %>% filter(is.na(rna_clonename)) %>% summarize(mean(nCount_RNA))
```
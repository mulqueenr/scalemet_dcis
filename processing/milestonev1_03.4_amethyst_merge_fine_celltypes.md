# Write out fine cell types into shared amethyst object.

```R

set.seed(111)
options(future.globals.maxSize= 500000*1024^2) #80gb limit for parallelizing
library(amethyst)
library(data.table)
library(rtracklayer)

#read in object from directory
task_cpus=300
processing_folder="03_fine_celltyping"
wd=paste(sep="/",project_data_directory,processing_folder)
setwd(wd)
obj<-readRDS(file=paste(project_data_directory,"02_copykit_cnv_calling","02_scaledcis.cnv_clones.amethyst.rds",sep="/"))

epithelial<-readRDS(file="03_scaledcis.fine_celltyping.epithelial.rds")
immune<-readRDS(file="03_scaledcis.fine_celltyping.immune.rds")
stromal<-readRDS(file="03_scaledcis.fine_celltyping.stromal.rds")

obj@metadata$fine_cluster_leiden<-NA
obj@metadata$fine_cluster_phenograph<-NA
obj@metadata$celltype<-NA
obj@metadata$fine_cluster_UMAP_X<-NA
obj@metadata$fine_cluster_UMAP_Y<-NA

obj@metadata[row.names(epithelial@metadata),]$fine_cluster_phenograph<-epithelial@metadata$fine_cluster_leidenclus
obj@metadata[row.names(immune@metadata),]$fine_cluster_phenograph<-immune@metadata$fine_cluster_leidenclus
obj@metadata[row.names(stromal@metadata),]$fine_cluster_phenograph<-stromal@metadata$fine_cluster_leidenclus

obj@metadata[row.names(epithelial@metadata),]$fine_cluster_phenograph<-epithelial@metadata$fine_cluster_phenograph
obj@metadata[row.names(immune@metadata),]$fine_cluster_phenograph<-immune@metadata$fine_cluster_phenograph
obj@metadata[row.names(stromal@metadata),]$fine_cluster_phenograph<-stromal@metadata$fine_cluster_phenograph

obj@metadata[row.names(epithelial@metadata),]$celltype<-epithelial@metadata$celltype
obj@metadata[row.names(immune@metadata),]$celltype<-immune@metadata$celltype
obj@metadata[row.names(stromal@metadata),]$celltype<-stromal@metadata$celltype

obj@metadata[row.names(epithelial@metadata),]$fine_cluster_UMAP_X<-epithelial@metadata$fine_cluster_UMAP_X
obj@metadata[row.names(immune@metadata),]$fine_cluster_UMAP_X<-immune@metadata$fine_cluster_UMAP_X
obj@metadata[row.names(stromal@metadata),]$fine_cluster_UMAP_X<-stromal@metadata$fine_cluster_UMAP_X

obj@metadata[row.names(epithelial@metadata),]$fine_cluster_UMAP_Y<-epithelial@metadata$fine_cluster_UMAP_Y
obj@metadata[row.names(immune@metadata),]$fine_cluster_UMAP_Y<-immune@metadata$fine_cluster_UMAP_Y
obj@metadata[row.names(stromal@metadata),]$fine_cluster_UMAP_Y<-stromal@metadata$fine_cluster_UMAP_Y

saveRDS(obj,file="03_scaledcis.final_celltypes.amethyst.rds")

obj<-generate_bigwig(obj=obj,
                        suffix="final_celltype",
                        groupBy="final_celltype",
                        outdir=getwd())

#assign nonlumhr and noncancer cell types as diploid, excluding any lumhr or cancer cells
table(obj@metadata[!(obj@metadata$cnv_ploidy_500kb=="aneuploid") & !(obj@metadata$celltype %in% c("lumhr","cancer")),]$celltype)

obj@metadata[!(obj@metadata$cnv_ploidy_500kb=="aneuploid") & !(obj@metadata$celltype %in% c("lumhr","cancer")),]$cnv_ploidy_500kb<-"diploid"

table(obj@metadata[obj@metadata$cnv_ploidy_500kb=="diploid",]$cnv_clonename)

#can now assign diploidy to NA samples with too low of read count (~1k cells)
obj@metadata[obj@metadata$cnv_ploidy_500kb=="diploid",]$cnv_clonename<-paste0(
                      obj@metadata[obj@metadata$cnv_ploidy_500kb=="diploid",]$Sample,"_diploid")

obj<-generate_bigwig(obj=obj,
                        suffix="cnv_clonename",
                        groupBy="cnv_clonename",
                        outdir=getwd())

saveRDS(obj,file="03_scaledcis.final_celltypes.amethyst.rds")
```

Generate UMAP of all cells 
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

obj@metadata$cnv_ploidy
obj<-generate_bigwig(obj=obj,
                        suffix="final_celltype",
                        groupBy="celltype",
                        outdir=getwd())

#renamed macro_mono to myeloid, since other cells present as well
obj@metadata[obj@metadata$celltype=="macro_mono",]$celltype<-"myeloid"
saveRDS(obj,file="03_scaledcis.final_celltypes.amethyst.rds")

#set colors
celltype_col=c(
'tcell'='#2e3fa3',
'bcell'='#00adea',
'myeloid'='#00a487',
'pericyte'='#c1d552',
'fibroblast'='#7f1911',
'endothelial'='#f0b243',
'basal'='#7200cc',
'lumsec'='#af00af',
'lumhr'='#d8007c',
'cancer'='#80FF80')

#final cluster of cell types

suffix="final_celltype_final_plots"
outdir="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/03_fine_celltyping/"
npc=15
reduction_name="irlba_final_celltype"
min_dist=1e-5
n_neighbors=10

plot_dir=paste0(outdir,"/","plot_",suffix)
system(paste0("mkdir -p ",plot_dir))

print(paste0("Saving plots to ",plot_dir))
print("Calculating PCA...")

pca <- t(obj@genomeMatrices[["vmr_matrix_cg_residuals"]]) %>%
    scale(center = T, scale = F) %>%
    prcomp_iterative(n = npc)  # this is number of principle components (npc in args)

pca_dims <- pca$x %>% 
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

  obj@metadata$final_celltype_UMAP_X<-NA
  obj@metadata$final_celltype_UMAP_Y<-NA

  obj@metadata$final_celltype_UMAP_X<-fine_cluster_UMAP1[row.names(obj@metadata)]
  obj@metadata$final_celltype_UMAP_Y<-fine_cluster_UMAP2[row.names(obj@metadata)]

  print("Plotting...")
  lapply(c("fine_cluster_leidenclus","fine_cluster_phenograph","Sample","mcg_pct","unique_reads","Group","celltype"), 
    function(groupby){
      if(groupby %in% colnames(obj@metadata)){
        if(groupby=="celltype"){
        plt_clus<-obj@metadata %>% 
          ggplot(aes(x = final_celltype_UMAP_X, y = final_celltype_UMAP_Y, color = obj@metadata[[groupby]])) +
            geom_point() +
            coord_fixed() + scale_color_manual(values=celltype_col)
          ggsave(plt_clus,file=paste0(plot_dir,"/","03.1.VMR_umap.",suffix,".",as.character(npc),".",as.character(min_dist),".","clus.",as.character(groupby),".pdf"),width=20,height=20)
 
        } else {
        plt_clus<-obj@metadata %>% 
          ggplot(aes(x = final_celltype_UMAP_X, y = final_celltype_UMAP_Y, color = obj@metadata[[groupby]])) +
            geom_point() +
            coord_fixed()
          ggsave(plt_clus,file=paste0(plot_dir,"/","03.1.VMR_umap.",suffix,".",as.character(npc),".",as.character(min_dist),".","clus.",as.character(groupby),".pdf"),width=20,height=20)
        }}
    })


#final celltype assignment by umap phenograph
#IF in a cluster, AND not cancer > assign
obj<-generate_bigwig(obj=obj,
                        suffix="final_celltype",
                        groupBy="celltype",
                        outdir=getwd())

#renamed macro_mono to myeloid, since other cells present as well
obj@metadata[obj@metadata$celltype=="macro_mono",]$celltype<-"myeloid"

saveRDS(obj,file="03_scaledcis.final_celltypes.amethyst.rds")


```

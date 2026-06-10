# Running final cell typing
1. check major lineages against reclustered data
2. subset out major cell lineages
3. perform reclustering on vmrs for just those cells
4. use multiple resolutions of clustering to combine ambiguous clusters
5. output bigwigs in IGV and use RNA markers for pairwise discrimination of cell states
6. use lineage to sublineage to cell type markers to define cells
# Prepare RNA marker genes per cell type

From processing/milestonev1_00.2_seurat_scrna_copykat.md seurat object.
```R

#compare dmr sites with rna marker genes overlap
library(Seurat)
library(GenomicRanges)

rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")
update_group<-c('BCMDCIS05T'='DCIS',
'BCMDCIS07T'='DCIS',
'BCMDCIS41T'='DCIS',
'BCMDCIS66T'='DCIS',
'BCMDCIS82T_24hTis'='DCIS',
'BCMDCIS99T'='DCIS',
'ECIS25T'='DCIS',
'ECIS26T'='DCIS',
'ECIS36T'='DCIS',
'BCMHBCA03R'='HBCA',
'BCMHBCA04R'='HBCA',
'BCMHBCA09R-3h'='HBCA',
'BCMHBCA12R-3h'='HBCA',
'BCMHBCA16R-3h'='HBCA',
'BCMHBCA17R-3h'='HBCA',
'BCMHBCA19R-4h'='HBCA',
'BCMHBCA22R-4h'='HBCA',
'BCMHBCA26L-24hTis-4h'='HBCA',
'BCMHBCA29L-2h'='HBCA',
'BCMHBCA38L-3h'='HBCA',
'BCMHBCA83L-3h'='HBCA',
'BCMHBCA85L-3h'='HBCA',
'BCMDCIS102T_24hTis'='IDC',
'BCMDCIS124T'='IDC',
'BCMDCIS22T'='IDC',
'BCMDCIS28T'='IDC',
'BCMDCIS49T'='IDC',
'BCMDCIS52T'='IDC',
'BCMDCIS65T'='IDC',
'BCMDCIS70T'='IDC',
'BCMDCIS74T'='IDC',
'BCMDCIS79T_24hTis_IDC'='IDC',
'BCMDCIS94T_24hTis'='IDC',
'BCMDCIS97T'='IDC',
'ECIS57T'='IDC',
'BCMDCIS32T'='Synchronous',
'BCMDCIS35T'='Synchronous',
'BCMDCIS79T_24hTis_DCIS'='Synchronous',
'BCMDCIS80T_24hTis'='Synchronous',
'BCMDCIS92T_24hTis'='Synchronous',
'ECIS48T'='Synchronous')
saveRDS(rna,"/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")

#just run on all cell types
table(rna$fine_celltype)
rna$celltype<-NA
rna$celltype<-rna$fine_celltype
Idents(rna)<-rna$celltype
table(Idents(rna))
rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)
rna <- ScaleData(rna)
rna <- JoinLayers(rna)

#treg vs cd4
Idents(rna)<-rna$fine_celltype
rna_markers<-FindMarkers(rna,assay="RNA",ident.1="tcell_nk",ident.2=NULL,only.pos=TRUE)
rna_markers %>% mutate(gene=row.names(rna_markers)) %>% filter(p_val_adj<0.05) %>% filter(avg_log2FC>2) %>% head(n=20) %>% select(gene)

Idents(rna)<-rna$fine_celltype

#example comparisons for methylome checking
rna_markers<-FindMarkers(rna,assay="RNA",ident.1="plasma",ident.2="bcell",only.pos=TRUE)
rna_markers %>% mutate(gene=row.names(rna_markers)) %>% filter(p_val_adj<0.05) %>% filter(avg_log2FC>2) %>% head(n=20) %>% select(gene)

rna_markers<-FindAllMarkers(rna,assay="RNA",only.pos=TRUE)
saveRDS(rna_markers,file="/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.rna_markers.rds")

#then tackle myeloid
rna_markers<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.rna_markers.rds")
rna_markers %>% filter(p_val_adj<0.05) %>% filter(avg_log2FC>1) %>% head(n=20) %>% select(gene)

rna_markers %>% filter(p_val_adj<0.05) %>% filter(avg_log2FC>1) %>% filter(cluster=="tcell_cd4") %>% head(n=20) %>% select(gene)

```


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

project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"

#read in object from directory
task_cpus=300
processing_folder="03_fine_celltyping"
wd=paste(sep="/",project_data_directory,processing_folder)
system(paste0("mkdir -p ",wd))
setwd(wd)
obj<-readRDS(file=paste(project_data_directory,"02_copykit_cnv_calling","02_scaledcis.cnv_clones.amethyst.rds",sep="/"))
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

```


#subcluster on epithelial cells

```R
obj<-readRDS(file=paste(project_data_directory,"02_copykit_cnv_calling","02_scaledcis.cnv_clones.amethyst.rds",sep="/"))

#filter to immune cells by fine_celltype defined lineages
obj<-subsetObject(obj,cells=row.names(obj@metadata[obj@metadata$celltype_lineage %in% c("epithelial"),]))
dim(obj@metadata)

suffix="epithelial"
outdir="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/03_fine_celltyping/"
npc=25
reduction_name="irlba_epithelial"
leiden_cluster_resolution=0.5e-5
pheno_k=200
min_dist=1e-6
n_neighbors=5

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
  lapply(c("fine_cluster_leidenclus","fine_cluster_phenograph","Sample","mcg_pct","unique_reads","Group"), 
    function(groupby){
      if(groupby %in% colnames(obj@metadata)){
        plt_clus<-obj@metadata %>% 
          ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = obj@metadata[[groupby]])) +
            geom_point() +
            coord_fixed()
          ggsave(plt_clus,file=paste0(plot_dir,"/","03.1.VMR_umap.",suffix,".",as.character(npc),".",as.character(min_dist),".","clus.",as.character(groupby),".pdf"),width=20,height=20)
        }
    })


epithelial<-generate_bigwig(obj=obj,
                        suffix="epithelial",
                        groupBy="fine_cluster_leidenclus",
                        outdir=getwd())

saveRDS(epithelial,file="03_scaledcis.fine_celltyping.epithelial.rds")




#now assign cancer to cells in the same cluster and CNV called cancer cells
hm<-table(epithelial@metadata$cnv_ploidy_500kb,epithelial@metadata$fine_cluster_leidenclus)
hm_count<-as.data.frame(hm) #convert it back into a data frame, now with the counts
hm_scale<-as.data.frame(scale(hm,center=F))
colnames(hm_count)<-c("group","cluster","count")
hm_count$scaled <- hm_scale$Freq #combine count and scaled data to plot both
plt<-ggplot(hm_count, aes(x=group, y=cluster, fill=scaled,label=count)) + 
  geom_tile() + 
  geom_text(color="white") + 
  theme(axis.text.x = element_text(angle = 90))
ggsave(plt,file="03_scaledcis.fine_celltyping.epithelial.cancer_cluster.heatmap.pdf")



#assign celltype as cancer to those with >90% aneuploid cells within cluster
#this is to pick up any cells with NA calls for aneuploidy
celltype_assignment=c(
  'epithelial_4'='lumsec',
  'epithelial_1'='lumsec',
  'epithelial_8'='basal',

  'epithelial_6'='cancer',
  'epithelial_10'='cancer',
  'epithelial_11'='cancer',
  'epithelial_9'='cancer',
  'epithelial_13'='cancer',
  'epithelial_3'='cancer',
  'epithelial_14'='cancer',
  'epithelial_12'='cancer',

  'epithelial_7'='lumhr',
  'epithelial_2'='lumhr',
  'epithelial_5'='lumhr')


#get clusters with cancer clones assigned, assign the whole cluster as cancer
epithelial@metadata$celltype<-"NA"
epithelial@metadata$celltype<-celltype_assignment[epithelial@metadata$fine_cluster_leidenclus]
epithelial@metadata[epithelial@metadata$cnv_ploidy_500kb=="aneuploid",]$celltype<-"cancer"
#and assign any aneuploidy cells as cancer


 plt_celltype<-epithelial@metadata %>% 
    ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = celltype)) +
    geom_point() +
    coord_fixed()
  ggsave(plt_celltype,file="03.1.epithelial.celltype.umap.pdf",width=20,height=20)

saveRDS(epithelial,file="03_scaledcis.fine_celltyping.epithelial.rds")




#assign cell types via IGV tracks


```

# Cluster on VMRs per stromal cells

```R
obj<-readRDS(file="03_scaledcis.final_celltypes.amethyst.rds")

#filter to stromal cells by fine_celltype defined lineages
stromal<-subsetObject(obj,cells=row.names(obj@metadata[obj@metadata$subcluster_group %in% c("stromal"),]))

#refining coarse cell types before subset and subcluster
stromal<-fine_celltype_vmr_cluster(obj=stromal,
                    suffix="stromal",
                    npc=15, #20 is decent, 30 is a bit better, trying 25dims
                    outdir=getwd(),
                    leiden_cluster_resolution=5e-6,
                    reduction_name="irlba_stromal",
                    pheno_k=150,
                    #feat_filt=10000,
                    min_dist=1e-4, #1e-5,
                    n_neighbors=8)

stromal<-generate_bigwig(obj=stromal,
                        suffix="stromal",
                        groupBy="fine_cluster_phenograph",
                        outdir=getwd())

#now assign cancer to cells in the same cluster and CNV called cancer cells
hm<-table(stromal@metadata$Group,stromal@metadata$fine_cluster_phenograph)
hm_count<-as.data.frame(hm) #convert it back into a data frame, now with the counts
hm_scale<-as.data.frame(scale(hm,center=F))
colnames(hm_count)<-c("group","cluster","count")
hm_count$scaled <- hm_scale$Freq #combine count and scaled data to plot both
plt<-ggplot(hm_count, aes(x=group, y=cluster, fill=scaled,label=count)) + 
  geom_tile() + 
  geom_text(color="white") + 
  theme(axis.text.x = element_text(angle = 90))
ggsave(plt,file="03_scaledcis.fine_celltyping.stromal.cancer_cluster.heatmap.pdf")

#assign cell types via IGV tracks
celltype_assignment=c(
  'stromal_1'='endothelial_1',
  'stromal_5'='endothelial_2',
  'stromal_2'='pericytes_1',
  'stromal_13'='pericytes_2',
  'stromal_8'='pericytes_2',
  'stromal_3'='fibroblast_1',
  'stromal_12'='fibroblast_1',
  'stromal_7'='fibroblast_1',
  'stromal_9'='fibroblast_2',
  'stromal_11'='fibroblast_2',
  'stromal_6'='fibroblast_3',
  'stromal_4'='fibroblast_3',
  'stromal_10'='fibroblast_4'
)

#get clusters with cancer clones assigned, assign the whole cluster as cancer
stromal@metadata$fine_celltype<-"NA"
stromal@metadata$fine_celltype<-celltype_assignment[stromal@metadata$fine_cluster_phenograph]

 plt_celltype<-stromal@metadata %>% 
    ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = fine_celltype)) +
    geom_point() +
    coord_fixed()
  ggsave(plt_celltype,file="03.1.stromal.celltype.umap.pdf",width=20,height=20)

saveRDS(stromal,file="03_scaledcis.fine_celltyping.stromal.rds")

```
# Cluster on VMRs per epithelial cells

```R
obj<-readRDS(file="03_scaledcis.final_celltypes.amethyst.rds")
obj@metadata[which(obj@metadata$coarse_celltype=="cancer"),]$celltype_lineage<-"epithelial"

#filter to immune cells by fine_celltype defined lineages
epitheial<-subsetObject(obj,cells=row.names(obj@metadata[obj@metadata$subcluster_group %in% c("epithelial"),]))

#refining coarse cell types before subset and subcluster
epithelial<-fine_celltype_vmr_cluster(obj=epithelial,
                    suffix="epithelial",
                    npc=20, #20 is decent, 30 is a bit better, trying 25dims
                    outdir=getwd(),
                    leiden_cluster_resolution=5e-6,
                    reduction_name="irlba_epithelial",
                    pheno_k=150,
                    #feat_filt=10000,
                    min_dist=1e-4, #1e-5,
                    n_neighbors=8)

epithelial<-generate_bigwig(obj=epithelial,
                        suffix="epithelial",
                        groupBy="fine_cluster_phenograph",
                        outdir=getwd())

#assign cell types via IGV tracks
celltype_assignment=c(
  'epithelial_8'='lumsec',
  'epithelial_11'='lumsec',
  'epithelial_21'='lumsec',
  'epithelial_15'='lumsec',
  'epithelial_12'='lumsec',
  'epithelial_3'='lumsec',
  'epithelial_19'='basal',
  'epithelial_16'='basal',
  'epithelial_18'='basal',
  'epithelial_13'='lumhr',
  'epithelial_22'='lumhr',
  'epithelial_4'='lumhr',
  'epithelial_5'='lumhr',
  'epithelial_1'='lumhr',
  'epithelial_26'='lumhr',
  'epithelial_7'='lumhr',
  'epithelial_17'='lumhr',
  'epithelial_14'='lumhr',
  'epithelial_9'='lumhr',
  'epithelial_24'='lumhr',
  'epithelial_20'='lumhr',
  'epithelial_6'='lumhr',
  'epithelial_29'='lumhr',
  'epithelial_23'='lumhr',
  'epithelial_27'='lumhr',
  'epithelial_30'='lumhr',
  'epithelial_25'='lumhr',
  'epithelial_28'='lumhr',
  'epithelial_2'='lumhr',
  'epithelial_10'='lumhr',
  'epithelial_32'='lumhr',
  'epithelial_31'='lumhr'
)

#get clusters with cancer clones assigned, assign the whole cluster as cancer
epithelial@metadata$fine_celltype<-"NA"
epithelial@metadata$fine_celltype<-celltype_assignment[epithelial@metadata$fine_cluster_phenograph]
epithelial@metadata[epithelial@metadata$cnv_clonename!="NA" & !endsWith(epithelial@metadata$cnv_clonename,suffix="diploid"),]$fine_celltype<-"cancer"

#now assign cancer to cells in the same cluster and CNV called cancer cells
hm<-table(epithelial@metadata$fine_celltype,epithelial@metadata$fine_cluster_phenograph)
hm_count<-as.data.frame(hm) #convert it back into a data frame, now with the counts
hm_scale<-as.data.frame(scale(hm,center=F))
colnames(hm_count)<-c("celltype","cluster","count")
hm_count$scaled <- hm_scale$Freq #combine count and scaled data to plot both
plt<-ggplot(hm_count, aes(x=celltype, y=cluster, fill=scaled,label=count)) + 
  geom_tile() + 
  geom_text(color="white") + 
  theme(axis.text.x = element_text(angle = 90))
ggsave(plt,file="03_scaledcis.fine_celltyping.epithelial.cancer_cluster.heatmap.pdf")

#cancer clusters have >50% cancer cells
cancer_assignment=c(
  'epithelial_7'='cancer',
  'epithelial_32'='cancer',
  'epithelial_31'='cancer',
  'epithelial_30'='cancer',
  'epithelial_29'='cancer',
  'epithelial_28'='cancer',
  'epithelial_27'='cancer',
  'epithelial_26'='cancer',
  'epithelial_24'='cancer',
  'epithelial_23'='cancer',
  'epithelial_20'='cancer',
  'epithelial_10'='cancer'
)
epithelial@metadata[epithelial@metadata$fine_cluster_phenograph %in% names(cancer_assignment) & epithelial@metadata$Group!="HBCA",]$fine_celltype<-"cancer"

 plt_celltype<-epithelial@metadata %>% 
    ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = fine_celltype)) +
    geom_point() +
    coord_fixed()
  ggsave(plt_celltype,file="03.1.epithelial.celltype.umap.pdf",width=20,height=20)

 plt_celltype<-epithelial@metadata %>% 
    ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = cnv_clonename)) +
    geom_point() +
    coord_fixed()
  ggsave(plt_celltype,file="03.1.epithelial.clonename.umap.pdf",width=20,height=20)

 plt_celltype<-epithelial@metadata %>% 
    ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = cnv_ploidy_500kb)) +
    geom_point() +
    coord_fixed()
  ggsave(plt_celltype,file="03.1.epithelial.ploidy.umap.pdf",width=20,height=20)


saveRDS(epithelial,file="03_scaledcis.fine_celltyping.epithelial.rds")
```

# Write out fine cell types into shared amethyst object.

```R

set.seed(111)
options(future.globals.maxSize= 500000*1024^2) #80gb limit for parallelizing
library(amethyst)
library(data.table)
library(rtracklayer)

#read in object from directory
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
processing_folder="03_fine_celltyping"
task_cpus=300
wd=paste(sep="/",project_data_directory,processing_folder)
setwd(wd)
obj<-readRDS(file="03_scaledcis.final_celltypes.amethyst.rds")

epithelial<-readRDS(file="03_scaledcis.fine_celltyping.epithelial.rds")
immune<-readRDS(file="03_scaledcis.fine_celltyping.immune.rds")
stromal<-readRDS(file="03_scaledcis.fine_celltyping.stromal.rds")

obj@metadata$fine_cluster_phenograph<-NA
obj@metadata$final_celltype<-NA

obj@metadata[row.names(epithelial@metadata),]$fine_cluster_phenograph<-epithelial@metadata$fine_cluster_phenograph
obj@metadata[row.names(immune@metadata),]$fine_cluster_phenograph<-immune@metadata$fine_cluster_phenograph
obj@metadata[row.names(stromal@metadata),]$fine_cluster_phenograph<-stromal@metadata$fine_cluster_phenograph

obj@metadata[row.names(epithelial@metadata),]$final_celltype<-epithelial@metadata$fine_celltype
obj@metadata[row.names(immune@metadata),]$final_celltype<-immune@metadata$fine_celltype
obj@metadata[row.names(stromal@metadata),]$final_celltype<-stromal@metadata$fine_celltype

obj<-fine_celltype_vmr_cluster(obj=obj,
                    suffix="final_celltype",
                    npc=8, #20 is decent, 30 is a bit better, trying 25dims
                    outdir=getwd(),
                    leiden_cluster_resolution=5e-6,
                    reduction_name="irlba_final_celltype",
                    pheno_k=100,
                    #feat_filt=10000,
                    min_dist=1e-5,
                    n_neighbors=8)

obj@metadata$final_celltype_UMAP_X<-obj@metadata$fine_cluster_UMAP_X
obj@metadata$final_celltype_UMAP_Y<-obj@metadata$fine_cluster_UMAP_Y

obj@metadata[row.names(epithelial@metadata),]$fine_cluster_UMAP_X<-epithelial@metadata$fine_cluster_UMAP_X
obj@metadata[row.names(immune@metadata),]$fine_cluster_UMAP_X<-immune@metadata$fine_cluster_UMAP_X
obj@metadata[row.names(stromal@metadata),]$fine_cluster_UMAP_X<-stromal@metadata$fine_cluster_UMAP_X

obj@metadata[row.names(epithelial@metadata),]$fine_cluster_UMAP_Y<-epithelial@metadata$fine_cluster_UMAP_Y
obj@metadata[row.names(immune@metadata),]$fine_cluster_UMAP_Y<-immune@metadata$fine_cluster_UMAP_Y
obj@metadata[row.names(stromal@metadata),]$fine_cluster_UMAP_Y<-stromal@metadata$fine_cluster_UMAP_Y

saveRDS(obj,file="03_scaledcis.final_celltypes.amethyst.rds")

plt_celltype<-obj@metadata %>% 
    ggplot(aes(x = final_celltype_UMAP_X, y = final_celltype_UMAP_Y, color = final_celltype)) +
    geom_point() +
    coord_fixed()
ggsave(plt_celltype,file=paste0("03.0.VMR_umap.","fine_celltype",".final_celltype.pdf"),width=20,height=20)

#355 NA cells being assigned by cluster cell type
#assign same cell types to shared phenograph cluster
tab<-table(obj@metadata$fine_celltype,obj@metadata$fine_cluster_phenograph)
max_tab <- setNames(nm=colnames(tab),row.names(tab)[apply(tab, 2, which.max)])
for(ind in which(is.na(obj@metadata$fine_celltype))){
  obj@metadata[ind,]$final_celltype<-max_tab[obj@metadata[ind,]$fine_cluster_phenograph]
}

plt_celltype<-obj@metadata %>% 
    ggplot(aes(x = final_celltype_UMAP_X, y = final_celltype_UMAP_Y, color = final_celltype)) +
    geom_point() +
    coord_fixed()
ggsave(plt_celltype,file=paste0("03.0.VMR_umap.","fine_celltype",".final_celltype.pdf"),width=20,height=20)
saveRDS(obj,file="03_scaledcis.final_celltypes.amethyst.rds")

obj<-generate_bigwig(obj=obj,
                        suffix="final_celltype",
                        groupBy="final_celltype",
                        outdir=getwd())

saveRDS(obj,file="03_scaledcis.final_celltypes.amethyst.rds")
```

For subclusters per celltype run DMRs to differentiate
-run fibroblasts one vs rest
  * fibroblast_1 calling all fibroblast
-run immune unknown one vs rest
  * immune_unknown_1 = AIF1+/PTPRClow/CST7+/CLEC9A+/FLT1+/LST1+ = dendritic_cell
  * immune_unknown_2 = AIF1-/PTPRClow/IRF7+/CD3D-/KIT+ = mast
  * immune_unknown_3 = AIF1-/PTPRChigh/CSF3R+/NAMPT+ = neutrophil
-run tcells one vs rest
  *tcell_1 = tcell
  *tcell_2 = tcell
  *tcell_3 = tcell
  *tcell_4 = tcell
  *tcell_5 = tcell
  *nk PTPRC+ CD3E- = nkcell

-run endothelial one vs rest
  * endothelial_1 = endothelial_lymphatic
  * endothelial_2 = endothelial_vascular
-run pericytes one v one
  *pericytes_1 = COL6A3 GUCY1A2 ADAMTS12 CCDC102B CD36 MYH11 MUSTN1 MYOCD = pericyte_VSMC 
  *pericytes_2 = COL6A3 GUCY1A2 ADAMTS12 CCDC102B CD36 MYH11 MUSTN1 MYOCD =  pericyte_VSMC

Check these via RNA markers first. Then DMR if necessary.

```R
table(Idents(rna))

rna_markers_all<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.rna_markers.rds")
rna_markers_all %>% filter(p_val_adj<0.05) %>% filter(avg_log2FC>1) %>% filter(cluster=="myeloid_mast") %>% head(n=20) %>% select(gene)

rna_markers<-FindMarkers(rna,assay="RNA",ident.1="tcell_cd4",ident.2="tcell_cd8",only.pos=TRUE)
rna_markers %>% mutate(gene=row.names(rna_markers)) %>% filter(p_val_adj<0.05) %>% filter(avg_log2FC>2) %>% head(n=20) %>% select(gene)
#pericytes_1 and pericytes_2: 
#pericytes COL6A3 GUCY1A2 ADAMTS12 CCDC102B CD36
#VSMC MYH11 MUSTN1 MYOCD
#they both express both markers, so calling pericyte_VSMC


#example of processing using immune_unknown_1 and checking IGV of top expressed genes
#using RNA markers for b cells vs plasma
#example comparisons for methylome checking
rna_markers<-FindMarkers(rna,assay="RNA",ident.1="plasma",ident.2="bcell",only.pos=TRUE)
rna_markers %>% mutate(gene=row.names(rna_markers)) %>% filter(p_val_adj<0.05) %>% filter(avg_log2FC>2) %>% head(n=20) %>% select(gene)
#bcell vs plasma
#bcell CD83 BANK1 HLA-DPB1
#plasma MZB1 DERL3 NUCB2 PRDM1
#"immune_unknown_1"="plasmacell"


celltype<-c(
  "fibroblast_1"="fibroblast",
  "fibroblast_2"="fibroblast",
  "fibroblast_3"="fibroblast",
  "fibroblast_4"="fibroblast",
  "fibroblast_5"="fibroblast",

  "immune_unknown_1"="myeloid_DC",
  "immune_unknown_2"="myeloid_mast",
  "immune_unknown_3"="myeloid_neutrophil",

  "tcell_1"="tcell",
  "tcell_2"="tcell",
  "tcell_3"="tcell",
  "tcell_4"="tcell",
  "tcell_5"="tcell",
  "nk"="nkcell",

  "endothelial_1"="endothelial_lymphatic",
  "endothelial_2"="endothelial_vascular",
  "pericytes_1"="pericyte_VSMC",
  "pericytes_2"="pericyte_VSMC"
)

#overwriting final cell type with this set
obj@metadata[obj@metadata$final_celltype %in% names(celltype),]$final_celltype<-celltype[obj@metadata[obj@metadata$final_celltype %in% names(celltype),]$final_celltype]


plt_celltype<-obj@metadata %>% 
    ggplot(aes(x = final_celltype_UMAP_X, y = final_celltype_UMAP_Y, color = final_celltype)) +
    geom_point() +
    coord_fixed()
ggsave(plt_celltype,file=paste0("03.0.VMR_umap.","fine_celltype",".final_celltype.pdf"),width=20,height=20)

obj<-generate_bigwig(obj=obj,
                        suffix="final_celltype",
                        groupBy="final_celltype",
                        outdir=getwd())

saveRDS(obj,file="03_scaledcis.final_celltypes.amethyst.rds")

#cluster a bit better

obj2<-fine_celltype_vmr_cluster(obj=obj,
                    suffix="fine_celltype_2",
                    npc=50, #20 is decent, 30 is a bit better, trying 25dims
                    outdir=getwd(),
                    leiden_cluster_resolution=5e-6,
                    reduction_name="irlba_fine_celltype2",
                    pheno_k=150,
                    #feat_filt=10000,
                    min_dist=1e-4, #1e-5,
                    n_neighbors=8)

plt_celltype<-obj2@metadata %>% 
    ggplot(aes(x = fine_cluster_UMAP_X, y = fine_cluster_UMAP_Y, color = final_celltype)) +
    geom_point() +
    coord_fixed()
ggsave(plt_celltype,file=paste0("03.0.VMR_umap.","fine_celltype",".final_celltype2.pdf"),width=20,height=20)
               
```



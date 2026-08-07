
NOTE: THIS DIDNT WORK WELL
PASSED DATA TO REEM TO TRY
Using CCA on matched genomic windows to try and integrate RNA and MET CNVs

```R
library(GenomicRanges)
library(copykit)
library(ComplexHeatmap)
library(parallel)
library(BiocParallel)
set.seed(111)
library(dplyr)
library(data.table)
library(scquantum)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dendextend)
library(amethyst)
library(irlba)
library(Seurat)
library(Signac)

task_cpus=150
seed=1234
register(MulticoreParam(progressbar = T, workers = task_cpus), default = T)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing


project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
resolution="500kb"

#read in object from directory
processing_folder="02_copykit_cnv_calling"
wd=paste(sep="/",project_data_directory,processing_folder)
setwd(wd)

#read in MET
obj<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/03_fine_celltyping/03_scaledcis.final_celltypes.amethyst.rds")
rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")

#read in cnv profiling (raw data, not integer copies)
scCNA<-readRDS(file="02_all_samples.all_cells.preintegerprocessing.copykit.Rds")

#read in RNA
scRNA<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.copykat_allcells.rds")
scRNA_meta<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.copykat_meta.rds")

#summarize scRNA over scCNA windows
scRNA_ranges<-data.frame(chr=unlist(lapply(strsplit(row.names(scRNA),"_"),"[",1)),
                start=as.numeric(unlist(lapply(strsplit(row.names(scRNA),"_"),"[",3))),
                end=as.numeric(unlist(lapply(strsplit(row.names(scRNA),"_"),"[",3)))+500 
                )

scRNA_ranges<-GRanges(scRNA_ranges)

overlap<-findOverlaps(query=scCNA@rowRanges,subject=scRNA_ranges)
overlap_list<-split(overlap@to,f=overlap@from)

#generating mean value for all bins overlapping 500kb filtered bins
binned_rna<-mclapply(1:length(overlap_list),function(i) {
    x=as.data.frame(colMeans(scRNA[overlap_list[[i]],],na.rm=T))
    colnames(x)<-names(overlap_list)[i]
    return(x)
},mc.cores=50)

rna_500kb <- do.call("cbind",binned_rna)
rna_500kb <-t(rna_500kb)

#filter any met CNV bins not present in the rna
met_500kb <- scCNA@assays@data$segment_ratios[as.numeric(row.names(rna_500kb)),]

#filter rna and met to only cancer or lumhr cells
met_500kb<-met_500kb[,colnames(met_500kb) %in% row.names(obj@metadata)[obj@metadata$celltype %in% c("cancer","lumhr")]]

rna_500kb<-rna_500kb[,colnames(rna_500kb) %in% row.names(rna@meta.data)[rna@meta.data$celltype %in% c("cancer","lumhr")]]


# Create standard Seurat objects
met_obj <- CreateSeuratObject(counts = met_500kb, 
                                data = met_500kb, 
                                assay = "cna")
met_obj$orig.ident<-"met"
met_obj<-AddMetaData(met_obj,
        setNames(nm=row.names(obj@metadata),
            obj@metadata$cnv_clonename_500kb),col.name="met_clonename")
met_obj<-AddMetaData(met_obj,
        setNames(nm=row.names(obj@metadata),
            obj@metadata$Sample),col.name="sample")
LayerData(object = met_obj, assay = "cna", layer = "scale.data")<-as.matrix(met_500kb)
met_obj <- RunPCA(met_obj,features=Features(met_obj),npcs=30)

rna_obj <- CreateSeuratObject(counts = rna_500kb, 
                                data=rna_500kb, 
                                assay = "cna")
rna_obj$orig.ident<-"rna"
rna_obj<-AddMetaData(rna_obj,
        setNames(nm=row.names(rna@meta.data),
                rna@meta.data$rna_clonename),col.name="rna_clonename")
rna_obj<-AddMetaData(rna_obj,
        setNames(nm=row.names(rna@meta.data),
                rna@meta.data$sample),col.name="sample")
LayerData(object = rna_obj, assay = "cna", layer = "scale.data")<-as.matrix(rna_500kb)
rna_obj <- RunPCA(rna_obj,features=Features(rna_obj),npcs=30)

anchors <- FindTransferAnchors(
  reference = met_obj,
  reference.assay="cna",
  query.assay="cna",
  reference.reduction="pca",
  query = rna_obj,
  reduction = "pcaproject",          
  features = Features(rna_obj),
  scale = FALSE
)

predicted_clones <- TransferData(
  anchorset = anchors,
  refdata = met_obj$met_clonename,k.weight=30)

rna_obj <- AddMetaData(rna_obj, metadata = predicted_clones)

rna_obj <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = met_obj,
  query = rna_obj,k.weight=30,
  new.reduction.name = "integrated_cna_pca"
)
met_obj[["integrated_cna_pca"]] <- met_obj[["pca"]]

merged_obj <- merge(met_obj, rna_obj,  merge.dr = c("integrated_cna_pca"))

# Run UMAP on the integrated coordinate reduction
merged_obj <- RunUMAP(merged_obj, reduction = "integrated_cna_pca", dims = 1:20)
plt1<-DimPlot(merged_obj, reduction="umap", group.by = c("orig.ident", "sample")) # View alignment of DNA and RNA cells
ggsave(plt1,file="test.pdf",width=50,height=50,limitsize=F)



merged_obj <- merge(met_obj, rna_obj)


merged_obj <- FindNeighbors(merged_obj, 
                            dims = 1:50, 
                            reduction = "pca")

merged_obj <- FindClusters(merged_obj, 
                            resolution = 2, 
                            cluster.name = "unintegrated_clusters")

merged_obj <- RunUMAP(merged_obj, 
                            dims = 1:50, 
                            reduction = "pca", 
                            reduction.name = "umap.unintegrated")

plt1<-DimPlot(merged_obj, 
                reduction = "umap.unintegrated", 
                group.by = c("orig.ident", "sample"))


merged_obj <- IntegrateLayers(object = merged_obj, 
                            method = CCAIntegration, 
                            orig.reduction = "pca", 
                            new.reduction = "integrated.cca", 
                            verbose = TRUE, 
                            features=Features(merged_obj))

merged_obj <- IntegrateLayers(object = merged_obj, 
                            method = RPCAIntegration, 
                            orig.reduction = "pca",
                            new.reduction = "integrated.rpca", 
                            verbose = TRUE,
                            features=Features(merged_obj))


plt1<-DimPlot(merged_obj, 
                reduction = "integrated.rpca", 
                group.by = c("orig.ident", "sample"))

ggsave(plt1,file="test.pdf",width=50,height=50,limitsize=F)

merged_obj <- IntegrateLayers(object = merged_obj, 
                                method = HarmonyIntegration, 
                                orig.reduction = "pca",
                                new.reduction = "harmony", 
                                verbose = TRUE,
                                features=Features(merged_obj))

plt1<-DimPlot(merged_obj, 
                reduction = "harmony", 
                group.by = c("orig.ident", "sample"))

ggsave(plt1,file="test.pdf",width=50,height=50,limitsize=F)

obj <- IntegrateLayers(object = obj, method = FastMNNIntegration, new.reduction = "integrated.mnn",
    verbose = FALSE)
obj <- IntegrateLayers(object = obj, method = scVIIntegration, new.reduction = "integrated.scvi",
    conda_env = "../miniconda3/envs/scvi-env", verbose = FALSE)

# Process DNA (CopyKit space)
met_obj <- ScaleData(met_obj)
met_obj <- RunPCA(met_obj, features=row.names(met_500kb),reduction.name = "dna_pca")

# Process DNA (CopyKat space)
rna_obj <- ScaleData(rna_obj)
rna_obj <- RunPCA(rna_obj, features=row.names(rna_500kb), reduction.name = "rna_pca")

anchors <- FindTransferAnchors(
  reference = met_obj,
  query = rna_obj,
  reference.reduction = "dna_pca",
  reduction = "pcaproject",
  features = row.names(rna_500kb)
)

rna_obj <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = met_obj,
  query = rna_obj,
  reductions = "pcaproject"
)

# Run a single joint UMAP

merged_obj <- RunUMAP(merged_obj, reduction = "dna_pca", dims = 1:30)
DimPlot(merged_obj, group.by = "orig.ident") # Colors by DNA vs RNA platform



library(CCA)
result <- cc(met_500kb, rna_500kb)

#now time to integrate with CCA
cross_cov <- t(scale(met_500kb)) %*% scale(rna_500kb)  # P_bins x P_bins
K=50
svd_cca <- irlba(cross_cov, nu = K, nv = K)

# Step 5: Extract Canonical Loadings (Feature Weights)
loadings_met <- svd_cca$u  # P_bins x K (CopyKit weights)
loadings_rna <- svd_cca$v  # P_bins x K (CopyKAT weights)

# Step 6: Compute the Low-Dimensional Aligned Cell Embeddings
# Project the scaled data matrices into the low-dimensional shared subspace
shared_embedding_X <- scale(met_500kb) %*% loadings_met  # N_cells x K
shared_embedding_Y <- scale(rna_500kb) %*% loadings_rna  # N_cells x K

# Step 7: L2-Normalize the cell embeddings
# This step is critical in single-cell integration to remove scale differences
shared_embedding_X <- shared_embedding_X / sqrt(rowSums(shared_embedding_X^2))
shared_embedding_Y <- shared_embedding_Y / sqrt(rowSums(shared_embedding_Y^2))

# Step 8: Create the integrated consensus matrix
integrated_cca_space <- (shared_embedding_X + shared_embedding_Y) / 2

umap_coords <- uwot::umap(integrated_cca_space)
pdf("test.pdf")
plot(umap_coords, col = "blue", pch = 16, main = "Integrated Space")
dev.off()
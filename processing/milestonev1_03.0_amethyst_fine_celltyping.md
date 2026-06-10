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

project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"

#read in object from directory
task_cpus=300
processing_folder="03_fine_celltyping"
wd=paste(sep="/",project_data_directory,processing_folder)
system(paste0("mkdir -p ",wd))
setwd(wd)
obj<-readRDS(file=paste(project_data_directory,"02_copykit_cnv_calling","02_scaledcis.cnv_clones.amethyst.rds",sep="/"))


#confirming cell type labelling by cnvs (cancer vs lumhr)
write.table(obj@metadata,col.names=T,row.names=T,paste0(project_data_directory,"/","metadata.csv"))

suffix="coarse_celltype"
# validation plots on clones and coarse cell types
plt_clonename<-obj@metadata %>% 
    ggplot(aes(x = coarse_cluster_UMAP_X, y = coarse_cluster_UMAP_Y, color = cnv_clonename)) +
    geom_point() +
    coord_fixed()
ggsave(plt_clonename,file=paste0("02.0.VMR_umap.",suffix,".clone.pdf"),width=20,height=20)

plt_celltype<-obj@metadata %>% 
    ggplot(aes(x = coarse_cluster_UMAP_X, y = coarse_cluster_UMAP_Y, color = coarse_celltype)) +
    geom_point() +
    coord_fixed()
ggsave(plt_celltype,file=paste0("02.0.VMR_umap.",suffix,".coarse_celltype.pdf"),width=20,height=20)

plt_celltype<-obj@metadata %>% 
    ggplot(aes(x = coarse_cluster_UMAP_X, y = coarse_cluster_UMAP_Y, color = celltype_lineage)) +
    geom_point() +
    coord_fixed()
ggsave(plt_celltype,file=paste0("02.0.VMR_umap.",suffix,".celltype_lineage.pdf"),width=20,height=20)

plt_celltype<-obj@metadata %>% 
    ggplot(aes(x = coarse_cluster_UMAP_X, y = coarse_cluster_UMAP_Y, color = cnv_ploidy_500kb)) +
    geom_point() +
    coord_fixed()
ggsave(plt_celltype,file=paste0("02.0.VMR_umap.",suffix,".cnv_ploidy_500kb.pdf"),width=20,height=20)

#theres one cluster that is a little suspect. has aneuploid cells but is labelled as basal
#coarse_cluster_phenograph 17 has strong basal epithelial signature despite aneuploidy 
#looks like its 97T c1 (which has big and clear CNV events), the diploid cells also seem to be from 97T which suggests some patient specific clustering?
#besides that, happy with clustering!
#ffpe looks like Sclerosing adenosis which is only known precurser of tnbc, maybe reflects basal like
#check for loss of er

```
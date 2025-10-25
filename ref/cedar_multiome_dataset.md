Took peaks called from CEDAR breast cancer multiome project to use with celltyping via methylation.

Performed on Exacloud and exported.
```R
dat<-readRDS("6_merged.celltyping.SeuratObject.rds")

markers <- presto:::wilcoxauc.Seurat(X = dat, group_by = "assigned_celltype", 
  groups_use=unname(unlist(unique(dat@meta.data$assigned_celltype))),
  y=unname(unlist(unique(dat@meta.data$assigned_celltype))), 
  assay = 'data', seurat_assay = "ATAC")

library(dplyr)
markers<-markers %>% filter(pval<0.05)
saveRDS(markers,file="celltypes.markers.peaks.df.rds")
write.table(markers,file="celltypes.markers.peaks.df.tsv",col.names=T,sep="\t",quote=F,row.names=F)
```

Stored on sbio at 
/data/rmulqueen/projects/scalebio_dcis/ref/CEDAR_multiome.celltypes.markers.peaks.df.tsv
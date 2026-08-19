```R

library(amethyst)

library(dplyr)
library(Matrix)

library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(cowplot)
library(patchwork)

library(parallel)
options(future.globals.maxSize= 500000*1024^2) #80gb limit for parallelizing
options(scipen = 0)
set.seed(111)
```

Set up environment
```R
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"



#set up output directory
output_directory=paste0(project_data_directory,"/","test_methyltreeplots")
system(paste0("mkdir -p ",output_directory))
setwd(output_directory)

#read methylation
obj<-readRDS(file=paste(project_data_directory,"03_fine_celltyping","03_scaledcis.final_celltypes.amethyst.rds",sep="/"))



library(ComplexHeatmap)
library(RColorBrewer)

dat<-read.table("/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/BCMDCIS41T_ALL_similarity.csv",sep=",",header=T,row.names=1)

cancerclone_col=setNames(
        nm=unique(as.character(obj@metadata[row.names(dat),]$cnv_clonename_500kb)),
        colorRampPalette(brewer.pal(8, "Spectral"))(length(unique(as.character(obj@metadata[row.names(dat),]$cnv_clonename_500kb)))))

#make column annotations
ha = rowAnnotation(
    cnv_clones=obj@metadata[row.names(dat),]$cnv_clonename_500kb,
    celltype=obj@metadata[row.names(dat),]$celltype,
    col= list(
        cnv_clones=cancerclone_col,
        celltype=celltype_col))


plt<-Heatmap(dat,show_row_names=FALSE,show_column_names=FALSE,right_annotation=ha)

pdf("test.pdf")
print(plt)
dev.off()

dat2<-read.table("/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/BCMDCIS41T_ALL_cnvfilt_similarity.csv",sep=",",header=T,row.names=1)

cancerclone_col=setNames(
        nm=unique(as.character(obj@metadata[row.names(dat2),]$cnv_clonename_500kb)),
        colorRampPalette(brewer.pal(8, "Spectral"))(length(unique(as.character(obj@metadata[row.names(dat2),]$cnv_clonename_500kb)))))

#make column annotations
ha = rowAnnotation(
    cnv_clones=obj@metadata[row.names(dat2),]$cnv_clonename_500kb,
    celltype=obj@metadata[row.names(dat2),]$celltype,
    col= list(
        cnv_clones=cancerclone_col,
        celltype=celltype_col))


plt<-Heatmap(dat2,show_row_names=FALSE,show_column_names=FALSE,right_annotation=ha)

pdf("test2_nocnv.pdf")
print(plt)
dev.off()

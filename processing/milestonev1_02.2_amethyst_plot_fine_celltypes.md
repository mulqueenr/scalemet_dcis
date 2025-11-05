```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
```


```R
source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)

obj<-readRDS(file="04_scaledcis.broad_celltype.amethyst.rds")
obj@metadata$fine_celltype<-obj@metadata$broad_celltype

#stromal
stromal<-"/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/fine_celltyping/stromal/04_scaledcis.stromal.fine_celltype.amethyst.rds"
stromal<-readRDS(stromal)
stromal_celltypes<-setNames(stromal@metadata[which(row.names(stromal@metadata) %in% row.names(obj@metadata)),]$fine_celltypes,nm=row.names(stromal@metadata))
obj@metadata[names(stromal_celltypes),]$fine_celltype<-stromal_celltypes

#immune
immune<-"/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/fine_celltyping/immune/04_scaledcis.immune.fine_celltype.amethyst.rds"
immune<-readRDS(immune)
immune_celltypes<-setNames(immune@metadata[which(row.names(immune@metadata) %in% row.names(obj@metadata)),]$fine_celltypes,nm=row.names(immune@metadata))
obj@metadata[names(immune_celltypes),]$fine_celltype<-immune_celltypes

#summarize celltype-dmr clusters over windows
obj<-dmr_and_1kb_window_gen(obj,
    prefix="fine_celltypes",
    groupBy="fine_celltype",
    threads=10,step=500)

saveRDS(obj,file="05_scaledcis.fine_celltype.amethyst.rds")
write.table(obj@metadata,file="05_scaledcis.fine_celltype.metadata.tsv",quote=F,col.names=T,row.names=T,sep="\t")

obj<-readRDS(file="01_celllines.amethyst.rds")
write.table(obj@metadata,file="01_celllines.metadata.tsv",quote=F,col.names=T,row.names=T,sep="\t")
```

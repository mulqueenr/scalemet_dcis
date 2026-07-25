library(amethyst)
library(dplyr)
library(ggplot2)
library(patchwork)

#source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"

#read in object from directory
processing_folder="10_fig_plots"
wd=paste(sep="/",project_data_directory,processing_folder)
outdir=paste0(wd,"/","figure_1")

system(paste0("mkdir -p ",wd))
setwd(wd)

obj<-readRDS(file=paste(project_data_directory,"03_fine_celltyping","03_scaledcis.final_celltypes.amethyst.rds",sep="/"))

#unique reads per cell type
t.test(unique_reads ~ method, data=obj@metadata, na.action=na.omit)

anova_model<-aov(mcg_pct ~ method*batch, data=obj@metadata, na.action=na.omit)
summary(anova_model)
TukeyHSD(anova_model)

t.test(mch_pct ~ method+celltype, data=obj@metadata, na.action=na.omit)


ggplot(obj@metadata,aes(x=method,y=log10(unique_reads)))+geom_violin()+geom_scatter()

# Plotting sample overview and cell typing (Figure 1)

```R
library(amethyst)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
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
outdir=paste0(wd,"/","figure_supp_celltyping")
system(paste0("mkdir -p ",wd))
setwd(wd)

obj<-readRDS(file=paste(project_data_directory,"03_fine_celltyping","03_scaledcis.final_celltypes.amethyst.rds",sep="/"))

rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")

####################################################
#           Supp Fig Celltype UMAPS                #
###################################################
#these are from 
#scalemet_dcis/processing/milestonev1_03.1_epithelial_fine_celltyping.md
#scalemet_dcis/processing/milestonev1_03.1_stromal_fine_celltyping.md
#scalemet_dcis/processing/milestonev1_03.1_immune_fine_celltyping.md
#plot_umap_panels_met and plot_umap_panels_rna function calls


####################################################
#           Supp Fig Celltype met tracks                #
###################################################
#these are from IGV tracks
#session saved here: https://www.dropbox.com/home/Ryan%20Mulqueen/navinlab/PAPERS/IN_PROGRESS/Methylation_DCIS/IGV_Tracks



####################################################
#           Supp Fig Celltype rna markers                #
###################################################

markers<-list()
markers[["fibroblast"]]<-c("DCN","LUM","COL1A1","COL1A2","APOD")
markers[["endothelial"]]<-c("PECAM1","CDH5","CD93","FLI1","BTNL9")
markers[["pericyte"]]<-c("RGS5","KCNJ8","PDGFRB","ACTA2","MCAM")
markers[["myeloid"]]<-c("AIF1","HLA-DRA","IL1B","MPO","LYZ")
markers[["bcell"]]<-c("CD19","MS4A1","CD22","CD79A","IGHM")
markers[["tcell"]]<-c("PTPRC","CD3E","CD4","CD6","IL7R")
markers[["basal"]]<-c("KRT14","KRT17","ACTA2-copy","KRT5","SFN")
markers[["lumsec"]]<-c("KIT","LTF","KRT15","MMP7","KLF5")
markers[["lumhr"]]<-c("AREG","AZGP1","KRT18","AGR2","ANKRD30A")
markers[["cancer"]]<-c("ESR1","PGR","ERBB2","GATA3","XBP1")

rna<-subset(rna,coarse_celltype %in% names(markers))
Idents(rna)<-factor(rna$coarse_celltype,levels=names(markers))

#add ACTA2 as second feature (ACTA2_copy) so i can use it twice
dat_insert <- t(as(as.matrix(rna@assays$RNA$data["ACTA2",]), "dgCMatrix"))
rownames(dat_insert) <- "ACTA2-copy"

counts_insert <- t(as(as.matrix(rna@assays$RNA$counts["ACTA2",]), "dgCMatrix"))
rownames(counts_insert) <- "ACTA2-copy"

new_assay <- CreateAssay5Object(counts = rbind2(rna@assays$RNA$counts, counts_insert), 
                                data =  rbind2(rna@assays$RNA$data, counts_insert))
rna[["RNA2"]]<-new_assay
plt<-DotPlot(rna, assay="RNA2",features = markers,cols=c("#c3c3c3","#e27428"))+  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(plt,file=paste0(outdir,"/supp_fig_celltyping_rna_markers.pdf"),height=30,width=10)



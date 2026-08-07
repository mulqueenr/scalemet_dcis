#Running scITD on both RNA and Methylation to assess factor weight across metadata

library(Seurat)
library(GenomicRanges)
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
library(scITD)

#read in object from directory
processing_folder="05_scITD"
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
wd=paste(sep="/",project_data_directory,processing_folder)
dir.create(wd)
setwd(wd)
obj<-readRDS(file=paste(project_data_directory,"03_fine_celltyping","03_scaledcis.final_celltypes.amethyst.rds",sep="/"))

#read in objects from fine cell typing step for now
rna<-readRDS(file=paste(project_data_directory,"03_fine_celltyping","03_rna.fine_celltyping.merged.rds",sep="/"))

#rna<-subset(rna,cells=row.names(rna@meta.data)[!is.na(rna$rna_clonename) & !endsWith(rna$rna_clonename,suffix="_diploid")])
rna_mat<-rna[["RNA"]]$counts

rna_meta<-rna@meta.data
rna_meta$ctypes<-rna_meta$coarse_celltype
rna_meta$donors<-rna_meta$sample

# set up project parameters
param_list <- initialize_params(ctypes_use = unique(rna_meta$ctypes),
                                ncores = 100, rand_seed = 123)

# create project container
rna_container <- make_new_container(count_data=rna_mat, 
                                     meta_data=rna_meta,
                                     params=param_list,
                                     metadata_cols=c("donors","ctypes",
                                     "Menopause","ER","PR","HER2",
                                     "Grade","Group","rna_clonename"),
                                     label_donor_sex = FALSE)

rna_container <- form_tensor(rna_container, donor_min_cells=0,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.1,
                              scale_var = TRUE, var_scale_power = 2)

print(length(rna_container[["all_vargenes"]]))

rna_container <- run_tucker_ica(rna_container, ranks=c(10,30),
                                 tucker_type = 'regular', rotation_type = 'hybrid')

# get donor scores-metadata associations
rna_container <- get_meta_associations(rna_container, vars_test=c('Group','PR','rna_clonename'), stat_use='pval')

# plot donor scores
rna_container <- plot_donor_matrix(rna_container, meta_vars=c('Group','PR','rna_clonename'),
                                    show_donor_ids = TRUE,
                                    add_meta_associations='pval')

pdf("test.factor.associations.pdf")
print(rna_container$plots$donor_matrix)
dev.off()

# get significant genes
rna_container <- get_lm_pvals(rna_container)

# generate the loadings plots
rna_container <- get_all_lds_factor_plots(rna_container, 
                                           use_sig_only=TRUE,
                                           nonsig_to_zero=TRUE,
                                           sig_thresh=.02,
                                           display_genes=FALSE,
                                           gene_callouts = TRUE,
                                           callout_n_gene_per_ctype=3,
                                           show_var_explained = TRUE)

# arrange the plots into a figure and show the figure
myfig <- render_multi_plots(rna_container,data_type='loadings')

pdf("test.factor.gene_loadings.pdf",width=20,height=20)
print(myfig)
dev.off()


# show the donor scores heatmap
pbmc_container$plots$donor_matrix
#not sure how to cut these windows, for now just test on rna
#met_mat<-obj@genomeMatrices
Integration of seurat object with methylation object.

1. Using cell type level DMRs 
2. assign DMRs to genes by overlap 
3. window per cell for those genes 
4. filter RNA to just those genes 
5. integrate with liger following:
https://welch-lab.github.io/liger/articles/rna_methylation.html
https://pmc.ncbi.nlm.nih.gov/articles/PMC8132955/

Use integration to finalize cell typing. Note some cell types were mislabelled based on small set of marker genes originally used. 
Using integration data for final classifications.

```R
BiocManager::install("HDF5Array")
install.packages('leidenAlg')
install.packages('rliger')
```


```R
#source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
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
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)

obj<-readRDS(file="07_scaledcis.integrated_celltyping.amethyst.rds")

#################################
#Get DMR per Integrated Celltypes
#################################
dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")
system(paste("mkdir -p", dmr_outdir))

#celltype
dmr_celltype_outdir=paste(sep="/",dmr_outdir,"integrated_celltype")
system(paste("mkdir -p", dmr_celltype_outdir))


celltype500bpwindows <- calcSmoothedWindows(obj, 
                                        type = "CG", 
                                        threads = 1,
                                        step = 500, # change to 500 for real data unless you have really low coverage
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "integrated_celltype",
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)
saveRDS(celltype500bpwindows,file=paste0(dmr_celltype_outdir,"/","dmr_analysis.integrated_celltype.500bp_windows.rds"))

pct_mat<-celltype500bpwindows[["pct_matrix"]] 
sum_mat<-celltype500bpwindows[["sum_matrix"]] 

#save clone object for future genome track plotting
obj@genomeMatrices[["cg_celltype_tracks"]] <- pct_mat #load it into amethyst object for plotting

dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
        eachVsAll = TRUE, # If TRUE, each group found in the sumMatrix will be tested against all others
        nminTotal = 3, # Min number observations across all groups to include the region in calculations
        nminGroup = 3) # Min number observations across either members or nonmembers to include the region

saveRDS(dmrs,file=paste0(dmr_celltype_outdir,"/","integrated_celltype_allcells",".dmr.rds"))

dmrs <- filterDMR(dmrs, 
            method = "bonferroni", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
            filter = FALSE, # If TRUE, removes insignificant results
            pThreshold = 0.05, # Maxmimum adjusted p value to allow if filter = TRUE
            logThreshold = 1.5) # Minimum absolute value of the log2FC to allow if filter = TRUE

collapsed_dmrs <- collapseDMR(obj, 
                        dmrs, 
                        maxDist = 2000, # Max allowable overlap between DMRs to be considered adjacent
                        minLength = 2000, # Min length of collapsed DMR window to include in the output
                        reduce = T, # Reduce results to unique observations (recommended)
                        annotate = T) # Add column with overlapping gene names

saveRDS(collapsed_dmrs,file=paste0(dmr_celltype_outdir,"/","integrated_celltype_allcells",".dmr_filt_collapse.rds"))

rename_dmr_output<-setNames(nm=1:length(unique(collapsed_dmrs$test)),gsub(colnames(sum_mat)[grepl(colnames(sum_mat),pattern="_t$")],pattern="_t",replacement=""))
collapsed_dmrs$type <- rename_dmr_output[collapsed_dmrs$test]
saveRDS(collapsed_dmrs,file=paste0(dmr_celltype_outdir,"/","integrated_celltype_allcells",".dmr_filt_collapse.rds"))

#plot dmr counts per clone
pal <- colorRampPalette(colors = c("#F9AB60", "#E7576E", "#630661", "#B5DCA5"))
COLS <- pal(length(unique(collapsed_dmrs$type)))

plt<-ggplot(collapsed_dmrs |> dplyr::group_by(type, direction) |> dplyr::summarise(n = n()), 
    aes(y = type, x = n, fill = type)) + 
    geom_col() + 
    facet_grid(vars(direction), scales = "free_y") + 
    scale_fill_manual(values = COLS) + 
    theme_classic()
ggsave(plt,file=paste0(dmr_celltype_outdir,"/","celltype_allcellsm",".dmr_counts.pdf"))


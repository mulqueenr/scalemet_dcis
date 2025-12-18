```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
source activate
#update data.table
```

<!-- 
```bash
export PATH="/home/rmulqueen/.local/bin:$PATH"
facet convert new_format.h5 old_format.h5
```
 -->
```R
#install.packages("data.table")
#ensure 0.9 version of amethyst
#system('wget https://github.com/lrylaarsdam/amethyst/releases/download/v0.0.0.9000/amethyst_0.0.0.9000.tar.gz')
#system('tar -xf amethyst_0.0.0.9000.tar.gz')
#devtools::install("./amethyst")
source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
set.seed(111)
library(amethyst)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)
obj<-readRDS(file="06_scaledcis.cnv_clones.amethyst.rds")



#update h5 files for cells
#h5_directory=paste0(project_data_directory,"/h5_files")
#system(paste("mkdir",h5_directory))
#for(i in unique(obj@h5paths$path)){
#   system(paste0("facet convert ",paste0(h5_directory,"/",basename(i))," ",i))
#}


```


Running through all the 500bp window generation split by:
- Celltype (fine_celltype)
- CNV based clones (cnv_clonenames)
- Celltype + Diagnosis (Group + fine_celltype)

Saving to folders for further processing.
For CNV based clones:
1. loop through cnv based clones samples to do DMR per sample per clone.
2. group together clones that have shared features (e.g. 1q amp) and compare to others
3. Plot DMR regions across genome per clone (see if clustered around major CNV changes or in trans)

For celltype by diagnosis
1. DMRs for TME cells per diagnosis (especially SYNCH vs non-SYNCH DCIS, and normal vs DCIS)

```R
#############################################
#Make 500bp windows once and slice for comparisons
#############################################
dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")
system(paste("mkdir -p", dmr_outdir))

#clones
dmr_clones_outdir=paste(sep="/",dmr_outdir,"cnv_clones")
system(paste("mkdir -p", dmr_clones_outdir))

clone500bpwindows <- amethyst::calcSmoothedWindows(obj, 
                                        type = "CG", 
                                        threads = 20,
                                        step = 500,
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "cnv_clonename",
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)
saveRDS(clone500bpwindows,file=paste0(dmr_clones_outdir,"/","dmr_analysis.cnv_clones.500bp_windows.rds"))

#celltype
dmr_celltype_outdir=paste(sep="/",dmr_outdir,"celltype")
system(paste("mkdir -p", dmr_celltype_outdir))
celltype500bpwindows <- calcSmoothedWindows(obj, 
                                        type = "CG", 
                                        threads = 1,
                                        step = 500, # change to 500 for real data unless you have really low coverage
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "fine_celltype",
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)
saveRDS(celltype500bpwindows,file=paste0(dmr_celltype_outdir,"/","dmr_analysis.celltype.500bp_windows.rds"))


#celltype by diagnosis
dmr_celltype_by_diag_outdir=paste(sep="/",dmr_outdir,"celltype_by_diagnosis")
system(paste("mkdir -p", dmr_celltype_by_diag_outdir))
obj@metadata$celltype_diag<-paste(sep="_",obj@metadata$Group,obj@metadata$fine_celltype)
celltypediag500bpwindows <- calcSmoothedWindows(obj, 
                                        type = "CG", 
                                        threads = 1,
                                        step = 500, # change to 500 for real data unless you have really low coverage
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "fine_celltype",
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)
saveRDS(celltypediag500bpwindows,file=paste0(dmr_celltype_by_diag_outdir,"/","dmr_analysis.celltype_by_diag.500bp_windows.rds"))

```


```R
#############################################
###### DMR Per CNV Clone Per Sample    ######
#############################################

#add plotting functions for DMR 
#https://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/pbmc_vignette/pbmc_vignette.html
#add dmr locations per clone over heatmap cnv plot

#output to CNV analysis folder
out_directory<-"/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/copykit/"

#run 500bp smoothed windows by clone split once, save and subset

dmr_per_clone<-function(obj,sample_name){
    obj_sample<-subsetObject(obj,cells=row.names(obj@metadata)[obj@metadata$Sample %in% c(sample_name)])
    if(length(unique(obj_sample@metadata$cnv_clones_split))>1){
        obj_sample<-subsetObject(obj_sample,cells=row.names(obj_sample@metadata)[!is.na(obj_sample@metadata$cnv_clones_split)])
        cluster500bpwindows <- calcSmoothedWindows(obj_sample, 
                                                type = "CG", 
                                                threads = 1,
                                                step = 500, # change to 500 for real data unless you have really low coverage
                                                smooth = 3,
                                                genome = "hg38",
                                                index = "chr_cg",
                                                groupBy = "cnv_clones_split",
                                                returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                                returnPctMatrix = TRUE)
        obj_sample@tracks[["cg_clone_tracks"]] <- cluster500bpwindows[["pct_matrix"]]
        saveRDS(obj_sample,file=paste0(out_directory,sample_name,"/",sample_name,".clones_amethyst.rds"))
        dmrs <- testDMR(cluster500bpwindows[["sum_matrix"]], # Sum of c and t observations in each genomic window per group
                eachVsAll = TRUE, # If TRUE, each group found in the sumMatrix will be tested against all others
                nminTotal = 5, # Min number observations across all groups to include the region in calculations
                nminGroup = 5) # Min number observations across either members or nonmembers to include the regio
        saveRDS(dmrs,file=paste0(out_directory,sample_name,"/",sample_name,".clones_dmr.rds"))
        dmrs <- filterDMR(dmrs, 
                  method = "bonferroni", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
                  filter = TRUE, # If TRUE, removes insignificant results
                  pThreshold = 0.01, # Maxmimum adjusted p value to allow if filter = TRUE
                  logThreshold = 1.5) # Minimum absolute value of the log2FC to allow if filter = TRUE
        collapsed_dmrs <- collapseDMR(obj_sample, 
                              dmrs, 
                              maxDist = 2000, # Max allowable overlap between DMRs to be considered adjacent
                              minLength = 2000, # Min length of collapsed DMR window to include in the output
                              reduce = T, # Reduce results to unique observations (recommended)
                              annotate = T) # Add column with overlapping gene names
        saveRDS(collapsed_dmrs,file=paste0(out_directory,sample_name,"/",sample_name,".clones_dmr_filt_collapse.rds"))
        rename_dmr_output<-setNames(nm=unique(collapsed_dmrs$test),unique(obj_sample@metadata$cnv_clones_split))
        collapsed_dmrs$type <- rename_dmr_output[collapsed_dmrs$test]
        #plot dmr counts per clone
        pal <- c("#F9AB60", "#E7576E", "#630661", "#B5DCA5") # makePalette(option = 7, n = 4) 
        plt<-ggplot(collapsed_dmrs |> dplyr::group_by(type, direction) |> dplyr::summarise(n = n()), 
            aes(y = type, x = n, fill = type)) + 
            geom_col() + 
            facet_grid(vars(direction), scales = "free_y") + 
            scale_fill_manual(values = pal) + 
            theme_classic()
        ggsave(plt,file=paste0(out_directory,sample_name,"/",sample_name,".clones_dmr_counts.pdf"))

    }else{
        print(paste(sample_name,"doesn't have more than one reported CNV clone."))
    }
}


dmr_per_clone(obj=obj,sample_name='BCMDCIS41T')


#run dmr analysis per clone


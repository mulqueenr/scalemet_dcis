#```bash
#singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
#source activate
#update data.table
#```

# To Do: Also try applying Bartlett's test for differential variability instead of just diff mean
see DNA methylation outliers in normal breast tissue  identify field defects that are enriched in cancer 2015

# Applying bartlett.test()
result = bartlett.test(len ~ interaction(supp, dose), 
                                  data = ToothGrowth)

```R
#install.packages("data.table")
#ensure 0.9 version of amethyst
#system('wget https://github.com/lrylaarsdam/amethyst/releases/download/v0.0.0.9000/amethyst_0.0.0.9000.tar.gz')
#system('tar -xf amethyst_0.0.0.9000.tar.gz')
#devtools::install("./amethyst")
#source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
#running locally

set.seed(111)
library(amethyst)
library(data.table)
library(dplyr)
library(msigdbr)
library(fgsea)
library(ggplot2)
library(data.table)

options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)
obj<-readRDS(file="06_scaledcis.cnv_clones.amethyst.rds")


```


Running through all the 500bp window generation split by:
- Celltype (fine_celltype)
- CNV based clones (cnv_clonenames)
- CNV based clones against all lumhr (cnv_clonename_alllumhr)
- Celltype + Diagnosis (Group + fine_celltype)

Saving to folders for further processing.
For CNV based clones:
1. loop through cnv based clones samples to do DMR per sample per clone.
2. group together clones that have shared features (e.g. 1q amp) and compare to others
3. Plot DMR regions across genome per clone (see if clustered around major CNV changes or in trans)

For celltype by diagnosis
1. DMRs for TME cells per diagnosis (especially SYNCH vs non-SYNCH DCIS, and normal vs DCIS)

```R
##################FUNCTIONS##################
testDMR <- function(sumMatrix, eachVsAll = TRUE, comparisons = NULL, nminTotal = 3,nminGroup = 3) {

  if (!eachVsAll && is.null(comparisons)) {
    stop("Please either specify eachVsAll = TRUE or provide a data frame of comparisons to make.")
  }

  # filter counts table
  data.table::setDT(sumMatrix)
  counts <- data.table::copy(sumMatrix)
  counts <- counts[rowSums(counts[, .SD, .SDcols = patterns("_c$|_t$")], na.rm = TRUE) >= nminTotal]

  # fast fisher's exact test developed by @zellerivo; see https://github.com/al2na/methylKit/issues/96
  fast.fisher <- function (
    cntg_table) {
    q <- cntg_table[1, 1]
    m <- cntg_table[1, 1] + cntg_table[2, 1]
    n <- cntg_table[1, 2] + cntg_table[2, 2]
    k <- cntg_table[1, 1] + cntg_table[1, 2]
    pval_right <- phyper(q = q, m = m, n = n, k = k, lower.tail = FALSE) +
      (0.5 * dhyper(q, m, n, k))
    pval_left <- phyper(q = q - 1, m = m, n = n, k = k, lower.tail = TRUE) +
      (0.5 * dhyper(q, m, n, k))
    return(ifelse(test = pval_right > pval_left, yes = pval_left *
                    2, no = pval_right * 2))
  }

  if (is.null(comparisons)) {
    # get unique groups
    groups <- as.list(sub("_c$", "", colnames(sumMatrix)[grep("_c$", colnames(sumMatrix))]))

    for (gr in groups) {
      m_c <- paste0(gr, "_c") # m = member
      m_t <- paste0(gr, "_t")

      nm_c <- setdiff(grep("_c$", colnames(counts), value = TRUE), m_c) # nm = nonmember
      nm_t <- setdiff(grep("_t$", colnames(counts), value = TRUE), m_t)

      counts <- counts[, `:=`(
        member_c = get(paste0(gr, "_c")),
        member_t = get(paste0(gr, "_t")),
        nonmember_c = rowSums(.SD[, mget(nm_c)]),
        nonmember_t = rowSums(.SD[, mget(nm_t)])
      )]

      # don't test where the minimum observations per group is not met
      counts <- counts[member_c + member_t <= nminGroup | nonmember_c + nonmember_t <= nminGroup, c("member_c", "member_t", "nonmember_c", "nonmember_t") := .(NA, NA, NA, NA)]

      # apply fast fishers exact test
      counts <- counts[, paste0(gr, "_all_pval") := apply(.SD, 1, function(x) fast.fisher(matrix(x, nrow = 2, byrow = TRUE))), .SDcols = c("member_c", "member_t", "nonmember_c", "nonmember_t")]
      counts <- counts[, paste0(gr, "_all_logFC") := round(log2((member_c / (member_c + member_t)) / (nonmember_c / (nonmember_c + nonmember_t))), 4)]
      counts <- counts[, c("member_c", "member_t", "nonmember_c", "nonmember_t") := NULL] # this line used to be outside the loop in < v1.0.2, causing nonmember variable buildup :(
      cat(paste0("Finished group ", gr, "\n"))
    }

  } else if (!is.null(comparisons)) {
    for (i in 1:nrow(comparisons)) {
      m <- unlist(strsplit(comparisons[i, "A"], ','))
      nm <- unlist(strsplit(comparisons[i, "B"], ',', fixed = FALSE))
      name <- comparisons[i, "name"]

      m_c <- paste0(m, "_c") # m = member
      m_t <- paste0(m, "_t")

      nm_c <- paste0(nm, "_c") # n = nonmember
      nm_t <- paste0(nm, "_t")

      counts <- counts[, `:=`(
        member_c = rowSums(.SD[, mget(m_c)]),
        member_t = rowSums(.SD[, mget(m_t)]),
        nonmember_c = rowSums(.SD[, mget(nm_c)]),
        nonmember_t = rowSums(.SD[, mget(nm_t)])
      )]

      # don't test where the minimum observations per group is not met
      counts <- counts[member_c + member_t <= nminGroup | nonmember_c + nonmember_t <= nminGroup, c("member_c", "member_t", "nonmember_c", "nonmember_t") := .(NA, NA, NA, NA)]

      # apply fast fishers exact test
      counts <- counts[, paste0(name, "_pval") := apply(.SD, 1, function(x) fast.fisher(matrix(x, nrow = 2, byrow = TRUE))), .SDcols = c("member_c", "member_t", "nonmember_c", "nonmember_t")]
      counts <- counts[, paste0(name, "_logFC") := round(log2((member_c / (member_c + member_t)) / (nonmember_c / (nonmember_c + nonmember_t))), 4)]
      counts <- counts[, c("member_c", "member_t", "nonmember_c", "nonmember_t") := NULL]
      cat(paste0("\nFinished testing ", name, ": ", paste0(m, collapse = ", "), " vs. ", paste0(nm, collapse = ", ")))
    }
  }
  return(counts)
}

dmr_celltype_by_all<-function(obj){
    #read 500bp windows
    celltype500bpwindows<-readRDS(paste0(dmr_celltype_outdir,"/","dmr_analysis.celltype.500bp_windows.rds"))
    
    pct_mat<-celltype500bpwindows[["pct_matrix"]] 
    sum_mat<-celltype500bpwindows[["sum_matrix"]] 

    #save clone object for future genome track plotting
    obj@genomeMatrices[["cg_celltype_tracks"]] <- pct_mat #load it into amethyst object for plotting

    dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
            eachVsAll = TRUE, # If TRUE, each group found in the sumMatrix will be tested against all others
            nminTotal = 3, # Min number observations across all groups to include the region in calculations
            nminGroup = 3) # Min number observations across either members or nonmembers to include the region
    saveRDS(dmrs,file=paste0(dmr_celltype_outdir,"/","celltype_allcells",".dmr.rds"))
   
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
    saveRDS(collapsed_dmrs,file=paste0(dmr_celltype_outdir,"/","celltype_allcells",".dmr_filt_collapse.rds"))

    rename_dmr_output<-setNames(nm=1:length(unique(collapsed_dmrs$test)),gsub(colnames(sum_mat)[grepl(colnames(sum_mat),pattern="_t$")],pattern="_t",replacement=""))
    collapsed_dmrs$type <- rename_dmr_output[collapsed_dmrs$test]
    saveRDS(collapsed_dmrs,file=paste0(dmr_celltype_outdir,"/","celltype_allcells",".dmr_filt_collapse.rds"))
    
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
}
dmr_clone_by_clone<-function(obj,sample_name){
    #set dir
    clone_outdir=paste0(dmr_clones_outdir,"/",sample_name)
    system(paste("mkdir -p",clone_outdir))
    
    print(paste("Subset object for sample:",sample_name))
    
    #subset object (sample name match AND fine cell type is lumhr or cancer (by original 500bp summaries above))
    cells_in=row.names(obj@metadata)[(obj@metadata$Sample %in% c(sample_name))]
    obj_sample<-subsetObject(obj,cells=cells_in)
    #subset cell types to just epithelial/cancer

    #read 500bp windows
    clone500bpwindows<-readRDS(paste0(dmr_clones_outdir,"/","dmr_analysis.cnv_clones.500bp_windows.rds"))
    
    #subset to just sample clones
    print(paste("Subset 500bp windows for sample:",sample_name))

    pct_columns_to_keep=c(1,2,3,which(grepl(pattern=sample_name,unlist(lapply(strsplit(colnames(clone500bpwindows[["pct_matrix"]]),"_"),"[",1)))))
    sum_columns_to_keep=c(1,2,3,which(grepl(pattern=sample_name,unlist(lapply(strsplit(colnames(clone500bpwindows[["sum_matrix"]]),"_"),"[",1)))))
    
    pct_mat<-clone500bpwindows[["pct_matrix"]] %>% select(all_of(pct_columns_to_keep))
    sum_mat<-clone500bpwindows[["sum_matrix"]] %>% select(all_of(sum_columns_to_keep))

    print(paste("Saving amethyst object for clones:",sample_name))
    #save clone object for future genome track plotting
    obj_sample@genomeMatrices[["cg_clone_tracks"]] <- pct_mat #load it into amethyst object for plotting
    saveRDS(obj_sample,file=paste0(clone_outdir,"/",sample_name,".clones_amethyst.rds"))

    print(paste("Running DMRs across clones for:",sample_name))
    dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
            eachVsAll = TRUE, # If TRUE, each group found in the sumMatrix will be tested against all others
            nminTotal = 3, # Min number observations across all groups to include the region in calculations
            nminGroup = 3) # Min number observations across either members or nonmembers to include the region
    saveRDS(dmrs,file=paste0(clone_outdir,"/",sample_name,".clones_dmr.rds"))
   
    print(paste("Filtering DMRs across clones for:",sample_name))
    dmrs <- filterDMR(dmrs, 
                method = "bonferroni", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
                filter = FALSE, # If TRUE, removes insignificant results
                pThreshold = 0.05, # Maxmimum adjusted p value to allow if filter = TRUE
                logThreshold = 1.5) # Minimum absolute value of the log2FC to allow if filter = TRUE
    collapsed_dmrs <- collapseDMR(obj_sample, 
                            dmrs, 
                            maxDist = 2000, # Max allowable overlap between DMRs to be considered adjacent
                            minLength = 2000, # Min length of collapsed DMR window to include in the output
                            reduce = T, # Reduce results to unique observations (recommended)
                            annotate = T) # Add column with overlapping gene names
    saveRDS(collapsed_dmrs,file=paste0(clone_outdir,"/",sample_name,".clones_dmr_filt_collapse.rds"))

    rename_dmr_output<-setNames(nm=1:length(unique(collapsed_dmrs$test)),gsub(colnames(sum_mat)[grepl(colnames(sum_mat),pattern="_t$")],pattern="_t",replacement=""))
    collapsed_dmrs$type <- rename_dmr_output[collapsed_dmrs$test]
    saveRDS(collapsed_dmrs,file=paste0(clone_outdir,"/",sample_name,".clones_dmr_filt_collapse.rds"))
    
    #plot dmr counts per clone
    pal <- colorRampPalette(colors = c("#F9AB60", "#E7576E", "#630661", "#B5DCA5"))
    COLS <- pal(length(unique(collapsed_dmrs$type)))

    plt<-ggplot(collapsed_dmrs |> dplyr::group_by(type, direction) |> dplyr::summarise(n = n()), 
        aes(y = type, x = n, fill = type)) + 
        geom_col() + 
        facet_grid(vars(direction), scales = "free_y") + 
        scale_fill_manual(values = COLS) + 
        theme_classic()
    ggsave(plt,file=paste0(clone_outdir,"/",sample_name,".clones_dmr_counts.pdf"))
}

dmr_clone_by_lumhr<-function(obj,sample_name){
    #set dir
    clone_outdir=paste0(dmr_clones_outdir,"/",sample_name)
    system(paste("mkdir -p",clone_outdir))
    
    print(paste("Subset object for sample:",sample_name))
    
    #subset object to all lumhr OR sample specific clones
    cells=row.names(obj@metadata)[obj@metadata$Sample %in% c(sample_name) & obj@metadata$cnv_clonename != "NA"]
    lumhr=row.names(obj@metadata)[obj@metadata$fine_celltype %in% c("lumhr")]
    cells<-unique(c(cells,lumhr))
    obj_sample<-subsetObject(obj,cells=cells)

    #read 500bp windows
    clone500bpwindows<-readRDS(paste0(dmr_clones_outdir,"/","dmr_analysis.cnv_clones_alllumhr.500bp_windows.rds"))
    
    #subset to just sample clones
    print(paste("Subset 500bp windows for sample:",sample_name))
    pct_columns_to_keep=c(1,2,3,which(grepl(pattern=sample_name,unlist(lapply(strsplit(colnames(clone500bpwindows[["pct_matrix"]]),"_"),"[",1)))))
    sum_columns_to_keep=c(1,2,3,4,5,which(grepl(pattern=sample_name,unlist(lapply(strsplit(colnames(clone500bpwindows[["sum_matrix"]]),"_"),"[",1)))))
    
    pct_mat<-clone500bpwindows[["pct_matrix"]] %>% select(all_of(pct_columns_to_keep))
    sum_mat<-clone500bpwindows[["sum_matrix"]] %>% select(all_of(sum_columns_to_keep))

    print(paste("Saving amethyst object for clones:",sample_name))
    #save clone object for future genome track plotting
    obj_sample@genomeMatrices[["cg_clone_tracks"]] <- pct_mat #load it into amethyst object for plotting
    saveRDS(obj_sample,file=paste0(clone_outdir,"/",sample_name,".clones_lumhrall_amethyst.rds"))

    #run tests of clones against lumhr, and all clones v lumhr
    comparisons=colnames(pct_mat)[5:ncol(pct_mat)]
    comparisons<-comparisons[!grepl(comparisons,pattern="_diploid")]
    comparisons<-comparisons[!is.na(comparisons)]

    compare_tests=data.frame(
      name=c(paste0(comparisons,"_lumhr"),paste0(sample_name,"_allclones_lumhr")),
      A=c(comparisons,paste(comparisons,collapse=",")),
      B=rep("lumhr",length(comparisons)+1))

    print(paste("Running DMRs across clones for:",sample_name))
    dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
            comparisons=compare_tests,
            nminTotal = 3, # Min number observations across all groups to include the region in calculations
            nminGroup = 3) # Min number observations across either members or nonmembers to include the region
    saveRDS(dmrs,file=paste0(clone_outdir,"/",sample_name,".clones_lumhrall_dmr.rds"))
   
    print(paste("Filtering DMRs across clones for:",sample_name))
    dmrs <- filterDMR(dmrs, 
                method = "bonferroni", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
                filter = FALSE, # If TRUE, removes insignificant results
                pThreshold = 0.05, # Maxmimum adjusted p value to allow if filter = TRUE
                logThreshold = 1.5) # Minimum absolute value of the log2FC to allow if filter = TRUE
    collapsed_dmrs <- collapseDMR(obj_sample, 
                            dmrs, 
                            maxDist = 2000, # Max allowable overlap between DMRs to be considered adjacent
                            minLength = 2000, # Min length of collapsed DMR window to include in the output
                            reduce = T, # Reduce results to unique observations (recommended)
                            annotate = T) # Add column with overlapping gene names
    saveRDS(collapsed_dmrs,file=paste0(clone_outdir,"/",sample_name,".clones_lumhrall_dmr_filt_collapse.rds"))

    rename_dmr_output<-setNames(nm=1:length(unique(collapsed_dmrs$test)),compare_tests$name)
    collapsed_dmrs$type <- rename_dmr_output[collapsed_dmrs$test]
    saveRDS(collapsed_dmrs,file=paste0(clone_outdir,"/",sample_name,".clones_lumhrall_dmr_filt_collapse.rds"))
    
    #plot dmr counts per clone
    pal <- colorRampPalette(colors = c("#F9AB60", "#E7576E", "#630661", "#B5DCA5"))
    COLS <- pal(length(unique(collapsed_dmrs$type)))

    plt<-ggplot(collapsed_dmrs |> dplyr::group_by(type, direction) |> dplyr::summarise(n = n()), 
        aes(y = type, x = n, fill = type)) + 
        geom_col() + 
        facet_grid(vars(direction), scales = "free_y") + 
        scale_fill_manual(values = COLS) + 
        theme_classic()
    ggsave(plt,file=paste0(clone_outdir,"/",sample_name,".clones_lumhrall_dmr_counts.pdf"))
}

```

```R
#############################################
#Make 500bp windows once and slice for comparisons
#############################################
dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")
system(paste("mkdir -p", dmr_outdir))

#clones (limited to lumhr and cancer)
dmr_clones_outdir=paste(sep="/",dmr_outdir,"cnv_clones")
system(paste("mkdir -p", dmr_clones_outdir))
obj_lumhr<-subsetObject(obj,cells=row.names(obj@metadata)[(obj@metadata$fine_celltype %in% c("lumhr","cancer"))])

clone500bpwindows <- amethyst::calcSmoothedWindows(obj_lumhr, 
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

#clones (all lumhr treated as same group)
dmr_clones_outdir=paste(sep="/",dmr_outdir,"cnv_clones_alllumhr")
system(paste("mkdir -p", dmr_clones_outdir))
obj_lumhr_all<-subsetObject(obj,cells=row.names(obj@metadata)[(obj@metadata$fine_celltype %in% c("lumhr","cancer"))])
obj_lumhr_all@metadata$cnv_clonename_alllumhr<-obj_lumhr_all@metadata$cnv_clonename
obj_lumhr_all@metadata[obj_lumhr_all@metadata$fine_celltype=="lumhr",]$cnv_clonename_alllumhr<-"lumhr"
clone500bpwindows <- amethyst::calcSmoothedWindows(obj_lumhr_all, 
                                        type = "CG", 
                                        threads = 20,
                                        step = 500,
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "cnv_clonename_alllumhr",
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)
saveRDS(clone500bpwindows,file=paste0(dmr_clones_outdir,"/","dmr_analysis.cnv_clones_alllumhr.500bp_windows.rds"))

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

#output to DMR analysis folder
dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")

#run 500bp smoothed windows by clone split once, save and subset

#samples with >= 2 clones and >=5% coverage per clone
#use dmr_analysis.cnv_clones.500bp_windows.rds matrix
clone_v_clone_analysis<-c('BCMDCIS05T',
'BCMDCIS124T',
'BCMDCIS28T',
'BCMDCIS35T',
'BCMDCIS41T',
'BCMDCIS52T',
'BCMDCIS66T',
'BCMDCIS74T',
'BCMDCIS79T',
'BCMDCIS82T',
'BCMDCIS97T',
'ECIS25T',
'ECIS26T',
'ECIS36T',
'ECIS48T')



lapply(clone_v_clone_analysis,function(x) dmr_clone_by_clone(obj,sample_name=x))
```

```R
#add plotting functions for DMR 
#https://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/pbmc_vignette/pbmc_vignette.html
#add dmr locations per clone over heatmap cnv plot

#output to DMR analysis folder
dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")

#run dmr analysis clone with all lumhr
#use dmr_analysis.cnv_clones_alllumhr.500bp_windows.rds matrix
#note some of these also have the coverage to do specific clone vs lumhr analysis, can merge matrices for this case
clone_v_lumhr_analysis<-c(
'BCMHBCA04R',
'BCMHBCA03R',
'BCMHBCA83L-3h',
'BCMDCIS05T',
'BCMDCIS82T',
'BCMDCIS99T',
'BCMHBCA26L',
'BCMDCIS22T',
'BCMDCIS79T',
'BCMDCIS28T',
'BCMDCIS52T',
'BCMDCIS80T',
'ECIS48T',
'BCMDCIS102T',
'BCMDCIS92T',
'BCMDCIS70T',
'BCMDCIS124T',
'ECIS36T',
'BCMDCIS74T',
'BCMDCIS97T',
'BCMDCIS94T',
'BCMDCIS41T',
'ECIS25T',
'ECIS26T',
'BCMDCIS66T',
'BCMDCIS35T'
)


lapply(clone_v_lumhr_analysis,function(x) dmr_clone_by_lumhr(obj,sample_name=x))


```

Run DMR results through gene ontology for interpretation (to be coded)


```R



#gsea of top DMR genes
gsea_enrichment_againstlumhr<-function(dmrs,species="human",
                          category="C3",
                          subcategory="TFT:GTRD",
                          out_setname="TFT",sample_name=sample_name,
                          outdir=outdir,
                          obj=obj){

  pathwaysDF <- msigdbr(species=species, 
                        collection=category, 
                        subcollection = subcategory)

  #limit pathways to genes in our data
  pathwaysDF<-pathwaysDF[pathwaysDF$ensembl_gene %in% unlist(lapply(strsplit(unique(obj@ref$gene_id),"[.]"),"[",1)),]
  
  pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

  #run plotting per group and per hyper and hypo? i think just using logFC is sufficient

  plt_list<-lapply(unique(dmrs$type),function(group1){
  #treat multiple gene overlaps as same logFC
  #set -Inf to -3 and Inf to 3
  group_features<-dmrs %>%
    dplyr::filter(type == group1) %>%
    dplyr::filter(dmr_padj<0.05) %>% 
    dplyr::filter(gene_names!="NA") %>% 
    dplyr::arrange(dmr_logFC) %>%
    dplyr::select(gene_names, dmr_logFC) %>%
    tidyr::separate_rows(gene_names) %>%
    dplyr::mutate(across(where(is.numeric), ~ replace(., .==-Inf, -3))) %>%
    dplyr::mutate(across(where(is.numeric), ~ replace(., .==Inf, 3)))

  ranks<-setNames(nm=group_features$gene_names,group_features$dmr_logFC)
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = ranks,
                    minSize  = 5,
                    nproc = 1)

  topPathwaysUp <- fgseaRes %>% filter(ES > 0) %>% slice_max(NES,n=10) %>% dplyr::select(pathway)
  topPathwaysDown <- fgseaRes %>% filter(ES < 0) %>% slice_max(abs(NES),n=10) %>% dplyr::select(pathway)
  topPathways <- unlist(c(topPathwaysUp, rev(topPathwaysDown)))
  saveRDS(fgseaRes,file=paste0(outdir,sample_name,".",out_setname,".GSEA_enrichment.clone_lumhrall.rds"))

  #not returned
  plt1<-plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam=0.5)+
        theme(axis.text.y = element_text( size = rel(0.2)),
        axis.text.x = element_text( size = rel(0.2)))

  # only plot the top 20 pathways NES scores
  nes_plt_dat<-rbind(
    fgseaRes  %>% slice_max(NES,n= 10),
    fgseaRes  %>% slice_min(NES,n= 10))

  plt2<-ggplot(nes_plt_dat, aes(reorder(pathway, NES), NES)) +
    geom_point(aes(color= padj,size=size)) + scale_color_gradient2(low="darkred",high="grey",mid="red",midpoint=0.05)+
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score \n Compared to normal LumHR \n (Hyper<->Hypo)",
        title="Hallmark pathways NES from GSEA") + 
    theme_minimal()+scale_fill_identity()+ggtitle(paste(group1,out_setname))+ylim(c(5,-5))
  return(plt2)
  })
  patchwork::wrap_plots(plt_list,ncol=1)

}

plot_gsea<-function(obj,
                    sample_name,
                    outdir){
  print(paste("Loading DMRS for sample:",sample_name))

  collapsed_dmrs<-readRDS(file=paste0(dmr_outdir,"/cnv_clones_alllumhr/",sample_name,"/",sample_name,".clones_lumhrall_dmr_filt_collapse.rds"))

  #run gsea enrichment on different sets
  print(paste("Calculating TF Binding Enrichment Compared to LumHR"))
  tft_plt<-gsea_enrichment_againstlumhr(species="human",
              category="C3",
              subcategory="TFT:GTRD",
              out_setname="TFT",sample_name=sample_name,
              dmrs=collapsed_dmrs,obj=obj,outdir=outdir)

  print(paste("Calculating Position Enrichment Compared to LumHR"))
  position_plt<-gsea_enrichment_againstlumhr(species="human",
              category="C1",
              subcategory=NULL,
              out_setname="position",sample_name=sample_name,
              dmrs=collapsed_dmrs,obj=obj,outdir=outdir)

  print(paste("Calculating Hallmark Enrichment Compared to LumHR"))
  hallmark_plt<-gsea_enrichment_againstlumhr(species="human",
              category="H",
              subcategory=NULL,
              out_setname="hallmark",sample_name=sample_name,
              dmrs=collapsed_dmrs,obj=obj,outdir=outdir)

  print(paste("Calculating Cancer Cell Atlas Enrichment Compared to LumHR"))
  hallmark_plt<-gsea_enrichment_againstlumhr(species="human",
              category="C4",
              subcategory="3CA",
              out_setname="3CA",sample_name=sample_name,
              dmrs=collapsed_dmrs,obj=obj,outdir=outdir)
  plt<-patchwork::wrap_plots(list((tft_plt),(position_plt),(hallmark_plt)),ncol=4,axes="collect")

  ggsave(plt,
          file=paste0(sample_name,".clones_lumhrall.GSEA_enrichment.pdf"),
          path=outdir,
          width=40,height=length(unique(collapsed_dmrs$type))*5,limitsize = FALSE)
  print(paste("Finished sample:",sample_name))
}

lapply(clone_v_lumhr_analysis, function(sample_name){
  plot_gsea(obj,
          sample_name=sample_name,
          outdir=paste0(dmr_outdir,"/cnv_clones_alllumhr/",sample_name,"/"))
})


```

#coregulation by gsea https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/geseca-tutorial.html


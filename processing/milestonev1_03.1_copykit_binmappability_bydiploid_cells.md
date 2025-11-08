```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
```
Use Sherman to simulate reads
https://github.com/FelixKrueger/Sherman
Use read simulations to define the 220kb bin sizes
Use diploid cells to blacklist regions beyond 3 SD from mean


Rework code below.


## Finding expected reads per diploid bin

Running copykit read counting diploid cells, then applying mappability correction for cnv clones.

```R
source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
library(Rsamtools)
library(copykit)
library(circlize)
detach("package:GeneNMF",unload=TRUE)
library(dendextend)
library(RColorBrewer)
library(ComplexHeatmap)

#download cytoband data
#system("wget -P /data/rmulqueen/projects/scalebio_dcis/ref https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz ")
#system("gzip -d /data/rmulqueen/projects/scalebio_dcis/ref/cytoBand.txt.gz")
cyto=read.table(file="/data/rmulqueen/projects/scalebio_dcis/ref/cytoBand.txt",sep="\t")
colnames(cyto)<-c("chr","start","end","band","stain")
table(cyto$stain) #set colors for these

#set colors
celltype_col=c(
'peri'='#c1d552',
'fibro1'='#7f1911',
'fibro2'='#e791f9',
'endo'='#f0b243',
'endo2'='#d0bd4a',
'tcell'='#2e3fa3',
'bcell'='#00adea',
'myeloid1'='#00a487',
'myeloid2'='#006455',
'basal'='#7200cc',
'lumsec'='#af00af',
'lumhr'='#d8007c')

#set environment and read in data
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)
obj<-readRDS(file="05_scaledcis.fine_celltype.amethyst.rds")

#make output directory
system(paste0("mkdir -p ",project_data_directory,"/copykit" ))

read_scalebio_bam<-function(obj_met,x,sample_name){
    #scalebio pipeline outputs bam files as Tn5 wells. so multiple cell IDs are in a bam. this function splits out the bam to the query cellid
    bam=obj_met$bam_path[x]
    cellid=strsplit(row.names(obj_met)[x],"[+]batch|[+]prelim")[[1]][1]

    what <- c("qname","rname", "pos")
    param <- ScanBamParam(what=what,
                            flag=scanBamFlag(isPaired=TRUE,
                                            isProperPair=TRUE,
                                            isSecondaryAlignment=FALSE,
                                            isDuplicate=FALSE,
                                            isSupplementaryAlignment=FALSE))

    input_bam<-Rsamtools::scanBam(bam,param=param)
    input_bam<-do.call("DataFrame", input_bam)
    input_bam$cellid<-gsub("^.*:", "", input_bam$qname)
    input_bam<-input_bam[input_bam$cellid==cellid,]
    input_bam$end<-input_bam$pos+1
    input_bam<-makeGRangesFromDataFrame(input_bam,seqnames.field="rname",start.field="pos",end.field="end")
    return(input_bam)
}

output_directory<-paste0(project_data_directory,"/copykit")
remove_Y = TRUE
min_bincount = 10
cores=100
genome = "hg38"
resolution="220kb"

# bindings for NSE and data
Chr <- chr <- strand <- GeneID <- NULL
reads_assigned_bins <- reads_duplicates <- reads_total <- NULL

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# genomic ranges (varbin scaffolds)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Reading hg38 VarBin ranges
hg38_grangeslist <- hg38_grangeslist

hg38_rg <- switch(resolution,
    "55kb" = hg38_grangeslist[["hg38_50kb"]],
    "110kb" = hg38_grangeslist[["hg38_100kb"]],
    "195kb" = hg38_grangeslist[["hg38_175kb"]],
    "220kb" = hg38_grangeslist[["hg38_200kb"]],
    "280kb" = hg38_grangeslist[["hg38_250kb"]],
    "500kb" = hg38_grangeslist[["hg38_500kb"]],
    "1Mb" = hg38_grangeslist[["hg38_1Mb"]],
    "2.8Mb" = hg38_grangeslist[["hg38_2Mb"]]
)

hg38_rg <- as.data.frame(hg38_rg)

rg <- hg38_rg %>%
    dplyr::rename(chr = "seqnames") %>%
    dplyr::mutate(GeneID = 1:nrow(hg38_rg))

if (remove_Y == TRUE) {
    rg <- dplyr::filter(rg,chr != "chrY")
}

message("Counting reads for genome ",genome," and resolution: ",resolution)

#filter to exclude any potential cancer cells (removing lumhr cells from DCIS and IDC)
#get list of bams and cellids
diploid_cells<-rbind(
    obj@metadata %>% filter(fine_celltype!="lumhr"),
    obj@metadata %>% filter(Group=="HBCA") %>% filter(fine_celltype=="lumhr"))

#return chr start position for reads filtered in bam to cell id
varbin_counts_list_all_fields<-mclapply(
                                    1:nrow(diploid_cells), 
                                    function(i) 
                                    read_scalebio_bam(obj_met=diploid_cells,x=i,sample_name=sample_name), 
                                    mc.cores=cores)

message("Read in all bam files.")

names(varbin_counts_list_all_fields)<- row.names(obj_met)
varbin_counts_list_all_fields<-as(varbin_counts_list_all_fields, "GRangesList")
ref<-as(rg,"GRanges")

varbin_counts_list <-mclapply(varbin_counts_list_all_fields,
                                function(x) 
                                GenomicRanges::countOverlaps(
                                query=ref,
                                subject=x,
                                type="any",
                                ignore.strand=TRUE))
message("Counted reads across all bins.")

varbin_counts_list <- lapply(varbin_counts_list,as.vector)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# filtering for minimal mean bin count
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# obtaining the index of the ones that FAIL to meet the min_bincount arg
min_bc <- which(vapply(varbin_counts_list, mean, numeric(1)) < min_bincount)
# subsetting counts list and main counts list

if (length(min_bc) > 0) {
    varbin_counts_list <- varbin_counts_list[-min_bc]
    varbin_counts_list_all_fields <- varbin_counts_list_all_fields[-min_bc]
    message(
        length(min_bc), " bam files had less than ", min_bincount,
        " mean bincounts and were removed."
    )
}

# LOWESS GC normalization

message("Performing GC correction.")

varbin_counts_list_gccor <-
    mclapply(varbin_counts_list, function(x) {
        gc_cor <- lowess(rg$gc_content, log(x + 1e-3), f = 0.05)
        gc_cor_z <- approx(gc_cor$x, gc_cor$y, rg$gc_content)
        exp(log(x) - gc_cor_z$y) * median(x)
    },mc.cores=cores
    )

varbin_counts_df <- round(dplyr::bind_cols(varbin_counts_list_gccor), 2)

# filtering low read counts where the sum of bins does not reach more than 0
good_cells <- names(varbin_counts_df[which(colSums(varbin_counts_df) != 0)])

varbin_counts_df <- varbin_counts_df[good_cells]

#add additional column to rg about mean coverage per window for diploid cells
met_mappability<-rowMeans(varbin_counts_df)
rg$met_mappability<-met_mappability

varbin_counts_df_map <-
    mclapply(1:ncol(varbin_counts_df), function(x) {
        met_mat <- lowess(rg$met_mappability, log(unlist(varbin_counts_df[,x]) + 1e-3), f = 0.05)
        map_cor_z <- approx(met_mat$x, met_mat$y, rg$met_mappability)
        exp(log(varbin_counts_df[,x]) - map_cor_z$y) * median(unlist(varbin_counts_df[,x]))
    },mc.cores=cores
    )

varbin_counts_df2 <- round(dplyr::bind_cols(varbin_counts_df_map), 2)

rg <- rg %>%
    dplyr::select(-strand, -GeneID)

rg_gr <- GenomicRanges::makeGRangesFromDataFrame(rg,
    ignore.strand = TRUE,
    keep.extra.columns = TRUE)

cna_obj <- CopyKit(
    assays = list(bincounts = varbin_counts_df2),
    rowRanges = rg_gr)

# Adding genome and resolution information to metadata
S4Vectors::metadata(cna_obj)$genome <- genome
S4Vectors::metadata(cna_obj)$resolution <- resolution

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# ADDING READS METRICS TO METADATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

# saving info and removing columns from list elements
bam_metrics <- obj_met[c("unique_reads","tss_enrich","mcg_pct","cg_cov","batch","plate_info","tgmt_well","i7_well","i5_well","fine_celltype")]

# making sure metrics match varbin_counts_df
bam_metrics <- bam_metrics[good_cells,]

bam_metrics$sample <- rownames(bam_metrics)
bam_metrics$sample_name="all_diploid_cells"
bam_metrics$reads_assigned_bins <- colSums(varbin_counts_df2)

# adding to metadata
SummarizedExperiment::colData(cna_obj) <-
    S4Vectors::DataFrame(bam_metrics)
colnames(cna_obj) <- names(varbin_counts_df2)

#runvarbin module
cna_obj <- runVst(cna_obj)
cna_obj <- runSegmentation(cna_obj)
cna_obj <- logNorm(cna_obj)

    # Mark euploid cells if they exist
    cna_obj <- findAneuploidCells(cna_obj)

    # Mark low-quality cells for filtering
    cna_obj <- findOutliers(cna_obj)

    # kNN smooth profiles
    cna_obj <- knnSmooth(cna_obj,k=10)

    # adds basic quality control information to colData
    cna_obj <- runMetrics(cna_obj)
    cna_obj <- runUmap(cna_obj)
    cna_obj <- findSuggestedK(cna_obj)
    S4Vectors::metadata(cna_obj)$suggestedK

    #define colors based on data
    log_col=colorRamp2(c(-2,-1,0,1,2),
                            c("darkblue","blue","white","red","darkred"))
    cg_perc_col=colorRamp2(c(40,60,80,100),
                            c("#4d2d18","#CABA9C","#4C6444","#102820"))
    reads_col=colorRamp2(c(min(log10(cna_obj@colData$unique_reads)),
                            max(log10(cna_obj@colData$unique_reads))),
                            c("white","black"))
    
    #plot heatmap
    ha = rowAnnotation(
        reads=log10(cna_obj@colData$unique_reads),
        cg_perc=cna_obj@colData$mcg_pct,
        celltype=cna_obj@colData$fine_celltype,
        col= list(
            celltype=celltype_col,
            reads=reads_col,
            cg_perc=cg_perc_col
        ))

    cyto_overlap<-GenomicRanges::findOverlaps(cna_obj@rowRanges,
                                                makeGRangesFromDataFrame(cyto,keep=TRUE),
                                                select="first")
    cna_obj@rowRanges$stain <- cyto[cyto_overlap,]$stain

    arm_col=c("p"="grey","q"="darkgrey")
    band_col=c("acen"="#99746F","gneg"="white","gpos100"="black","gpos25"="lightgrey","gpos50"="grey","gpos75"="darkgrey","gvar"="#446879")

    column_ha = HeatmapAnnotation(
        arm = cna_obj@rowRanges$arm,
        band = cna_obj@rowRanges$stain,
        col=list(arm=arm_col,band=band_col))

    plt<-Heatmap(
        t(cna_obj@assays@data$logr),
        left_annotation=ha,col=log_col,
        show_column_names=FALSE,show_row_names=FALSE,
        top_annotation=column_ha,cluster_columns=FALSE,cluster_column_slices=FALSE,column_split=seqnames(cna_obj@rowRanges),
        name="logr")

    pdf(paste0(output_directory,"/copykit.",sample_name,".",resolution,".pdf"),width=20)
    print(plt)
    dev.off()

    saveRDS(cna_obj,file=paste0(output_directory,"/copykit",".",sample_name[1],".",resolution,".rds"))

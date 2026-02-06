
# Finding expected reads per diploid bin

From running the standard bins for copykit (both 220kb and 500kb) it is clear that some bins have mappability issues due to the bisulfite conversion. Using diploid cells (HBCA samples and non epithelial DCIS and IDC) to calculate coverage across bins.

Tried multiple ways to account for this, but finally just ended up doing 500kb bins (less cell dropout due to lower read counts) and remove bins with excessively high or low mapping based on the diploid counts.

Things I tried (and decided against):

### Resize windows by coverage
Running copykit read counting diploid cells, then applying mappability correction for cnv clones.
Plan is:
- Run on standard 220kb windows
- Count median reads per window across all cells
- Tile (subdivide/split 200kb windows by 10)
- Count up to median reads to merge window fractions into new windows
This controls for mappability, then just add into granges object as new ref
(This ended up decreasing CNV calls substantially for some reason.)

### Normalize by mappability per bin
- use 220kb bins then matrix multiply by bin mappability weighting (scale bin counts 0-1 then divide) (This created more noisy coverage, mappability is not linear per read count per cell)

### Final decision:

*Ended up using the diploid coverage per bin to remove outlier mapping bins (those over mean +/- 1.5 SD away). Didn't use diploid mapping for bin-level correction. But did plot it as annotation bar.*


## Read in diploid cells from amethyst metadata
```R
library(Rsamtools)
library(copykit)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(dendextend)
library(amethyst)
library(dplyr)
#read in cyto info
cyto=read.table(file="/data/rmulqueen/projects/scalebio_dcis/ref/cytoBand.txt",sep="\t")
colnames(cyto)<-c("chr","start","end","band","stain")
cyto$arm<-substring(cyto$band, 1, 1)
cyto<-cyto[!is.na(cyto$band),]
cyto<-cyto[cyto$chr %in% c(paste0("chr",1:22),"chrX"),]
table(cyto$stain) #set colors for these

#set environment and read in data
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)
obj<-readRDS(file="05_scaledcis.fine_celltype.amethyst.rds")
#write.table(as.data.frame(obj@metadata),file="05_scaledcis.fine_celltype.amethyst.metadata.csv",header=T,sep=",")

#make output directory
system(paste0("mkdir -p ",project_data_directory,"/copykit" ))

#count reads in 220kb, then use varbin_count data to split reads 
read_scalebio_bam<-function(obj_met,x,sample_name){
    #scalebio pipeline outputs bam files as Tn5 wells. so multiple cell IDs are in a bam. this function splits out the bam to the query cellid
    bam=obj_met$bam_path[x]
    cellid=strsplit(row.names(obj_met)[x],"[+]batch|[+]prelim")[[1]][1]
    print(paste("Running sample",cellid))
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

#count nonlumhr cells in high resolution bins, redefine bin sizes to roughly match 220kb coverage
output_directory<-paste0(project_data_directory,"/copykit")
remove_Y = TRUE
min_bincount = 10
cores=100
genome = "hg38"
#resolution="220kb"
#run at 500kb and 280 as well
#resolution="500kb"
resolution="280kb"

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

#dim(diploid_cells)
#[1] 16632    39

####SUBSET FOR SPEED####
#diploid_cells<-diploid_cells[sample(row.names(diploid_cells),size=2000),]

#return chr start position for reads filtered in bam to cell id
varbin_counts_list_all_fields<-mclapply(
                                    1:nrow(diploid_cells), 
                                    function(i) 
                                    read_scalebio_bam(obj_met=diploid_cells,x=i,sample_name=sample_name), 
                                    mc.cores=cores)
                            
message("Read in all bam files.")

names(varbin_counts_list_all_fields)<- row.names(diploid_cells)
varbin_counts_list_all_fields<-as(varbin_counts_list_all_fields, "GRangesList")
ref<-as(rg,"GRanges")
varbin_counts_list <-mclapply(varbin_counts_list_all_fields,
                                function(x) 
                                GenomicRanges::countOverlaps(
                                query=ref,
                                subject=x,
                                type="any",
                                ignore.strand=TRUE),
                                mc.cores=cores)


saveRDS(varbin_counts_list,file=paste0("/data/rmulqueen/projects/scalebio_dcis/ref/copykit.met_windows.",resolution,".diploidcount.rds"))

varbin_counts_list<-readRDS(file=paste0("/data/rmulqueen/projects/scalebio_dcis/ref/copykit.met_windows.",resolution,".diploidcount.rds"))

#TRY MAPPABILITY CORRECTION INSTEAD OF VARIABLE BIN SIZES, DISTRIBUTION TO BIG AND I THINK IT MASKS CNVS
dip_cov<-do.call("cbind",varbin_counts_list)
#normalize by read counts per cell
dip_cov<-as.data.frame(scale(dip_cov, center=FALSE, scale=colSums(dip_cov)))
#get average percentage coverage per bin
varbin_perc_mean<-scale(rowMeans(dip_cov),center=FALSE)
ref$diploid_cov<-varbin_perc_mean

cyto_overlap<-as.data.frame(findOverlaps(ref,makeGRangesFromDataFrame(cyto)))
cyto_overlap<-cyto_overlap[unique(cyto_overlap$queryHits),]

ref$band<-cyto[cyto_overlap$subjectHits,]$band
ref$arm<-cyto[cyto_overlap$subjectHits,]$arm
ref$stain<-cyto[cyto_overlap$subjectHits,]$stain

saveRDS(ref,file=paste0("/data/rmulqueen/projects/scalebio_dcis/ref/copykit.met_windows.",resolution,".diploidcorrected.ref.rds"))
```

## Count original 220kb original windows for metrics

```R
varbin_counts_list<-readRDS(file=paste0("/data/rmulqueen/projects/scalebio_dcis/ref/copykit.met_windows.",resolution,".diploidcount.rds"))

#hg38_grangeslist[["hg38_200kb"]]<-readRDS(file=paste0("/data/rmulqueen/projects/scalebio_dcis/ref/copykit.met_windows.",resolution,".diploidcorrected.ref.rds"))
hg38_grangeslist[["hg38_500kb"]]<-readRDS(file=paste0("/data/rmulqueen/projects/scalebio_dcis/ref/copykit.met_windows.",resolution,".diploidcorrected.ref.rds"))
hg38_grangeslist[["hg38_250kb"]]<-readRDS(file=paste0("/data/rmulqueen/projects/scalebio_dcis/ref/copykit.met_windows.",resolution,".diploidcorrected.ref.rds"))

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

ref<-as(rg,"GRanges")


varbin_counts_list <-mclapply(varbin_counts_list,
                                function(x) 
                                GenomicRanges::countOverlaps(
                                query=ref,
                                subject=x,
                                type="any",
                                ignore.strand=TRUE),
                                mc.cores=cores)
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

# Correct by mappability first, then do GC correction
# LOWESS GC normalization
message("Performing GC correction.")

#runs per cell
varbin_counts_list_gccor <-
    mclapply(varbin_counts_list, function(x) {
        x <- unlist(x) + 1e-3
        x<-base::round(x*(x/rg$diploid_cov),digits=4) #added for coverage correction, multiply by expected bin coverage based on diploids
        gc_cor <- lowess(rg$gc_content, log(x + 1e-3), f = 0.05)
        gc_cor_z <- approx(gc_cor$x, gc_cor$y, rg$gc_content)
        exp(log(x) - gc_cor_z$y) * median(x)
    },mc.cores=cores)

varbin_counts_df <- round(dplyr::bind_cols(varbin_counts_list_gccor), 2)
if(length(min_bc>0)){
names(varbin_counts_df)<-row.names(diploid_cells)[-min_bc]
} else {
    names(varbin_counts_df)<-row.names(diploid_cells)
}
# filtering low read counts where the sum of bins does not reach more than 0
good_cells <- names(varbin_counts_df[which(base::colSums(varbin_counts_df,na.rm=T) != 0)])
varbin_counts_df <- varbin_counts_df[good_cells]

cna_obj <- CopyKit(
    assays = list(bincounts = varbin_counts_df),
    rowRanges = ref)

# Adding genome and resolution information to metadata
S4Vectors::metadata(cna_obj)$genome <- genome
S4Vectors::metadata(cna_obj)$resolution <- resolution #note this is different windows than regular copykit, easiest way to fit it in the functions though

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# ADDING READS METRICS TO METADATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

# saving info and removing columns from list elements
bam_metrics <- diploid_cells[c("unique_reads","tss_enrich","mcg_pct","cg_cov","batch","plate_info","tgmt_well","i7_well","i5_well","fine_celltype")]

# making sure metrics match varbin_counts_df
bam_metrics <- bam_metrics[good_cells,]

bam_metrics$sample <- rownames(bam_metrics)
bam_metrics$sample_name="all_diploid_cells"
bam_metrics$reads_assigned_bins <- colSums(varbin_counts_df)

# adding to metadata
SummarizedExperiment::colData(cna_obj) <-
    S4Vectors::DataFrame(bam_metrics)
colnames(cna_obj) <- names(varbin_counts_df)

sample_name="diploid_2200kb_notcorrected"

saveRDS(cna_obj,file=paste0(output_directory,"/copykit",".",sample_name,".",resolution,".rds"))

#runvarbin module
cna_obj <- runVst(cna_obj)
cna_obj <- runSegmentation(cna_obj)
cna_obj <- logNorm(cna_obj)

# kNN smooth profiles
cna_obj <- knnSmooth(cna_obj)

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
    paste0(output_directory,"/copykit.",sample_name,".",resolution,".pdf")

    saveRDS(cna_obj,file=paste0(output_directory,"/copykit",".",sample_name,".",resolution,".rds"))

    #replot different resolutions on -2 to 2 scale  scale
    #define colors based on data
    log_col=colorRamp2(c(-2,-1,0,1,2),c("darkblue","blue","white","red","darkred"))
    lapply(c("200kb","250kb","500kb"),function(res){
        cna_obj<-readRDS(file=paste0(output_directory,"/copykit",".",sample_name,".",resolution,".rds"))

    plt<-Heatmap(
        t(cna_obj@assays@data$logr),
        left_annotation=ha,col=log_col,
        show_column_names=FALSE,show_row_names=FALSE,
        top_annotation=column_ha,cluster_columns=FALSE,cluster_column_slices=FALSE,column_split=seqnames(cna_obj@rowRanges),
        name="logr")

    pdf(paste0(output_directory,"/copykit.",sample_name,".",resolution,".pdf"),width=20)
    print(plt)
    dev.off()
    paste0(output_directory,"/copykit.",sample_name,".",resolution,".pdf")})

```



## Changing window sizes to correct for mappability


Decimate every range to add more modular scaling of size for met coverage

Now rerun on ~22kb bins and resize to match mean bincounts of passing bins

```R
varbin_counts_list<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/ref/copykit.met_windows.220kb.diploidcount.rds")

hg38_grangeslist[["hg38_200kb"]]<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/ref/copykit.met_windows.220kb.diploidcorrected.ref.rds")

rg<-hg38_grangeslist[["hg38_200kb"]]

ref<-as(rg,"GRanges")
#scaled diploid cov per window
ref$diploid_cov
tiles <- unlist(tile(ref, n = 10))


 varbin_counts_list <-mclapply(varbin_counts_list_all_fields,
                                 function(x) 
                                 GenomicRanges::countOverlaps(
                                 query=tiles,
                                 subject=x,
                                 type="any",
                                 ignore.strand=TRUE),
                                 mc.cores=cores)
 message("Counted reads across all bins.")
 varbin_counts_list <- lapply(varbin_counts_list,as.vector)
 saveRDS(varbin_counts_list,file="/data/rmulqueen/projects/scalebio_dcis/ref/copykit.met_windows.diploidcount.decimated_win.rds")
varbin_counts_list<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/ref/copykit.met_windows.diploidcount.decimated_win.rds")

#make into data frame
decibin_counts<-do.call("cbind",varbin_counts_list)
#normalize by total counts per cell
decibin_counts<-as.data.frame(
    scale(decibin_counts,
    center=FALSE,
    scale=colSums(decibin_counts)))

decibin_mean<-rowMeans(decibin_counts)
summary(decibin_mean)
#       Min.     1st Qu.      Median        Mean     3rd Qu.        Max. 
#0.000000000 0.000007387 0.000008920 0.000008875 0.000010347 0.000019900 

#split by chr
chr_bin_means<-GenomicRanges::split(decibin_mean,f=as.character(seqnames(tiles)))
bin_cov_split<-mean(decibin_mean)*10

#split by mean bin count (of ~22kb bins* 10)
 win_split=data.frame()
 for(chr in unique(names(chr_bin_means))){
     j=1
     for(i in 1:length(chr_bin_means[[chr]])){
         if(sum(chr_bin_means[[chr]][seq(j,i)],na.rm=T)>=bin_cov_split){
             win_split<-rbind(win_split,c(chr,j,i))
             j<-i
         } else if(i==length(chr_bin_means[[chr]])){
             win_split<-rbind(win_split,c(chr,j,i)) #cap off windows
         }
     }
 }

 colnames(win_split)<-c("chr","start_idx","end_idx")
 win_split<-win_split[win_split$chr %in% c(paste0("chr",1:22),"chrX"),]
 ref_df<-as.data.frame(tiles)

 gr_met_resize<-lapply(1:nrow(win_split), function(i){
     chr_in<-win_split[i,]$chr
     ref_chr<-ref_df[ref_df$seqnames==chr_in,] #subset to chr for matched indexing
     row.names(ref_chr)<-as.numeric(1:nrow(ref_chr))
     ref_win_resize<-data.frame(chr=chr_in,
                                 start=ref_chr[as.integer(win_split[i,]$start_idx),]$start,
                                 end=ref_chr[as.integer(win_split[i,]$end_idx),]$end)
     return(ref_win_resize)
 })

 gr_met<-makeGRangesFromDataFrame(do.call("rbind", gr_met_resize),)

 summary(width(gr_met))
 #Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
 #  45693   197580   223438   267482   263000 30348946 


#from https://www.biostars.org/p/478444/ user ATpoint
 GetGC <- function(i){
   seqs <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr_met[i])
   return(as.numeric(Biostrings::letterFrequency(x = seqs, letters = "GC", as.prob = TRUE)))

 }

gr_met$gc_content<-unlist(mclapply(1:length(gr_met),GetGC,mc.cores=100)) #add gc content
gr_met$abspos<-start(gr_met) #add abspos
gr_met$GeneID<-1:length(gr_met) #add GeneID

#add cyto info
cyto_overlap<-as.data.frame(findOverlaps(gr_met,makeGRangesFromDataFrame(cyto)))
cyto_overlap<-cyto_overlap[unique(cyto_overlap$queryHits),]
gr_met$band<-cyto[cyto_overlap$subjectHits,]$band
gr_met$arm<-cyto[cyto_overlap$subjectHits,]$arm
gr_met$stain<-cyto[cyto_overlap$subjectHits,]$stain

#filter out gr_met windows that are >2x mean size
gr_met2<-gr_met[(width(gr_met) <= 2*mean(width(gr_met)))]

#plot bin size of gr_met
library(patchwork)
library(ggplot2)
plt1<-ggplot()+geom_histogram(aes(x=width(gr_met)),bins=100)+theme_minimal()+xlim(c(0,1000000))
plt2<-ggplot()+geom_histogram(aes(x=width(gr_met2)),bins=100)+theme_minimal()+xlim(c(0,1000000))
ggsave(plt1/plt2,file="diploid_resized_window.width.histo.pdf")

summary(width(gr_met))
summary(width(gr_met2))

length(gr_met)

saveRDS(gr_met2,file="/data/rmulqueen/projects/scalebio_dcis/ref/copykit.220kb.met_windows.rds")
```

Perform CNV calling on diploid with new resized windows

```R
hg38_gr<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/ref/copykit.220kb.met_windows.rds")

hg38_grangeslist[["hg38_200kb"]]<-hg38_gr
hg38_rg <- as.data.frame(hg38_gr)

rg <- hg38_rg %>%
    dplyr::rename(chr = "seqnames")

#count reads per bin
varbin_counts_list <-mclapply(varbin_counts_list_all_fields,
                                function(x) 
                                GenomicRanges::countOverlaps(
                                query=hg38_gr,
                                subject=x,
                                type="any",
                                ignore.strand=TRUE))

message("Counted reads across all bins.")
varbin_counts_list <- lapply(varbin_counts_list,as.vector)

# LOWESS GC normalization
message("Performing GC correction.")
#runs per cell
varbin_counts_list_gccor <-
    mclapply(varbin_counts_list, function(x) {
        gc_cor <- lowess(rg$gc_content, log(x + 1e-3), f = 0.05)
        gc_cor_z <- approx(gc_cor$x, gc_cor$y, rg$gc_content)
        exp(log(x) - gc_cor_z$y) * median(x)
    },mc.cores=cores)

varbin_counts_df <- round(dplyr::bind_cols(varbin_counts_list_gccor), 2)
names(varbin_counts_df)<-row.names(diploid_cells)[-min_bc]

# filtering low read counts where the sum of bins does not reach more than 0
good_cells <- names(varbin_counts_df[which(colSums(varbin_counts_df) != 0)])
varbin_counts_df <- varbin_counts_df[good_cells]

cna_obj <- CopyKit(
    assays = list(bincounts = varbin_counts_df),
    rowRanges = hg38_gr)

# Adding genome and resolution information to metadata
S4Vectors::metadata(cna_obj)$genome <- genome
S4Vectors::metadata(cna_obj)$resolution <- "220kb" #note this is different windows than regular copykit, easiest way to fit it in the functions though

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# ADDING READS METRICS TO METADATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

# saving info and removing columns from list elements
bam_metrics <- diploid_cells[c("unique_reads","tss_enrich","mcg_pct","cg_cov","batch","plate_info","tgmt_well","i7_well","i5_well","fine_celltype")]

# making sure metrics match varbin_counts_df
bam_metrics <- bam_metrics[good_cells,]

bam_metrics$sample <- rownames(bam_metrics)
bam_metrics$sample_name="all_diploid_cells"
bam_metrics$reads_assigned_bins <- colSums(varbin_counts_df)

# adding to metadata
SummarizedExperiment::colData(cna_obj) <-S4Vectors::DataFrame(bam_metrics)
colnames(cna_obj) <- names(varbin_counts_df)

#runvarbin module
cna_obj <- runVst(cna_obj)
cna_obj <- runSegmentation(cna_obj,method="multipcf") #CBS fails at merging levels?
cna_obj <- logNorm(cna_obj)

#check S4Vectors::metadata(cna_obj)$vst
sample_name="all_diploid_binresized"
saveRDS(cna_obj,file=paste0(output_directory,"/copykit",".",sample_name,".",resolution,".rds"))

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

sample_name="all_diploid_binresized"
pdf(paste0(output_directory,"/copykit.",sample_name,".",resolution,".pdf"),width=20)

print(plt)
dev.off()

saveRDS(cna_obj,file=paste0(output_directory,"/copykit",".",sample_name,".",resolution,".rds"))

```



# Test run on 41T (lots of clonal structure)
1. Regular 220kb bins
2. Regular 500kb bins
3. 220kb resized bins
3. 220kb bins with mappability correction

```R
library(Rsamtools)
library(copykit)
library(circlize)
library(dendextend)
library(RColorBrewer)
library(ComplexHeatmap)
library(parallel)

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
'lumhr'='#d8007c',
'cancer'="#DFFF00")



#read in cyto info
cyto=read.table(file="/data/rmulqueen/projects/scalebio_dcis/ref/cytoBand.txt",sep="\t")
colnames(cyto)<-c("chr","start","end","band","stain")
cyto$arm<-substring(cyto$band, 1, 1)
cyto<-cyto[!is.na(cyto$band),]
cyto<-cyto[cyto$chr %in% c(paste0("chr",1:22),"chrX"),]
table(cyto$stain) #set colors for these



runCountReads_amethyst <- function(obj,
                        sample_name,
                        genome = "hg38",
                        resolution = c(
                                        "220kb",
                                        "55kb",
                                        "110kb",
                                        "195kb",
                                        "280kb",
                                        "500kb",
                                        "1Mb",
                                        "2.8Mb"),
                        remove_Y = TRUE,
                        min_bincount = 10,
                        cores=100,
                        subclone_addition=5,
                        superclone_addition=2,
                        clus_distance="euclidean") {
    output_directory=paste0(project_data_directory,"/copykit/",sample_name[1])
    system(paste0("mkdir -p ",project_data_directory,"/copykit/",sample_name[1]))

    resolution <- match.arg(resolution)
    #resolution="220kb"

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

    #get list of bams and cellids
    obj_met<-obj@metadata[obj@metadata$Sample %in% sample_name,]
    #return chr start position for reads filtered in bam to cell id
    varbin_counts_list_all_fields<-mclapply(
                                        1:nrow(obj_met), 
                                        function(i) 
                                        read_scalebio_bam(obj_met=obj_met,x=i,sample_name=sample_name), 
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
                                    ignore.strand=TRUE),
                                    mc.cores=cores)
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

    rg <- rg %>%
        dplyr::select(-strand, -GeneID)

    rg_gr <- GenomicRanges::makeGRangesFromDataFrame(rg,
        ignore.strand = TRUE,
        keep.extra.columns = TRUE)

    cna_obj <- CopyKit(
        assays = list(bincounts = varbin_counts_df),
        rowRanges = rg_gr)

    # Adding genome and resolution information to metadata
    S4Vectors::metadata(cna_obj)$genome <- genome
    S4Vectors::metadata(cna_obj)$resolution <- resolution

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sun Feb 14 20:55:01 2021
    # ADDING READS METRICS TO METADATA
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sun Feb 14 20:55:24 2021

    # saving info and removing columns from list elements
    bam_metrics <- obj_met[c("unique_reads","tss_enrich","mcg_pct","cg_cov","batch","plate_info","tgmt_well","i7_well","i5_well","fine_celltype")]

    # making sure metrics match varbin_counts_df
    bam_metrics <- bam_metrics[good_cells,]

    bam_metrics$sample <- rownames(bam_metrics)
    bam_metrics$sample_name=sample_name[1]
    bam_metrics$reads_assigned_bins <- colSums(varbin_counts_df)

    # adding to metadata
    SummarizedExperiment::colData(cna_obj) <-
        S4Vectors::DataFrame(bam_metrics)
    colnames(cna_obj) <- names(varbin_counts_df)
    #runvarbin module
    cna_obj <- runVst(cna_obj)
    cna_obj <- runSegmentation(cna_obj)
    cna_obj <- logNorm(cna_obj)

    # Mark euploid cells if they exist
    #cna_obj <- findAneuploidCells(cna_obj)

    # Mark low-quality cells for filtering
    #cna_obj <- findOutliers(cna_obj)

    # kNN smooth profiles
    cna_obj <- knnSmooth(cna_obj)

    # adds basic quality control information to colData
    cna_obj <- runMetrics(cna_obj)
    cna_obj <- if(nrow(obj_met)<50){
        copykit::runUmap(cna_obj,n_neighbors=nrow(obj_met)-10)
    }else{
        copykit::runUmap(cna_obj)
    }

    dend <- t(cna_obj@assays@data$logr) %>% 
            dist(method=clus_distance) %>% hclust(method="ward.D2") %>% as.dendrogram
    k_optimal=find_k(dend, krange = 2:10)
    print(paste("optimal k value for cutting hclust:", k_optimal$k))
    superclones=dendextend::cutree(dend,k=k_optimal$k+superclone_addition)
    subclones=dendextend::cutree(dend,k=k_optimal$k+subclone_addition)
    cna_obj@colData$subclones<-subclones[row.names(cna_obj@colData)]
    cna_obj@colData$superclones<-superclones[row.names(cna_obj@colData)]

    #define colors based on data
    #updated to be -4 to 4 instead of -2 to 2
    log_col=colorRamp2(c(-4,-2,0,2,4), 
                            c("darkblue","blue","white","red","darkred"))
    cg_perc_col=colorRamp2(c(40,60,80,100),
                            c("#4d2d18","#CABA9C","#4C6444","#102820"))
    reads_col=colorRamp2(c(min(log10(cna_obj@colData$unique_reads)),
                            max(log10(cna_obj@colData$unique_reads))),
                            c("white","black"))

    superclone_col=setNames(nm=unique(as.character(cna_obj@colData$superclones)),
                            colorRampPalette(brewer.pal(9, "Pastel1"))(length(unique(as.character(cna_obj@colData$superclones)))))
    subclone_col=setNames(nm=unique(as.character(cna_obj@colData$subclones)),
                            colorRampPalette(brewer.pal(9, "Spectral"))(length(unique(as.character(cna_obj@colData$subclones)))))
    
    #plot heatmap
    ha = rowAnnotation(
        reads=log10(cna_obj@colData$unique_reads),
        cg_perc=cna_obj@colData$mcg_pct,
        celltype=cna_obj@colData$fine_celltype,
        superclones=as.character(cna_obj@colData$superclones),
        subclones=as.character(cna_obj@colData$subclones),
        col= list(
            celltype=celltype_col,
            reads=reads_col,
            cg_perc=cg_perc_col,
            superclones=superclone_col,
            subclones=subclone_col
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
        row_split=as.character(cna_obj@colData$superclones),
        show_column_names=FALSE,show_row_names=FALSE,
        top_annotation=column_ha,cluster_columns=FALSE,cluster_column_slices=FALSE,column_split=seqnames(cna_obj@rowRanges),
        name="logr")

    pdf(paste0(output_directory,"/copykit.",sample_name[1],".",resolution,".pdf"),width=20)
    print(plt)
    dev.off()

    print(paste("Plotted... ",paste0(output_directory,"/copykit.",sample_name[1],".",resolution,".pdf")))
    cna_obj <- calcConsensus(cna_obj)
    cna_obj <- runConsensusPhylo(cna_obj)
    plt_umap<-plotUmap(cna_obj,label="subclones")
    
    pdf(paste0(output_directory,"/copykit.",sample_name[1],".",resolution,".umap.pdf"))
    print(plt_umap)
    dev.off()

    saveRDS(cna_obj,file=paste0(output_directory,"/copykit",".",sample_name[1],".",resolution,".rds"))
    return(cna_obj)
}


runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS41T'),resolution='220kb',superclone_addition=15) 
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS41T'),resolution='500kb',superclone_addition=15) 


runCountReads_amethyst <- function(obj,
                        sample_name,
                        genome = "hg38",
                        resolution = c(
                                        "220kb",
                                        "55kb",
                                        "110kb",
                                        "195kb",
                                        "280kb",
                                        "500kb",
                                        "1Mb",
                                        "2.8Mb"),
                        remove_Y = TRUE,
                        min_bincount = 10,
                        cores=100,
                        subclone_addition=5,
                        superclone_addition=2,
                        clus_distance="euclidean") {
    output_directory=paste0(project_data_directory,"/copykit/",sample_name[1])
    system(paste0("mkdir -p ",project_data_directory,"/copykit/",sample_name[1]))

    resolution <- match.arg(resolution)
    #resolution="220kb"

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

    #get list of bams and cellids
    obj_met<-obj@metadata[obj@metadata$Sample %in% sample_name,]
    #return chr start position for reads filtered in bam to cell id
    varbin_counts_list_all_fields<-mclapply(
                                        1:nrow(obj_met), 
                                        function(i) 
                                        read_scalebio_bam(obj_met=obj_met,x=i,sample_name=sample_name), 
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
                                    ignore.strand=TRUE),
                                    mc.cores=cores)
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

    rg <- rg %>%
        dplyr::select(-strand, -GeneID)

    rg_gr <- GenomicRanges::makeGRangesFromDataFrame(rg,
        ignore.strand = TRUE,
        keep.extra.columns = TRUE)

    cna_obj <- CopyKit(
        assays = list(bincounts = varbin_counts_df),
        rowRanges = rg_gr)

    # Adding genome and resolution information to metadata
    S4Vectors::metadata(cna_obj)$genome <- genome
    S4Vectors::metadata(cna_obj)$resolution <- resolution

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sun Feb 14 20:55:01 2021
    # ADDING READS METRICS TO METADATA
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sun Feb 14 20:55:24 2021

    # saving info and removing columns from list elements
    bam_metrics <- obj_met[c("unique_reads","tss_enrich","mcg_pct","cg_cov","batch","plate_info","tgmt_well","i7_well","i5_well","fine_celltype")]

    # making sure metrics match varbin_counts_df
    bam_metrics <- bam_metrics[good_cells,]

    bam_metrics$sample <- rownames(bam_metrics)
    bam_metrics$sample_name=sample_name[1]
    bam_metrics$reads_assigned_bins <- colSums(varbin_counts_df)

    # adding to metadata
    SummarizedExperiment::colData(cna_obj) <-
        S4Vectors::DataFrame(bam_metrics)
    colnames(cna_obj) <- names(varbin_counts_df)
    #runvarbin module
    cna_obj <- runVst(cna_obj)
    cna_obj <- runSegmentation(cna_obj)
    cna_obj <- logNorm(cna_obj)

    # Mark euploid cells if they exist
    #cna_obj <- findAneuploidCells(cna_obj)

    # Mark low-quality cells for filtering
    #cna_obj <- findOutliers(cna_obj)

    # kNN smooth profiles
    cna_obj <- knnSmooth(cna_obj)

    # adds basic quality control information to colData
    cna_obj <- runMetrics(cna_obj)
    cna_obj <- if(nrow(obj_met)<50){
        copykit::runUmap(cna_obj,n_neighbors=nrow(obj_met)-10)
    }else{
        copykit::runUmap(cna_obj)
    }

    dend <- t(cna_obj@assays@data$logr) %>% 
            dist(method=clus_distance) %>% hclust(method="ward.D2") %>% as.dendrogram
    k_optimal=find_k(dend, krange = 2:10)
    print(paste("optimal k value for cutting hclust:", k_optimal$k))
    superclones=dendextend::cutree(dend,k=k_optimal$k+superclone_addition)
    subclones=dendextend::cutree(dend,k=k_optimal$k+subclone_addition)
    cna_obj@colData$subclones<-subclones[row.names(cna_obj@colData)]
    cna_obj@colData$superclones<-superclones[row.names(cna_obj@colData)]

    #define colors based on data
    #updated to be -4 to 4 instead of -2 to 2
    log_col=colorRamp2(c(-4,-2,0,2,4), 
                            c("darkblue","blue","white","red","darkred"))
    cg_perc_col=colorRamp2(c(40,60,80,100),
                            c("#4d2d18","#CABA9C","#4C6444","#102820"))
    reads_col=colorRamp2(c(min(log10(cna_obj@colData$unique_reads)),
                            max(log10(cna_obj@colData$unique_reads))),
                            c("white","black"))

    superclone_col=setNames(nm=unique(as.character(cna_obj@colData$superclones)),
                            colorRampPalette(brewer.pal(9, "Pastel1"))(length(unique(as.character(cna_obj@colData$superclones)))))
    subclone_col=setNames(nm=unique(as.character(cna_obj@colData$subclones)),
                            colorRampPalette(brewer.pal(9, "Spectral"))(length(unique(as.character(cna_obj@colData$subclones)))))
    
    #plot heatmap
    ha = rowAnnotation(
        reads=log10(cna_obj@colData$unique_reads),
        cg_perc=cna_obj@colData$mcg_pct,
        celltype=cna_obj@colData$fine_celltype,
        superclones=as.character(cna_obj@colData$superclones),
        subclones=as.character(cna_obj@colData$subclones),
        col= list(
            celltype=celltype_col,
            reads=reads_col,
            cg_perc=cg_perc_col,
            superclones=superclone_col,
            subclones=subclone_col
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
        row_split=as.character(cna_obj@colData$superclones),
        show_column_names=FALSE,show_row_names=FALSE,
        top_annotation=column_ha,cluster_columns=FALSE,cluster_column_slices=FALSE,column_split=seqnames(cna_obj@rowRanges),
        name="logr")

    pdf(paste0(output_directory,"/copykit.",sample_name[1],".",resolution,".pdf"),width=20)
    print(plt)
    dev.off()

    print(paste("Plotted... ",paste0(output_directory,"/copykit.",sample_name[1],".",resolution,".pdf")))
    cna_obj <- calcConsensus(cna_obj)
    cna_obj <- runConsensusPhylo(cna_obj)
    plt_umap<-plotUmap(cna_obj,label="subclones")
    
    pdf(paste0(output_directory,"/copykit.",sample_name[1],".",resolution,".umap.pdf"))
    print(plt_umap)
    dev.off()

    saveRDS(cna_obj,file=paste0(output_directory,"/copykit",".",sample_name[1],".",resolution,".rds"))
    return(cna_obj)
}

```

# old version where i resize windows

<!---

#make output directory
system(paste0("mkdir -p ",project_data_directory,"/copykit" ))

read_scalebio_bam<-function(obj_met,x,sample_name){
    #scalebio pipeline outputs bam files as Tn5 wells. so multiple cell IDs are in a bam. this function splits out the bam to the query cellid
    bam=obj_met$bam_path[x]
    cellid=strsplit(row.names(obj_met)[x],"[+]batch|[+]prelim")[[1]][1]
    print(paste("Running sample",cellid))
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

#count nonlumhr cells in high resolution bins, redefine bin sizes to roughly match 220kb coverage
output_directory<-paste0(project_data_directory,"/copykit")
remove_Y = TRUE
min_bincount = 10
cores=300
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

```

Count original 220kb original windows for metrics
```R

names(varbin_counts_list_all_fields)<- row.names(diploid_cells)
varbin_counts_list_all_fields<-as(varbin_counts_list_all_fields, "GRangesList")
ref<-as(rg,"GRanges")

varbin_counts_list <-mclapply(varbin_counts_list_all_fields,
                                function(x) 
                                GenomicRanges::countOverlaps(
                                query=ref,
                                subject=x,
                                type="any",
                                ignore.strand=TRUE),
                                mc.cores=cores)

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

#4627 bam files had less than 10 mean bincounts and were removed.
length(varbin_counts_list)
mean(unlist(varbin_counts_list))
#33.41627
median(unlist(varbin_counts_list))
#22

```
Now rerun on 20kb bins and resize to match mean bincounts of passing bins

decimate every range to add more modular scaling of size for met coverage

```R

ref<-as(rg,"GRanges")
tiles <- unlist(tile(ref, n = 10))

varbin_counts_list <-mclapply(varbin_counts_list_all_fields,
                                function(x) 
                                GenomicRanges::countOverlaps(
                                query=tiles,
                                subject=x,
                                type="any",
                                ignore.strand=TRUE),
                                mc.cores=cores)
message("Counted reads across all bins.")
varbin_counts_list <- lapply(varbin_counts_list,as.vector)

#mean per index
varbin_counts_mean<-rowMeans(do.call("cbind",varbin_counts_list))

#split by chr
chr_bin_means<-split(varbin_counts_mean,f=as.character(seqnames(tiles)))

#split by 33.4
win_split=data.frame()
for(chr in unique(names(chr_bin_means))){
    j=1
    for(i in 1:length(chr_bin_means[[chr]])){
        if(sum(chr_bin_means[[chr]][seq(j,i)],na.rm=T)>=33.4){
            win_split<-rbind(win_split,c(chr,j,i))
            j<-i
        } else if(i==length(chr_bin_means[[chr]])){
            win_split<-rbind(win_split,c(chr,j,i)) #cap off windows
        }
    }
}



colnames(win_split)<-c("chr","start_idx","end_idx")
win_split<-win_split[win_split$chr %in% c(paste0("chr",1:22),"chrX"),]
ref_df<-as.data.frame(tiles)

gr_met_resize<-lapply(1:nrow(win_split), function(i){
    chr_in<-win_split[i,]$chr
    ref_chr<-ref_df[ref_df$seqnames==chr_in,] #subset to chr for matched indexing
    row.names(ref_chr)<-as.numeric(1:nrow(ref_chr))
    ref_win_resize<-data.frame(chr=chr_in,
                                start=ref_chr[as.integer(win_split[i,]$start_idx),]$start,
                                end=ref_chr[as.integer(win_split[i,]$end_idx),]$end)
    return(ref_win_resize)
})

gr_met<-makeGRangesFromDataFrame(do.call("rbind", gr_met_resize))

#from https://www.biostars.org/p/478444/ user ATpoint
GetGC <- function(i){
  seqs <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr_met[i])
  return(as.numeric(Biostrings::letterFrequency(x = seqs, letters = "GC", as.prob = TRUE)))

}

gr_met$gc_content<-unlist(mclapply(1:length(gr_met),GetGC,mc.cores=100)) #add gc content
gr_met$abspos<-start(gr_met) #add abspos
gr_met$GeneID<-1:length(gr_met) #add GeneID
cyto_overlap<-as.data.frame(findOverlaps(gr_met,makeGRangesFromDataFrame(cyto)))
cyto_overlap<-cyto_overlap[unique(cyto_overlap$queryHits),]

gr_met$band<-cyto[cyto_overlap$subjectHits,]$band
gr_met$arm<-cyto[cyto_overlap$subjectHits,]$arm
gr_met$stain<-cyto[cyto_overlap$subjectHits,]$stain

gr_met
saveRDS(gr_met,file="/data/rmulqueen/projects/scalebio_dcis/ref/copykit.220kb.met_windows.rds")

```
Perform CNV calling with new sized windows on diploid population looking for an improvement and correction for mappability noise.

```R
hg38_gr<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/ref/copykit.220kb.met_windows.rds")

hg38_grangeslist[["hg38_200kb"]]<-hg38_gr
hg38_rg <- as.data.frame(hg38_gr)

rg <- hg38_rg %>%
    dplyr::rename(chr = "seqnames")

#count reads per bin
varbin_counts_list <-mclapply(varbin_counts_list_all_fields,
                                function(x) 
                                GenomicRanges::countOverlaps(
                                query=hg38_gr,
                                subject=x,
                                type="any",
                                ignore.strand=TRUE))

message("Counted reads across all bins.")
varbin_counts_list <- lapply(varbin_counts_list,as.vector)

# LOWESS GC normalization
message("Performing GC correction.")
#runs per cell
varbin_counts_list_gccor <-
    mclapply(varbin_counts_list, function(x) {
        gc_cor <- lowess(rg$gc_content, log(x + 1e-3), f = 0.05)
        gc_cor_z <- approx(gc_cor$x, gc_cor$y, rg$gc_content)
        exp(log(x) - gc_cor_z$y) * median(x)
    },mc.cores=cores)

varbin_counts_df <- round(dplyr::bind_cols(varbin_counts_list_gccor), 2)
names(varbin_counts_df)<-row.names(diploid_cells)[-min_bc]

# filtering low read counts where the sum of bins does not reach more than 0
good_cells <- names(varbin_counts_df[which(colSums(varbin_counts_df) != 0)])
varbin_counts_df <- varbin_counts_df[good_cells]

cna_obj <- CopyKit(
    assays = list(bincounts = varbin_counts_df),
    rowRanges = hg38_gr)

# Adding genome and resolution information to metadata
S4Vectors::metadata(cna_obj)$genome <- genome
S4Vectors::metadata(cna_obj)$resolution <- "220kb" #note this is different windows than regular copykit, easiest way to fit it in the functions though

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# ADDING READS METRICS TO METADATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

# saving info and removing columns from list elements
bam_metrics <- diploid_cells[c("unique_reads","tss_enrich","mcg_pct","cg_cov","batch","plate_info","tgmt_well","i7_well","i5_well","fine_celltype")]

# making sure metrics match varbin_counts_df
bam_metrics <- bam_metrics[good_cells,]

bam_metrics$sample <- rownames(bam_metrics)
bam_metrics$sample_name="all_diploid_cells"
bam_metrics$reads_assigned_bins <- colSums(varbin_counts_df)

# adding to metadata
SummarizedExperiment::colData(cna_obj) <-
    S4Vectors::DataFrame(bam_metrics)
colnames(cna_obj) <- names(varbin_counts_df)

#runvarbin module
cna_obj <- runVst(cna_obj)
cna_obj <- runSegmentation(cna_obj)
cna_obj <- logNorm(cna_obj)

sample_name="all_diploid_binresized"
saveRDS(cna_obj,file=paste0(output_directory,"/copykit",".",sample_name,".",resolution,".rds"))

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

    sample_name="all_diploid_binresized"
    pdf(paste0(output_directory,"/copykit.",sample_name,".",resolution,".pdf"),width=20)

    print(plt)
    dev.off()

    saveRDS(cna_obj,file=paste0(output_directory,"/copykit",".",sample_name,".",resolution,".rds"))

> -->


<!-- 

Use BSBolt to simulate reads for genome, then use that read depth to define window sizes

```bash
mkdir -p /home/rmulqueen/ref/hg38_bsbolt_simulated
cd /home/rmulqueen/ref/hg38_bsbolt_simulated

bsbolt Simulate -G /home/rmulqueen/ref/refdata-gex-GRCm39-2024-A/fasta/genome.fa \
-O /home/rmulqueen/ref/hg38_bsbolt_simulated/hg38_sim_1 \
-NS -PE -CH \
-RD 1 -RL 125 

gzip *fq

bsbolt Align -F1 hg38_sim_1_1.fq.gz -F2 hg38_sim_1_2.fq.gz \
-t 100 -OT 100 -O hg38_sim_1 \
-DB /home/rmulqueen/ref/hg38_bsbolt \
-j >> hg38_sim.bsbolt.log 2>> hg38_sim.bsbolt.log &

samtools view -f 3 -F 3852 hg38_sim_1.bam
samtools markdup -r -S -m s -@ 100 hg38_sim_1.bam -O BAM > hg38_sim_1.dedup.bam
#105105299 reads
samtools view -f 3 -F 3852 hg38_sim_1.bam -O BAM > hg38_sim_1.filt.bam

#431448 reads
samtools sort -@ 100 -m 2G \
--write-index -O BAM \
-o hg38_sim_1.sorted.filt.bam hg38_sim_1.filt.bam

#start with 220kb windows of genome
#make ref genome by chr sizes
bedtools makewindows -g /home/rmulqueen/ref/refdata-gex-GRCm39-2024-A/star/chrNameLength.txt -w 20000 > hg38.20kb.chr.bed
grep "^chr" hg38.20kb.chr.bed > hg38.20kb.chr.filt.bed
awk '{print $1}' hg38.20kb.chr.filt.bed | uniq
#136442 hg38.20kb.chr.bed

#count reads over that and get average (exclude scaffolds and poorly aligned (multi chromosome))
samtools view -f 3 -F 3852 hg38_sim_1.sorted.bam | \
awk -v chr=$chr_in -F'\t' '$3 ~ chr {print $0}' | \
awk -F'\t' '$7 ~ /=/ {print $0}' | \
awk 'OFS="\t" {if ($8<$4) {start=$8; end=$4} else {start=$4; end=$8} print $3,start,end}' | 
sort --parallel=20 -k1,1 -k2,2n - | \
bedtools intersect \
-C -a hg38.20kb.chr.filt.bed -b - \
> 20kb_unadjusted.bed.coverage

136442 20kb_unadjusted.bed.coverage

#new method, use 20kb windows, count and once count exceeds an amount output as 220kb
#get avg read per 11 windows
grep "^chr" 20kb_unadjusted.bed.coverage | grep -v "chrY" | awk '{if (NR % 6 == 0) {a+=$4; print a; a=0} else {a+=$4}}' | awk '{ sum += $1}  END { print sum, NR, sum / NR }'

#17381784 21934 792.458

#then determine number to split at per chr
grep "^chr" 220kb_unadjusted.bed.coverage | grep -v "chrY" | awk '{ sum += $4}  END { print sum, NR, sum / NR }'
#median 790
#mean 820

#split windows by number
count_per_chr() {
    chr_in=$1
    samtools view -f 3 -F 3852 hg38_sim_1.sorted.bam | \
    awk -v chr=$chr_in -F'\t' '$3 ~ chr {print $0}' | \
    awk -F'\t' '$7 ~ /=/ {print $0}' | \
    sort --parallel=20 -k4,4n - | \
    awk 'NR % 820 == 0 {print $3,$4}' | \
    sort -k2,2n #| \
    #awk 'OFS="\t" { if (NR == 1) {a=1} print $1,a,$2; a=$2 }' 
}

export -f count_per_chr 

count_per_chr("chr2")

parallel count_per_chr ::: $(awk '{print $1}' hg38.chr.bed | grep '^chr' | uniq) > 220kb_adjusted.bed


#get gc content 
bedtools nuc -fi /home/rmulqueen/ref/refdata-gex-GRCm39-2024-A/fasta/genome.fa -bed 220kb_adjusted.bed  > 220kb_adjusted.gc.bed 

#for chromosome, split into windows with 79
#simulated reads are aligned to different chr despite being paired?? just using real reads for now
 -->
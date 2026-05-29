
# Generate CopyKit for each sample
Running copykit with filtered bins see:

```
scalemet_dcis/processing/milestonev1_02.1_copykit_binmappability_bydiploid_cells.md
```


Cell types with CNVs are assigned to the broad_celltype metdata column as "cancer".

for calculating filter. Went with 500kb bins after running both 220kb and 500kb. 500kb had less cells filtered out due to coverage issues, and captured the same clones (despite some losses of focal amp/del at lower resolution).

for integer copy number, used scquantum, modified segmentation to not require 3 bins a row (since our bins are bigger)
for scquantum, it failed on some cells (reported as haploid) in those instances i used the segment ratio mean+1 to bring to near diploid

```R

library(Rsamtools)
library(GenomicRanges)
library(copykit)
library(circlize)
library(dendextend)
library(RColorBrewer)
library(ComplexHeatmap)
library(parallel)
library(BiocParallel)
library(Rsamtools)
library(copykit)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(dendextend)
library(amethyst)
library(dplyr)

#set environment and read in data
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=150

project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
obj<-readRDS(file=paste(project_data_directory,"01_amethyst_initial_object","01.0.patient.filt.amethyst.rds",sep="/"))
resolution="500kb"

#make output directory (should already exist from 02.1 script)
processing_folder="02_copykit_cnv_calling"
wd=paste(sep="/",project_data_directory,processing_folder)
system(paste0("mkdir -p ",wd))
setwd(wd)

#set colors
celltype_col=c(
'tcell'='#2e3fa3',
'bcell'='#00adea',
'myeloid'='#00a487',
'peri_VSMC'='#c1d552',
'fibroblast'='#7f1911',
'endothelial'='#f0b243',
'adipocyte'='#d0bd4a',
'basal'='#7200cc',
'lumsec'='#af00af',
'lumhr'='#d8007c')

#read in cyto info
cyto=read.table(file="/data/rmulqueen/projects/scalebio_dcis/ref/cytoBand.txt",sep="\t")
colnames(cyto)<-c("chr","start","end","band","stain")
cyto$arm<-substring(cyto$band, 1, 1)
cyto<-cyto[!is.na(cyto$band),]
cyto<-cyto[cyto$chr %in% c(paste0("chr",1:22),"chrX"),]
table(cyto$stain) #set colors for these

#using my own granges list with the coverage and cyto information added
#setting it and updating it here because both bin counting and running segmentation use it 
#hg38_grangeslist[["hg38_200kb"]]<-readRDS(file=paste0("/data/rmulqueen/projects/scalebio_dcis/ref/copykit.met_windows.220kb.diploidcorrected.ref.rds")) #11268
#hg38_grangeslist[["hg38_250kb"]]<-readRDS(file=paste0("/data/rmulqueen/projects/scalebio_dcis/ref/copykit.met_windows.280kb.diploidcorrected.ref.rds")) #8747
hg38_grangeslist[["hg38_500kb"]]<-readRDS(file=paste0("copykit.met_windows.",resolution,".diploidcorrected.ref.rds"))
length(hg38_grangeslist[["hg38_500kb"]])
#4107

#filtering genomic bins by coverage, takes about 10% of bins
#hg38_grangeslist[["hg38_500kb"]]<-hg38_grangeslist[["hg38_500kb"]][
#    which(
#        hg38_grangeslist[["hg38_500kb"]]$diploid_cov < mean(hg38_grangeslist[["hg38_500kb"]]$diploid_cov)+(1.5*sd(hg38_grangeslist#[["hg38_500kb"]]$diploid_cov)) &
#        hg38_grangeslist[["hg38_500kb"]]$diploid_cov > mean(hg38_grangeslist[["hg38_500kb"]]$diploid_cov)-(1.5*sd(hg38_grangeslist#[["hg38_500kb"]]$diploid_cov))),]
#9691

#hg38_grangeslist[["hg38_250kb"]]<-hg38_grangeslist[["hg38_250kb"]][
#    which(
#        hg38_grangeslist[["hg38_250kb"]]$diploid_cov < mean(hg38_grangeslist[["hg38_250kb"]]$diploid_cov)+(1.5*sd(hg38_grangeslist[["hg38_250kb"]]$diploid_cov)) &
#        hg38_grangeslist[["hg38_250kb"]]$diploid_cov > mean(hg38_grangeslist[["hg38_250kb"]]$diploid_cov)-(1.5*sd(hg38_grangeslist[["hg38_250kb"]]$diploid_cov))),]
#7548

#limit our cnv calling to those with less than 1.5SD in diploid
hg38_grangeslist[["hg38_500kb"]]<-hg38_grangeslist[["hg38_500kb"]][
    which(
        hg38_grangeslist[["hg38_500kb"]]$diploid_cov < mean(hg38_grangeslist[["hg38_500kb"]]$diploid_cov)+(1.5*sd(hg38_grangeslist[["hg38_500kb"]]$diploid_cov)) &
        hg38_grangeslist[["hg38_500kb"]]$diploid_cov > mean(hg38_grangeslist[["hg38_500kb"]]$diploid_cov)-(1.5*sd(hg38_grangeslist[["hg38_500kb"]]$diploid_cov))),]
#3559

#initial counting and object creation. after this we will perform integer estimation via scquantum

read_scalebio_bam<-function(obj_met,x,sample_name){
    #scalebio pipeline outputs bam files as Tn5 wells. so multiple cell IDs are in a bam. this function splits out the bam to the query cellid
    bam=obj_met[obj_met$Sample %in% c(sample_name),]$bam_path[x]
    cellid=strsplit(row.names(obj_met)[x],"[+]batch|[+]prelim")[[1]][1]

    what <- c("qname","rname", "pos")
    param <- Rsamtools::ScanBamParam(what=what,
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
    input_bam<-GenomicRanges::makeGRangesFromDataFrame(input_bam,seqnames.field="rname",start.field="pos",end.field="end")
    message(paste("Finished counting bins for:",cellid))
    return(input_bam)
}

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
                        k_smooth=10,
                        min_bincount = 10,
                        cores=100,
                        subclone_addition=5,
                        superclone_addition=2,
                        clus_distance="euclidean",
                        correct_mappability=FALSE) {
    output_directory=paste0(sample_name[1]) #make directory per sample as subdir
    system(paste0("mkdir -p ",output_directory))
    
    # bindings for NSE and data
    Chr <- chr <- strand <- GeneID <- NULL
    reads_assigned_bins <- reads_duplicates <- reads_total <- NULL

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # genomic ranges (varbin scaffolds)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Reading hg38 VarBin ranges
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

    print("Read in diploid corrected mappability bins reference.")

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
                                        mc.cores=task_cpus)

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
                                    mc.cores=task_cpus)
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
    if(correct_mappability){
        message("Performing GC AND mappability correction.")
        #add correction similar to hmmcopy and copykit
        #https://github.com/shahcompbio/HMMcopy/blob/master/R/correction.R
        #just divide by mappability (1 normalized) per bin
        varbin_counts_list_mapgccor <-
            mclapply(varbin_counts_list_gccor, function(x) {
                x<-unlist(x)
                cov_cor <- lowess(rg$diploid_cov, log(x + 1e-3), f = 0.05)
                cov_cor_z <- approx(cov_cor$x, cov_cor$y, rg$diploid_cov)
                (exp(log(x) - gc_cor_z$y) * median(x))/ref$diploid_cov #added /ref$diploid_cov for mappability coverage
            },mc.cores=task_cpus)
        varbin_counts_df <- round(dplyr::bind_cols(varbin_counts_list_mapgccor), 2)
    } else {
        message("Performing GC correction.")
        varbin_counts_list_gccor <-
            mclapply(varbin_counts_list, function(x) {
                gc_cor <- lowess(rg$gc_content, log(x + 1e-3), f = 0.05)
                gc_cor_z <- approx(gc_cor$x, gc_cor$y, rg$gc_content)
                (exp(log(x) - gc_cor_z$y) * median(x)) 
            },mc.cores=task_cpus)
        varbin_counts_df <- round(dplyr::bind_cols(varbin_counts_list_gccor), 2)
    }

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
    S4Vectors::metadata(cna_obj)$genome <- "hg38"
    S4Vectors::metadata(cna_obj)$resolution <- resolution

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sun Feb 14 20:55:01 2021
    # ADDING READS METRICS TO METADATA
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sun Feb 14 20:55:24 2021

    # saving info and removing columns from list elements
    bam_metrics <- obj_met[c("unique_reads","tss_enrich","mcg_pct","cg_cov","batch","plate_info","tgmt_well","i7_well","i5_well","coarse_celltype")]

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

    # kNN smooth profiles
    cna_obj <- knnSmooth(cna_obj,k=k_smooth) #knn smoothing reruns segmentation

    # adds basic quality control information to colData
    cna_obj <- runMetrics(cna_obj)
    cna_obj <- if(nrow(obj_met)<50){
        copykit::runUmap(cna_obj,n_neighbors=nrow(obj_met)-10)
    }else{
        copykit::runUmap(cna_obj)
    }

    clus_distance="euclidean"
    dend <- t(cna_obj@assays@data$logr) %>% 
            dist(method=clus_distance) %>% hclust(method="ward.D2") %>% as.dendrogram
    k_optimal=find_k(dend, krange = 2:10)
    saveRDS(dend,file=paste0(output_directory,"/copykit.",sample_name[1],".",resolution,".dendrogram.rds"))

    print(paste("optimal k value for cutting hclust:", k_optimal$k))
    superclones=dendextend::cutree(dend,k=k_optimal$k+superclone_addition)
    subclones=dendextend::cutree(dend,k=k_optimal$k+subclone_addition)
    cna_obj@colData$subclones<-subclones[row.names(cna_obj@colData)]
    cna_obj@colData$superclones<-superclones[row.names(cna_obj@colData)]

    #define colors based on data
    #updated to be -4 to 4 instead of -2 to 2
    log_col=colorRamp2(c(-2,-1.5,0,1.5,2), 
                            c("darkblue","blue","white","red","darkred"))
    cg_perc_col=colorRamp2(c(40,60,80,100),
                            c("#4d2d18","#CABA9C","#4C6444","#102820"))
    reads_col=colorRamp2(c(min(log10(cna_obj@colData$unique_reads)),
                            max(log10(cna_obj@colData$unique_reads))),
                            c("white","black"))
    integer_col=colorRamp2(c(-4,-2,0,2,4), 
                            c("darkblue","blue","white","red","darkred"))
    superclone_col=setNames(nm=unique(as.character(cna_obj@colData$superclones)),
                            colorRampPalette(brewer.pal(9, "Pastel1"))(length(unique(as.character(cna_obj@colData$superclones)))))
    subclone_col=setNames(nm=unique(as.character(cna_obj@colData$subclones)),
                            colorRampPalette(brewer.pal(9, "Spectral"))(length(unique(as.character(cna_obj@colData$subclones)))))
    
    #plot heatmap
    ha = rowAnnotation(
        reads=log10(cna_obj@colData$unique_reads),
        cg_perc=cna_obj@colData$mcg_pct,
        celltype=cna_obj@colData$coarse_celltype,
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
    dip_cov=colorRamp2(c(0.5,1,1.5), 
                            c("white","grey","black"))
    column_ha = HeatmapAnnotation(
        mappability=cna_obj@rowRanges$diploid_cov,
        arm = cna_obj@rowRanges$arm,
        band = cna_obj@rowRanges$stain,
        col=list(mappability=dip_cov,arm=arm_col,band=band_col))

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

    print(paste("Plotted... ",paste0(output_directory,"/copykit.",sample_name[1],".",resolution,".pdf"), "Log Segment Ratios"))

    #cna_obj <- calcInteger(cna_obj, method = 'scquantum', assay = 'smoothed_bincounts')

    #cna_obj <- calcConsensus(cna_obj)
    #cna_obj <- runConsensusPhylo(cna_obj)
    #plt_umap <- plotUmap(cna_obj,label="subclones")
    
    #pdf(paste0(output_directory,"/copykit.",sample_name[1],".",resolution,".umap.pdf"))
    #print(plt_umap)
    #dev.off()

    saveRDS(cna_obj,file=paste0(output_directory,"/copykit",".",sample_name[1],".",resolution,".rds"))
    return(cna_obj)
}

#running with no bin filter
register(MulticoreParam(progressbar = T, workers = 125), default = T)

#rerun all of this at 500kb as well (since it filters out less cells)
res='500kb' #'500kb'
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS05T'),resolution=res,superclone_addition=15) #done #rerun
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS07T'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS102T_24hTis'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS124T'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS22T'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS28T'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS32T'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS35T'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS41T'),resolution=res,superclone_addition=15)  #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS49T'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS52T'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS65T'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS66T'),resolution=res,superclone_addition=15) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS70T'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS74T'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS79T_24hTis_DCIS','BCMDCIS79T_24hTis_IDC'),resolution=res,superclone_addition=15) #done #rerun
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS80T_24hTis'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS82T_24hTis'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS92T_24hTis'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS94T_24hTis'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS97T'),resolution=res,superclone_addition=15) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS99T'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA03R'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA04R'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA09R-3h'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA12R-3h'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA16R-3h'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA17R-3h'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA19R-4h'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA22R-4h'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA26L-24hTis-4h'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA29L-2h'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA38L-3h'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA83L-3h'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA85L-3h'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('ECIS25T'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('ECIS26T'),resolution=res,superclone_addition=15) #done
runCountReads_amethyst(obj=obj,sample_name=c('ECIS36T'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('ECIS48T'),resolution=res) #done
runCountReads_amethyst(obj=obj,sample_name=c('ECIS57T'),resolution=res) #done
```



# Assign clone names based on log ratios
Assign aneuploid and diploid clones, and subclones per sample.

```R
library(copykit)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(amethyst)


project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
obj<-readRDS(file=paste(project_data_directory,"01_amethyst_initial_object","01.0.patient.filt.amethyst.rds",sep="/"))
resolution="500kb"

#make output directory (should already exist from 02.1 script)
processing_folder="02_copykit_cnv_calling"
wd=paste(sep="/",project_data_directory,processing_folder)
system(paste0("mkdir -p ",wd))
setwd(wd)

#read in cyto info
cyto=read.table(file="/data/rmulqueen/projects/scalebio_dcis/ref/cytoBand.txt",sep="\t")
colnames(cyto)<-c("chr","start","end","band","stain")
cyto$arm<-substring(cyto$band, 1, 1)
cyto<-cyto[!is.na(cyto$band),]
cyto<-cyto[cyto$chr %in% c(paste0("chr",1:22),"chrX"),]
table(cyto$stain) #set colors for these


#set colors
celltype_col=c(
'tcell'='#2e3fa3',
'bcell'='#00adea',
'myeloid'='#00a487',
'peri_VSMC'='#c1d552',
'fibroblast'='#7f1911',
'endothelial'='#f0b243',
'adipocyte'='#d0bd4a',
'basal'='#7200cc',
'lumsec'='#af00af',
'lumhr'='#d8007c',
'cancer'='#80FF80')


copykit_output<-list.files(path=paste0(project_data_directory,"/",processing_folder),recursive=TRUE,full.names=TRUE,pattern="kb.rds")
#remove diploid cell call rds used for bin correction
copykit_output<-copykit_output[!grepl(copykit_output,pattern="diploid")]

#make clones a named list to collapse overclustering or low cell counts/cluster
assign_copykit_aneuploid_clonename<-function(sample_name,cancer_clones,split_on="superclones",resolution='220kb',dat=dat){
    tmp<-readRDS(paste0(project_data_directory,"/",processing_folder,"/",sample_name,"/copykit.",sample_name,".",resolution,".rds"))

    tmp@colData$ploidy<-"diploid"
    tmp@colData$coarse_celltype<-dat@metadata[row.names(tmp@colData),]$coarse_celltype

    tmp@rowRanges<-sort(tmp@rowRanges)
    if(length(cancer_clones)>0){
        if(split_on=="subclones"){
        tmp@colData$clones_split<-"subclones"
        tmp@colData[tmp@colData$subclones %in% cancer_clones,]$ploidy<-"aneuploid"
        tmp@colData$clonename<-unlist(paste(sample_name,names(cancer_clones[match(tmp@colData$subclones,cancer_clones)]),sep="_"))
        tmp@colData$clonename<-gsub("_NA", replacement = "_diploid", x = tmp@colData$clonename)
        tmp@colData[tmp@colData$ploidy=="aneuploid",]$coarse_celltype<-"cancer"
        }else{
        tmp@colData$clones_split<-"superclones"
        tmp@colData[tmp@colData$superclones %in% cancer_clones,]$ploidy<-"aneuploid"
        tmp@colData$clonename<-unlist(paste(sample_name,names(cancer_clones[match(tmp@colData$superclones,cancer_clones)]),sep="_"))
        tmp@colData$clonename<-gsub("_NA", replacement = "_diploid", x = tmp@colData$clonename)
        tmp@colData[tmp@colData$ploidy=="aneuploid",]$coarse_celltype<-"cancer"
        }} else {
        tmp@colData$clones_split<-"all_diploid"
        tmp@colData$clonename<-paste(sample_name,"diploid",sep="_")
        }
    #define colors based on data
    #updated to be -4 to 4 instead of -2 to 2
    log_col=colorRamp2(c(-3,-2,-1,0,1,2,3), 
                            c("#053061","#2166ac","#4393c3","white","#d6604d","#b2182b","#67001f"))
    cg_perc_col=colorRamp2(c(40,60,80,100),
                            c("#4d2d18","#CABA9C","#4C6444","#102820"))
    reads_col=colorRamp2(c(min(log10(tmp@colData$unique_reads)),
                            max(log10(tmp@colData$unique_reads))),
                            c("white","black"))

    superclone_col=setNames(nm=unique(as.character(tmp@colData$superclones)),
                            colorRampPalette(brewer.pal(9, "Pastel1"))(length(unique(as.character(tmp@colData$superclones)))))
    subclone_col=setNames(nm=unique(as.character(tmp@colData$subclones)),
                            colorRampPalette(brewer.pal(9, "Spectral"))(length(unique(as.character(tmp@colData$subclones)))))
    cancerclone_col=setNames(nm=unique(as.character(tmp@colData$clonename)),
                            colorRampPalette(brewer.pal(8, "Pastel2"))(length(unique(as.character(tmp@colData$clonename)))))
    
    #plot heatmap
    ha = rowAnnotation(
        reads=log10(tmp@colData$unique_reads),
        cg_perc=tmp@colData$mcg_pct,
        celltype=tmp@colData$coarse_celltype,
        superclones=as.character(tmp@colData$superclones),
        subclones=as.character(tmp@colData$subclones),
        cancerclone=as.character(tmp@colData$clonename),
        col= list(
            celltype=celltype_col,
            reads=reads_col,
            cg_perc=cg_perc_col,
            superclones=superclone_col,
            subclones=subclone_col,
            cancerclone_col
        ))

    cyto_overlap<-GenomicRanges::findOverlaps(tmp@rowRanges,
                                                makeGRangesFromDataFrame(cyto,keep=TRUE),
                                                select="first")
    tmp@rowRanges$stain <- cyto[cyto_overlap,]$stain
    tmp@rowRanges$arm <- cyto[cyto_overlap,]$arm

    arm_col=c("p"="grey","q"="darkgrey")
    band_col=c("acen"="#99746F","gneg"="white","gpos100"="black","gpos25"="lightgrey","gpos50"="grey","gpos75"="darkgrey","gvar"="#446879")
    dip_cov=colorRamp2(c(0.5,1,1.5), 
                            c("white","grey","black"))
    column_ha = HeatmapAnnotation(
        mappability=tmp@rowRanges$diploid_cov,
        arm = tmp@rowRanges$arm,
        band = tmp@rowRanges$stain,
        col=list(mappability=dip_cov,arm=arm_col,band=band_col))

    logdat<-t(tmp@assays@data$logr)
    plt<-Heatmap(
        logdat,
        clustering_distance_rows = "manhattan", 
        left_annotation=ha, col=log_col,
        row_split=as.character(tmp@colData$clonename),
        show_column_names=TRUE,show_row_names=FALSE,
        top_annotation=column_ha,
        cluster_columns=FALSE,
        cluster_column_slices=FALSE,
        column_split=seqnames(tmp@rowRanges),
        name="logr")

    pdf(paste0(project_data_directory,"/",processing_folder,"/",sample_name,"/copykit.",sample_name,".",resolution,".cancerclone.pdf"),width=20)
    print(plt)
    dev.off()

    saveRDS(tmp,file=paste0(project_data_directory,"/",processing_folder,"/",sample_name,"/copykit.",sample_name,".",resolution,".rds"))
}
```

# for 500kb
```R
res="500kb"

assign_copykit_aneuploid_clonename(dat=obj,resolution=res,sample_name="BCMDCIS05T",
                                cancer_clones=c("c1"='10',"c1"='4',"c1"='6',
                                'c2'='16','c2'='5',
                                'c3'='14','c3'='15','c3'='13',
                                'c4'='8')) #done

assign_copykit_aneuploid_clonename(dat=obj,resolution=res,sample_name='BCMDCIS07T',cancer_clones=c()) #done

assign_copykit_aneuploid_clonename(dat=obj,resolution=res,sample_name='BCMDCIS102T_24hTis',split_on='subclones',cancer_clones=c("c1"='5','c1'='7',
                                                                                                                            'c2'='4',
                                                                                                                            'c3'='3',
                                                                                                                            'c4'='2')) #done
assign_copykit_aneuploid_clonename(dat=obj,resolution=res,sample_name='BCMDCIS124T',cancer_clones=c("c1"='3',
                                                                                                'c2'='4',
                                                                                                "c3"='5')) #done

assign_copykit_aneuploid_clonename(dat=obj,resolution=res,sample_name='BCMDCIS22T',cancer_clones=c("c1"='2',
                                                                                                    'c2'='3')) #done 


assign_copykit_aneuploid_clonename(dat=obj,resolution=res,sample_name='BCMDCIS28T',cancer_clones=c("c1"='4',
                                                                                                    "c2"='2')) #done

assign_copykit_aneuploid_clonename(dat=obj,resolution=res,sample_name='BCMDCIS32T',cancer_clones=c()) #done all diploid

assign_copykit_aneuploid_clonename(dat=obj,resolution=res,sample_name='BCMDCIS35T',cancer_clones=c(
                                                            "c1"='1',"c1"='2',"c1"='3',"c1"='4',"c1"='6')) #done 

assign_copykit_aneuploid_clonename(dat=obj,resolution=res,sample_name='BCMDCIS41T',cancer_clones=c('c1'='2',
                                                                                        'c2'='5',
                                                                                        'c3'='3','c3'='12','c3'='18','c3'='4',
                                                                                        'c4'='10','c4'='8','c4'='14',
                                                                                        'c5'='7','c5'='17',
                                                                                        'c6'='9','c6'='13','c6'='16')) #done

assign_copykit_aneuploid_clonename(dat=obj,resolution=res,sample_name='BCMDCIS49T',cancer_clones=c()) #diploid
assign_copykit_aneuploid_clonename(dat=obj,resolution=res,sample_name='BCMDCIS52T',cancer_clones=c('c1'='6','c1'='2','c1'='1')) #done
assign_copykit_aneuploid_clonename(dat=obj,resolution=res,sample_name='BCMDCIS65T',cancer_clones=c('c1'='7','c1'='5')) #done


assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMDCIS66T',cancer_clones=c('c1'='6','c1'='5',
'c2'='1',
'c3'='15','c3'='16','c3'='13','c3'='10','c3'='9','c3'='11','c3'='14','c3'='12'))#done

assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMDCIS70T',split_on='subclones',cancer_clones=c('c1'='6',
                                    'c2'='5','c2'='2','c2'='4')) #done

assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMDCIS74T',split_on='subclones',cancer_clones=c('c1'='5',
                                                                                'c2'='3','c2'='12','c2'='9',
                                                                                'c3'='13',
                                                                                'c4'='7',
                                                                                'c5'='10','c5'='8',
                                                                                'c6'='6')) #done

assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMDCIS79T_24hTis_DCIS',
                                    cancer_clones=c('c1'='7','c1'='9','c1'='5',
                                                    'c2'='3',
                                                    'c3'='8','c3'='11','c3'='5','c3'='6','c3'='14','c3'='15','c3'='16','c3'='13')) #done


assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMDCIS80T_24hTis',cancer_clones=c('c1'='5',
                                                                                    'c2'='4',
                                                                                    'c3'='6')) #done

assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMDCIS82T_24hTis',cancer_clones=c('c1'='3','c1'='4')) #done

assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMDCIS92T_24hTis',cancer_clones=c('c1'='3','c1'='4')) #done

assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMDCIS94T_24hTis',split_on='subclones',cancer_clones=c('c1'='7',
                                                                                                        'c2'='5','c2'='3','c2'='6')) #done

assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMDCIS97T',cancer_clones=c(
                            'c1'='3','c1'='7','c1'='14','c1'='5','c1'='10','c1'='12','c1'='6','c1'='9',
                                                                            'c2'='13',
                                                                            'c3'='16')) #done

assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMDCIS99T',cancer_clones=c('c1'='5')) #done

assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMHBCA03R',cancer_clones=c('c1'='4')) #looks real to me

assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMHBCA04R',cancer_clones=c()) #done 
assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMHBCA09R-3h',cancer_clones=c()) #done
assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMHBCA12R-3h',cancer_clones=c()) #done
assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMHBCA16R-3h',cancer_clones=c()) #done
assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMHBCA17R-3h',cancer_clones=c()) #done
assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMHBCA19R-4h',cancer_clones=c()) #maybe an aneuploid cell?
assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMHBCA22R-4h',cancer_clones=c()) #done
assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMHBCA26L-24hTis-4h',cancer_clones=c()) #done
assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMHBCA29L-2h',cancer_clones=c()) #done
assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMHBCA38L-3h',cancer_clones=c()) #done
assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMHBCA83L-3h',cancer_clones=c('c1'='6')) #done
assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='BCMHBCA85L-3h',cancer_clones=c()) #done
assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='ECIS25T',split_on="subclones",cancer_clones=c('c1'='7','c1'='4',
                                                                                            'c2'='1','c2'='5',
                                                                                            'c3'='3',
                                                                                            'c4'='6')) #done

assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='ECIS26T',cancer_clones=c('c1'='15','c1'='14','c1'='20','c1'='18',
                                                                    'c2'='16',
                                                                    'c3'='3','c3'='11','c3'='5','c3'='17','c3'='4','c3'='6','c3'='9','c3'='19',
                                                                    'c4'='7',
                                                                    'c5'='1','c5'='10','c5'='21','c5'='13','c5'='2'
                                                                    )) #done

assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='ECIS36T',split_on="subclones",cancer_clones=c('c1'='2','c1'='7','c3'='1','c2'='5','c3'='6')) #done
assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='ECIS48T',cancer_clones=c()) #maybe one cell showing aneuploidy?
assign_copykit_aneuploid_clonename(dat=dat,resolution=res,sample_name='ECIS57T',split_on="subclones",cancer_clones=c('c1'='5','c1'='4','c1'='6','c1'='2','c2'='1','c1'='7')) #done

```




<!--
## for 220kb

```R
assign_copykit_aneuploid_clonename(dat=dat,sample_name="BCMDCIS05T",cancer_clones=c("c1"='2',"c2"='3',"c2"='4')) #done
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS07T',cancer_clones=c()) #done all diploid #x loss?
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS102T_24hTis',cancer_clones=c("c1"='3','c2'='4')) #done some x loss?
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS124T',cancer_clones=c("c1="='2',"c2"='4')) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS22T',cancer_clones=c("c1"='2','c2'='4')) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS28T',,cancer_clones=c("c1"='5',"c2"='3','c2'='6','c2'='7')) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS32T',cancer_clones=c()) #all diploid
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS35T',cancer_clones=c("c1"='1',"c1"='2',"c1"='3',"c1"='5',"c1"='6',"c1"='7',"c1"='8'))
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS41T',split_on='subclones',cancer_clones=c('c1'='4','c1'='6',
                                                                                        'c2'='1',
                                                                                        'c3'='3','c3'='7',
                                                                                        'c4'='2')) 
#assign_copykit_aneuploid_clonename(sample_name='BCMDCIS49T',cancer_clones=c()) #not run
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS52T',cancer_clones=c('c1'='3','c2'='2','c2'='4')) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS65T',cancer_clones=c('c1'='5')) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS66T',cancer_clones=c('c1'='2',
                                                                            'c2'='4','c2'='11',
                                                                            'c3'='14','c3'='12','c3'='7','c3'='15','c3'='13','c3'='8','c3'='5'))                                   
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS70T',cancer_clones=c('c1'='2','c1'='3','c1'='4')) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS74T',cancer_clones=c('c1'='11','c1'='8','c1'='9',
                                                                            'c2'='6',
                                                                            'c3'='10',
                                                                            'c4'='3',
                                                                            'c5'='2',
                                                                            'c6'='4',
                                                                            'c8'='12',
                                                                            'c9'='7')) #lotsa clear evolution in this one
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS79T_24hTis_DCIS',cancer_clones=c('c1'='2','c2'='4')) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS80T_24hTis',cancer_clones=c('c1'='3',
                                                                                    'c1'='5','c1'='4')) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS82T_24hTis',cancer_clones=c('c1'='3','c1'='4','c1'='2')) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS92T_24hTis',cancer_clones=c('c1'='3','c1'='4')) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS94T_24hTis',split_on='subclones',cancer_clones=c('c1'='5','c1'='4',
                                                                                                        'c2'='3')) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS97T',cancer_clones=c('c1'='11','c1'='5','c1'='10','c1'='12','c1'='17','c1'='7',
                                                                            'c2'='6','c2'='2','c2'='8',
                                                                            'c4'='15',
                                                                            'c3'='9')) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMDCIS99T',cancer_clones=c('c1'='4')) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMHBCA03R',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMHBCA04R',split_on='subclones',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMHBCA09R-3h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMHBCA12R-3h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMHBCA16R-3h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMHBCA17R-3h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMHBCA19R-4h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMHBCA22R-4h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMHBCA26L-24hTis-4h',cancer_clones=c()) #i dont see any cancer clones previously split_on="subclones",cancer_clones=c('c1'='6')
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMHBCA29L-2h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMHBCA38L-3h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMHBCA83L-3h',cancer_clones=c('c1'='4')) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='BCMHBCA85L-3h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='ECIS25T',split_on="subclones",cancer_clones=c('c1'='7','c1'='4','c1'='6',
                                                                                            'c2'='5','c2'='1','c2'='8',
                                                                                            'c3'='9','c3'='2')) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='ECIS26T',cancer_clones=c('c1'='3','c1'='17','c1'='15','c1'='9','c1'='10','c1'='5','c1'='11',
                                                                        'c1'='13','c1'='4','c1'='1','c1'='6','c1'='16','c1'='7','c1'='8',
                                                                        'c2'='14')) 

assign_copykit_aneuploid_clonename(dat=dat,sample_name='ECIS36T',,cancer_clones=c('c1'='1','c2'='4','c3'='2')) 
assign_copykit_aneuploid_clonename(dat=dat,sample_name='ECIS48T',cancer_clones=c('c1'='5')) #might have a cancer precursor in the diploid pop chr16 loss in some lumhr, 1q gain in lumsec?
assign_copykit_aneuploid_clonename(dat=dat,sample_name='ECIS57T',cancer_clones=c('c1'='2','c1'='3','c1'='4')) 
```
-->

Assign both 220kb and 500kb clone names into amethyst metadata. Plot alluvial plot to show consistency in names.

```R
library(Rsamtools)
library(GenomicRanges)
library(copykit)
library(circlize)
library(dendextend)
library(RColorBrewer)
library(ComplexHeatmap)
library(parallel)
library(BiocParallel)
library(amethyst)
set.seed(111)

project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
obj<-readRDS(file=paste(project_data_directory,"01_amethyst_initial_object","01.0.patient.filt.amethyst.rds",sep="/"))
resolution="500kb"

#make output directory (should already exist from 02.1 script)
processing_folder="02_copykit_cnv_calling"
wd=paste(sep="/",project_data_directory,processing_folder)
system(paste0("mkdir -p ",wd))
setwd(wd)

#read in cyto info
cyto=read.table(file="/data/rmulqueen/projects/scalebio_dcis/ref/cytoBand.txt",sep="\t")
colnames(cyto)<-c("chr","start","end","band","stain")
cyto$arm<-substring(cyto$band, 1, 1)
cyto<-cyto[!is.na(cyto$band),]
cyto<-cyto[cyto$chr %in% c(paste0("chr",1:22),"chrX"),]
table(cyto$stain) #set colors for these


#set colors
celltype_col=c(
'tcell'='#2e3fa3',
'bcell'='#00adea',
'myeloid'='#00a487',
'peri_VSMC'='#c1d552',
'fibroblast'='#7f1911',
'endothelial'='#f0b243',
'adipocyte'='#d0bd4a',
'basal'='#7200cc',
'lumsec'='#af00af',
'lumhr'='#d8007c',
'cancer'='#80FF80')


copykit_output<-list.files(path=paste0(project_data_directory,"/",processing_folder),recursive=TRUE,full.names=TRUE,pattern="kb.rds")
#remove diploid cell call rds used for bin correction
copykit_output<-copykit_output[!grepl(copykit_output,pattern="diploid")]

#remove diploid cell call rds used for bin correction
copykit_output_500kb<-copykit_output[!grepl(copykit_output,pattern="diploid")]

#read in all meta data from copykit, append to amethyst object
read_meta_copykit<-function(x){
    tmp<-readRDS(x)
    meta<-as.data.frame(tmp@colData[c("sample_name","reads_assigned_bins","plate_info","superclones","subclones","ploidy","clonename","coarse_celltype","clones_split")])
    return(meta)
}

cnv_meta_500kb<-do.call("rbind",lapply(copykit_output_500kb,read_meta_copykit))

obj_backup<-obj
#obj@metadata$cnv_ploidy_220kb<-"NA" #na values are too low read count to call
#obj@metadata[row.names(cnv_meta_220kb),]$cnv_ploidy_220kb<-cnv_meta_220kb$ploidy
#obj@metadata$cnv_superclones_220kb<-"NA"
#obj@metadata[row.names(cnv_meta_220kb),]$cnv_superclones_220kb<-cnv_meta_220kb$superclones
#obj@metadata$cnv_subclones_220kb<-"NA"
#obj@metadata[row.names(cnv_meta_220kb),]$cnv_subclones_220kb<-cnv_meta_220kb$subclones
#obj@metadata$cnv_clonename_220kb<-"NA"
#obj@metadata[row.names(cnv_meta_220kb),]$cnv_clonename_220kb<-cnv_meta_220kb$clonename
#obj@metadata$cnv_clones_split_220kb<-"NA"
#obj@metadata[row.names(cnv_meta_220kb),]$cnv_clones_split_220kb<-cnv_meta_220kb$clones_split

obj@metadata$cnv_ploidy_500kb<-"NA" #na values are too low read count to call
obj@metadata[row.names(cnv_meta_500kb),]$cnv_ploidy_500kb<-cnv_meta_500kb$ploidy
obj@metadata$cnv_superclones_500kb<-"NA"
obj@metadata[row.names(cnv_meta_500kb),]$cnv_superclones_500kb<-cnv_meta_500kb$superclones
obj@metadata$cnv_subclones_500kb<-"NA"
obj@metadata[row.names(cnv_meta_500kb),]$cnv_subclones_500kb<-cnv_meta_500kb$subclones
obj@metadata$cnv_clonename_500kb<-"NA"
obj@metadata[row.names(cnv_meta_500kb),]$cnv_clonename_500kb<-cnv_meta_500kb$clonename
obj@metadata$cnv_clones_split_500kb<-"NA"
obj@metadata[row.names(cnv_meta_500kb),]$cnv_clones_split_500kb<-cnv_meta_500kb$clones_split


#final cnv clones based on 500kb calling
obj@metadata$cnv_clonename<-obj@metadata$cnv_clonename_500kb

#also looks like cluster 17 is mislabelled (shows basal and lumhr markers, see KRT5 and ANKRD30A resp. So calling lumhr since its the expected source for ER+ cancer, but noted to investigate further, patient 97T)
obj@metadata[which(obj@metadata$coarse_cluster_phenograph=="17"),]$coarse_celltype<-"lumhr"

obj@metadata[which(obj@metadata$cnv_ploidy_500kb=="aneuploid"),]$coarse_celltype<-"cancer"
saveRDS(obj,file="02_scaledcis.cnv_clones.amethyst.rds")


```

```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
```

# Generate CopyKit for each sample

```R
source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
library(Rsamtools)
library(copykit)
library(circlize)
detach("package:GeneNMF",unload=TRUE)
library(dendextend)
library(RColorBrewer)
library(ComplexHeatmap)

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

#set environment and read in data
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)
obj<-readRDS(file="05_scaledcis.fine_celltype.amethyst.rds")

#using custom resized windows for methylation mapping
#hg38_grangeslist[["hg38_200kb"]]<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/ref/copykit.220kb.met_windows.rds")


#read in cyto info
cyto=read.table(file="/data/rmulqueen/projects/scalebio_dcis/ref/cytoBand.txt",sep="\t")
colnames(cyto)<-c("chr","start","end","band","stain")
cyto$arm<-substring(cyto$band, 1, 1)
cyto<-cyto[!is.na(cyto$band),]
cyto<-cyto[cyto$chr %in% c(paste0("chr",1:22),"chrX"),]
table(cyto$stain) #set colors for these

#make output directory
system(paste0("mkdir -p ",project_data_directory,"/copykit" ))

read_scalebio_bam<-function(obj_met,x,sample_name){
    #scalebio pipeline outputs bam files as Tn5 wells. so multiple cell IDs are in a bam. this function splits out the bam to the query cellid
    bam=obj_met[obj_met$Sample %in% c(sample_name),]$bam_path[x]
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

runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS05T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS07T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS102T_24hTis'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS124T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS22T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS28T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS32T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS35T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS41T'),resolution='220kb',superclone_addition=15) 
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS49T'),resolution='220kb') #59 cells
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS52T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS65T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS66T'),resolution='220kb',superclone_addition=15)
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS70T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS74T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS79T_24hTis_DCIS','BCMDCIS79T_24hTis_IDC'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS80T_24hTis'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS82T_24hTis'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS92T_24hTis'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS94T_24hTis'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS97T'),resolution='220kb',superclone_addition=15)
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS99T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA03R'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA04R'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA09R-3h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA12R-3h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA16R-3h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA17R-3h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA19R-4h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA22R-4h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA26L-24hTis-4h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA29L-2h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA38L-3h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA83L-3h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA85L-3h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('ECIS25T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('ECIS26T'),resolution='220kb',superclone_addition=15)
runCountReads_amethyst(obj=obj,sample_name=c('ECIS36T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('ECIS48T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('ECIS57T'),resolution='220kb') #41 cells
```

# Read all CopyKit RDS objects and plot together
Assign aneuploid and diploid clones, and subclones per sample

```R

copykit_output<-list.files(path=paste0(project_data_directory,"/copykit/"),recursive=TRUE,full.names=TRUE,pattern=".220kb.rds")
#remove diploid cell call rds used for bin correction
copykit_output<-copykit_output[!grepl(copykit_output,pattern="diploid")]

#make clones a named list to collapse overclustering or low cell counts/cluster

assign_copykit_aneuploid_clonename<-function(sample_name,cancer_clones,split_on="superclones"){
    tmp<-readRDS(paste0("/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/copykit/",sample_name,"/copykit.",sample_name,".220kb.rds"))
    tmp@colData$ploidy<-"diploid"

    if(length(cancer_clones)>0){
        if(split_on=="subclones"){
        tmp@colData$clones_split<-"subclones"
        tmp@colData[tmp@colData$subclones %in% cancer_clones,]$ploidy<-"aneuploid"
        tmp@colData$clonename<-unlist(paste(sample_name,names(cancer_clones[match(tmp@colData$subclones,cancer_clones)]),sep="_"))
        tmp@colData$clonename<-gsub("_NA", replacement = "_diploid", x = tmp@colData$clonename)
        tmp@colData[tmp@colData$ploidy=="aneuploid",]$fine_celltype<-"cancer"
        }else{
        tmp@colData$clones_split<-"superclones"
        tmp@colData[tmp@colData$superclones %in% cancer_clones,]$ploidy<-"aneuploid"
        tmp@colData$clonename<-unlist(paste(sample_name,names(cancer_clones[match(tmp@colData$superclones,cancer_clones)]),sep="_"))
        tmp@colData$clonename<-gsub("_NA", replacement = "_diploid", x = tmp@colData$clonename)
        tmp@colData[tmp@colData$ploidy=="aneuploid",]$fine_celltype<-"cancer"
        }} else {
        tmp@colData$clones_split<-"all_diploid"
        tmp@colData$clonename<-paste(sample_name,"diploid",sep="_")

        }
    #define colors based on data
    #updated to be -3 to 3 instead of -2 to 2
    log_col=colorRamp2(c(-1,-0.5,0,0.5,1), 
                            c("darkblue","plum","white","tomato","darkred"))
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
        celltype=tmp@colData$fine_celltype,
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

    arm_col=c("p"="grey","q"="darkgrey")
    band_col=c("acen"="#99746F","gneg"="white","gpos100"="black","gpos25"="lightgrey","gpos50"="grey","gpos75"="darkgrey","gvar"="#446879")

    column_ha = HeatmapAnnotation(
        arm = tmp@rowRanges$arm,
        band = tmp@rowRanges$stain,
        col=list(arm=arm_col,band=band_col))

    plt<-Heatmap(
        t(tmp@assays@data$logr),
        left_annotation=ha,col=log_col,
        row_split=as.character(tmp@colData$clonename),
        show_column_names=FALSE,show_row_names=FALSE,
        top_annotation=column_ha,cluster_columns=FALSE,cluster_column_slices=FALSE,column_split=seqnames(tmp@rowRanges),
        name="logr")

    pdf(paste0("/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/copykit/",sample_name,"/copykit.",sample_name,".220kb.","cancerclone.pdf"),width=20)
    print(plt)
    dev.off()

    saveRDS(tmp,file=paste0("/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/copykit/",sample_name,"/copykit.",sample_name,".220kb.rds"))
}

assign_copykit_aneuploid_clonename(sample_name="BCMDCIS05T",cancer_clones=c("c1"='3',"c2"='4'))
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS07T',cancer_clones=c()) #all diploid
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS102T_24hTis',cancer_clones=c("c1"='3'))
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS124T',cancer_clones=c("c1="='2',"c2"='4'))
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS22T',cancer_clones=c("c1"='2'))
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS28T',split_on="subclones",cancer_clones=c("c1"='4',"c2"='5'))
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS32T',cancer_clones=c()) #all diploid
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS35T',cancer_clones=c("c1"='1',"c2"='2',"c3"='3',"c3"='5',"c3"='6',"c3"='7',"c3"='8'))
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS41T',cancer_clones=c('c1'='5','c1'='3','c2'='12',
                                                                            'c3'='11','c3'='15',
                                                                            'c4'='9',
                                                                            'c5'='8','c5'='1','c5'='2',
                                                                            'c6'='4','c6'='14',
                                                                            'c7'='13'))
#assign_copykit_aneuploid_clonename(sample_name='BCMDCIS49T',cancer_clones=c()) #not run
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS52T',cancer_clones=c('c1'='3','c2'='4'))
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS65T',cancer_clones=c('c1'='5')) 
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS66T',cancer_clones=c('c1'='4',
                                                                            'c2'='2',
                                                                            'c3'='10',
                                                                            'c4'='5',
                                                                            'c5'='7',
                                                                            'c6'='9',
                                                                            'c7'='8'))                                      
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS70T',split_on='subclones',cancer_clones=c('c1'='3','c1'='4','c2'='6'))

assign_copykit_aneuploid_clonename(sample_name='BCMDCIS74T',cancer_clones=c(
                                                                            'c1'='9','c1'='5',
                                                                            'c2'='7','c2'='7','c2'='3',
                                                                            'c3'='12','c3'='11',
                                                                            'c4'='8',
                                                                            'c5'='1',
                                                                            'c6'='10'))

assign_copykit_aneuploid_clonename(sample_name='BCMDCIS79T_24hTis_DCIS',split_on='superclones',cancer_clones=c('c1'='6','c1'='4','c2'='2'))
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS80T_24hTis',cancer_clones=c('c1'='4'))
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS82T_24hTis',cancer_clones=c('c1'='3','c2'='4','c2'='2'))
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS92T_24hTis',cancer_clones=c('c1'='2')) 
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS94T_24hTis',split_on="subclones",cancer_clones=c('c1'='4','c1'='3'))
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS97T',split_on="superclones",cancer_clones=c(
                                                                            'c1'='14','c1'='7','c1'='16','c1'='4','c1'='15',
                                                                            'c2'='9','c2'='5',
                                                                            'c3'='10','c3'='2',
                                                                            'c4'='8','c4'='13','c4'='12')) 
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS99T',cancer_clones=c('c1'='4')) 

assign_copykit_aneuploid_clonename(sample_name='BCMHBCA03R',split_on="subclones",cancer_clones=c('c1'='5'))
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA04R',cancer_clones=c('c1'='5'))
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA09R-3h',cancer_clones=c())
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA12R-3h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA16R-3h',cancer_clones=c())
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA17R-3h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA19R-4h',cancer_clones=c())
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA22R-4h',cancer_clones=c())
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA26L-24hTis-4h',split_on="subclones",cancer_clones=c('c1'='6')) 
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA29L-2h',cancer_clones=c())
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA38L-3h',cancer_clones=c())
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA83L-3h',cancer_clones=c('c1'='4')) 
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA85L-3h',cancer_clones=c())

assign_copykit_aneuploid_clonename(sample_name='ECIS25T',split_on="subclones",cancer_clones=c('c1'='5',
                                                                                            'c2'='1','c2'='4',
                                                                                            'c3'='6','c3'='7','c3'='2'))


assign_copykit_aneuploid_clonename(sample_name='ECIS36T',split_on="subclones",cancer_clones=c('c1'='1','c2'='3','c3'='4','c1'='2','c3'='5')) #might be some more cancer cells with low read count in diploid pop
assign_copykit_aneuploid_clonename(sample_name='ECIS48T',cancer_clones=c('c1'='6','c2'='7')) #might have a cancer precursor in the diploid pop chr16 loss in some lumhr, 1q gain in lumsec?
#assign_copykit_aneuploid_clonename(sample_name='ECIS57T',cancer_clones=c()) #not run cell count too low

assign_copykit_aneuploid_clonename(sample_name='ECIS26T',split_on="superclones",cancer_clones=c('c1'='18','c1'='16','c1'='5','c1'='6','c1'='9','c1'='13','c1'='7','c1'='17','c1'='10','c1'='3','c1'='1','c1'='8',
'c2'='11','c2'='4','c2'='12',
'c3'='14'))

#read in all meta data from copykit, append to amethyst object
read_meta_copykit<-function(x){
    tmp<-readRDS(x)
    meta<-as.data.frame(tmp@colData[c("sample_name","reads_assigned_bins","plate_info","superclones","subclones","ploidy","clonename","fine_celltype","clones_split")])
    return(meta)
}
cnv_meta<-do.call("rbind",lapply(copykit_output,read_meta_copykit))

obj_backup<-obj
obj@metadata$cnv_ploidy<-"NA" #na values are too low read count to call
obj@metadata[row.names(cnv_meta),]$cnv_ploidy<-cnv_meta$ploidy
obj@metadata$cnv_superclones<-"NA"
obj@metadata[row.names(cnv_meta),]$cnv_superclones<-cnv_meta$superclones
obj@metadata$cnv_subclones<-"NA"
obj@metadata[row.names(cnv_meta),]$cnv_subclones<-cnv_meta$subclones
obj@metadata$cnv_clonename<-"NA"
obj@metadata[row.names(cnv_meta),]$cnv_clonename<-cnv_meta$clonename
obj@metadata$cnv_clones_split<-"NA"
obj@metadata[row.names(cnv_meta),]$cnv_clones_split<-cnv_meta$clones_split

obj@metadata[which(obj@metadata$cnv_ploidy=="aneuploid"),]$fine_celltype<-"cancer"
saveRDS(obj,file="06_scaledcis.cnv_clones.amethyst.rds")

```

Plot all clones together (with clustering for shared cross-patient cnvs)

```R
library(ComplexHeatmap)

copykit_output<-list.files(path=paste0(project_data_directory,"/copykit"),recursive=TRUE,full.names=TRUE,pattern="*rds")
copykit_output<-copykit_output[!grepl(copykit_output,pattern="diploid")]
cna_obj<-readRDS(copykit_output[1]) #just to grab row ranges
output_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/copykit/"
#read in all meta data from copykit
read_meta_copykit<-function(x){
    tmp<-readRDS(x)
    meta<-as.data.frame(tmp@colData[c("sample_name","reads_assigned_bins","plate_info","subclones","fine_celltype","clonename","ploidy")])
    return(meta)
}
cnv_meta<-do.call("rbind",lapply(copykit_output,read_meta_copykit))

#read in all logr from copykit
read_logr_copykit<-function(x){
    tmp<-readRDS(x)
    logr<-tmp@assays@data$logr
    return(logr)
}
cnv_logr<-do.call("cbind",lapply(copykit_output,read_logr_copykit))

#get 220kb windows ranges
copykit<-readRDS(copykit_output[1])
windows<-copykit@rowRanges

#relevant CNV genes from curtis work
#from https://www.nature.com/articles/s41416-024-02804-6#Sec20
#change RAB7L1 to RAB29
#lost RAB7L1
cnv_genes<-c('ESR1','PGR','DLEU2L', 'TRIM46', 'FASLG', 'KDM5B', 'RAB7L1', 'PFN2', 'PIK3CA', 'EREG', 'AIM1', 'EGFR', 'ZNF703', 'MYC', 'SEPHS1', 'ZMIZ1', 'EHF', 'POLD4', 'CCND1', 'P2RY2', 'NDUFC2-KCTD14', 'FOXM1', 'MDM2', 'STOML3', 'NEMF', 'IGF1R', 'TP53I13', 'ERBB2', 'SGCA', 'RPS6KB1', 'BIRC5', 'NOTCH3', 'CCNE1', 'RCN3', 'SEMG1', 'ZNF217', 'TPD52L2', 'PCNT', 'CDKN2AIP', 'LZTS1', 'PPP2R2A', 'CDKN2A', 'PTEN', 'RB1', 'CAPN3', 'CDH1', 'MAP2K4', 'GJC2', 'TERT', 'RAD21', 'ST3GAL1', 'SOCS1')
cnv_genes_class<-c('amp','amp','amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'del', 'del', 'del', 'del', 'del', 'del', 'del', 'del', 'del', 'amp', 'amp', 'amp', 'amp', 'amp')
cnv_genes<-setNames(cnv_genes_class,nm=cnv_genes)

#use gtf file to get gene locations
gtf_file="/container_ref/gencode.v43.annotation.gtf.gz"

gtf <- rtracklayer::readGFF(gtf_file)
gtf<- gtf %>% 
    filter(type=="gene" & gene_type %in% c("protein_coding")) %>% 
    filter(gene_name %in% names(cnv_genes))

cnv_genes_windows<-gtf[gtf$gene_name %in% names(cnv_genes),] #filter annotation to genes we want
cnv_genes_windows<-cnv_genes_windows[!duplicated(cnv_genes_windows$gene_name),] #remove duplicates
cnv_genes_windows$cnv_gene_class<-unname(cnv_genes[match(cnv_genes_windows$gene_name, names(cnv_genes))]) #add amp/del

cnv_genes_windows<-makeGRangesFromDataFrame(cnv_genes_windows,keep.extra.columns=TRUE) #make granges
wind<-GenomicRanges::findOverlaps(windows,cnv_genes_windows) #do overlap to get window indexes
wind<-as.data.frame(wind)
wind<-wind[!duplicated(wind$subjectHits),] #take first subject hit

annot<-data.frame(
  window_loc=wind$queryHits,
  gene=cnv_genes_windows$gene_name[wind$subjectHits],
  cnv_class=cnv_genes_windows$cnv_gene_class[wind$subjectHits])

annot$col<-ifelse(annot$cnv_class=="amp","red","blue")

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

hc = columnAnnotation(common_cnv = anno_mark(at = annot$window_loc, 
                        labels = annot$gene,
                        which="column",side="bottom",
                        labels_gp=gpar(col=annot$col)))

#define colors based on data
log_col=colorRamp2(c(-2,-1,0,1,2),
                        c("darkblue","blue","white","red","darkred"))

```

Only cancer

```R
aneu<-obj@metadata %>% filter(cnv_ploidy=="aneuploid")
aneu_logr<-cnv_logr[row.names(aneu)]
aneu_meta<-cnv_meta[row.names(aneu),]

#plot heatmap
ha = rowAnnotation(
    sample=aneu_meta$sample,
    celltype=aneu_meta$fine_celltype,
    ploidy=aneu_meta$ploidy,
    clones=aneu_meta$clonenames,
    col= list(
        #sample=sample_col,
        celltype=celltype_col
        #ploidy=ploidy_col,
        #clones=clone_names_col
    ))

pdf(paste0(output_directory,"/","aneuploid.cnv.heatmap.pdf"),height=90,width=40)
Heatmap(t(aneu_logr),
  col=log_col,
  cluster_columns=FALSE,
  cluster_rows=TRUE,
  show_row_names = FALSE, row_title_rot = 0,
  show_column_names = FALSE,
  cluster_row_slices = TRUE,
  bottom_annotation=hc,
  top_annotation=column_ha,
  left_annotation=ha,
  row_split=aneu_meta$clonename,
  column_split=seqnames(windows),
  border = TRUE)
dev.off()
print(paste0(output_directory,"/","aneuploid.cnv.heatmap.pdf"))

```

```R
diploid<-obj@metadata %>% filter(cnv_ploidy=="diploid") 
diploid_logr<-cnv_logr[row.names(diploid)]
diploid_meta<-cnv_meta[row.names(diploid),]

#plot heatmap
ha = rowAnnotation(
    sample=diploid_meta$sample,
    celltype=diploid_meta$fine_celltype,
    ploidy=diploid_meta$ploidy,
    clones=diploid_meta$clonenames,
    col= list(
        #sample=sample_col,
        celltype=celltype_col
        #ploidy=ploidy_col,
        #clones=clone_names_col
    ))

pdf(paste0(output_directory,"/","diploid.cnv.heatmap.pdf"),height=90,width=40)
Heatmap(t(diploid_logr),
  col=log_col,
  cluster_columns=FALSE,
  cluster_rows=TRUE,
  show_row_names = FALSE, row_title_rot = 0,
  show_column_names = FALSE,
  cluster_row_slices = TRUE,
  bottom_annotation=hc,
  top_annotation=column_ha,
  left_annotation=ha,
  row_split=diploid_meta$clonename,
  column_split=seqnames(windows),
  border = TRUE)
dev.off()
print(paste0(output_directory,"/","diploid.cnv.heatmap.pdf"))

```


#Run PCA on aneuploid cells plot feature loadings over genome

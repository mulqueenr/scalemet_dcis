
Estimate integer copy number for clones

Using scquantum, if scquantum fails (estimates ploidy below 2, use segment ratio mean)

```R
library(GenomicRanges)
library(copykit)
library(ComplexHeatmap)
library(parallel)
library(BiocParallel)
set.seed(111)
library(dplyr)
library(data.table)
library(scquantum)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dendextend)
library(amethyst)
task_cpus=150
seed=1234
register(MulticoreParam(progressbar = T, workers = task_cpus), default = T)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing

#set environment and read in data
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)


hg38_grangeslist[["hg38_500kb"]]<-readRDS(file=paste0("/data/rmulqueen/projects/scalebio_dcis/ref/copykit.met_windows.500kb.diploidcorrected.ref.rds")) #4107
hg38_grangeslist[["hg38_500kb"]]<-hg38_grangeslist[["hg38_500kb"]][
    which(
        hg38_grangeslist[["hg38_500kb"]]$diploid_cov < mean(hg38_grangeslist[["hg38_500kb"]]$diploid_cov)+(1.5*sd(hg38_grangeslist[["hg38_500kb"]]$diploid_cov)) &
        hg38_grangeslist[["hg38_500kb"]]$diploid_cov > mean(hg38_grangeslist[["hg38_500kb"]]$diploid_cov)-(1.5*sd(hg38_grangeslist[["hg38_500kb"]]$diploid_cov))),]
copykit_output_500kb <- list.files(path=paste0(project_data_directory,"/copykit/"),recursive=TRUE,full.names=TRUE,pattern=".500kb.rds")

#remove diploid cell call rds used for bin correction
copykit_output_500kb <- copykit_output_500kb[!grepl(copykit_output_500kb,pattern="diploid")]

#merge all copykit and run at same time (to help with integer estimation)
merged_copykit <- do.call("cbind",lapply(copykit_output_500kb,function(x){
    obj<-readRDS(x)
    obj@colData<-obj@colData[c("unique_reads","tss_enrich","mcg_pct","cg_cov","batch","plate_info",
                            "tgmt_well","i7_well",
                            "i5_well","fine_celltype",
                            "sample","sample_name",
                            "reads_assigned_bins","overdispersion",
                            "breakpoint_count","subclones",
                            "superclones","ploidy",
                            "clones_split","clonename",
                            "celltype")]
    obj@assays@data<-obj@assays@data[c("bincounts","smoothed_bincounts","ratios","ft","logr")]
return(obj)}))


#filter to just aneuploid, or not
merged_copykit<-merged_copykit[,!endsWith(colData(merged_copykit)$clonename,"_diploid")]

############ Segmentation ############

    #modified to require less regions width, since im using bigger windows
    runSegmentation <- function(scCNA,
                                method = c("CBS", "multipcf"),
                                seed = 17,
                                alpha = 1e-5,
                                merge_levels_alpha = 1e-5,
                                gamma = 40,
                                undo.splits = "prune",
                                name = "segment_ratios",
                                BPPARAM = bpparam()) {
        # Args
        method <- match.arg(method)

        # bindings for NSE and data
        chr <- arm <- chrarm <- NULL

        # checks
        if (!is.numeric(alpha)) stop("Argument alpha must be numeric.")
        if (!is.numeric(gamma)) stop("Argument gamma must be numeric.")
        if (!is.numeric(seed)) stop("Argument seed must be numeric.")

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Genome Assembly
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # adding genome assembly info to metadata
        if (S4Vectors::metadata(scCNA)$genome == "hg19") {
            genome <- "hg19"
        }

        if (S4Vectors::metadata(scCNA)$genome == "hg38") {
            genome <- "hg38"
        }

        # Reading hg38 VarBin ranges
        if (genome == "hg38") {
            resolution <- S4Vectors::metadata(scCNA)$resolution

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

            hg38_rg <- as.data.frame(hg38_rg) %>%
                dplyr::rename(chr = "seqnames")

            hg38_rg_mod <- hg38_rg
            # match for chrY presence
            chr_sccna <-
                as.character(as.data.frame(rowRanges(scCNA))$seqnames)
            hg38_rg_mod <-
                hg38_rg_mod[which(hg38_rg_mod$chr %in% chr_sccna), ]

            hg38_rg_mod <- hg38_rg_mod %>%
                dplyr::mutate(
                    chr = gsub("X", "23", chr),
                    chr = gsub("Y", "24", chr)
                )

            chr_info <-
                as.numeric(gsub("chr", "", hg38_rg_mod$chr))

            ref <- hg38_rg_mod
        }

        # reading hg19 varbin ranges
        if (genome == "hg19") {
            hg19_rg_mod <- hg19_rg
            # match for chrY presence
            chr_sccna <-
                as.character(as.data.frame(SummarizedExperiment::rowRanges(scCNA))$seqnames)
            hg19_rg_mod <-
                hg19_rg_mod[which(hg19_rg_mod$chr %in% chr_sccna), ]

            hg19_rg_mod <- hg19_rg_mod %>%
                dplyr::mutate(
                    chr = gsub("X", "23", chr),
                    chr = gsub("Y", "24", chr)
                )

            chr_info <-
                as.numeric(gsub("chr", "", hg19_rg_mod$chr))

            ref <- hg19_rg_mod
        }

        ref_chrarm <- ref %>%
            dplyr::mutate(chrarm = paste0(gsub("chr", "", chr), arm))

        levels_chrarm <- gtools::mixedsort(unique(ref_chrarm$chrarm))

        ref_chrarm <- ref_chrarm %>%
            dplyr::mutate(chrarm = as.factor(chrarm)) %>%
            dplyr::mutate(chrarm = forcats::fct_relevel(chrarm, levels_chrarm))

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Data Setup
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (S4Vectors::metadata(scCNA)$vst == "ft") {
            counts_df <- SummarizedExperiment::assay(scCNA, "ft")
        }

        if (S4Vectors::metadata(scCNA)$vst == "log") {
            counts_df <- SummarizedExperiment::assay(scCNA, "log")
        }

        # smoothing data
        message("Smoothing outlier bins.")
        smooth_counts <-
            BiocParallel::bplapply(as.data.frame(counts_df), function(x) {
                CNA_object <-
                    DNAcopy::CNA(x,
                        ref_chrarm$chrarm,
                        ref$start,
                        data.type = "logratio",
                        sampleid = names(x)
                    )
                withr::with_seed(seed,
                                smoothed_CNA_counts <- DNAcopy::smooth.CNA(CNA_object)[, 3]

                )

            }, BPPARAM = BPPARAM)

        smooth_counts_df <- dplyr::bind_cols(smooth_counts) %>%
            as.data.frame() %>%
            round(2)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Segmentation methods
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        message(
            "Running segmentation algorithm: ",
            method,
            " for genome ",
            genome
        )

        if (method == "CBS") {
            seg_list <-
                BiocParallel::bplapply(
                    smooth_counts_df,
                    FUN = function(x) {
                        CNA_object <-
                            DNAcopy::CNA(x,
                                ref_chrarm$chrarm,
                                ref$start,
                                data.type = "logratio",
                                sampleid = names(x)
                            )

                        withr::with_seed(seed,
                                        segment_smoothed_CNA_object <-
                                            .quiet(
                                                DNAcopy::segment(
                                                    CNA_object,
                                                    alpha = alpha,
                                                    min.width = 2,
                                                    undo.splits = undo.splits
                                                )
                                            )
                                        )


                        short_cbs <- segment_smoothed_CNA_object[[2]]
                        log_seg_mean_LOWESS <-
                            rep(short_cbs$seg.mean, short_cbs$num.mark)
                    },
                    BPPARAM = BPPARAM
                )

            seg_df <- dplyr::bind_cols(seg_list) %>%
                as.data.frame() %>%
                round(2)
        }



        if (method == "multipcf") {
            smooth_multipcf <- cbind(
                as.numeric(gsub("chr", "", ref_chrarm$chr)),
                ref_chrarm$start,
                smooth_counts_df
            )

            mpcf <- .multipcf(smooth_multipcf,
                                        gamma = gamma,
                                        arms = vapply(
                                            regmatches(ref_chrarm$chrarm,
                                                        regexec("[pq]",
                                                                ref_chrarm$chrarm)),
                                            FUN =  "[",
                                            1,
                                            FUN.VALUE = character(1)
                                        )
            )

            seg_df <- apply(mpcf[, 6:ncol(mpcf)], 2, function(x) {
                rep.int(x, mpcf$n.probes)
            })

            seg_df <- round(as.data.frame(seg_df), 2)
        }

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # merge levels
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        message("Merging levels.")

        if (S4Vectors::metadata(scCNA)$vst == "ft") {
            smooth_counts_df[smooth_counts_df == 0] <- 1e-4
            seg_df[seg_df == 0] <- 1e-4
            smooth_counts_df <- log2(smooth_counts_df)
            seg_df <- log2(seg_df)
        }

        seg_ml_list <- BiocParallel::bplapply(seq_along(seg_df), function(i) {
            cell_name <- names(seg_df)[i]
            smoothed_cell_ct <- smooth_counts_df[, i]
            seg_means_cell <- seg_df[, i]
            seg_means_ml <- mergeLevels(smoothed_cell_ct,
                seg_means_cell,
                verbose = 0,
                pv.thres = merge_levels_alpha
            )$vecMerged
        })

        names(seg_ml_list) <- names(seg_df)
        seg_ml_df <- dplyr::bind_cols(seg_ml_list)

        if (S4Vectors::metadata(scCNA)$vst == "ft") {
            seg_ml_df <- round(2^seg_ml_df, 2)
        }

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # reverting the transformation back to ratios
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (S4Vectors::metadata(scCNA)$vst == "ft") {
            SummarizedExperiment::assay(scCNA, "smoothed_bincounts") <-
                .invft(2^smooth_counts_df)

            seg_ratio_df <- .invft(seg_ml_df)
        }

        if (S4Vectors::metadata(scCNA)$vst == "log") {
            SummarizedExperiment::assay(scCNA, "smoothed_bincounts") <-
                2^smooth_counts_df

            seg_ratio_df <- 2^seg_ml_df
        }

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # saving information to the scCNA object
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # saving as segment ratios
        seg_ratios <- sweep(seg_ratio_df, 2, apply(seg_ratio_df, 2, mean), "/")
        SummarizedExperiment::assay(scCNA, name) <- round(seg_ratios, 2)


        # calculating ratios from the bincounts, used for ratio plots
        scCNA <- calcRatios(scCNA, assay = "smoothed_bincounts")

        message("Done.")

        return(scCNA)

    }

    prefix="all_samples.aneuploid_only"
    output_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/copykit"

    merged_copykit<-runSegmentation(merged_copykit,
        method = "CBS",
        undo.splits = "none",
        name = "segment_ratios")

    saveRDS(merged_copykit,file=paste0(output_directory,"/","all_samples.aneuploid.preintegerprocessing.copykit.Rds"))
merged_copykit<-readRDS(file=paste0(output_directory,"/","all_samples.aneuploid.preintegerprocessing.copykit.Rds"))
    scCNA<-merged_copykit

    scCNA <- calcConsensus(scCNA,
                            consensus_by="clonename",
                            assay='segment_ratios',
                            fun="median")

############ Estimating ploidy by segment ratios ############

    #assign ploidy by mean ratios, old way to do it
    #changing this ploidy assignment to scquantum to account for larger deviation from diploid
    est_ploidy_clone<-setNames(nm=names(colMeans(scCNA@consensus)),colMeans(scCNA@consensus)*2)
    SummarizedExperiment::colData(scCNA)$estimated_ratio_ploidy<-unname(est_ploidy_clone[SummarizedExperiment::colData(scCNA)$clonename])
    SummarizedExperiment::colData(scCNA)$ploidy<-SummarizedExperiment::colData(scCNA)$estimated_ratio_ploidy

    #calc integer using segment ratio values
    scCNA<-calcInteger(scCNA,
                        assay = "segment_ratios",
                        method = "metadata",
                        name = "integer",
                        penalty = 15)
    SummarizedExperiment::assay(scCNA, 'singlecell_integer_by_estimated_ratio')<-SummarizedExperiment::assay(scCNA, 'integer')

    #consensus to per clone slot (integer)
    scCNA <- calcConsensus(scCNA,consensus_by="clonename",assay='singlecell_integer_by_estimated_ratio',fun="median")

    consensus_integer<-mclapply(1:ncol(scCNA@assays@data$singlecell_integer_by_estimated_ratio),function(i){
        cellid<-colnames(scCNA@assays@data$singlecell_integer_by_estimated_ratio)[i]
        cloneid<-colData(scCNA)[cellid,]$clonename
        print(paste(i,cellid,cloneid))
        return(scCNA@consensus[[cloneid]])
    },mc.cores=100)

    consensus_integer<-as.data.frame(do.call("cbind",consensus_integer))
    colnames(consensus_integer)<-colnames(scCNA@assays@data$singlecell_integer_by_estimated_ratio)
    SummarizedExperiment::assay(scCNA, 'consensus_integer_by_estimated_ratio') <- consensus_integer

    SummarizedExperiment::assay(scCNA, 'consensus_integer_by_estimated_ratio_discrete')<-SummarizedExperiment::assay(scCNA, 'consensus_integer_by_estimated_ratio')
    SummarizedExperiment::assay(scCNA, 'consensus_integer_by_estimated_ratio_discrete') <- as.data.frame(apply(SummarizedExperiment::assay(scCNA, 'consensus_integer_by_estimated_ratio_discrete'), 2, as.integer))
    SummarizedExperiment::assay(scCNA, 'consensus_integer_by_estimated_ratio_discrete')[which(SummarizedExperiment::assay(scCNA, 'consensus_integer_by_estimated_ratio_discrete')>=6,arr.ind=T)]<-6
    table(unlist(SummarizedExperiment::assay(scCNA, 'consensus_integer_by_estimated_ratio_discrete')),useNA="ifany")

    SummarizedExperiment::assay(scCNA, 'singlecell_integer_by_estimated_ratio_discrete')<-SummarizedExperiment::assay(scCNA, 'singlecell_integer_by_estimated_ratio')
    SummarizedExperiment::assay(scCNA, 'singlecell_integer_by_estimated_ratio_discrete')[which(SummarizedExperiment::assay(scCNA, 'singlecell_integer_by_estimated_ratio')>=6,arr.ind=T)]<-6
    table(t(SummarizedExperiment::assay(scCNA, "singlecell_integer_by_estimated_ratio_discrete")),useNA="ifany")

    saveRDS(scCNA,file=paste0(output_directory,"/",prefix,".copykit.Rds"))

############################################################
###############Estimating ploidy by scQuantum ##############
############################################################

    #check to ensure ploidy
    rg <- as.data.frame(SummarizedExperiment::rowRanges(scCNA))
    seg <- SummarizedExperiment::assay(scCNA, "segment_ratios")
    bin <- SummarizedExperiment::assay(scCNA, 'smoothed_bincounts')

    sc_quants <- mclapply(seq_along(seg), function(z) {
        # extracting segments rle id and lengths
        segnums <- cumsum(c(TRUE, abs(diff(seg[, z])) > 0.00001))
        seg_length <- rle(seg[, z])$lengths

        # extracting segment-wise means and index of dispersion
        seg_bins_mean <- tapply(bin[, z], segnums, mean)

        if (any(seg_length <= 3)) {
            iod.est <- timeseries.iod(bin[,z])
        } else {
            iod.est <- tapply(bin[, z], segnums, timeseries.iod)
        }
        # bincount mean estimate
        mean.est <- mean(bin[, z])
        estimates <- scquantum::ploidy.inference(
            x = seg_bins_mean,
            chrom = NULL,
            start = NULL,
            end = NULL,
            seg_length = seg_length,
            iod = iod.est,
            mean_bincount = mean.est,
            do_segmentation = FALSE)

        },mc.cores=100)

    sc_ploidies <- vapply(sc_quants, function(x) x$ploidy, numeric(1)) # extracting ploidies from sc_quantum object
    sc_confidence <- vapply(sc_quants, function(x) x$confidence_ratio, numeric(1)) # extracting ploidies from sc_quantum object
    ploidy_score <- abs(1-sc_confidence) # calculating ploidy score from scquantum confidence ratio

    SummarizedExperiment::colData(scCNA)$scquantum_ploidy <- sc_ploidies
    SummarizedExperiment::colData(scCNA)$scquantum_confidence_ratio <- sc_confidence
    SummarizedExperiment::colData(scCNA)$scquantum_ploidy_score <- ploidy_score

#consensus sc_quantum ploidy per clone
scquant_consensus <- SummarizedExperiment::colData(scCNA) %>% 
                    as.data.frame() %>% 
                    group_by(clonename) %>% 
                    summarize(consensus_scquantum_ploidy=median(scquantum_ploidy))
                    
scquant_consensus <- setNames(nm=scquant_consensus$clonename,scquant_consensus$consensus_scquantum_ploidy)

#set a floor, so if ploidy < 2, then use estimated ratio ploidy
scquant_consensus[scquant_consensus<2] <- est_ploidy_clone[names(scquant_consensus[scquant_consensus<2])]

SummarizedExperiment::colData(scCNA)$scquantum_ploidy_cloneconsensus <- scquant_consensus[SummarizedExperiment::colData(scCNA)$clonename]

SummarizedExperiment::colData(scCNA)$ploidy<-SummarizedExperiment::colData(scCNA)$scquantum_ploidy_cloneconsensus #calcinteger by metadata requires it be called "ploidy"

scCNA<-calcInteger(scCNA,
                    assay = "segment_ratios",
                    method = "metadata",
                    name = "scquantum_integer",
                    penalty = 15)

#consensus to per clone slot (integer)
scCNA <- calcConsensus(scCNA,consensus_by="clonename",assay='scquantum_integer',fun="median")

consensus_integer<-mclapply(1:ncol(scCNA@assays@data$scquantum_integer),function(i){
    cellid<-colnames(scCNA@assays@data$scquantum_integer)[i]
    cloneid<-colData(scCNA)[cellid,]$clonename
    print(paste(i,cellid,cloneid))
    return(scCNA@consensus[[cloneid]])
},mc.cores=100)

consensus_integer<-as.data.frame(do.call("cbind",consensus_integer))
colnames(consensus_integer)<-colnames(scCNA@assays@data$scquantum_integer)

SummarizedExperiment::assay(scCNA, 'scquantum_singlecell_integer')<-SummarizedExperiment::assay(scCNA, 'scquantum_integer')
SummarizedExperiment::assay(scCNA, 'scquantum_singlecell_integer_discrete')<-SummarizedExperiment::assay(scCNA, 'scquantum_singlecell_integer')
SummarizedExperiment::assay(scCNA, 'scquantum_singlecell_integer_discrete')[which(SummarizedExperiment::assay(scCNA, 'scquantum_singlecell_integer')>=6,arr.ind=T)]<-6

SummarizedExperiment::assay(scCNA, 'scquantum_consensus_integer')<-consensus_integer
SummarizedExperiment::assay(scCNA, 'scquantum_consensus_integer_discrete')<-SummarizedExperiment::assay(scCNA, 'scquantum_consensus_integer')
SummarizedExperiment::assay(scCNA, 'scquantum_consensus_integer_discrete') <- as.data.frame(apply(SummarizedExperiment::assay(scCNA, 'scquantum_consensus_integer_discrete'), 2, as.integer))
SummarizedExperiment::assay(scCNA, 'scquantum_consensus_integer_discrete')[which(SummarizedExperiment::assay(scCNA, 'scquantum_consensus_integer_discrete')>=6,arr.ind=T)]<-6
table(unlist(SummarizedExperiment::assay(scCNA, 'scquantum_consensus_integer_discrete')),useNA="ifany")
#######################################
############# Plotting #############
#######################################

#read in cyto info
cyto=read.table(file="/data/rmulqueen/projects/scalebio_dcis/ref/cytoBand.txt",sep="\t")
colnames(cyto)<-c("chr","start","end","band","stain")
cyto$arm<-substring(cyto$band, 1, 1)
cyto<-cyto[!is.na(cyto$band),]
cyto<-cyto[cyto$chr %in% c(paste0("chr",1:22),"chrX"),]
table(cyto$stain) #set colors for these
cyto_overlap<-GenomicRanges::findOverlaps(scCNA@rowRanges,
                                            makeGRangesFromDataFrame(cyto,keep=TRUE),
                                            select="first")
scCNA@rowRanges$stain <- cyto[cyto_overlap,]$stain
scCNA@rowRanges$arm <- cyto[cyto_overlap,]$arm

arm_col=c("p"="grey","q"="darkgrey")
band_col=c("acen"="#99746F","gneg"="white","gpos100"="black","gpos25"="lightgrey","gpos50"="grey","gpos75"="darkgrey","gvar"="#446879")
dip_cov=colorRamp2(c(0.5,1,1.5), 
                        c("white","grey","black"))

column_ha = HeatmapAnnotation(
    mappability=scCNA@rowRanges$diploid_cov,
    arm = scCNA@rowRanges$arm,
    band = scCNA@rowRanges$stain,
    col=list(mappability=dip_cov,arm=arm_col,band=band_col))

#updated to be -4 to 4 instead of -2 to 2
seg_ratio_col=colorRamp2(c(0,0.5,1,1.5,2),
                        c("darkblue","blue","white","red","darkred"))
log_col=colorRamp2(c(-3,-2,-1,0,1,2,3), 
                        c("#053061","#2166ac","#4393c3","white","#d6604d","#b2182b","#67001f"))
int_col=c("0"="#053061","1"="#4393c3","2"="#f7f7f7","3"="#f4a582","4"="#b2182b","5"="#67001f","6"="#3d0229")
                        
cg_perc_col=colorRamp2(c(40,60,80,100),
                        c("#4d2d18","#CABA9C","#4C6444","#102820"))
reads_col=colorRamp2(c(min(log10(scCNA@colData$unique_reads)),
                        max(log10(scCNA@colData$unique_reads))),
                        c("white","black"))

superclone_col=setNames(nm=unique(as.character(scCNA@colData$superclones)),
                        colorRampPalette(brewer.pal(9, "Pastel1"))(length(unique(as.character(scCNA@colData$superclones)))))
subclone_col=setNames(nm=unique(as.character(scCNA@colData$subclones)),
                        colorRampPalette(brewer.pal(9, "Spectral"))(length(unique(as.character(scCNA@colData$subclones)))))
cancerclone_col=setNames(nm=unique(as.character(scCNA@colData$clonename)),
                        colorRampPalette(brewer.pal(8, "Pastel2"))(length(unique(as.character(scCNA@colData$clonename)))))

#set colors
celltype_col=c(
    "basal"="#87529a",
    "lumsec"="#e0b0ff",
    "lumhr"="#c500e8",
    "cancer"="#ff00ff",
    "pericyte_VSMC"="#ff6666",
    "fibroblast"="#9b1c31",
    "CAF"="#ff2222",
    "endothelial"="#ffab5f",
    "TEC"="#ffe922",
    "monocyte"="#98d3b9",
    "macrophage"="#00af5f",
    "DC"="#008080",
    "TAM"="#ccff00",
    "nk_tnk"="#00ffff",
    "tcell_cd4"="#00bae5",
    "tcell_cd8"="#1800ff",
    "tcell_cd8_2"="#0016b7",
    "bcell"="#87ceeb",
    "plasma"="#73abdb")

#make column annotations
ha = rowAnnotation(
    reads=log10(scCNA@colData$unique_reads),
    cg_perc=scCNA@colData$mcg_pct,
    celltype=scCNA@colData$celltype,
    superclones=as.character(scCNA@colData$superclones),
    subclones=as.character(scCNA@colData$subclones),
    cancerclone=as.character(scCNA@colData$clonename),
    col= list(
        celltype=celltype_col,
        reads=reads_col,
        cg_perc=cg_perc_col,
        superclones=superclone_col,
        subclones=subclone_col,
        cancerclone_col))

dend <- t(scCNA@assays@data$scquantum_consensus_integer_discrete) %>% 
        dist(method="euclidean") %>% 
        hclust(method="ward.D2") %>% 
        as.dendrogram
saveRDS(dend,file=paste0(output_directory,"/all_cells.integer.heatmap.segment_ratios.dendrogram.rds"))

#use shared clustering for both plots
print("Plotting heatmap...")
plt1<-Heatmap(t(SummarizedExperiment::assay(scCNA, "logr")),
            col=log_col,
            cluster_columns=FALSE,
            cluster_rows=dend,
            show_row_names = FALSE, row_title_rot = 0,
            show_column_names = FALSE,
            cluster_row_slices = TRUE,
            top_annotation = column_ha, left_annotation = ha,
            name="logr",
            #row_split=scCNA@colData$clonename,
            column_split=scCNA@rowRanges@seqnames,
            border = FALSE)

plt2<-Heatmap(t(SummarizedExperiment::assay(scCNA, "segment_ratios")),
            col=seg_ratio_col,
            cluster_columns=FALSE,
            cluster_rows=dend,
            show_row_names = FALSE, row_title_rot = 0,
            show_column_names = FALSE,
            cluster_row_slices = TRUE,
            top_annotation = column_ha, left_annotation = ha,
            name="Segment Ratios",
            #row_split=scCNA@colData$clonename,
            column_split=scCNA@rowRanges@seqnames,
            border = FALSE)

plt3<-Heatmap(t(SummarizedExperiment::assay(scCNA, "singlecell_integer_by_estimated_ratio_discrete")),
            col=int_col,
            cluster_columns=FALSE,
            cluster_rows=dend,
            show_row_names = FALSE, row_title_rot = 0,
            show_column_names = FALSE,
            cluster_row_slices = TRUE,
            top_annotation = column_ha, left_annotation = ha,
            name="Single cell copy number",
            #row_split=scCNA@colData$clonename,
            column_split=scCNA@rowRanges@seqnames,
            border = FALSE)

plt4<-Heatmap(t(SummarizedExperiment::assay(scCNA, "consensus_integer_by_estimated_ratio_discrete")),
            col=int_col,
            cluster_columns=FALSE,
            cluster_rows=dend,
            show_row_names = FALSE, row_title_rot = 0,
            show_column_names = FALSE,
            cluster_row_slices = TRUE,
            top_annotation = column_ha, left_annotation = ha,
            name="Consensus copy number",
            #row_split=scCNA@colData$clonename,
            column_split=scCNA@rowRanges@seqnames,
            border = FALSE)

plt5<-Heatmap(t(SummarizedExperiment::assay(scCNA, "scquantum_singlecell_integer_discrete")),
            col=int_col,
            cluster_columns=FALSE,
            cluster_rows=dend,
            show_row_names = FALSE, row_title_rot = 0,
            show_column_names = FALSE,
            cluster_row_slices = TRUE,
            top_annotation = column_ha, left_annotation = ha,
            name="Single cell copy number",
            #row_split=scCNA@colData$clonename,
            column_split=scCNA@rowRanges@seqnames,
            border = FALSE)

plt6<-Heatmap(t(SummarizedExperiment::assay(scCNA, "scquantum_consensus_integer_discrete")),
            col=int_col,
            cluster_columns=FALSE,
            cluster_rows=dend,
            show_row_names = FALSE, row_title_rot = 0,
            show_column_names = FALSE,
            cluster_row_slices = TRUE,
            top_annotation = column_ha, left_annotation = ha,
            name="Consensus copy number",
            #row_split=scCNA@colData$clonename,
            column_split=scCNA@rowRanges@seqnames,
            border = FALSE)

pdf_outname=paste0(output_directory,"/","all_cells.aneuploid.integer.heatmap.pdf")
pdf(pdf_outname,width=40,height=20)
print(plt1+plt2)
print(plt3+plt4)
print(plt5+plt6)
dev.off()

saveRDS(scCNA,file=paste0(output_directory,"/","all_samples.aneuploid.copykit.Rds"))


#plot consensus integer by clone and cell count
library(ggrepel)
obj <- readRDS(file="08_scaledcis.final_celltype.amethyst.rds")

plotting_metadata <- SummarizedExperiment::colData(scCNA) %>% as.data.frame() %>% select(scquantum_ploidy_cloneconsensus,clonename)
plotting_metadata <- plotting_metadata[!duplicated(plotting_metadata$clonename),]
row.names(plotting_metadata) <- plotting_metadata$clonename
obj_meta<-obj@metadata[!duplicated(obj@metadata$cnv_clonename),]

plt_data <- obj@metadata %>% group_by(cnv_clonename) %>% summarize(count=n(),group=first(Group)) %>% as.data.frame()
row.names(plt_data) <- plt_data$cnv_clonename

plt_data <- cbind(plt_data,plotting_metadata[row.names(plt_data),])


plt<-ggplot(plt_data, aes(x = scquantum_ploidy_cloneconsensus, y = log10(count), label = clonename, fill=group,color=group)) +
  geom_point(size=0.2,alpha=0.5)+geom_text_repel(dat=plt_data,max.overlaps=Inf,size=2) + theme_minimal()
  #geom_text_repel(size = 0.1,max.overlaps=Inf)

ggsave(plt,file=paste0(output_directory,"/","all_samples.scquantum.ploidy.pdf"),width=10,height=10)

#now add cnv integer copy as assay in amethyst object
obj <- readRDS(file="08_scaledcis.final_celltype.amethyst.rds")
obj@genomeMatrices[["aneuploid_cnv"]] <- SummarizedExperiment::assay(scCNA, "scquantum_consensus_integer_discrete")
row.names(obj@genomeMatrices[["aneuploid_cnv"]]) <- paste(as.character(seqnames(scCNA@rowRanges)),start(scCNA@rowRanges),end(scCNA@rowRanges),sep="_")
scquantum_ploidy<-setNames(nm=row.names(colData(scCNA)),colData(scCNA)$ploidy)
obj@metadata$scquantum_ploidy <- NA
obj@metadata[names(scquantum_ploidy),]$scquantum_ploidy <- scquantum_ploidy
saveRDS(obj,file="09_scaledcis.final_ploidy.amethyst.rds")



#And summarize methylation over cnv windows
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing

bed<-data.frame(chr=unlist(lapply(strsplit(row.names(obj@genomeMatrices[["aneuploid_cnv"]]),"_"),"[",1)),
                start=unlist(lapply(strsplit(row.names(obj@genomeMatrices[["aneuploid_cnv"]]),"_"),"[",2)),
                end=unlist(lapply(strsplit(row.names(obj@genomeMatrices[["aneuploid_cnv"]]),"_"),"[",3)))

win_out <- makeWindows(obj, 
                    bed = bed,
                    type = "CG", 
                    metric = "percent", 
                    threads = 50, 
                    index = "chr_cg", 
                    nmin = 2) 
obj@genomeMatrices[["cg_cnv_segments"]]<-win_out
saveRDS(obj,file="09_scaledcis.final_ploidy.amethyst.rds")

        
```




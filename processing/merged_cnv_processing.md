
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis/ /data/rmulqueen/projects/scalebio_dcis/singularity/copykit.sif

#maybe just move sc_bams to new folder per sample instead?
```R
library(copykit)
library(BiocParallel)
library(optparse)
library(circlize)
library(dendextend,lib.loc="/home/rmulqueen/R/x86_64-pc-linux-gnu-library/4.3") #add this
library(cluster)
library(ComplexHeatmap)
library(ggplot2)
BiocParallel::bpparam()
library(plyranges,lib.loc="/home/rmulqueen/R/x86_64-pc-linux-gnu-library/4.3")


#new cnv calling > put bed files into a tmp folder to run

option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=".", 
              help="List of single-cell bam files", metavar="character"),
  make_option(c("-p", "--output_prefix"), type="character", default="allcells", 
              help="Prefix of output"),
  make_option(c("-c", "--task_cpus"), type="integer", default=300, 
              help="Integer number of cpus")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
cpu_count=opt$task_cpus
prefix=opt$output_prefix
system("mkdir -p merged_cnv")
prefix="./merged_cnv/scaledcis"
register(MulticoreParam(progressbar = T, workers = cpu_count), default = T)

#mappability of bsconverted DNA (PE 100bp)
if(!exists(k100.umap.bedgraph.gz)){
    system("wget https://bismap.hoffmanlab.org/raw/hg38/k100.umap.bedgraph.gz")
}

#ref setup
genome <- "hg38"
resolution <- "200kb"
min_bincount = 10
hg38_grangeslist <- hg38_grangeslist
hg38_rg<-hg38_grangeslist[[paste0("hg38_",resolution)]]

#summarize mappability over ref bincounts
mappability<-read_bed("k100.umap.bedgraph.gz")
mappability$mapp<-as.numeric(mappability$name)

map_out<- hg38_rg %>%
  group_by_overlaps(mappability) %>%
  summarise(mapp = mean(mapp,na.rm=T))
hg38_rg$mapp<-NA
hg38_rg[map_out$query,]$mapp<-map_out$mapp
#fill in mapp with NA values of average scores
hg38_rg[is.na(hg38_rg$mapp),]$mapp<-mean(hg38_rg$mapp,na.rm=T)

hg38_rg <- as.data.frame(hg38_rg)
rg<-NULL
rg <- hg38_rg %>%
dplyr::rename(chr = "seqnames") %>%
dplyr::mutate(GeneID = 1:nrow(hg38_rg))
rg <- dplyr::filter(rg,chr != "chrY")

setwd("/data/rmulqueen/projects/scalebio_dcis/data")
merged_met<-read.csv("scaledcis.metadata.tsv",sep="\t")
col_fun = colorRamp2(c(-2,-1.5,1, 0, 1, 1.5,2), c("darkblue","blue","white", "white", "white","red","darkred"))


cnv_per_sample<-function(samp_name){

    samp<-merged_met[merged_met$sample==samp_name,]
    samp<-samp[!is.na(samp$sc_bam),]
    # bindings for NSE and data
    Chr <- chr <- strand <- GeneID <- NULL
    reads_assigned_bins <- reads_duplicates <- reads_total <- NULL

    files <-samp$sc_bam
    files_names <- row.names(samp)

    varbin_counts_list_all_fields <-BiocParallel::bplapply(
                files,Rsubread::featureCounts,
                ignoreDup = TRUE,
                countMultiMappingReads = FALSE,
                annot.ext = rg,
                useMetaFeatures = FALSE,
                verbose = TRUE,
                isPairedEnd = TRUE,
                BPPARAM = bpparam())
                #maxMOp=100, #this is added for methylation alignments
                #primaryOnly=TRUE, #this is added for methylation alignments

    names(varbin_counts_list_all_fields) <- files_names
    varbin_counts_list <- lapply(varbin_counts_list_all_fields,"[[",1)
    varbin_counts_list <- lapply(varbin_counts_list,as.vector)
    min_bc <- which(vapply(varbin_counts_list, mean, numeric(1)) < min_bincount)
    if (length(min_bc) > 0) {
    varbin_counts_list <- varbin_counts_list[-min_bc]
    varbin_counts_list_all_fields <- varbin_counts_list_all_fields[-min_bc]
    message(length(min_bc), " bam files had less than ", min_bincount," mean bincounts and were removed.")}

    # LOWESS GC normalization
    message("Performing GC correction.")
    varbin_counts_list_gccor <-BiocParallel::bplapply(varbin_counts_list, function(x) {
        gc_cor <- lowess(rg$gc_content, log(x + 1e-3), f = 0.05)
        gc_cor_z <- approx(gc_cor$x, gc_cor$y, rg$gc_content)
        exp(log(x) - gc_cor_z$y) * median(x)
        },
    BPPARAM = bpparam())


    pdf(paste0(prefix,".",samp_name,".subclone.loessGC.pdf"))
    gplots::plotLowess(gc_cor)
    dev.off()
    # LOWESS MAPPABILITY normalization (also on 0-1 scale)
    message("Performing MAPPABILITY correction.")
    varbin_counts_list_mappcor <-BiocParallel::bplapply(varbin_counts_list, function(x) {
        map_cor <- lowess(rg$mapp, log(x + 1e-3), f = 0.05)
        map_cor_z <- approx(map_cor$x, map_cor$y, rg$mapp)
        exp(log(x) - map_cor_z$y) * median(x)
        },
    BPPARAM = bpparam())

    varbin_counts_df <- round(dplyr::bind_cols(varbin_counts_list_mappcor), 2)

    # filtering low read counts where the sum of bins does not reach more than 0
    good_cells <- names(varbin_counts_df[which(colSums(varbin_counts_df) != 0)])
    varbin_counts_df <- varbin_counts_df[good_cells]
    rg <- rg %>%dplyr::select(-strand, -GeneID)
    rg_gr <- GenomicRanges::makeGRangesFromDataFrame(rg,ignore.strand = TRUE,keep.extra.columns = TRUE)
    cna_obj <- CopyKit(assays = list(bincounts = varbin_counts_df),rowRanges = rg_gr)

    # Adding genome and resolution information to metadata
    S4Vectors::metadata(cna_obj)$genome <- genome
    S4Vectors::metadata(cna_obj)$resolution <- resolution

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sun Feb 14 20:55:01 2021
    # ADDING READS METRICS TO METADATA
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sun Feb 14 20:55:24 2021

    varbin_reads_list <- lapply(varbin_counts_list_all_fields,"[[",4)

    # saving info and removing columns from list elements
    metadata_info_names <- varbin_reads_list[[1]][c(1, 2, 8, 9, 12, 14), 1]
    metadata_info_names <-
    c(
    "reads_assigned_bins",
    "reads_unmapped",
    "reads_duplicates",
    "reads_multimapped",
    "reads_unassigned",
    "reads_ambiguous"
    )

    varbin_reads_info <-lapply(seq_along(varbin_reads_list), function(x) {
    # RSubread seems to change underlines to dot on some cases
    # Have to make more complicated lapply to extract the name of the list
    # and guarantee that the cell is properly named
    name <- names(varbin_reads_list)[[x]]
    df <- varbin_reads_list[[x]][c(1, 2, 8, 9, 12, 14), -1, drop = FALSE]
    names(df) <- name
    df
    })

    names(varbin_reads_list) <- names(varbin_counts_list_all_fields)
    bam_metrics <- dplyr::bind_cols(varbin_reads_info)

    # making sure metrics match varbin_counts_df
    bam_metrics <- bam_metrics[good_cells]
    rownames(bam_metrics) <- metadata_info_names
    bam_metrics <- as.data.frame(t(bam_metrics))

    # adding total
    reads_tot <- rowSums(bam_metrics)
    bam_metrics$sample <- rownames(bam_metrics)
    bam_metrics <-
    dplyr::relocate(bam_metrics, sample, .before = reads_assigned_bins)
    bam_metrics <- bam_metrics %>%
    dplyr::mutate(reads_total = reads_tot,percentage_duplicates = round(reads_duplicates / reads_total, 3))

    # adding to metadata
    SummarizedExperiment::colData(cna_obj) <-S4Vectors::DataFrame(bam_metrics)
    colnames(cna_obj) <- names(varbin_counts_df)

    dat<-cna_obj
    dat@metadata$resolution <- "220kb"

    #dat <- dat[,colData(dat)$reads_assigned_bins >= 100000]
    colData(dat)$sample_name<-samp_name
    colData(dat)$celltype<-NA
    colData(dat)$cg_cov<-NA
    colData(dat)$mcg_pct<-NA

    colData(dat)[row.names(samp),]$celltype<-samp$celltype
    colData(dat)[row.names(samp),]$cg_cov<-samp$cg_cov
    colData(dat)[row.names(samp),]$mcg_pct<-samp$mcg_pct

    #to denoise some more increase penalties
    #alpha = 1e-20, merge_levels_alpha = 1e-20,gamma = 60
    #STD: alpha = 1e-5,merge_levels_alpha = 1e-5,gamma = 40,
    dat<-runVst(dat)
    dat<-runSegmentation(dat,method="CBS",merge_levels_alpha=1e-5,alpha=1e-5,undo.splits="prune")

    pdf(paste0(prefix,".",samp_name,".subclone.heatmap.outlier.pdf"))
       plt<-plotHeatmap(dat, label = c('reads_assigned_bins','sample_name',"mcg_pct","cg_cov","celltype"),order='hclust',n_threads=cpu_count)
       draw(plt)
    dev.off()

#col=col_fun
    pdf(paste0(prefix,".",samp_name,".subclone.heatmap.outlier.annot.pdf"))
    plt<-HeatmapAnnotation(gc = dat@rowRanges$gc_content, map =dat@rowRanges$mapp, 
    simple_anno_size = unit(2, "cm"))
    draw(plt)
    dev.off()

    dat  <- runMetrics(dat)

    pdf(paste0(prefix,".",samp_name,".qc_metrics.pdf"))
    plotMetrics(dat, metric = c("overdispersion", 
                                "breakpoint_count",
                                "reads_total",
                                "reads_duplicates",
                                "reads_assigned_bins",
                                "percentage_duplicates"),label = "reads_total")
    dev.off()

    dat <- findAneuploidCells(dat)
    dat <- findOutliers(dat)

    pdf(paste0(prefix,".",samp_name,".subclone.heatmap.outlier.pdf"))
    print(plotHeatmap(dat, row_split='outlier',n_threads=cpu_count,col=col_fun)) 
    dev.off()

    # kNN smooth profiles
    dat <- knnSmooth(dat,k=10)

    #distMat(dat) <- segment_ratios(dat) %>% t() %>% dist(method="manhattan") 
    #dend<-distMat(dat) %>% stats::hclust("ward.D2") %>% as.dendrogram()
    #dend_order<-dend %>% order()
    #find_k_opt<-find_k(dend, krange = 3:10)
    #k_clus <- cutree(dend,k = find_k_opt$nc, order_clusters_as_data = FALSE) 
    #colData(dat)$subclones<-NA
    #colData(dat)[names(k_clus),]$subclones<-factor(unname(k_clus),levels=unique(k_clus))
    
    # Create a umap embedding 
    dat <- runUmap(dat)
    dat<-findSuggestedK(dat,k_range=2:50) #10
    dat <- findClusters(dat,k_subclones=30)
    
    pdf(paste0(prefix,".",samp_name,".subclone.findk.pdf"),width=20)
    print(plotSuggestedK(dat,geom="tile"))
    dev.off()

    pdf(paste0(prefix,".",samp_name,".subclone.umap.pdf"))
    print(plotUmap(dat, label = 'subclones'))
    dev.off()

    # Calculate consensus profiles for each subclone, 
    # and order cells by cluster for visualization with plotHeatmap
    dat <- calcConsensus(dat,consensus_by="subclones")
    dat <- runConsensusPhylo(dat)
    dat <- runPhylo(dat, metric = 'manhattan')
    dat <- calcInteger(dat,assay="segment_ratios",method="scquantum")

    best_cell<-row.names(dat@colData[which(dat@colData$reads_assigned_bins==max(dat@colData$reads_assigned_bins)),])[1]

    worst_cell<-row.names(dat@colData[which(dat@colData$reads_assigned_bins==min(dat@colData$reads_assigned_bins)),])[1]

    best_cell_plot<-plotRatio(dat, best_cell)
    worse_cell_plot<-plotRatio(dat, worst_cell)

    pdf(paste0(prefix,".",samp_name,".cell.ratioplots.pdf"))
    print(best_cell_plot)
    print(worse_cell_plot)
    dev.off()

    #col_fun = colorRamp2(c(-0.5, 0, 5), c("blue", "white", "red"))
    # Plot a copy number heatmap with clustering annotation
    pdf(paste0(prefix,".",samp_name,".subclone.heatmap.segment_ratios.pdf"))
    draw(plotHeatmap(dat, 
        label = c("subclones",'reads_assigned_bins','sample_name',"celltype","mcg_pct","cg_cov"),
        assay="segment_ratios",order='hclust',
        n_threads=cpu_count,col=col_fun))
    dev.off()

    col_int_fun = colorRamp2(c(0,1,2,3), c("darkblue","white","red","darkred"))

    pdf(paste0(prefix,".",samp_name,".subclone.heatmap.integer.pdf"))
        draw(plotHeatmap(dat, label = c('reads_assigned_bins','sample_name',"celltype"),assay="integer",order='hclust',n_threads=cpu_count,col=col_int_fun))
    dev.off()

    pdf(paste0(prefix,".",samp_name,".subclone.phylo.pdf"))
    print(plotPhylo(dat, label = 'subclones'))
    dev.off()

    saveRDS(dat,file=paste0(prefix,".",samp_name,".scCNA.rds"))
    write.table(as.data.frame(dat@colData),file=paste0(prefix,".",samp_name,".scCNA.tsv"),sep="\t",col.names=T,row.names=T)
}

lapply(unique(merged_met$sample),cnv_per_sample)
```
#run line by line
#singularity shell --bind /home/rmulqueen/projects/scalebio_dcis/ ~/singularity/copykit.sif

#example run submission
#singularity exec \
#--bind /home/rmulqueen/projects/scalebio_dcis/ \
#~/singularity/copykit.sif \
#Rscript ~/projects/scalebio_dcis/tools/scalemet_dcis/src/copykit_cnvcalling.R \
#--input_dir /home/rmulqueen/projects/scalebio_dcis/data/250329_RM_scalebio_batch1_initseq/homebrew_dat/sc_bams \
#--output_prefix homemade \
#--task_cpus 150

library(copykit)
library(BiocParallel)
library(optparse)
library(circlize)
BiocParallel::bpparam()

option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=".", 
              help="Directory with single-cell bam files", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=".", 
              help="Output directory for plots and RDS file"),
  make_option(c("-p", "--output_prefix"), type="character", default="homemade", 
              help="Prefix of output"),
  make_option(c("-c", "--task_cpus"), type="integer", default=200, 
              help="Integer number of cpus")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
outdir=ifelse(endsWith(opt$output_dir,suffix="/"),opt$output_dir, paste0(opt$output_dir,"/"))
cpu_count=opt$task_cpus
prefix=opt$output_prefix

register(MulticoreParam(progressbar = T, workers = cpu_count), default = T)

dat <- runVarbin(opt$input_dir, is_paired_end = TRUE, remove_Y = TRUE)
dat  <- runMetrics(dat)

#note that "sample" is a default metadata column for the cell label, so we need to use "sample_name"
colData(dat)$sample_name<-unlist(lapply(strsplit(row.names(colData(dat)),"[.]"),"[",1))
colData(dat)$method<-prefix

copykit_per_sample<-function(dat,sample_name){
    pdf(paste0(outdir,prefix,".qc_metrics.pdf"))
    plotMetrics(dat, metric = c("overdispersion", 
                                "breakpoint_count",
                                "reads_total",
                                "reads_duplicates",
                                "reads_assigned_bins",
                                "percentage_duplicates"),label = "reads_total")
    dev.off()

    dat <- findAneuploidCells(dat)
    dat <- findOutliers(dat)

    pdf(paste0(outdir,prefix,".subclone.heatmap.outlier.pdf"))
    plotHeatmap(dat, row_split='outlier',n_threads=cpu_count) 
    dev.off()

    # kNN smooth profiles
    dat <- knnSmooth(dat)

    # Create a umap embedding 
    dat <- runUmap(dat)
    k_clones<-findSuggestedK(dat) #10
    dat  <- findClusters(dat, k_superclones=k_clones@metadata$suggestedK-3, k_subclones=k_clones@metadata$suggestedK+3)#output from k_clones

    #pdf(paste0(prefix,".subclone.umap.pdf"))
    #plotUmap(dat, label = 'subclones')
    #dev.off()

    pdf(paste0(outdir,prefix,".samples.umap.pdf"))
    plotUmap(dat, label = 'sample_name')
    dev.off()

    # Calculate consensus profiles for each subclone, 
    # and order cells by cluster for visualization with plotHeatmap
    dat <- calcConsensus(dat)
    dat <- runConsensusPhylo(dat)
    dat <- runPhylo(dat, metric = 'manhattan')
    dat <- calcInteger(dat, method = 'scquantum',assay="segment_ratios")
    col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    # Plot a copy number heatmap with clustering annotation

    #pdf(paste0(prefix,".subclone.phylo.pdf"))
    #plotPhylo(dat, label = 'subclones')
    #dev.off()

    pdf(paste0(outdir,prefix,".subclone.heatmap.segment_ratios.pdf"))
    plotHeatmap(dat, label = c('reads_total','sample_name','method'),  
        order_cells = 'consensus_tree',
        assay="segment_ratios",
        n_threads=cpu_count, 
        col=col_fun)
    dev.off()

    pdf(paste0(outdir,prefix,".subclone.heatmap.smoothed_bincounts.pdf"))
    plotHeatmap(dat, label = c('reads_total','sample_name','method'),  
        order_cells = 'consensus_tree',
        assay="smoothed_bincounts",
        n_threads=cpu_count)
    dev.off()

    pdf(paste0(outdir,prefix,".subclone.heatmap.integer.pdf"))
    plotHeatmap(dat, label = c('reads_total','sample_name','method'), 
        order_cells = 'consensus_tree', 
        assay="integer",
        n_threads=cpu_count)
    dev.off()

    saveRDS(dat,file=paste0(outdir,prefix,".scCNA.rds"))
    write.table(as.data.frame(dat@colData),file=paste0(outdir,prefix,".scCNA.tsv"),sep="\t",col.names=T,row.names=T)

}
#run line by line
#singularity shell --bind /data/rmulqueen/projects/scalebio_dcis/ ~/singularity/copykit.sif


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
library(dendextend,lib.loc="/home/rmulqueen/R/x86_64-pc-linux-gnu-library/4.3") #add this
library(cluster)
BiocParallel::bpparam()

option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=".", 
              help="Directory with single-cell bam files", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=".", 
              help="Output directory for plots and RDS file"),
  make_option(c("-p", "--output_prefix"), type="character", default="scale", 
              help="Prefix of output"),
  make_option(c("-t", "--task_cpus"), type="integer", default=125, 
              help="Integer number of cpus")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
outdir=ifelse(endsWith(opt$output_dir,suffix="/"),opt$output_dir, paste0(opt$output_dir,"/"))
cpu_count=opt$task_cpus
prefix=opt$output_prefix

register(MulticoreParam(progressbar = T, workers = cpu_count), default = T)

dat <- runVarbin(opt$input_dir, is_paired_end = TRUE, remove_Y = TRUE,resolution="500kb")
dat  <- runMetrics(dat)

#note that "sample" is a default metadata column for the cell label, so we need to use "sample_name"
colData(dat)$sample_name<-unlist(lapply(strsplit(row.names(colData(dat)),"[.]"),"[",1))
colData(dat)$method<-prefix

copykit_per_sample<-function(dat_in,sample_name_in="BCMDCIS07T"){
    print(paste("Processing sample:",sample_name_in))
    dat_in <- dat_in[,colData(dat_in)$sample_name == sample_name_in]
    if(ncol(dat_in)<50){
        print(paste("Sample",sample_name_in,"too few cells, skipping."))
    } else {

    pdf(paste0(outdir,prefix,".",sample_name_in,".qc_metrics.pdf"))
    plt<-plotMetrics(dat_in, metric = c("overdispersion", 
                                "breakpoint_count",
                                "reads_total",
                                "reads_duplicates",
                                "reads_assigned_bins",
                                "percentage_duplicates"),label = "reads_total")
    print(plt)
    dev.off()

    dat_in <- findAneuploidCells(dat_in,simul=FALSE)
    dat_in <- findOutliers(dat_in)

    pdf(paste0(outdir,prefix,".",sample_name_in,".subclone.heatmap.outlier.pdf"))
    plt<-plotHeatmap(dat_in, row_split='outlier',n_threads=cpu_count) 
    print(plt)
    dev.off()

    # kNN smooth profiles
    dat_in <- knnSmooth(dat_in,k=10)
    col_fun = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

    #cluster with ward d2
    dend<-dist(t(dat_in@assays@data$segment_ratios)) %>%
                            hclust(method="ward.D2") %>%
                            as.dendrogram()
    k_clus<-find_k(dend, krange = 3:10)

    pdf(paste0(outdir,prefix,".",sample_name_in,".cluster_estimate.pdf"))    
    plot(k_clus,
        xlab = "Number of clusters (k)",
        ylab = "Average silhouette width",
        main = "Estimating the number of clusters using\n average silhouette width")
    plot(silhouette(k_clus$pamobject))
    dev.off()
    memb <- cutree(dend, k = k_clus$nc)

    dat_in@distMat<-dist(t(dat_in@assays@data$segment_ratios))
    colData(dat_in)$cnv_clus<-as.character(unname(memb[row.names(dat_in@colData)]))

    pdf(paste0(outdir,prefix,".",sample_name_in,".subclone.heatmap.smoothed_bincounts.pdf"))
    plt<-plotHeatmap(dat_in, 
        label = c('reads_total','sample_name','method','cnv_clus'),  
        assay="segment_ratios",
        n_threads=cpu_count,
        row_split="cnv_clus",
        col=col_fun,
        order_cells='hclust')
    print(plt)
    dev.off()


    # Create a umap embedding 
    dat_in <- runUmap(dat_in,assay="segment_ratios")

    pdf(paste0(outdir,prefix,".",sample_name_in,".clus.umap.pdf"))
    plt<-plotUmap(dat_in, label = 'cnv_clus')
    print(plt)
    dev.off()

    pdf(paste0(outdir,prefix,".",sample_name_in,".superclones.umap.pdf"))
    plotUmap(dat_in, label = 'cnv_clus')
    dev.off()

    # Calculate consensus profiles for each subclone, 
    # and order cells by cluster for visualization with plotHeatmap
    dat_in <- calcConsensus(dat_in,consensus_by="cnv_clus")
    dat_in <- runConsensusPhylo(dat_in)
    dat_in <- runPhylo(dat_in, metric = 'manhattan')
    dat_in <- calcInteger(dat_in, method = 'scquantum',assay="segment_ratios")
    col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    # Plot a copy number heatmap with clustering annotation

    #pdf(paste0(prefix,".subclone.phylo.pdf"))
    #plotPhylo(dat, label = 'subclones')
    #dev.off()

    pdf(paste0(outdir,prefix,".",sample_name_in,".subclone.heatmap.segment_ratios.pdf"))
    plt<-plotHeatmap(dat_in, label = c('reads_total','sample_name','method'),  
        order_cells = 'consensus_tree',
        assay="segment_ratios",
        n_threads=cpu_count, 
        col=col_fun)
    print(plt)
    dev.off()

    pdf(paste0(outdir,prefix,".",sample_name_in,".subclone.heatmap.smoothed_bincounts.pdf"))
    plt<-plotHeatmap(dat_in, label = c('reads_total','sample_name','method'),  
        order_cells = 'consensus_tree',
        assay="smoothed_bincounts",
        n_threads=cpu_count)
    print(plt)
    dev.off()

    pdf(paste0(outdir,prefix,".",sample_name_in,".subclone.heatmap.integer.pdf"))
    plt<-plotHeatmap(dat_in, label = c('reads_total','sample_name','method'), 
        order_cells = 'consensus_tree', 
        assay="integer",
        n_threads=cpu_count)
    print(plt)
    dev.off()

    saveRDS(dat_in,file=paste0(outdir,prefix,".",sample_name_in,".scCNA.rds"))
    write.table(as.data.frame(dat_in@colData)[,c("sample","sample_name","overdispersion","cnv_clus")],file=paste0(outdir,prefix,".",sample_name_in,".scCNA.tsv"),sep="\t",col.names=F,row.names=F)
    }
}

sample_list<-unique(colData(dat)$sample_name)

lapply(sample_list,function(x) {copykit_per_sample(dat,x)})
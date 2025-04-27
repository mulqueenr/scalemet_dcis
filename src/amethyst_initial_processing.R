#run per line 
#singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif

#singularity exec \
#--bind /data/rmulqueen/projects/scalebio_dcis \
#~/singularity/amethyst.sif
#Rscript /data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_initial_processing.R \
#--input_dir ${runDir}/scale_dat \
#--task_cpus 150

library(amethyst)
library(rhdf5)
library(data.table)
library(ggplot2)
library(patchwork)
library(tibble)
library(tidyr)
library(plyr); library(dplyr)
library(future)
library(furrr)
library(purrr)
library(cowplot)
library(pheatmap)
library(optparse,lib.loc="/home/rmulqueen/R/x86_64-conda-linux-gnu-library/4.4") #add this


option_list = list(
    make_option(c("-i", "--input_dir"), type="character", default="/data/rmulqueen/projects/scalebio_dcis/data/240202_prelim1/scale_dat", 
                help="Run Directory, output from ScaleMethyl pipeline", metavar="character"),
    make_option(c("-p", "--output_prefix"), type="character", default="scale", 
                help="Prefix of output for all samples merged amethyst output."),
    make_option(c("-c","--copykit_input"),type="character",default="scale.merged.cnv.tsv",
                help="Merged output from COPYKIT function call, tsv."),
    make_option(c("-t", "--task_cpus"), type="integer", default=125, 
                help="Integer number of cpus")
);



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
cpu_count=opt$task_cpus
prefix=opt$output_prefix
cnv<-read.table(opt$copykit_input,col.names=c("cell_id","sample_name","overdispersion","superclones","subclones"))
row.names(cnv)<-cnv$sample

#read in all sample/csv/sample.passingCellsMapMethylStats.csv files into data frame
#make a dataframe of all h5 files also <sample>\t<h5location>
in_dir=opt$input_dir
setwd(in_dir)

samples_list_meta<-list.files("./samples",pattern="*allCells.csv",full.names=T)


prepare_amethyst_obj<-function(sample_meta="./samples/BCMDCIS07T.allCells.csv",cnv_in=cnv){
    sample_name<-gsub(sample_meta,pattern=".allCells.csv",replace="")
    sample_name<-gsub(sample_name,pattern="./samples/",replace="")
    sample_meta<-read.csv(sample_meta)
    sample_meta<-sample_meta[which(sample_meta$pass=="pass"),]
    sample_meta$h5_path<-unlist(lapply(1:nrow(sample_meta),function(x){
        well=sample_meta[x,]$tgmt_well
        h5_file=list.files(paste0(getwd(),"/samples/methylation_coverage/amethyst/",sample_name),pattern=well,full.names=T,include.dirs=T)
        return(h5_file)
        }))
    print(paste("Generating amethyst object for :",as.character(sample_name)))

    obj <- createObject()

    #metadata MUST have a column called mcg_pct for score calculation
    #metadata MUST have a column called cov to regress coverage mias
    obj@metadata<-sample_meta
    row.names(obj@metadata)<-obj@metadata$cell_id

    #add in CNV information
    obj@metadata$overdispersion<-cnv_in[row.names(obj@metadata),"overdispersion"]
    obj@metadata$superclones<-cnv_in[row.names(obj@metadata),"superclones"]
    obj@metadata$subclones<-cnv_in[row.names(obj@metadata),"subclones"]
    cg_per_read<-mean(obj@metadata$cg_cov/obj@metadata$unique_reads)
    print(paste("Plotting QC for :",as.character(sample_name)))

    head(obj@metadata)
    plt1<-ggplot(obj@metadata, aes(x=unique_reads/total_reads,y=log10(unique_reads)))+geom_point(size=0.2)+theme_minimal()+xlim(c(0,1))+ylim(c(0,7))
    plt2<-ggplot(obj@metadata, aes(x=sample, y = log10(unique_reads))) + geom_violin() + geom_jitter()+ylim(c(0,7))+theme_minimal()
    plt3<-ggplot(obj@metadata, aes(x=sample, y = log10(cg_cov))) + geom_violin() + geom_jitter()+ylim(c(0,7))+theme_minimal()
    plt4<-ggplot(obj@metadata, aes(x=unique_reads, y = cg_cov)) +geom_point() +  geom_smooth(method = "lm", se = FALSE)+theme_minimal()+ggtitle(paste("CG covered per read: ",cg_per_read))
    plt5<-ggplot(obj@metadata, aes(x=sample, y = mcg_pct))+ geom_violin() + geom_jitter() +ylim(c(0,100))+theme_minimal()

    ggsave((plt1|plt2|plt3)/(plt4|plt5),file=paste0("./amethyst/",sample_name,".cov_plot.pdf"))

    obj@h5paths <- data.frame(row.names = c(rownames(obj@metadata)), paths = obj@metadata$h5_path)

    # index files
    obj@index[["chr_cg"]] <- indexChr(obj, type = "CG", threads = cpu_count)
    print(paste("Saving amethyst object for :",as.character(sample_name)))
    saveRDS(obj,file=paste0(sample_name,".amethyst.rds"))
    return(obj)
}

obj_list<-lapply(samples_list_meta,prepare_amethyst_obj)

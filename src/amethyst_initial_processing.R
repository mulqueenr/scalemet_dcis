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
    make_option(c("-i", "--input_dir"), type="character", default="/data/rmulqueen/projects/scalebio_dcis/data/240523_prelim2/scale_dat", 
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
row.names(cnv)<-unlist(lapply(strsplit(cnv$cell_id,"[.]"),"[",3))

#set up ref
gtf <- rtracklayer::readGFF("/container_ref/gencode.v43.annotation.gtf.gz")
for (i in c("gene_name", "exon_number")) {
    gtf$i <- unlist(lapply(gtf$attributes, extractAttributes, i))}
gtf <- dplyr::mutate(gtf, location = paste0(seqid, "_", start, "_", end))

#read in all sample/csv/sample.passingCellsMapMethylStats.csv files into data frame
#make a dataframe of all h5 files also <sample>\t<h5location>
in_dir=opt$input_dir

samples_list_meta<-list.files(paste0(in_dir,"/samples"),pattern="*allCells.csv",full.names=T)

correct_h5_cellnames<-function(h5,run_id){
    h5list = h5ls(h5)
    print(paste("Correcting",run_id,"metadata cell names in",basename(h5)))

    for (i in 1:nrow(h5list)){
        tryCatch({if(endsWith(h5list[i,"group"],"/CG")){ #just run on CG
            if( !(paste0(h5list[i,"name"],"+",run_id) %in% h5list$name) && !endsWith(h5list[i,"name"],run_id)){
                celldat<-h5read(h5,paste0("CG/",h5list[i,"name"]))
                h5write(celldat, file=h5, name=paste0("CG/",h5list[i,"name"],"+",run_id))
            }
        }}, error =function(e) { cat("Proceeding past line",i,"for",basename(h5),"\n")} )      
    }
}

prepare_amethyst_obj<-function(sample_meta="./samples/BCMDCIS07T.allCells.csv",cnv_in=cnv){
    sample_name<-gsub(sample_meta,pattern=".allCells.csv",replace="")
    run_id=strsplit(in_dir,"/")[[1]][length(strsplit(in_dir,"/")[[1]])-1]
    sample_name<-basename(sample_name)
    sample_meta<-read.csv(sample_meta)
    sample_meta<-sample_meta[which(sample_meta$pass=="pass"),]
    sample_meta$h5_path<-unlist(lapply(1:nrow(sample_meta),function(x){
        well=sample_meta[x,]$tgmt_well
        h5_file=list.files(paste0(in_dir,"/samples/methylation_coverage/amethyst/",sample_name),pattern=well,full.names=T,include.dirs=T)
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

    ggsave((plt1|plt2|plt3)/(plt4|plt5),file=paste0(sample_name,".cov_plot.pdf"))

    obj@h5paths <- data.frame(row.names = c(rownames(obj@metadata)), paths = obj@metadata$h5_path)

    #correct h5 names
    print(paste("Appended",run_id,"to names in h5 files."))
    lapply(unique(obj@h5paths$paths),function(i){correct_h5_cellnames(h5=i,run_id)})
    h5closeAll()

    #correct cell id names
    print(paste("Corrected",run_id,"metadata cell names."))
    row.names(obj@h5paths)<-paste0(row.names(obj@h5paths),"+",run_id)
    row.names(obj@metadata)<-paste0(row.names(obj@metadata),"+",run_id)

    #add ref
    obj@ref<-gtf

    # index files
    obj@index[["chr_cg"]] <- indexChr(obj, type = "CG", threads = cpu_count)

    #initial window genomeMatrix
    obj@genomeMatrices[["cg_100k_score"]] <- makeWindows(obj,
                                                      stepsize = 100000, type = "CG", 
                                                      metric = "score", index = "chr_cg", nmin = 2,
                                                      threads = as.numeric(cpu_count)) 

    print(paste("Saving amethyst object for :",as.character(sample_name)))
    saveRDS(obj,file=paste0(sample_name,".amethyst.rds"))
    return(obj)
}

obj_list<-lapply(samples_list_meta,prepare_amethyst_obj)

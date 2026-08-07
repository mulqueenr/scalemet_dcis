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
    make_option(c("-i", "--input_dir"), type="character", default="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/scalemethyl_pipeline_out", 
                help="Run Directory, output from ScaleMethyl pipeline", metavar="character"),
    make_option(c("-p", "--output_prefix"), type="character", default="scaledcis", 
                help="Prefix of output for all samples merged amethyst output."),
    make_option(c("-t", "--task_cpus"), type="integer", default=300, 
                help="Integer number of cpus")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
cpu_count=opt$task_cpus
prefix=opt$output_prefix
#cnv<-read.table(opt$copykit_input,col.names=c("cell_id","sample_name","overdispersion","superclones","subclones"))
#row.names(cnv)<-unlist(lapply(strsplit(cnv$cell_id,"[.]"),"[",3))

#set up ref
gtf <- rtracklayer::readGFF("/container_ref/gencode.v43.annotation.gtf.gz")
for (i in c("gene_name", "exon_number")) {
    gtf$i <- unlist(lapply(gtf$attributes, extractAttributes, i))}
gtf <- dplyr::mutate(gtf, location = paste0(seqid, "_", start, "_", end))

#read in all sample/csv/sample.passingCellsMapMethylStats.csv files into data frame
#make a dataframe of all h5 files also <sample>\t<h5location>
in_dir=opt$input_dir
system(paste0("mkdir -p ",in_dir,"/","amethyst_plate_obj"))

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

read_plate_sample_meta<-function(sample_meta){
    sample_name<-gsub(basename(sample_meta),pattern=".allCells.csv",replace="")
    plate_info=strsplit(gsub(sample_meta,pattern=in_dir,replace=""),"/")[[1]][2]
    batch<-strsplit(plate_info,"_")[[1]][1]
    method<-strsplit(plate_info,"_")[[1]][2]
    plate<-strsplit(plate_info,"_")[[1]][3]

    sample_meta<-read.csv(file=sample_meta,colClasses=c("i7_well"="character","i5_well"="character","tgmt_well"="character"))
    sample_meta<-sample_meta[which(sample_meta$pass=="pass"),]
    sample_meta$batch<-batch
    sample_meta$method<-method
    sample_meta$plate<-plate
    sample_meta$plate_info<-plate_info
    sample_meta$sample<-sample_name

    sample_meta$h5_path<-unlist(lapply(1:nrow(sample_meta),function(x){
        well=sample_meta[x,]$tgmt_well
        h5_file=list.files(paste(in_dir,plate_info,"samples/methylation_coverage/amethyst",sample_name,sep="/"),pattern=well,full.names=T,include.dirs=T)
        return(h5_file)
        }))
    sample_meta$bam_path<-unlist(lapply(1:nrow(sample_meta),function(x){
        well=sample_meta[x,]$tgmt_well
        sample_name=sample_meta[x,]$sample
        bam_file=list.files(paste(in_dir,plate_info,"alignments/dedup",paste0(sample_name,".",well),sep="/"),pattern="dedup.bam",full.names=T,include.dirs=T)
    }))

    print(paste("Generating amethyst object for plate:",as.character(plate_info), "Sample:", as.character(sample_name)))
    return(sample_meta)
}

prepare_amethyst_obj<-function(sample_plate_meta){
    #sample_meta is a list of all samples per plate.

    plate_meta<-do.call("rbind",lapply(sample_plate_meta,read_plate_sample_meta))
    plate_info<-plate_meta$plate_info[1]
    obj <- createObject()

    #metadata MUST have a column called mcg_pct for score calculation
    #metadata MUST have a column called cov to regress coverage mias
    obj@metadata<-plate_meta
    row.names(obj@metadata)<-obj@metadata$cell_id

    #add in CNV information
    #obj@metadata$overdispersion<-cnv_in[row.names(obj@metadata),"overdispersion"]
    #obj@metadata$superclones<-cnv_in[row.names(obj@metadata),"superclones"]
    #obj@metadata$subclones<-cnv_in[row.names(obj@metadata),"subclones"]
    
    cg_per_read<-mean(obj@metadata$cg_cov/obj@metadata$unique_reads)
    obj@h5paths <- data.frame(row.names = c(rownames(obj@metadata)), paths = obj@metadata$h5_path)

    #correct h5 names
    print(paste("Appended",plate_info,"to names in h5 files."))
    lapply(unique(obj@h5paths$paths),function(i){correct_h5_cellnames(h5=i,plate_info)})
    h5closeAll()

    #correct cell id names
    print(paste("Corrected",plate_info,"metadata cell names."))
    row.names(obj@h5paths)<-paste0(row.names(obj@h5paths),"+",plate_info)
    row.names(obj@metadata)<-paste0(row.names(obj@metadata),"+",plate_info)

    #add ref
    obj@ref<-gtf

    # index files
    obj@index[["chr_cg"]] <- indexChr(obj, type = "CG", threads = cpu_count)

    #initial window genomeMatrix
    obj@genomeMatrices[["cg_100k_score"]] <- makeWindows(obj,
                                                      stepsize = 100000, type = "CG", 
                                                      metric = "score", index = "chr_cg", nmin = 2,
                                                      threads = as.numeric(cpu_count)) 

    print(paste("Saving amethyst object for :",as.character(plate_info)))
    saveRDS(obj,file=paste0(in_dir,"/amethyst_plate_obj/",plate_info,".amethyst.rds"))
    return(obj)
}

plates<-list.files(in_dir)
plates<-plates[grep(plates,pattern="_homebrew_|_scalebio_")]


plates<-plates[grep(plates,invert=T,pattern="batch1|batch2")]

plates_meta<-lapply(plates,function(plate) list.files(paste(in_dir,plate,"samples",sep="/"),pattern="*allCells.csv",full.names=T))

plates_objs<-lapply(plates_meta,prepare_amethyst_obj)
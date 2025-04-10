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
library(dplyr)
library(tibble)
library(tidyr)
library(plyr)
library(future)
library(furrr)
library(purrr)
library(cowplot)
library(pheatmap)
library(plyr)
library(optparse) #add this


option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default="/data/rmulqueen/projects/scalebio_dcis/data/250329_RM_scalebio_batch1_initseq/scale_dat", 
              help="Run Directory, output from ScaleMethyl pipeline", metavar="character"),
  make_option(c("-p", "--output_prefix"), type="character", default="scale", 
              help="Prefix of output for all samples merged amethyst output."),
  make_option(c("-c", "--task_cpus"), type="integer", default=125, 
              help="Integer number of cpus")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
cpu_count=opt$task_cpus
prefix=opt$output_prefix

#read in all sample/csv/sample.passingCellsMapMethylStats.csv files into data frame
#make a dataframe of all h5 files also <sample>\t<h5location>
in_dir=opt$input_dir
setwd(in_dir)
system("mkdir -p ./amethyst")

samples_list_meta<-list.files("./samples",pattern="*allCells.csv",full.names=T)


methyltree_output<-function(obj=obj,prefix="DCIS-41T",sample="DCIS-41T",filt_min_pct=10,filt_max_pct=80,threads=1){
        obj@metadata$methyltree_group<-"all"
        #make 500bp windows with methylation percentages
        methyltreewindows <- calcSmoothedWindows(obj, 
                                            type = "CG", 
                                            threads = threads,
                                            step = 500,
                                            smooth = 1,
                                            index = "chr_cg",
                                            groupBy = "methyltree_group", 
                                            returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                            returnPctMatrix = TRUE)
        print(paste("Starting number of windows:",as.character(nrow(methyltreewindows[["pct_matrix"]]))))
        methyltreewindows[["pct_matrix"]]<-methyltreewindows[["pct_matrix"]][methyltreewindows[["pct_matrix"]]$all>=filt_min_pct & methyltreewindows[["pct_matrix"]]$all<=filt_max_pct,]
        #filter to windows to middling methylation values
        print(paste("Filtering by m0 >=", as.character(filt_min_pct), "m1 <=", as.character(filt_max_pct),as.character(nrow(methyltreewindows[["pct_matrix"]]))))
        #merge windows that are touching
        methyltreewindows<-reduce(GenomicRanges::makeGRangesFromDataFrame(methyltreewindows[["pct_matrix"]]))
        print(paste("Filtered window count:",as.character(nrow(as.data.frame((methyltreewindows))))))
        print(paste("Filtered window average width:",as.character(mean(width(methyltreewindows)))))
        print(paste("Total genome covered:",as.character(sum(width(methyltreewindows))/1000000),"Mbp"))
        #make a merged windows percentile matrix per cell for methyltree

        methyltreeoutput<-makeWindows(obj,
                    type = "CG", 
                    metric = "percent", 
                    bed = as.data.frame(methyltreewindows,col.names=NULL)[,1:3],
                    threads = threads, 
                    index = "chr_cg", 
                    nmin = 1) 
      print(paste("Mean percent cells covered per window:",
            mean((rowSums(!is.na(methyltreeoutput))/ncol(methyltreeoutput))*100)))
      print("Filtering to windows with >10% of cells with coverage")
      methyltreeoutput<-methyltreeoutput[(rowSums(!is.na(methyltreeoutput))/ncol(methyltreeoutput)*100)>=10,]
      methyltreewindows<-data.frame(do.call("rbind",strsplit(row.names(methyltreeoutput),"_")))
      colnames(methyltreewindows)<-c("chr","start","end")
      methyltreewindows<-GenomicRanges::makeGRangesFromDataFrame(methyltreewindows)
      print(paste("Final Filtered window count:",as.character(nrow(as.data.frame((methyltreewindows))))))
      print(paste("Final Filtered window average width:",as.character(mean(width(methyltreewindows)))))
      print(paste("Final Total genome covered:",as.character(sum(width(methyltreewindows))/1000000),"Mbp"))
      methyltreeoutput<-makeWindows(obj,
                  type = "CG", 
                  metric = "percent", 
                 bed = as.data.frame(methyltreewindows,col.names=NULL)[,1:3],
                  threads = threads, 
                  index = "chr_cg", 
                  nmin = 1) 

      methyltreeoutput$genomic_region_id<-row.names(methyltreeoutput)
      methyltreeoutput <- methyltreeoutput |> 
          pivot_longer(
          cols = !genomic_region_id, 
          names_to = "cell_id",
          values_to = "value",
          values_drop_na = TRUE)
    #make a metadata sheet with cluster info
    out_metadata<-obj@metadata[,c("pass","cluster_id","cg_cov","mcg_pct","subclones")]
    colnames(obj_metadata)<-c("HQ","celltype","nCG","met_rate","large_clone_id") #match names
    out_metadata$sample<-row.names(out_metadata) #sample (cell) names
    out_metadata$met_rate<-out_metadata$met_rate/100 #percentage to rate

    if(file.exists(paste0(prefix,"_methyltree_input.h5"))){
        system(paste0("rm -rf ",prefix,"_methyltree_input.h5"))
    }
      h5createFile(file=paste0(prefix,"_methyltree_input.h5"))
      h5write(methyltreeoutput,file=paste0(prefix,"_methyltree_input.h5"),name="data")
      h5write(out_metadata,file=paste0(prefix,"_methyltree_input.h5"),name="metadata")
}



prepare_amethyst_obj<-function(sample_meta="./samples/BCMDCIS07T.allCells.csv"){
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
    saveRDS(obj,file=paste0("./amethyst/",sample_name,".amethyst.rds"))

    print(paste("Generating Methyltree input for :",as.character(sample_name)))

    methyltree_output(obj,prefix=sample_name,sample=sample_name,filt_min_pct=10,filt_max_pct=80,threads=cpu_count)
    return(obj)
}

obj_list<-lapply(samples_list_meta,prepare_amethyst_obj)

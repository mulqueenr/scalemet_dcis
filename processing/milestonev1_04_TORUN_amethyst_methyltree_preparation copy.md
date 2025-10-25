
# Output data for methyltree input as well.


```R

source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
library(Rsamtools)
library(copykit)

#set environment and read in data
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)
obj<-readRDS(file="04_scaledcis.amethyst.broad_celltype.rds")

methyltree_output<-function(obj=obj,output_directory=".",
                            sample="DCIS-41T",
                            filt_min_pct=10,
                            filt_max_pct=80,
                            threads=1){
        obj_met<-subsetObject(obj, cells = row.names(obj@metadata[obj@metadata$sample %in% sample,]))
        obj_met@metadata$methyltree_group<-"all"
        #make 500bp windows with methylation percentages
        methyltreewindows <- calcSmoothedWindows(obj_met, 
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

        methyltreeoutput<-makeWindows(obj_met,
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

      methyltreeoutput<-makeWindows(obj_met,
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
    out_metadata<-obj_met@metadata[,c("pass","broad_celltype","cg_cov","mcg_pct","subclones")]
    colnames(out_metadata)<-c("HQ","celltype","nCG","met_rate","large_clone_id") #match names
    out_metadata$sample<-row.names(out_metadata) #sample (cell) names
    out_metadata$met_rate<-out_metadata$met_rate/100 #percentage to rate

    methyltree_input_file=paste(output_directory,paste("/methyltree",sample_name,"_methyltree_input.h5",sep="."),sep="/")
    if(file.exists(methyltree_input_file)){
        system(paste0("rm -rf ",methyltree_input_file))
    }
      h5createFile(file=methyltree_input_file)
      h5write(methyltreeoutput,file=methyltree_input_file,name="data")
      h5write(out_metadata,file=methyltree_input_file,name="metadata")
    }

output_directory=paste0(dirname(getwd()),"/methyltree")
system(paste("mkdir -p",output_directory))
lapply(unique(obj@metadata$sample),function(x) 
methyltree_output(obj=obj,
                    output_directory=output_directory,
                    sample=x,
                    filt_min_pct=10,
                    filt_max_pct=80,
                    threads=1))

```
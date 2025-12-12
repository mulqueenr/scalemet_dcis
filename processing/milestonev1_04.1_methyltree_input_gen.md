```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
```

# Generate methyltree formatted input for each sample

```R
source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in

#set environment and read in data
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)
obj<-readRDS(file="06_scaledcis.cnv_clones.amethyst.rds")


#make output directory
system(paste0("mkdir -p ",project_data_directory,"/methyltree"))


# Output files for methyltree format
```R

methyltree_output<-function(obj=obj,
                            sample_name="BCMDCIS41T",
                            filt_min_pct=20,
                            filt_max_pct=70,
                            threads=1){
        
        output_directory=paste0(project_data_directory,"/methyltree/",sample_name[1])
        system(paste0("mkdir -p ",output_directory))

        obj_met<-subsetObject(obj, cells = row.names(obj@metadata[obj@metadata$sample %in% sample_name,]))
        obj_met@metadata$methyltree_group<-"all"
        obj_met@metadata$pass<-"TRUE"
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
    out_metadata<-obj_met@metadata[,c("pass","fine_celltype","cg_cov","mcg_pct","cnv_clonename")]
    colnames(out_metadata)<-c("HQ","celltype","nCG","met_rate","cnv_clonename") #match names
    out_metadata$sample<-row.names(out_metadata) #sample (cell) names
    out_metadata$met_rate<-out_metadata$met_rate/100 #percentage to rate
    out_metadata$large_clone_id<-as.numeric(as.factor(as.character(out_metadata$cnv_clonename))) #factorize clones (to run on methyltree)

    methyltree_input_file=paste(output_directory,
        paste("methyltree",sample_name[1],"methyltree_input.h5",sep="."),sep="/")
    if(file.exists(methyltree_input_file)){
        system(paste0("rm -rf ",methyltree_input_file))
    }
      h5createFile(file=methyltree_input_file)
      h5write(methyltreeoutput,file=methyltree_input_file,name="data")
      h5write(out_metadata,file=methyltree_input_file,name="metadata")
    }

methyltree_output(obj=obj,sample_name=c('BCMDCIS05T'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS07T'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS102T_24hTis'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS124T'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS22T'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS28T'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS32T'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS35T'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS41T'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS49T'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS52T'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS65T'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS66T'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS70T'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS74T'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS79T_24hTis_DCIS'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS80T_24hTis'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS82T_24hTis'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS92T_24hTis'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS94T_24hTis'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS97T'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMDCIS99T '),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMHBCA03R'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMHBCA04R'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMHBCA09R-3h'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMHBCA12R-3h'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMHBCA16R-3h'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMHBCA17R-3h'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMHBCA19R-4h'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMHBCA22R-4h'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMHBCA26L-24hTis-4h'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMHBCA29L-2h'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMHBCA38L-3h'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMHBCA83L-3h'),threads=1)
methyltree_output(obj=obj,sample_name=c('BCMHBCA85L-3h'),threads=1)
methyltree_output(obj=obj,sample_name=c('ECIS25T'),threads=1)
methyltree_output(obj=obj,sample_name=c('ECIS26T'),threads=1)
methyltree_output(obj=obj,sample_name=c('ECIS36T'),threads=1)
methyltree_output(obj=obj,sample_name=c('ECIS48T'),threads=1)
methyltree_output(obj=obj,sample_name=c('ECIS57T'),threads=1)
```
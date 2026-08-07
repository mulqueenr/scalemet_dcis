#```bash
#singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
#```

# Generate methyltree formatted input for each sample

```R
set.seed(111)
#options(future.globals.maxSize= 50000*1024^2) #80gb limit for parallelizing
library(amethyst)
library(data.table)
library(dplyr)
library(GenomicRanges)
library(parallel)
library(rhdf5)
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"

#read in object from directory
task_cpus=300
processing_folder="04_methyltree_input"
wd=paste(sep="/",project_data_directory,processing_folder)
system(paste0("mkdir -p ",wd))
setwd(wd)

#make output directory
output_directory=paste0(wd,"/methyltree_input")
system(paste0("mkdir -p ",output_directory))

obj<-readRDS(file=paste(project_data_directory,"03_fine_celltyping","03_scaledcis.final_celltypes.amethyst.rds",sep="/"))


# Output files for methyltree format
methyltree_output<-function(obj=obj,
                            sample_name="BCMDCIS35T",
                            filt_min_pct=30,
                            filt_max_pct=70,
                            output_dir=output_directory,
                            threads=1){
        
        output_directory=paste0(output_dir,"/",sample_name[1])
        system(paste0("mkdir -p ",output_directory))

        obj_met<-subsetObject(obj, cells = row.names(obj@metadata[obj@metadata$Sample %in% sample_name,]))

        obj_met@metadata$methyltree_group<-"all"
        obj_met@metadata$pass<-"TRUE"
        #make 500bp windows with methylation percentages, this is just used to filter
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
        
        #save unmerged data as well
        methyltreewindows_unmerged<-GenomicRanges::makeGRangesFromDataFrame(methyltreewindows[["pct_matrix"]])

        #merge windows that are touching
        methyltreewindows<-IRanges::reduce(methyltreewindows_unmerged)
        print(paste("Filtered window count:",as.character(nrow(as.data.frame((methyltreewindows))))))
        print(paste("Filtered window average width:",as.character(mean(width(methyltreewindows)))))
        print(paste("Total genome covered:",as.character(sum(width(methyltreewindows))/1000000),"Mbp"))
        #make a merged windows percentile matrix per cell for methyltree
        options(future.globals.maxSize= 500000*1024^2) #80gb limit for parallelizing

        #merged output
        methyltreeoutput<-makeWindows(obj_met,
                    type = "CG", 
                    metric = "percent", 
                    bed = as.data.frame(methyltreewindows,col.names=NULL)[,1:3],
                    threads = threads, 
                    index = "chr_cg", 
                    nmin = 1) 

        #cnv filtered output
        #first regenerate consensus integer copy per clone
        #get cnv data per cell
        cnv<-obj_met@genomeMatrices$scquantum_cnv[,grep(colnames(obj_met@genomeMatrices$scquantum_cnv),pattern="metadata",invert=T)]
        
        #get clone assignment in metadata
        clones<-setNames(obj_met@metadata$cnv_clonename,nm=row.names(obj_met@metadata))
        #filter and order clones
        clones<-clones[colnames(cnv)]

        #get windows where clone has 95% cells reporting diploid
        clone_consensus <- cnv %>% 
            t() %>% as.data.frame() %>% 
            mutate(clones=clones) %>% 
            group_by(clones) %>% 
            summarise(across(where(is.numeric), function(x){
                ifelse(sum(as.numeric(x)==2)>(length(x)*0.95),TRUE,FALSE)
                }
            )) %>% 
            tibble::column_to_rownames(var = "clones") %>% 
            as.data.frame()

        #grab all windows with colMeans < 1 (so any clone has less that 95% diploid)
        filter_windows<-colnames(clone_consensus)[which(colMeans(clone_consensus)<1)]
        print(paste("Percent of cnv windows with all clones >95% diploid:",
        as.character((ncol(clone_consensus)-length(filter_windows))/ncol(clone_consensus)*100)))
        
        #get genomic regions with cnvs
        #filter out regions from unmerged
        #make windows for output
        print(paste("Mean percent cells covered per window:",
            mean((rowSums(!is.na(methyltreeoutput))/ncol(methyltreeoutput))*100)))

      print("Filtering to windows with >5% of cells with coverage")
      methyltreeoutput<-methyltreeoutput[(rowSums(!is.na(methyltreeoutput))/ncol(methyltreeoutput)*100)>=5,]
      methyltreewindows<-data.frame(do.call("rbind",strsplit(row.names(methyltreeoutput),"_")))

      colnames(methyltreewindows)<-c("chr","start","end")
      methyltreewindows<-GenomicRanges::makeGRangesFromDataFrame(methyltreewindows)
      print(paste("Final Filtered window count:",as.character(nrow(as.data.frame((methyltreewindows))))))
      print(paste("Final Filtered window average width:",as.character(mean(width(methyltreewindows)))))
      print(paste("Final Total genome covered:",as.character(sum(width(methyltreewindows))/1000000),"Mbp"))
      
      #methyltreeoutput<-makeWindows(obj_met,
      #                              type = "CG", 
      #                              metric = "percent", 
      #                              bed = as.data.frame(methyltreewindows,col.names=NULL)[,1:3],
      #                              threads = threads, 
      #                              index = "chr_cg", 
      #                              nmin = 1) 

      methyltreeoutput$genomic_region_id<-row.names(methyltreeoutput)

      methyltreeoutput <- methyltreeoutput |> 
          tidyr::pivot_longer(
          cols = !genomic_region_id, 
          names_to = "cell_id",
          values_to = "value",
          values_drop_na = TRUE)

    #make a metadata sheet with cluster info
    out_metadata<-obj_met@metadata[,c("pass","celltype","cg_cov","mcg_pct","cnv_clonename_500kb")]
    colnames(out_metadata)<-c("HQ","celltype","nCG","met_rate","cnv_clonename") #match names
    out_metadata$sample<-row.names(out_metadata) #sample (cell) names
    out_metadata$met_rate<-out_metadata$met_rate/100 #percentage to rate
    out_metadata$large_clone_id<-as.numeric(as.factor(as.character(out_metadata$cnv_clonename))) #factorize clones (to run on methyltree)

    #cnv windows to remove
    filter_windows_granges<-data.frame(do.call("rbind",strsplit(filter_windows,"_")))

    colnames(filter_windows_granges)<-c("chr","start","end")
    filter_windows_granges<-GenomicRanges::makeGRangesFromDataFrame(filter_windows_granges)
    overlaps<-GenomicRanges::findOverlaps(methyltreewindows,filter_windows_granges)
    non_cnv_overlap_windows<-methyltreewindows[!(1:length(methyltreewindows) %in% unique(overlaps@from)),]

    print(paste("Percent of methyltree windows in diploid regions for all clones:",length(non_cnv_overlap_windows)/length(methyltreewindows)*100))

    win_filt<-paste(
        seqnames(non_cnv_overlap_windows),
        start(non_cnv_overlap_windows),
        end(non_cnv_overlap_windows),
        sep="_")

    
    methyltreewindows_cnv_filtered<-methyltreeoutput[methyltreeoutput$genomic_region_id %in% win_filt,]

    methyltree_input_file=paste(output_directory,
        paste("methyltree",sample_name[1],"methyltree_input.h5",sep="."),sep="/")
    
    if(file.exists(methyltree_input_file)){file.remove(methyltree_input_file)}

    #initiate data frame and append in chunks (was hitting a seg fault error doing it all at once)
        methyltreeoutput<-as.matrix(methyltreeoutput)
        total_rows<-nrow(methyltreeoutput)
        total_cols<-ncol(methyltreeoutput)

        rhdf5::h5createFile(file=methyltree_input_file)
        h5createDataset(file=methyltree_input_file, dataset="data", 
                        dims=c(total_rows, total_cols), 
                        maxdims = c(total_rows, total_cols), 
                        storage.mode = storage.mode(methyltreeoutput), 
                        chunk = c(min(5000, total_rows), total_cols))
        chunk_size <- 1e6
        start_row <- 1

    while (start_row <= total_rows) {
        end_row <- min(start_row + chunk_size - 1, total_rows)
        chunk_data <- methyltreeoutput[start_row:end_row,]
        h5write(obj = chunk_data, 
            file = methyltree_input_file, 
            name = "data", 
            index = list(start_row:end_row, 1:total_cols))
        start_row <- end_row + 1
        }

    #do the same for cnv filtered data
    methyltreewindows_cnv_filtered<-as.matrix(methyltreewindows_cnv_filtered)
    total_rows<-nrow(methyltreewindows_cnv_filtered)
    total_cols<-ncol(methyltreewindows_cnv_filtered)

    h5createDataset(file=methyltree_input_file, dataset="data_cnv_filt", 
                        dims=c(total_rows, total_cols), 
                        maxdims = c(total_rows, total_cols), 
                        storage.mode = storage.mode(methyltreewindows_cnv_filtered), 
                        chunk = c(min(5000, total_rows), total_cols))
    chunk_size <- 1e6
    start_row <- 1

    while (start_row <= total_rows) {
        end_row <- min(start_row + chunk_size - 1, total_rows)
        chunk_data <- methyltreewindows_cnv_filtered[start_row:end_row,]
        h5write(obj = chunk_data, 
            file = methyltree_input_file, 
            name = "data_cnv_filt", 
            index = list(start_row:end_row, 1:total_cols))
        start_row <- end_row + 1
        }
      rhdf5::h5write(as.data.frame(out_metadata),file=methyltree_input_file,name="metadata")
      h5closeAll()
}

sample_names_priority<-table(obj@metadata$Sample,obj@metadata$cnv_ploidy_500kb) %>% 
            as.data.frame() %>% 
            filter(Var2=="aneuploid" & Freq >=100) %>% 
            select(Var1) %>% unlist() %>% droplevels() %>% levels()

#segfault of 66T for some reason
sample_names_priority<-sample_names_priority[!(sample_names_priority) %in% c('BCMDCIS79T_24hTis_DCIS','BCMDCIS79T_24hTis_IDC')]

for(i in sample_names_priority){
    print(strrep("-", 30))
    print(paste("Running ",i))
    print(strrep("-", 30))

    methyltree_output(obj=obj,sample_name=c(i),threads=1) 
    print(strrep("-", 30))
    print(paste("Done with ",i))
    print(strrep("-", 30))
}

#run 79T together since from same sample
methyltree_output(obj=obj,sample_name=c('BCMDCIS79T_24hTis_DCIS','BCMDCIS79T_24hTis_IDC'),threads=1) 

#now running lower priority samples
#segfault of 66T for some reason
sample_names<-unique(obj@metadata$Sample)
sample_names_lowpriority<-sample_names[!(sample_names %in% sample_names_priority) & !(sample_names %in% c('BCMDCIS79T_24hTis_DCIS','BCMDCIS79T_24hTis_IDC')) & !is.na(sample_names)]


for(i in sample_names_lowpriority){
    print(strrep("-", 30))
    print(paste("Running ",i))
    print(strrep("-", 30))

    methyltree_output(obj=obj,sample_name=c(i),threads=1) 
    print(strrep("-", 30))
    print(paste("Done with ",i))
    print(strrep("-", 30))
}

#prioritizing samples with >= 100 cancer cells for riley. but going to run all samples to help with celltyping lineage distinctions


```


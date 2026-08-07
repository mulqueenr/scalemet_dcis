
Collapsing cell type and clone bam reads into single files.

- Run BSbolt variant calling to assess VAF
- Run BISCUIT to assess structural variants
- Use cell type information for long read deconvolution

#script to make a merged bam per clone, and somatic lineage

```R
set.seed(111)
library(amethyst)
library(data.table)
library(dplyr)
library(GenomicRanges)
library(parallel)
library(Rsamtools)

options(future.globals.maxSize= 500000*1024^2) #80gb limit for parallelizing

project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"

#read in object from directory
task_cpus=300
processing_folder="04_collapsed_bam_and_dmr"
wd=paste(sep="/",project_data_directory,processing_folder)
system(paste0("mkdir -p ",wd))
setwd(wd)

system(paste0("mkdir -p ",wd,"/tmp"))


filter_bam_tmp<-function(cellid_list,idx){
    cellid<-cellid_list[idx,"cellid"]
    cellid_nobatch=strsplit(cellid,"[+]batch|[+]prelim")[[1]][1] #split down to bam index (remove batch info)
    bam<-cellid_list[idx,"bam_path"]
    
    tmp_bam=paste0(wd,"/tmp","/tmp.",idx,".bam")
    tmp_bam_sorted=paste0(wd,"/tmp","/tmp.",idx,".sorted")

    print(paste("Running cellid",cellid, "in",basename(bam)))
    if(!file.exists(paste0(bam,".bai"))){
          indexBam(bam)
    }
    param <- ScanBamParam(what=c("qname"),
                            flag=scanBamFlag(isPaired=TRUE,
                                            isProperPair=TRUE,
                                            isSecondaryAlignment=FALSE,
                                            isDuplicate=FALSE,
                                            isSupplementaryAlignment=FALSE))
    filt=FilterRules(function(x){gsub("^.*:", "", x$qname)==cellid_nobatch})
    filterBam(file=bam,destination=tmp_bam,filter=filt,param=param)
    sortBam(file=tmp_bam,destination=tmp_bam_sorted)
    indexBam(file=paste0(tmp_bam_sorted,".bam"))
}

output_merged_bam<-function(cellid_list,output_directory,output_sample,output_name,task_cpus=100){
    #scalebio pipeline outputs bam files as Tn5 wells. so multiple cell IDs are in a bam. this function splits out the bam to the query cellid
    #cellid_list is data frame with columns cellid and bam_path from metadata
    #clear tmp directory

    if(length(list.files(paste0(wd,"/tmp")))>0){
        file.remove(list.files(paste0(wd,"/tmp"),full.names=T))
    }

    #make tmp bams from filtered cell ids
    mclapply(1:nrow(cellid_list), function(i) filter_bam_tmp(cellid_list=cellid_list,idx=i),mc.cores=task_cpus)
    
    #merge bams in tmp directory then remove them
    system(paste0("mkdir -p ",output_directory,"/",output_sample))
    bam_merge=paste0(wd,"/tmp/merged.tmp.bam")
    bam_out=paste0(output_directory,"/",output_sample,"/",output_name,".bam")

    # Convert the list to a flat character vector and write to file
    bam_merge_list=paste0(wd,"/tmp/tmp.bam.merge.txt")
    writeLines(unlist(list.files(path=paste0(wd,"/tmp"),pattern="sorted.bam$",full.names=T)), bam_merge_list)
    system(paste0("samtools cat -b ",bam_merge_list," -o ",bam_merge))
    print(paste("Wrote out merged bam for",output_name))

    #samtools sort
    print(paste("Sorting and indexing...",output_name))
    system(paste0("samtools sort -m 2G --write-index -@ ",task_cpus," -o ",bam_out, " ",bam_merge))
    print(paste("Done! Final bam here:", bam_out))
    
    #clean up tmp bam files
    if(length(list.files(paste0(wd,"/tmp")))>0){
        file.remove(list.files(paste0(wd,"/tmp"),full.names=T))
    }
   
}


obj<-readRDS(file=paste(project_data_directory,"03_fine_celltyping","03_scaledcis.final_celltypes.amethyst.rds",sep="/"))
obj<-subsetObject(obj,cells=row.names(obj@metadata[!(obj@metadata$celltype %in% c("stromal_unknown")),]))
obj<-subsetObject(obj,cells=row.names(obj@metadata[!(obj@metadata$cnv_clonename_500kb %in% c("NA")),]))

#filter to samples with at least one clone

clones <- table(obj@metadata$cnv_clonename_500kb)
count_clones_per_sample <-table(unlist(lapply(strsplit(names(clones),"_"),"[",1)))
samples_passing_filter <- names(count_clones_per_sample[which(count_clones_per_sample>1)])

clones_passing_filter <- names(clones)[unlist(lapply(strsplit(names(clones),"_"),"[",1)) %in% samples_passing_filter]

obj<-subsetObject(obj,cells=row.names(obj@metadata[obj@metadata$cnv_clonename_500kb %in% clones_passing_filter,]))

cellid_list_clones<-obj@metadata %>% group_by(cnv_clonename_500kb) %>% select(cellid,bam_path,cnv_clonename_500kb) %>% group_split()

for(i in 1:length(cellid_list_clones)){
    cellid_list<-as.data.frame(cellid_list_clones[[i]])
    clone_name=cellid_list_clones[[i]]$cnv_clonename_500kb[1]
    sample_name=unlist(lapply(strsplit(clone_name,"_"),"[",1))
    print(paste("Running...",sample_name,clone_name))
    output_merged_bam(cellid_list=cellid_list,
                        output_name=clone_name,
                        output_sample=sample_name,
                        output_directory=paste0(wd,"/bam_merged_clones"),
                        task_cpus=300)
}

#note that diploid includes all cell types. 
#if the goal is to look for somatic events within the lumhr preceeding clonal changes, we can also filter to lumhr diploid

````
bam="BCMDCIS41T.3B04.dedup.bam"
bsbolt BamIndex \
-I $bam

bsbolt CallVariation \
-t 10 \
-I $bam \
-DB /data/rmulqueen/projects/scalebio_dcis/ref/hg38_bsbolt \
-O ${bam::-4}

#make temp bams per clone (and somatic by non-epithelial cells) and call variants
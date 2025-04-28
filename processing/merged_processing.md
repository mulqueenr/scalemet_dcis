```bash
#run per line 
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif

cd /data/rmulqueen/projects/scalebio_dcis
```

```R
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
library(parallel)
library(optparse,lib.loc="/home/rmulqueen/R/x86_64-conda-linux-gnu-library/4.4") #add this

task_cpus=50
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data"
amethyst_files=list.files(path=project_data_directory,pattern="*.amethyst.rds",recursive=TRUE,full.names=TRUE)

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

read_amethyst_obj<-function(x){
    dat_tmp<-readRDS(x)
    
    #set up ref
    gtf <- rtracklayer::readGFF("/container_ref/gencode.v43.annotation.gtf.gz")
    for (i in c("gene_name", "exon_number")) {
        gtf$i <- unlist(lapply(gtf$attributes, extractAttributes, i))}
    gtf <- dplyr::mutate(gtf, location = paste0(seqid, "_", start, "_", end))
    dat_tmp@ref<-gtf

    run_id<-gsub(x,pattern=project_data_directory,replacement="")
    run_id<-strsplit(dirname(run_id),"/")[[1]][2]

    #correct h5 names
    print(paste("Appended",run_id,"to names in h5 files in",basename(x)))
    lapply(unique(dat_tmp@h5paths$paths),function(i){correct_h5_cellnames(h5=i,run_id)})
    
    #correct cell id names
    print(paste("Corrected",run_id,"metadata cell names in",basename(x)))
    row.names(dat_tmp@h5paths)<-paste0(row.names(dat_tmp@h5paths),"+",run_id)
    row.names(dat_tmp@metadata)<-paste0(row.names(dat_tmp@metadata),"+",run_id)

    #initial window genomeMatrix
    dat_tmp@genomeMatrices[["cg_100k_score"]] <- makeWindows(dat_tmp,
                                                      stepsize = 100000, type = "CG", 
                                                      metric = "score", index = "chr_cg", nmin = 2,
                                                      threads = as.numeric(task_cpus)) 
    return(dat_tmp)
}

dat_list<-for(x in amethyst_files){read_amethyst_obj(x)}


dat <- combineObject(objList = dat_list, genomeMatrices="cg_100k_score")
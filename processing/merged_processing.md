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


read_amethyst_obj<-function(x){
    dat_tmp<-readRDS(x)
    

    return(dat_tmp)
}

dat_list<-for(x in amethyst_files){read_amethyst_obj(x)}


dat <- combineObject(objList = dat_list, genomeMatrices="cg_100k_score")
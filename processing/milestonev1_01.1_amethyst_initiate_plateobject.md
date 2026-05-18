# Generate amethyst object and perform initial clustering.

Initialize each sequenced plate in amethyst.

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
#set up ref
gtf <- rtracklayer::readGFF("/data/rmulqueen/projects/scalebio_dcis/ref/gencode.v43.annotation.gtf.gz")
for (i in c("gene_name", "exon_number")) {
    gtf$i <- unlist(lapply(gtf$attributes, extractAttributes, i))}

gtf <- dplyr::mutate(gtf, location = paste0(seqid, "_", start, "_", end))
gtf<-gtf %>% filter(seqid %in% c(paste0("chr",c(1:22,"X")))) #filter to autosomes + X

#read in all sample/csv/sample.passingCellsMapMethylStats.csv files into data frame
#make a dataframe of all h5 files also <sample>\t<h5location>
in_dir="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/scalemethyl_pipeline_out"
system(paste0("mkdir -p ",in_dir,"/","amethyst_plate_obj"))

#updated indexChr for chrlist from https://github.com/lrylaarsdam/amethyst/blob/main/R/index.R
#modified h5 read in to be old format
indexChr <- function(obj,
                     type,
                     chrList = NULL,
                     threads = 1) {

  # get barcodes and paths from amethyst object
    barcodes <- row.names(obj@h5paths)
    paths <- obj@h5paths$path

  if (threads > 1) {
    future::plan(future::multicore, workers = threads)
  }

  if (is.null(chrList)) {
    cat("Notice: If a chrList is not specified, chromosomes with underscores will be removed. This is to avoid indexing unwanted data like alternative contigs.\n")
  }

  output <- list()
  output <- furrr::future_pmap(.l = list(paths, barcodes), .f = function(path, bar) {
    tryCatch({
      unique_id <- paste0(bar)
      h5 <- data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", bar)))
      h5[, index := 1:.N]
      # for each gene, filter the reads based on corresponding chromosome number and make sure the position is in between the gene's start and end.
      # The minimum row number (index) is the start and the number of rows is the count. This portion of the hdf5 file can now be quickly read in for future gene-specific functions.

      #define chr groups
      if (is.null(chrList)) {
        chrs <- as.list(unique(h5$chr[!grepl("_|EBV|M", h5$chr)]))
      } else {
        chrs <- chrList
      }

      index <- furrr::future_map(chrs,
                                 .f = function(x) {
                                   ind <- h5[chr == x]
                                   ind <- ind[, .(cell_id = unique_id, chr = x, start = min(index), count = .N)]
                                   ind
                                 }, .progress = FALSE)
      index <- do.call(rbind, index) # bind to make one data frame
    }, error = function(e) {
      cat("Error processing data for barcode", bar, ":", conditionMessage(e), "\n")
      return(NULL)
    })
  }, .progress = TRUE) # show a progress bar

  if (threads > 1) {
    future::plan(future::sequential)
    gc()
  }
  output <- do.call(rbind, output)
  output <- split(output, output$chr)
  output
}

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
        h5_file=list.files(paste(in_dir,plate_info,"samples/methylation_coverage/amethyst",sample_name,sep="/"),
                            pattern=well,
                            full.names=T,
                            include.dirs=T)
        return(h5_file)
        }))
    sample_meta$bam_path<-unlist(lapply(1:nrow(sample_meta),function(x){
        well=sample_meta[x,]$tgmt_well
        sample_name=sample_meta[x,]$sample
        bam_file=list.files(paste(in_dir,plate_info,"alignments/dedup",paste0(sample_name,".",well),sep="/"),pattern="dedup.bam",full.names=T,include.dirs=T)
    }))

    print(paste("Generated amethyst object for plate:",as.character(plate_info), "Sample:", as.character(sample_name)))
    return(sample_meta)
}

prepare_amethyst_obj<-function(sample_plate_meta,cpu_count=200){
    #sample_meta is a list of all samples per plate.

    plate_meta<-do.call("rbind",lapply(sample_plate_meta,read_plate_sample_meta))
    plate_info<-plate_meta$plate_info[1]
    obj <- createObject()

    #metadata MUST have a column called mcg_pct for score calculation
    #metadata MUST have a column called cov to regress coverage mias
    obj@metadata<-plate_meta
    row.names(obj@metadata)<-obj@metadata$cell_id

    cg_per_read<-mean(obj@metadata$cg_cov/obj@metadata$unique_reads)
    obj@h5paths <- data.frame(row.names = c(rownames(obj@metadata)), paths = obj@metadata$h5_path)

    #correct h5 names
    print(paste("Appended",plate_info,"to names in h5 files."))
    mclapply(unique(obj@h5paths$paths),function(i){correct_h5_cellnames(h5=i,plate_info)},mc.cores=cpu_count)
    h5closeAll()

    #correct cell id names
    print(paste("Corrected",plate_info,"metadata cell names."))
    row.names(obj@h5paths)<-paste0(row.names(obj@h5paths),"+",plate_info)
    row.names(obj@metadata)<-paste0(row.names(obj@metadata),"+",plate_info)

    #add ref
    obj@ref<-gtf

    # index files
    obj@index[["chr_cg"]] <- indexChr(obj, 
                                        type = "CG", 
                                        chrList=factor(paste0("chr",c(1:22,"X"))),
                                        threads = cpu_count)

    #initial window genomeMatrix
    obj@genomeMatrices[["cg_5mb_score"]] <- makeWindows(obj,
                                                      stepsize = 5000000, 
                                                      type = "CG", 
                                                      metric = "score", 
                                                      index = "chr_cg", 
                                                      nmin = 2,
                                                      threads = cpu_count) 

    print(paste("Saving amethyst object for :",as.character(plate_info)))
    saveRDS(obj,file=paste0(in_dir,"/amethyst_plate_obj/",plate_info,".amethyst.rds"))
    print(sample_plate_meta)
}

plates<-list.files(in_dir)
plates<-plates[grep(plates,pattern="_homebrew_|_scalebio_")]
plates_meta<-lapply(plates,function(plate) list.files(paste(in_dir,plate,"samples",sep="/"),pattern="*allCells.csv",full.names=T))
plates_objs<-lapply(plates_meta,prepare_amethyst_obj)

```


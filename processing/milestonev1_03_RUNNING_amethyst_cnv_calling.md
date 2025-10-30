```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
```

# Generate CopyKit for each sample

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
obj<-readRDS(file="05_scaledcis.fine_celltype.amethyst.rds")

#make output directory
system(paste0("mkdir -p ",project_data_directory,"/copykit" ))
system(paste0("mkdir -p ",project_data_directory,"/methyltree" ))

read_scalebio_bam<-function(obj_met,x,sample_name){
    #scalebio pipeline outputs bam files as Tn5 wells. so multiple cell IDs are in a bam. this function splits out the bam to the query cellid
    bam=obj_met[obj_met$sample==sample_name,]$bam_path[x]
    cellid=strsplit(row.names(obj_met)[x],"[+]batch|[+]prelim")[[1]][1]

    what <- c("qname","rname", "pos")
    param <- ScanBamParam(what=what,
                            flag=scanBamFlag(isPaired=TRUE,
                                            isProperPair=TRUE,
                                            isSecondaryAlignment=FALSE,
                                            isDuplicate=FALSE,
                                            isSupplementaryAlignment=FALSE))

    input_bam<-Rsamtools::scanBam(bam,param=param)
    input_bam<-do.call("DataFrame", input_bam)
    input_bam$cellid<-gsub("^.*:", "", input_bam$qname)
    input_bam<-input_bam[input_bam$cellid==cellid,]
    input_bam$end<-input_bam$pos+1
    input_bam<-makeGRangesFromDataFrame(input_bam,seqnames.field="rname",start.field="pos",end.field="end")
    return(input_bam)
}

runCountReads_amethyst <- function(obj,
                        sample_name,
                        genome = "hg38",
                        resolution = c(
                                        "220kb",
                                        "55kb",
                                        "110kb",
                                        "195kb",
                                        "280kb",
                                        "500kb",
                                        "1Mb",
                                        "2.8Mb"),
                        remove_Y = TRUE,
                        min_bincount = 10,
                        cores=100) {
    output_directory=paste0(project_data_directory,"/copykit/",sample_name)
    system(paste0("mkdir -p ",project_data_directory,"/copykit/",sample_name))

    resolution <- match.arg(resolution)
    #resolution="220kb"

    # bindings for NSE and data
    Chr <- chr <- strand <- GeneID <- NULL
    reads_assigned_bins <- reads_duplicates <- reads_total <- NULL

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # genomic ranges (varbin scaffolds)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Reading hg38 VarBin ranges
    hg38_grangeslist <- hg38_grangeslist

    hg38_rg <- switch(resolution,
        "55kb" = hg38_grangeslist[["hg38_50kb"]],
        "110kb" = hg38_grangeslist[["hg38_100kb"]],
        "195kb" = hg38_grangeslist[["hg38_175kb"]],
        "220kb" = hg38_grangeslist[["hg38_200kb"]],
        "280kb" = hg38_grangeslist[["hg38_250kb"]],
        "500kb" = hg38_grangeslist[["hg38_500kb"]],
        "1Mb" = hg38_grangeslist[["hg38_1Mb"]],
        "2.8Mb" = hg38_grangeslist[["hg38_2Mb"]]
    )

    hg38_rg <- as.data.frame(hg38_rg)

    rg <- hg38_rg %>%
        dplyr::rename(chr = "seqnames") %>%
        dplyr::mutate(GeneID = 1:nrow(hg38_rg))

    if (remove_Y == TRUE) {
        rg <- dplyr::filter(rg,chr != "chrY")
    }

    message("Counting reads for genome ",genome," and resolution: ",resolution)

    #get list of bams and cellids
    obj_met<-obj@metadata[obj@metadata$Sample==sample_name,]
    #return chr start position for reads filtered in bam to cell id
    varbin_counts_list_all_fields<-mclapply(
                                        1:nrow(obj_met), 
                                        function(i) 
                                        read_scalebio_bam(obj_met=obj_met,x=i,sample_name=sample_name), 
                                        mc.cores=cores)

    message("Read in all bam files.")

    names(varbin_counts_list_all_fields)<- row.names(obj_met)
    varbin_counts_list_all_fields<-as(varbin_counts_list_all_fields, "GRangesList")
    ref<-as(rg,"GRanges")

    varbin_counts_list <-mclapply(varbin_counts_list_all_fields,
                                    function(x) 
                                    GenomicRanges::countOverlaps(
                                    query=ref,
                                    subject=x,
                                    type="any",
                                    ignore.strand=TRUE))
    message("Counted reads across all bins.")

    varbin_counts_list <- lapply(varbin_counts_list,as.vector)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # filtering for minimal mean bin count
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # obtaining the index of the ones that FAIL to meet the min_bincount arg
    min_bc <- which(vapply(varbin_counts_list, mean, numeric(1)) < min_bincount)
    # subsetting counts list and main counts list

    if (length(min_bc) > 0) {
        varbin_counts_list <- varbin_counts_list[-min_bc]
        varbin_counts_list_all_fields <- varbin_counts_list_all_fields[-min_bc]
        message(
            length(min_bc), " bam files had less than ", min_bincount,
            " mean bincounts and were removed."
        )
    }

    # LOWESS GC normalization

    message("Performing GC correction.")

    varbin_counts_list_gccor <-
        mclapply(varbin_counts_list, function(x) {
            gc_cor <- lowess(rg$gc_content, log(x + 1e-3), f = 0.05)
            gc_cor_z <- approx(gc_cor$x, gc_cor$y, rg$gc_content)
            exp(log(x) - gc_cor_z$y) * median(x)
        },mc.cores=cores
        )

    varbin_counts_df <- round(dplyr::bind_cols(varbin_counts_list_gccor), 2)

    # filtering low read counts where the sum of bins does not reach more than 0
    good_cells <- names(varbin_counts_df[which(colSums(varbin_counts_df) != 0)])

    varbin_counts_df <- varbin_counts_df[good_cells]

    rg <- rg %>%
        dplyr::select(-strand, -GeneID)

    rg_gr <- GenomicRanges::makeGRangesFromDataFrame(rg,
        ignore.strand = TRUE,
        keep.extra.columns = TRUE)


    cna_obj <- CopyKit(
        assays = list(bincounts = varbin_counts_df),
        rowRanges = rg_gr)

    # Adding genome and resolution information to metadata
    S4Vectors::metadata(cna_obj)$genome <- genome
    S4Vectors::metadata(cna_obj)$resolution <- resolution

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sun Feb 14 20:55:01 2021
    # ADDING READS METRICS TO METADATA
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sun Feb 14 20:55:24 2021


    # saving info and removing columns from list elements
    bam_metrics <- obj_met[c("unique_reads","tss_enrich","mcg_pct","cg_cov","batch","plate_info","tgmt_well","i7_well","i5_well","broad_celltype")]

    # making sure metrics match varbin_counts_df
    bam_metrics <- bam_metrics[good_cells,]

    bam_metrics$sample <- rownames(bam_metrics)
    bam_metrics$sample_name=sample_name
    bam_metrics$reads_assigned_bins <- colSums(varbin_counts_df)

    # adding to metadata
    SummarizedExperiment::colData(cna_obj) <-
        S4Vectors::DataFrame(bam_metrics)
    colnames(cna_obj) <- names(varbin_counts_df)
    #runvarbin module
    cna_obj <- runVst(cna_obj)
    cna_obj <- runSegmentation(cna_obj)
    cna_obj <- logNorm(cna_obj)

    # Mark euploid cells if they exist
    cna_obj <- findAneuploidCells(cna_obj)

    # Mark low-quality cells for filtering
    cna_obj <- findOutliers(cna_obj)

    # kNN smooth profiles
    cna_obj <- knnSmooth(cna_obj,k=10)

    # adds basic quality control information to colData
    cna_obj <- runMetrics(cna_obj)
    cna_obj <- runUmap(cna_obj)
    cna_obj <- findSuggestedK(cna_obj)
    S4Vectors::metadata(cna_obj)$suggestedK
    cna_obj <- findClusters(cna_obj)
    cna_obj <- calcConsensus(cna_obj)
    cna_obj <- runConsensusPhylo(cna_obj)

    plt_umap<-plotUmap(cna_obj,label="subclones")
    logr_heatmap<-plotHeatmap(
        cna_obj,
        assay = "logr",
        label = c('subclones','mcg_pct','plate_info','reads_assigned_bins','broad_celltype'),
        genes = c("TP53", "BRAF", "MYC"),
        order_cells = 'consensus_tree',
        n_threads=cores)
    segratio_heatmap<-plotHeatmap(
        cna_obj,
        assay = "segment_ratios",
        label = c('subclones','mcg_pct','plate_info','reads_assigned_bins'),
        genes = c("TP53", "BRAF", "MYC"),
        order_cells = 'consensus_tree',
        n_threads=cores)

    pdf(paste0(output_directory,"/copykit.",sample_name,".",resolution,".pdf"))
    print(plt_umap)
    print(logr_heatmap)
    print(segratio_heatmap)
    dev.off()
    print(paste("Plotted... ",paste0(output_directory,"/copykit.",sample_name,".",resolution,".pdf")))

    saveRDS(cna_obj,file=paste0(output_directory,"/copykit",sample_name,".",resolution,".rds"))
    return(cna_obj)
}


sample_list<-sort(unique(obj@metadata$Sample))

lapply(sample_list[4:length(sample_list)],function(x) runCountReads_amethyst(obj=obj,sample_name=x,resolution="220kb"))



list.files(paste0(project_data_directory,"/copykit"))
```
# Read all CopyKit RDS objects and plot together

```R
library(ComplexHeatmap)

copykit_output<-list.files(path=paste0(project_data_directory,"/copykit"),
    recursive=TRUE,full.names=TRUE,pattern="*rds")

#read in all meta data from copykit
read_meta_copykit<-function(x){
    tmp<-readRDS(x)
    meta<-as.data.frame(tmp@colData[c("sample_name","reads_assigned_bins","plate_info","subclones")])
    return(meta)
}
cnv_meta<-do.call("rbind",lapply(copykit_output,read_meta_copykit))

#read in all logr from copykit
read_logr_copykit<-function(x){
    tmp<-readRDS(x)
    logr<-tmp@assays@data$logr
    return(logr)
}
cnv_logr<-do.call("cbind",lapply(copykit_output,read_logr_copykit))

#get 220kb windows ranges
copykit<-readRDS(copykit_output[1])
windows<-copykit@rowRanges

cnv_col<-c("0"="#002C3E", "0.5"="#78BCC4", "1"="#F7F8F3", "1.5"="#F7444E", "2"="#aa1407", "3"="#440803")

#relevant CNV genes from curtis work
#from https://www.nature.com/articles/s41416-024-02804-6#Sec20
#change RAB7L1 to RAB29
#lost RAB7L1
cnv_genes<-c('ESR1','PGR','DLEU2L', 'TRIM46', 'FASLG', 'KDM5B', 'RAB7L1', 'PFN2', 'PIK3CA', 'EREG', 'AIM1', 'EGFR', 'ZNF703', 'MYC', 'SEPHS1', 'ZMIZ1', 'EHF', 'POLD4', 'CCND1', 'P2RY2', 'NDUFC2-KCTD14', 'FOXM1', 'MDM2', 'STOML3', 'NEMF', 'IGF1R', 'TP53I13', 'ERBB2', 'SGCA', 'RPS6KB1', 'BIRC5', 'NOTCH3', 'CCNE1', 'RCN3', 'SEMG1', 'ZNF217', 'TPD52L2', 'PCNT', 'CDKN2AIP', 'LZTS1', 'PPP2R2A', 'CDKN2A', 'PTEN', 'RB1', 'CAPN3', 'CDH1', 'MAP2K4', 'GJC2', 'TERT', 'RAD21', 'ST3GAL1', 'SOCS1')
cnv_genes_class<-c('amp','amp','amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'amp', 'del', 'del', 'del', 'del', 'del', 'del', 'del', 'del', 'del', 'amp', 'amp', 'amp', 'amp', 'amp')
cnv_genes<-setNames(cnv_genes_class,nm=cnv_genes)

#use gtf file to get gene locations
gtf_file="/container_ref/gencode.v43.annotation.gtf.gz"

gtf <- rtracklayer::readGFF(gtf_file)
gtf<- gtf %>% 
    filter(type=="gene" & gene_type %in% c("protein_coding")) %>% 
    filter(gene_name %in% names(cnv_genes))

cnv_genes_windows<-gtf[gtf$gene_name %in% names(cnv_genes),] #filter annotation to genes we want
cnv_genes_windows<-cnv_genes_windows[!duplicated(cnv_genes_windows$gene_name),] #remove duplicates
cnv_genes_windows$cnv_gene_class<-unname(cnv_genes[match(cnv_genes_windows$gene_name, names(cnv_genes))]) #add amp/del

cnv_genes_windows<-makeGRangesFromDataFrame(cnv_genes_windows,keep.extra.columns=TRUE) #make granges
wind<-GenomicRanges::findOverlaps(windows,cnv_genes_windows,select="first") #do overlap to get window indexes
wind<-as.data.frame(wind)
wind<-wind[!duplicated(wind$subjectHits),] #take first subject hit

annot<-data.frame(
  window_loc=wind$queryHits,
  gene=cnv_genes_windows$gene_name[wind$subjectHits],
  cnv_class=cnv_genes_windows$cnv_gene_class[wind$subjectHits])

annot$col<-ifelse(annot$cnv_class=="amp","red","blue")

hc = columnAnnotation(common_cnv = anno_mark(at = annot$window_loc, 
                        labels = annot$gene,
                        which="column",side="bottom",
                        labels_gp=gpar(col=annot$col)))

hr = rowAnnotation(sample=cnv_meta$sample,
                    subclones=cnv_meta$subclones,
                    reads=anno_barplot(log10(cnv_meta$reads_assigned_bins)))

#### CNV call output folder
output_directory=paste0(dirname(getwd()),"/copykit")

pdf(paste0(output_directory,"/","all_met.cnv.heatmap.pdf"),height=90,width=40)
Heatmap(t(cnv_logr),
  #col=cnv_col,
  cluster_columns=FALSE,
  cluster_rows=TRUE,
  show_row_names = FALSE, row_title_rot = 0,
  show_column_names = FALSE,
  cluster_row_slices = TRUE,
  bottom_annotation=hc,
  left_annotation=hr,
  row_split=paste(cnv_meta$sample_name,cnv_meta$subclones),
  column_split=seqnames(windows),
  border = TRUE)
dev.off()
print(paste0(output_directory,"/","all_met.cnv.heatmap.pdf"))

```


# Output files for methyltree format

```R
methyltree_output<-function(obj=obj,
                            sample_name="DCIS-41T",
                            filt_min_pct=10,
                            filt_max_pct=80,
                            threads=1){
        
        output_directory=paste0(project_data_directory,"/methyltree/",sample_name)
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
    out_metadata<-obj_met@metadata[,c("pass","fine_celltype","cg_cov","mcg_pct","subclones")]
    colnames(out_metadata)<-c("HQ","celltype","nCG","met_rate","large_clone_id") #match names
    out_metadata$sample<-row.names(out_metadata) #sample (cell) names
    out_metadata$met_rate<-out_metadata$met_rate/100 #percentage to rate

    methyltree_input_file=paste(output_directory,
        paste("methyltree",sample_name,"_methyltree_input.h5",sep="."),sep="/")
    if(file.exists(methyltree_input_file)){
        system(paste0("rm -rf ",methyltree_input_file))
    }
      h5createFile(file=methyltree_input_file)
      h5write(methyltreeoutput,file=methyltree_input_file,name="data")
      h5write(out_metadata,file=methyltree_input_file,name="metadata")
    }

lapply(unique(obj@metadata$sample),function(x) 
methyltree_output(obj=obj,
                   sample_name=x,
                    filt_min_pct=10,
                    filt_max_pct=80,
                    threads=1))

```

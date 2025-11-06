```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
```

# Generate CopyKit for each sample

```R
source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
library(Rsamtools)
library(copykit)
library(circlize)
detach("package:GeneNMF",unload=TRUE)
library(dendextend)
library(RColorBrewer)
library(ComplexHeatmap)

#download cytoband data
#system("wget -P /data/rmulqueen/projects/scalebio_dcis/ref https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz ")
#system("gzip -d /data/rmulqueen/projects/scalebio_dcis/ref/cytoBand.txt.gz")
cyto=read.table(file="/data/rmulqueen/projects/scalebio_dcis/ref/cytoBand.txt",sep="\t")
colnames(cyto)<-c("chr","start","end","band","stain")
table(cyto$stain) #set colors for these

#set colors
celltype_col=c(
'peri'='#c1d552',
'fibro1'='#7f1911',
'fibro2'='#e791f9',
'endo'='#f0b243',
'endo2'='#d0bd4a',
'tcell'='#2e3fa3',
'bcell'='#00adea',
'myeloid1'='#00a487',
'myeloid2'='#006455',
'basal'='#7200cc',
'lumsec'='#af00af',
'lumhr'='#d8007c')

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
    bam=obj_met[obj_met$Sample %in% c(sample_name),]$bam_path[x]
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
    output_directory=paste0(project_data_directory,"/copykit/",sample_name[1])
    system(paste0("mkdir -p ",project_data_directory,"/copykit/",sample_name[1]))

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
    obj_met<-obj@metadata[obj@metadata$Sample %in% sample_name,]
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
    bam_metrics <- obj_met[c("unique_reads","tss_enrich","mcg_pct","cg_cov","batch","plate_info","tgmt_well","i7_well","i5_well","fine_celltype")]

    # making sure metrics match varbin_counts_df
    bam_metrics <- bam_metrics[good_cells,]

    bam_metrics$sample <- rownames(bam_metrics)
    bam_metrics$sample_name=sample_name[1]
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

    #add own clusters based on logr calling (based on hclust of logr heatmap)
    #filters to top 20% of variable sites then clusters on those
    # var_df <- t(cna_obj@assays@data$logr) %>% as.data.frame() %>%
    #         summarise(across(.cols = everything(), .fns = var)) %>% 
    #         pivot_longer(cols = everything()) %>% 
    #         slice_max(n = round(nrow(cna_obj@assays@data$logr)/10), order_by = value) %>%
    #         mutate(name=gsub("V","",name))

    #dend <- t(cna_obj@assays@data$logr[row.names(cna_obj@assays@data$logr) %in% var_df$name,]) %>% 
    #        dist(method="manhattan") %>% hclust(method="ward.D2") %>% as.dendrogram

    dend <- t(cna_obj@assays@data$logr) %>% 
            dist(method="manhattan") %>% hclust(method="ward.D2") %>% as.dendrogram
    k_optimal=find_k(dend, krange = 2:10)
    print(paste("optimal k value for cutting hclust:", k_optimal$k))
    superclones=dendextend::cutree(dend,k=k_optimal$k+2)
    subclones=dendextend::cutree(dend,k=k_optimal$k+5)
    cna_obj@colData$subclones<-subclones[row.names(cna_obj@colData)]
    cna_obj@colData$superclones<-superclones[row.names(cna_obj@colData)]

    #define colors based on data
    log_col=colorRamp2(c(-2,-1,0,1,2),
                            c("darkblue","blue","white","red","darkred"))
    cg_perc_col=colorRamp2(c(40,60,80,100),
                            c("#4d2d18","#CABA9C","#4C6444","#102820"))
    reads_col=colorRamp2(c(min(log10(cna_obj@colData$unique_reads)),
                            max(log10(cna_obj@colData$unique_reads))),
                            c("white","black"))

    superclone_col=setNames(nm=unique(as.character(cna_obj@colData$superclones)),
                            colorRampPalette(brewer.pal(9, "Pastel1"))(length(unique(as.character(cna_obj@colData$superclones)))))
    subclone_col=setNames(nm=unique(as.character(cna_obj@colData$subclones)),
                            colorRampPalette(brewer.pal(9, "Spectral"))(length(unique(as.character(cna_obj@colData$subclones)))))
    
    #plot heatmap
    ha = rowAnnotation(
        reads=log10(cna_obj@colData$unique_reads),
        cg_perc=cna_obj@colData$mcg_pct,
        celltype=cna_obj@colData$fine_celltype,
        superclones=as.character(cna_obj@colData$superclones),
        subclones=as.character(cna_obj@colData$subclones),
        col= list(
            celltype=celltype_col,
            reads=reads_col,
            cg_perc=cg_perc_col,
            superclones=superclone_col,
            subclones=subclone_col
        ))

    cyto_overlap<-GenomicRanges::findOverlaps(cna_obj@rowRanges,
                                                makeGRangesFromDataFrame(cyto,keep=TRUE),
                                                select="first")
    cna_obj@rowRanges$stain <- cyto[cyto_overlap,]$stain

    arm_col=c("p"="grey","q"="darkgrey")
    band_col=c("acen"="#99746F","gneg"="white","gpos100"="black","gpos25"="lightgrey","gpos50"="grey","gpos75"="darkgrey","gvar"="#446879")

    column_ha = HeatmapAnnotation(
        arm = cna_obj@rowRanges$arm,
        band = cna_obj@rowRanges$stain,
        col=list(arm=arm_col,band=band_col))

    plt<-Heatmap(
        t(cna_obj@assays@data$logr),
        left_annotation=ha,col=log_col,
        row_split=as.character(cna_obj@colData$superclones),
        show_column_names=FALSE,show_row_names=FALSE,
        top_annotation=column_ha,cluster_columns=FALSE,cluster_column_slices=FALSE,column_split=seqnames(cna_obj@rowRanges),
        name="logr")

    pdf(paste0(output_directory,"/copykit.",sample_name[1],".",resolution,".pdf"),width=20)
    print(plt)
    dev.off()

    print(paste("Plotted... ",paste0(output_directory,"/copykit.",sample_name[1],".",resolution,".pdf")))
    cna_obj <- calcConsensus(cna_obj)
    cna_obj <- runConsensusPhylo(cna_obj)
    plt_umap<-plotUmap(cna_obj,label="subclones")
    saveRDS(cna_obj,file=paste0(output_directory,"/copykit",".",sample_name[1],".",resolution,".rds"))
    return(cna_obj)
}

runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS05T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS07T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS102T_24hTis'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS124T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS22T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS28T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS32T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS35T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS41T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS49T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS52T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS65T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS66T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS70T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS74T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS79T_24hTis_DCIS','BCMDCIS79T_24hTis_IDC'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS80T_24hTis'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS82T_24hTis'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS92T_24hTis'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS94T_24hTis'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS97T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMDCIS99T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA03R'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA04R'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA09R-3h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA12R-3h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA16R-3h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA17R-3h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA19R-4h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA22R-4h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA26L-24hTis-4h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA29L-2h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA38L-3h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA83L-3h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('BCMHBCA85L-3h'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('ECIS25T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('ECIS26T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('ECIS36T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('ECIS48T'),resolution='220kb')
runCountReads_amethyst(obj=obj,sample_name=c('ECIS57T'),resolution='220kb')

#assign aneuploidy based on CNV profiles.

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
    meta<-as.data.frame(tmp@colData[c("sample_name","reads_assigned_bins","plate_info","subclones","broad_celltype")])
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
wind<-GenomicRanges::findOverlaps(windows,cnv_genes_windows) #do overlap to get window indexes
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
                    celltype=cnv_meta$broad_celltype,
                    subclones=cnv_meta$subclones,
                    reads=anno_barplot(log10(cnv_meta$reads_assigned_bins)))

#### CNV call output folder
setwd(project_data_directory)
output_directory=paste0(project_data_directory,"/copykit")

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
    out_metadata<-obj_met@metadata[,c("pass","fine_celltype","cg_cov","mcg_pct","subclones")]
    colnames(out_metadata)<-c("HQ","celltype","nCG","met_rate","large_clone_id") #match names
    out_metadata$sample<-row.names(out_metadata) #sample (cell) names
    out_metadata$met_rate<-out_metadata$met_rate/100 #percentage to rate

    methyltree_input_file=paste(output_directory,
        paste("methyltree.",sample_name[1],"_methyltree_input.h5",sep="."),sep="/")
    if(file.exists(methyltree_input_file)){
        system(paste0("rm -rf ",methyltree_input_file))
    }
      h5createFile(file=methyltree_input_file)
      h5write(methyltreeoutput,file=methyltree_input_file,name="data")
      h5write(out_metadata,file=methyltree_input_file,name="metadata")
    }

#making methyltree output
#note i think some samples could use improvements on subclones, but for now just running

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
methyltree_output(obj=obj,sample_name=c('BCMDCIS79T_24hTis_DCIS','BCMDCIS79T_24hTis_IDC'),threads=1)
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
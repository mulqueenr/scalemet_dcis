Not quite sure what to do here. Running PCA on CNV events to group common events into a single continuous metric (the PC loading).

- Make barplot per PC loading for aneuploid cells (see if sample specific)
- Bin cells by PC loading and do differential methylation?

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


#set environment and read in data
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)
obj<-readRDS(file="06_scaledcis.cnv_clones.amethyst.rds")

#read in cyto info
cyto=read.table(file="/data/rmulqueen/projects/scalebio_dcis/ref/cytoBand.txt",sep="\t")
colnames(cyto)<-c("chr","start","end","band","stain")
cyto$arm<-substring(cyto$band, 1, 1)
cyto<-cyto[!is.na(cyto$band),]
cyto<-cyto[cyto$chr %in% c(paste0("chr",1:22),"chrX"),]
table(cyto$stain) #set colors for these

copykit_output<-list.files(path=paste0(project_data_directory,"/copykit"),recursive=TRUE,full.names=TRUE,pattern="*rds")
copykit_output<-copykit_output[!grepl(copykit_output,pattern="diploid")]
cna_obj<-readRDS(copykit_output[1]) #just to grab row ranges
output_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/copykit/"
#read in all meta data from copykit
read_meta_copykit<-function(x){
    tmp<-readRDS(x)
    meta<-as.data.frame(tmp@colData[c("sample_name","reads_assigned_bins","plate_info","subclones","fine_celltype","clonename","ploidy")])
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

hc = columnAnnotation(common_cnv = anno_mark(at = annot$window_loc, 
                        labels = annot$gene,
                        which="column",side="bottom",
                        labels_gp=gpar(col=annot$col)))

#define colors based on data
col=colorRamp2(log(c(0,0.02,0.04)),
                        c("white","red","darkred"))


#run PCA on cnv_logr
library(umap)

#PCA on windows to determine correlations and shared events
pca_result <- prcomp(t(cnv_logr), center = FALSE, scale = FALSE)

#define PC cutoff by elbow
plt<-ggplot()+geom_line(aes(y=pca_result$sdev[1:100],x=1:100))+theme_minimal()
ggsave(paste0(output_directory,"/","cnv_loadings.elbowplot.pdf"))
#i'm going to say 1:25

pdf(paste0(output_directory,"/","cnv_loadings.heatmap.pdf"),height=20,width=40)
Heatmap(t(pca_result$rotation[,1:25]),
  col=col,
  cluster_columns=FALSE,
  cluster_rows=TRUE,
  show_row_names = TRUE, row_title_rot = 0,
  show_column_names = FALSE,
  cluster_row_slices = TRUE,
  bottom_annotation=hc,
  top_annotation=column_ha,
  column_split=seqnames(windows),
  border = TRUE)
dev.off()
print(paste0(output_directory,"/","cnv_loadings.heatmap.pdf"))


#add CNV_PC loadings as metadata to amethyst object
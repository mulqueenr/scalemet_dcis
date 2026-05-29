Not quite sure what to do here. Running PCA on CNV events to group common events into a single continuous metric (the PC loading).

- Make barplot per PC loading for aneuploid cells (see if sample specific)
- Bin cells by PC loading and do differential methylation?

# Generate CopyKit for each sample

```R
library(Rsamtools)
library(copykit)
library(circlize)
library(dendextend)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggplot2)
library(seriation)
library(dplyr)
library(umap)

#set environment and read in data
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)
obj<-readRDS(file="07_scaledcis.cnv_clones.amethyst.rds")

#read in cyto info
cyto=read.table(file="/data/rmulqueen/projects/scalebio_dcis/ref/cytoBand.txt",sep="\t")
colnames(cyto)<-c("chr","start","end","band","stain")
cyto$arm<-substring(cyto$band, 1, 1)
cyto<-cyto[!is.na(cyto$band),]
cyto<-cyto[cyto$chr %in% c(paste0("chr",1:22),"chrX"),]
table(cyto$stain) #set colors for these

copykit_output<-list.files(path=paste0(project_data_directory,"/copykit"),recursive=TRUE,full.names=TRUE,pattern="*.500kb.rds")
copykit_output<-copykit_output[!grepl(copykit_output,pattern="diploid")]
cna_obj<-readRDS(copykit_output[1]) #just to grab row ranges
output_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/copykit/"

#read in methylation per logr data set
met_cnv<-readRDS(file=paste0(output_directory,"/","07_scaledcis.cnv_windows_methylation.correlation.rds"))
met_cnv<-met_cnv %>% group_by(window) %>% summarize(cor(met,cnv)) %>% as.data.frame()
colnames(met_cnv)<-c("window","cor")

cnv_meta<-readRDS(file=paste0(output_directory,"/","07_scaledcis.cnv_windows_meta.matrix.rds"))

#get 500kb windows ranges
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
gtf_file="/data/rmulqueen/projects/scalebio_dcis/ref/gencode.v43.annotation.gtf.gz"

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

#pHaplo and pTriplo
gene_dosage<-read.table("/data/rmulqueen/projects/scalebio_dcis/ref/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz")
colnames(gene_dosage)<-c("gene","pHaplo","pTriplo")
#assign window scores by gene overlap
gene_dosage<-rtracklayer::readGFF(gtf_file) %>% 
    filter(type=="gene" & gene_type %in% c("protein_coding")) %>% 
    mutate(gene=gene_name) %>%
    filter(gene %in% gene_dosage$gene) %>%
    inner_join(gene_dosage,by="gene")
wind_dosage<-GenomicRanges::findOverlaps(windows,GRanges(gene_dosage)) #do overlap to get window indexes
wind_dosage<-as.data.frame(wind_dosage)
wind_dosage<-wind_dosage[!duplicated(wind_dosage$subjectHits),] #take first subject hit
gene_dosage$window<-NA
gene_dosage[wind_dosage$subjectHits,]$window<-wind_dosage$queryHits
gene_dosage<-gene_dosage[!is.na(gene_dosage$window),]
gene_dosage<-gene_dosage %>% group_by(window) %>% summarize(pHaplo=mean(pHaplo,na.rm=T),pTriplo=mean(pTriplo,na.rm=T))

cna_obj@rowRanges$pHaplo<-NA
cna_obj@rowRanges[gene_dosage$window,]$pHaplo<-gene_dosage$pHaplo

cna_obj@rowRanges$pTriplo<-NA
cna_obj@rowRanges[gene_dosage$window,]$pTriplo<-gene_dosage$pTriplo

arm_col=c("p"="grey","q"="darkgrey")
band_col=c("acen"="#99746F","gneg"="white","gpos100"="black","gpos25"="lightgrey","gpos50"="grey","gpos75"="darkgrey","gvar"="#446879")

cor_col=colorRamp2(c(-0.5,0,0.5),c("darkblue","white","darkred"))

pHaplo_col=colorRamp2(c(0,0.5,0.8,1),c("white","white","blue","darkblue"))
pTriplo_col=colorRamp2(c(0,0.5,0.8,1),c("white","white","red","darkred"))

column_ha = HeatmapAnnotation(
                            arm = cna_obj@rowRanges$arm,
                            band = cna_obj@rowRanges$stain,
                            pHaplo=cna_obj@rowRanges$pHaplo,
                            pTriplo=cna_obj@rowRanges$pTriplo,
                            cnv_met_cor=met_cnv$cor,
                            col=list(arm=arm_col,band=band_col,pHaplo=pHaplo_col,pTriplo=pTriplo_col,cnv_met_cor=cor_col))

hc = columnAnnotation(common_cnv = anno_mark(at = annot$window_loc, 
                        labels = annot$gene,
                        which="column",side="bottom",
                        labels_gp=gpar(col=annot$col)))

#run PCA on cnv_logr
cna_obj@rowRanges
cnv_logr<-readRDS(file=paste0(output_directory,"/","07_scaledcis.cnv_windows_logr.matrix.rds"))
row.names(cnv_logr)<-paste(seqnames(cna_obj@rowRanges),start(cna_obj@rowRanges),end(cna_obj@rowRanges),sep="_")
cnv_logr<-cnv_logr[colnames(cnv_logr) %in% row.names(obj@metadata[!endsWith(obj@metadata$cnv_clonename,suffix="diploid") & (obj@metadata$cnv_clonename!="NA"),])]

dat<-t(cnv_logr)
dat <- dat %>% mutate_if(is.character, as.numeric)

#PCA on windows to determine correlations and shared events
pca_result <- irlba::prcomp_irlba(as.matrix(dat),n=50)

#define PC cutoff by elbow
plt<-ggplot()+geom_line(aes(y=pca_result$sdev[1:50],x=1:50))+geom_vline(xintercept=10,color="red")+theme_minimal()
ggsave(paste0(output_directory,"/","cnv_loadings.elbowplot.pdf"))

#define colors based on data
col=colorRamp2(c(-2,0,2),c("blue","white","red"))
mat<-scale(t(pca_result$rotation[,1:10]))
#o1 = seriate(dist(mat), method = "TSP")


loading_col2=colorRamp2(breaks=seq(from = -3, to = 3, length.out = 11),colors=rev(c("#67001FFF", "#B2182BFF", "#D6604DFF", "#F4A582FF", "#FDDBC7FF", "#F7F7F7FF", "#D1E5F0FF", "#92C5DEFF", "#4393C3FF", "#2166ACFF", "#053061FF")))

plt1<-Heatmap(mat[c("PC1","PC2"),],
  cluster_columns=FALSE,
  col=loading_col2,
  clustering_distance_rows = "maximum",
  show_row_names = TRUE, row_title_rot = 0,
  show_column_names = FALSE,
  bottom_annotation=hc,
  top_annotation=column_ha,
  column_split=seqnames(windows),
  border = TRUE,
  width = 10)
  

#add CNV_PC loadings as metadata to amethyst object
met<-obj@metadata
pca_res<-as.data.frame(pca_result$x[,1:10])
row.names(pca_res)<-row.names(dat)
met<-cbind(met[row.names(pca_res),],pca_res)

met<-met %>% group_by(cnv_clonename)

clone_pc<-as.data.frame(met %>% select(starts_with("PC")) %>% summarize_all(mean, na.rm = TRUE))
row.names(clone_pc)<-clone_pc[,1]
clone_pc<-clone_pc[,2:ncol(clone_pc)]
loading_col=colorRamp2(c(-5,0,5),c("blue","white","red"))


#https://emilhvitfeldt.github.io/r-color-palettes/discrete/NatParksPalettes/Acadia/
loading_col=colorRamp2(breaks=seq(from = -10, to = 10, length.out = 11),colors=rev(c("#7F3B08FF", "#B35806FF", "#E08214FF", "#FDB863FF", "#FEE0B6FF", "#F7F7F7FF", "#D8DAEBFF", "#B2ABD2FF", "#8073ACFF", "#542788FF", "#2D004BFF")))

plt2<-Heatmap(t(clone_pc),
  cluster_columns=TRUE,
  col=loading_col,
  cluster_rows=FALSE,
  show_row_names = TRUE, row_title_rot = 0,
  show_column_names = TRUE,
  clustering_distance_columns = "maximum",
  border = TRUE, column_names_gp = gpar(fontsize = 5),
  width=5)

pdf(paste0(output_directory,"/","clone_cnv_loadings.heatmap.pdf"),height=10,width=20)
plt2
plt1
dev.off()
print(paste0(output_directory,"/","clone_cnv_loadings.heatmap.pdf"))

umap<-umap(clone_pc)
umap<-as.data.frame(umap$layout)
umap$label<-row.names(umap)
plt<-ggplot(umap,aes(x=V1,y=V2,label=row.names(umap)))+geom_text()+theme_void()

ggsave(plt,file=paste0(output_directory,"/","cnv_loadings.umap.pdf"))

#to be added to cell level metadata to associate with other factors
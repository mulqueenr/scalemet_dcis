Set up infercnv with singularity
```bash
cd /home/rmulqueen/singularity
wget https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/InferCNV/infercnv-1.20.0.simg
```

Run locally
Note, went forward with copykat for output instead.

```R

setwd("/data/rmulqueen/projects/scalebio_dcis/rna")
#read in transcript gene locations from cellranger tool
annotation<-rtracklayer::readGFF("/geo_seq/code/3rd_party/10X/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf")
saveRDS(as.data.frame(annotation),file="GRCh38.genes.annot.rds")
```



Run infercnv/copykat in singularity
```bash
singularity shell \
--writable-tmpfs \
-e \
-B /data/rmulqueen/projects/scalebio_dcis/rna \
-B /home/rmulqueen/R/x86_64-pc-linux-gnu-library/4.3 \
/home/rmulqueen/singularity/infercnv-1.20.0.simg
```
```R
library(Seurat) #local
library(ggplot2) #local
library(patchwork) #local
library(ComplexHeatmap,lib.loc="/home/rmulqueen/R/x86_64-pc-linux-gnu-library/4.3/") #local
library(infercnv)

setwd("/data/rmulqueen/projects/scalebio_dcis/rna")
obj<-readRDS("tenx_dcis.pf.rds")

#create gene order annotation
gene_order<-readRDS(file="GRCh38.genes.annot.rds")
gene_order<-gene_order[gene_order$type=="gene",]
gene_order<-gene_order[!duplicated(gene_order$gene_name),]
gene_order<-gene_order[gene_order$gene_name %in% Features(obj),]
gene_order<-gene_order[gene_order$seqid %in% c(1:22,"X"),]
gene_order$seqid<-factor(gene_order$seqid,levels=c(1:22,"X")) # set chr order
gene_order<-with(gene_order, gene_order[order(seqid, start),]) #order by chr and start position
#write out gene order list
write.table(gene_order[c("gene_name","seqid","start","end")],file="gene_order.infercnv.txt",sep="\t",col.names=F,row.names=F,quote=F)
gene_order<-read.csv(file="gene_order.infercnv.txt",sep="\t",header=F,row.names=1)
#infercnv
#make directory
system(paste0("mkdir ","./infercnv/"))



infercnv_per_sample<-function(dat,sample_name="BCMDCIS80T"){
  system(paste0("mkdir -p ",paste0("./infercnv/",sample_name)))
  dat<-subset(obj,sample==sample_name) #subset data to sample specified by x and outname
  DefaultAssay(dat)<-"RNA" #using raw counts, and not SOUPX corrected counts for this

  #set all other celltypes as ref
  dat$cnv_ref<-"FALSE"
  dat@meta.data[!(dat$coarse_celltype %in% c("lumhr")),]$cnv_ref<-"TRUE" #set cnv ref by cell type


  if(sum(dat$cnv_ref=="FALSE")>50 && sum(dat$cnv_ref=="TRUE")>50){
    annot<-as.data.frame(dat@meta.data$cnv_ref)
    row.names(annot)<-row.names(dat@meta.data)

    raw_data<-t(FetchData(object = dat, vars=row.names(gene_order) ,layer = "counts"))

  infercnv_obj = CreateInfercnvObject(
    raw_counts_matrix=raw_data,
    annotations_file=annot,
    gene_order_file="gene_order.infercnv.txt",
    delim="\t",
    ref_group_names=c("TRUE"))

  infercnv_obj = infercnv::run(infercnv_obj,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir=paste0("./infercnv/",sample_name), 
    cluster_by_groups=TRUE, 
    denoise=TRUE,
    HMM=TRUE,
    HMM_report_by="cell",
    resume_mode=F,
    analysis_mode='subclusters',
    HMM_type='i6',
    num_threads=200)
  saveRDS(infercnv_obj,paste0("./infercnv/",sample_name,"/",sample_name,".inferCNV.rds"))
  }
}


lapply(unique(obj$sample), function(x) infercnv_per_sample(dat=obj,sample_name=x))


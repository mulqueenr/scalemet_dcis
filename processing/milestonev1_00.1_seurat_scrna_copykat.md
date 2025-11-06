# Run copykat 
```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
```
```R
library(copykat)
library(circlize)
library(Seurat)
library(ComplexHeatmap)
library(dendextend)
library(RColorBrewer)
setwd("/data/rmulqueen/projects/scalebio_dcis/rna")
obj<-readRDS("tenx_dcis.pf.rds")
output_dir="/data/rmulqueen/projects/scalebio_dcis/rna"



#copykat
#make directory
system(paste0("mkdir -p ",paste0(output_dir,"/copykat/")))
copykat_per_sample<-function(dat,sample_name="BCMDCIS80T"){
  system(paste0("mkdir -p ",paste0(output_dir,"/copykat/",sample_name)))
  dat<-subset(obj,sample==sample_name) 
  DefaultAssay(dat)<-"RNA" #using raw counts, and not SOUPX corrected counts for this
  norm_cells=setNames(dat@meta.data$fine_celltype,nm=row.names(dat@meta.data))
  norm_cells<-names(norm_cells[!(norm_cells %in% c("lumhr"))]) #all nonlumhr are considered normal
  exp.rawdata <- LayerData(dat, assay = "RNA", layer = "counts")
  copykat.test <- copykat(rawmat=exp.rawdata, 
                          sam.name=sample_name, 
                          distance="euclidean", 
                          norm.cell.names=norm_cells,
                          output.seg="FLASE", 
                          plot.genes="FALSE", 
                          genome="hg20",n.cores=20) #hg20 is same as hg38
  CNA.test <- data.frame(copykat.test$CNAmat)
  pred.test <- data.frame(copykat.test$prediction)
  cnv_col<-circlize::colorRamp2(colors=c("#002C3E", "#78BCC4", "#F7F8F3",  "#aa1407", "#440803"),breaks=c(-0.5,-0.25,0,0.25,0.5))

  broad_celltype_annot=setNames(dat@meta.data$broad_celltype,nm=gsub(pattern="-",replacement=".",row.names(dat@meta.data)))


  dend <- t(CNA.test[,4:ncol(CNA.test)]) %>% 
            dist(method="manhattan") %>% hclust(method="ward.D2") %>% as.dendrogram
    k_optimal=find_k(dend, krange = 2:10)

    print(paste("optimal k value for cutting hclust:", k_optimal$k))
    superclones=dendextend::cutree(dend,k=k_optimal$k+2)
    subclones=dendextend::cutree(dend,k=k_optimal$k+5)
    names(subclones)<-gsub("[.]","-",names(subclones))
    names(superclones)<-gsub("[.]","-",names(superclones))
    copykat.test$prediction$subclones<-subclones[copykat.test$prediction$cell.names]
    copykat.test$prediction$superclones<-superclones[copykat.test$prediction$cell.names]

    superclone_col=setNames(nm=unique(as.character(copykat.test$prediction$superclones)),
                            brewer.pal(length(unique(as.character(copykat.test$prediction$superclones))), "Pastel1"))
    subclone_col=setNames(nm=unique(as.character(copykat.test$prediction$subclones)),
                            brewer.pal(length(unique(as.character(copykat.test$prediction$subclones))), "Spectral"))
    
  celltype_col=c(
  'perivasc'='#c1d552',
  'vascular'='#d5a852',
  'lymphatic'='#80d552',

  'fibro'='#7f1911',
  'endo'='#f0b243',

  'tcell'='#2e3fa3',
  'bcell'='#00adea',
  'myeloid'='#00a487',
  'plasma'='#006455',

  'basal'='#7200cc',
  'lumsec'='#af00af',
  'lumhr'='#d8007c')

  row_annot = rowAnnotation(
    broad_celltype = broad_celltype_annot[names(broad_celltype_annot) %in% row.names(t(CNA.test))],
    col=list(broad_celltype=celltype_col))

  plt<-Heatmap(t(CNA.test[,4:ncol(CNA.test)]),
          cluster_columns=FALSE,
          left_annotation=row_annot,
          col=cnv_col,
          column_split=CNA.test$chrom, 
          show_row_names=FALSE,
          cluster_rows=TRUE,
          show_column_names=FALSE)

  pdf(copykat.test,paste0(output_dir,"/copykat/",sample_name,"/",sample_name,".copykat.heatmap.pdf"))
  print(plt)
  dev.off()

  saveRDS(copykat.test,paste0(output_dir,"/copykat/",sample_name,"/",sample_name,".copykat.rds"))

}

copykat_per_sample(dat=obj,sample_name='BCMDCIS05T')
copykat_per_sample(dat=obj,sample_name='BCMDCIS07T')
copykat_per_sample(dat=obj,sample_name='BCMDCIS102T-4h')
copykat_per_sample(dat=obj,sample_name='BCMDCIS124T')
copykat_per_sample(dat=obj,sample_name='BCMDCIS22T')
copykat_per_sample(dat=obj,sample_name='BCMDCIS28T')
copykat_per_sample(dat=obj,sample_name='BCMDCIS32T')
copykat_per_sample(dat=obj,sample_name='BCMDCIS35T-3h')
copykat_per_sample(dat=obj,sample_name='BCMDCIS41T')
copykat_per_sample(dat=obj,sample_name='BCMDCIS49T-24hTis')
copykat_per_sample(dat=obj,sample_name='BCMDCIS52T')
copykat_per_sample(dat=obj,sample_name='BCMDCIS65T')
copykat_per_sample(dat=obj,sample_name='BCMDCIS66T')
copykat_per_sample(dat=obj,sample_name='BCMDCIS70T')
copykat_per_sample(dat=obj,sample_name='BCMDCIS74T')
copykat_per_sample(dat=obj,sample_name='BCMDCIS80T')
copykat_per_sample(dat=obj,sample_name='BCMDCIS82T-24hTis')
copykat_per_sample(dat=obj,sample_name='BCMDCIS94T-24hTis')
copykat_per_sample(dat=obj,sample_name='BCMDCIS97T')
copykat_per_sample(dat=obj,sample_name='BCMDCIS99T')
copykat_per_sample(dat=obj,sample_name='BCMHBCA03R')
copykat_per_sample(dat=obj,sample_name='BCMHBCA04R')
copykat_per_sample(dat=obj,sample_name='BCMHBCA09R-3h')
copykat_per_sample(dat=obj,sample_name='BCMHBCA12R-3h')
copykat_per_sample(dat=obj,sample_name='BCMHBCA22R-4h')
copykat_per_sample(dat=obj,sample_name='BCMHBCA26L-24hTis-4h')
copykat_per_sample(dat=obj,sample_name='BCMHBCA29L-2h')
copykat_per_sample(dat=obj,sample_name='BCMHBCA38L-3h')
copykat_per_sample(dat=obj,sample_name='BCMHBCA85L-3h')
copykat_per_sample(dat=obj,sample_name='DCIS-79T')
copykat_per_sample(dat=obj,sample_name='DCIS-92T')
copykat_per_sample(dat=obj,sample_name='ECIS25T')
copykat_per_sample(dat=obj,sample_name='ECIS26T')
copykat_per_sample(dat=obj,sample_name='ECIS36T')
copykat_per_sample(dat=obj,sample_name='ECIS48T')
copykat_per_sample(dat=obj,sample_name='ECIS57T')
copykat_per_sample(dat=obj,sample_name='HBCA-16R')
copykat_per_sample(dat=obj,sample_name='HBCA-19T')
copykat_per_sample(dat=obj,sample_name='HBCA-83L')
copykat_per_sample(dat=obj,sample_name='HBCA17T')
copykat_per_sample(dat=obj,sample_name='IDC-79T')
```



Set up infercnv with singularity
```bash
cd /home/rmulqueen/singularity
wget https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/InferCNV/infercnv-1.20.0.simg
```

Run locally
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
/home/rmulqueen/singularity/infercnv-1.20.0.simg
```
```R
library(Seurat) #local
library(ggplot2) #local
library(patchwork) #local
library(infercnv)
library(ComplexHeatmap)

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
  system(paste0("mkdir ",paste0("./infercnv/",sample_name)))
  dat<-subset(obj,sample==sample_name) #subset data to sample specified by x and outname
  DefaultAssay(dat)<-"RNA" #using raw counts, and not SOUPX corrected counts for this

  #set all other celltypes as ref
  dat$cnv_ref<-"FALSE"
  dat@meta.data[!(dat$broad_celltype %in% c("lumhr")),]$cnv_ref<-"TRUE" #set cnv ref by cell type


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
    HMM_type='i6',
    num_threads=50)
  saveRDS(infercnv_obj,paste0("./infercnv/",sample_name,"/",sample_name,".inferCNV.rds"))
  }
}


lapply(unique(obj$sample), function(x) infercnv_per_sample(dat=obj,sample_name=x))


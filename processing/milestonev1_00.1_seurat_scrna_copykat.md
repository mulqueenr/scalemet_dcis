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
library(GenomicRanges)
setwd("/data/rmulqueen/projects/scalebio_dcis/rna")
obj<-readRDS("tenx_dcis.pf.rds")
output_dir="/data/rmulqueen/projects/scalebio_dcis/rna"

#read in cyto info
cyto=read.table(file="/data/rmulqueen/projects/scalebio_dcis/ref/cytoBand.txt",sep="\t")
colnames(cyto)<-c("chr","start","end","band","stain")
table(cyto$stain) #set colors for these

#copykat
#make directory
system(paste0("mkdir -p ",paste0(output_dir,"/copykat/")))


copykat_per_sample<-function(dat,sample_name="BCMDCIS80T"){
  system(paste0("mkdir -p ",paste0(output_dir,"/copykat/",sample_name[1])))
  dat<-subset(obj,sample %in% c(sample_name)) 
  DefaultAssay(dat)<-"RNA" #using raw counts, and not SOUPX corrected counts for this
  norm_cells=setNames(dat@meta.data$fine_celltype,nm=row.names(dat@meta.data))
  norm_cells<-names(norm_cells[!(norm_cells %in% c("lumhr"))]) #all nonlumhr are considered normal
  read_counts<-dat$nCount_RNA
  exp.rawdata <- LayerData(dat, assay = "RNA", layer = "counts")
  copykat.test <- copykat(rawmat=exp.rawdata, 
                          sam.name=sample_name[1], 
                          distance="euclidean", 
                          norm.cell.names=norm_cells,
                          output.seg="FLASE", 
                          plot.genes="FALSE", 
                          genome="hg20",n.cores=20) #hg20 is same as hg38
  CNA.test <- data.frame(copykat.test$CNAmat)
  #filter CNA.test[] to where there are predictions.

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

    #pred.test is for plotting
    pred.test <- data.frame(copykat.test$prediction)
    cells_in=row.names(t(CNA.test)[4:ncol(CNA.test),])
    pred.test <- pred.test[cells_in,]
    
    #define colors based on data
    log_col=colorRamp2(c(-0.5,-0.25,0,0.25,0.5),
                          c("darkblue","blue","white","red","darkred"))

    reads_col=colorRamp2(c(min(log10(read_counts)),
                            max(log10(read_counts))),
                            c("white","black"))

    superclone_col=setNames(nm=unique(as.character(pred.test$superclones[!is.na(pred.test$superclones)])),
                            colorRampPalette(brewer.pal(9, "Pastel1"))(length(unique(as.character(pred.test$superclones[!is.na(pred.test$superclones)])))))
    subclone_col=setNames(nm=unique(as.character(pred.test$subclones[!is.na(pred.test$subclones)])),
                            colorRampPalette(brewer.pal(9, "Spectral"))(length(unique(as.character(pred.test$subclones[!is.na(pred.test$subclones)])))))
    #names(superclone_col[is.na(names(superclone_col))])<-"Not.Assigned"
    #names(subclone_col[is.na(names(subclone_col))])<-"Not.Assigned"

    celltype_col=c(
    'perivasc'='#c1d552',
    'fibro'='#7f1911',
    'lymphatic'='#f0b243',
    'vascular'='#d0bd4a',
    'tcell'='#2e3fa3',
    'bcell'='#00adea',
    'plasma'='#006455',
    'basal'='#7200cc',
    'lumsec'='#af00af',
    'lumhr'='#d8007c')

    #plot heatmap
    ha = rowAnnotation(
        reads_annot=log10(read_counts[gsub("[.]","-",cells_in)]),
        celltype_annot=dat@meta.data[gsub("[.]","-",cells_in),]$broad_celltype,
        superclones_annot=pred.test$superclones,
        subclones_annot=pred.test$subclones,
        col=list(
            celltype_annot=celltype_col,
            reads_annot=reads_col,
            superclones_annot=superclone_col,
            subclones_annot=subclone_col
        ))

    copykat.test$CNAmat[copykat.test$CNAmat$chrom=="23",]$chrom<-"X"
    grange=makeGRangesFromDataFrame(data.frame(chr=paste0("chr",copykat.test$CNAmat$chrom),
                              start=copykat.test$CNAmat$chrompos,
                              end=copykat.test$CNAmat$chrompos))
    cyto_overlap<-GenomicRanges::findOverlaps(grange,
                                                makeGRangesFromDataFrame(cyto,keep=TRUE),
                                                select="first")
    grange$stain <- cyto[cyto_overlap,]$stain
    grange$arm <- substr(cyto[cyto_overlap,]$band,1,1)

    arm_col=c("p"="grey","q"="darkgrey")
    band_col=c("acen"="#99746F","gneg"="white","gpos100"="black","gpos25"="lightgrey","gpos50"="grey","gpos75"="darkgrey","gvar"="#446879")

    column_ha = HeatmapAnnotation(
        arm = grange$arm,
        band = grange$stain,
        col=list(arm=arm_col,band=band_col))

  plt<-Heatmap(
          t(CNA.test[,4:ncol(CNA.test)]),
          col=log_col,
          clustering_distance_rows = "manhattan",
          clustering_method_rows="ward.D2",
          cluster_columns=FALSE,
          left_annotation=ha,
          top_annotation=column_ha,
          column_split=CNA.test$chrom, 
          cluster_column_slices=FALSE,
          show_row_names=FALSE,
          cluster_rows=TRUE,
          row_split=pred.test$superclones,
          show_column_names=FALSE)

  pdf(file=paste0(output_dir,"/copykat/",sample_name[1],"/copykat.",sample_name[1],".heatmap.pdf"),width=20)
  print(plt)
  dev.off()
  saveRDS(copykat.test,paste0(output_dir,"/copykat/",sample_name[1],"/copykat.",sample_name[1],".rds"))

}

copykat_per_sample(dat=obj,sample_name=c('BCMDCIS05T'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS07T'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS102T-4h'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS124T'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS22T'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS28T'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS32T'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS35T'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS41T'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS49T'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS52T'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS65T'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS66T'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS70T'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS74T'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS80T'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS79T_24hTis_DCIS','BCMDCIS79T_24hTis_IDC'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS80T_24hTis'))

copykat_per_sample(dat=obj,sample_name=c('BCMDCIS82T-24hTis'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS92T_24hTis'))

copykat_per_sample(dat=obj,sample_name=c('BCMDCIS94T-24hTis'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS97T'))
copykat_per_sample(dat=obj,sample_name=c('BCMDCIS99T'))
copykat_per_sample(dat=obj,sample_name=c('BCMHBCA03R'))
copykat_per_sample(dat=obj,sample_name=c('BCMHBCA04R'))
copykat_per_sample(dat=obj,sample_name=c('BCMHBCA09R-3h'))
copykat_per_sample(dat=obj,sample_name=c('BCMHBCA12R-3h'))
copykat_per_sample(dat=obj,sample_name=c('BCMHBCA16R-3h'))
copykat_per_sample(dat=obj,sample_name=c('BCMHBCA17R-3h'))
copykat_per_sample(dat=obj,sample_name=c('BCMHBCA22R-4h'))
copykat_per_sample(dat=obj,sample_name=c('BCMHBCA26L-24hTis-4h'))
copykat_per_sample(dat=obj,sample_name=c('BCMHBCA29L-2h'))
copykat_per_sample(dat=obj,sample_name=c('BCMHBCA38L-3h'))
copykat_per_sample(dat=obj,sample_name=c('BCMHBCA83L-3h'))
copykat_per_sample(dat=obj,sample_name=c('BCMHBCA85L-3h'))
copykat_per_sample(dat=obj,sample_name=c('ECIS25T'))
copykat_per_sample(dat=obj,sample_name=c('ECIS26T'))
copykat_per_sample(dat=obj,sample_name=c('ECIS36T'))
copykat_per_sample(dat=obj,sample_name=c('ECIS48T'))
copykat_per_sample(dat=obj,sample_name=c('ECIS57T'))
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


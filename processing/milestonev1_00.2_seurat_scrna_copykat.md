# Run copykat 
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


copykat_per_sample<-function(obj,sample_name="BCMDCIS41T",winsize=25,cores=100){
#winsize 25 is default
  setwd(output_dir)
  system(paste0("mkdir -p ",paste0(output_dir,"/copykat/",sample_name[1])))
  setwd(paste0(output_dir,"/copykat/",sample_name[1]))
  dat<-subset(obj,sample %in% c(sample_name)) 
  DefaultAssay(dat)<-"RNA" #using raw counts, and not SOUPX corrected counts for this
  norm_cells=setNames(dat@meta.data$fine_celltype,nm=row.names(dat@meta.data))
  norm_cells<-names(norm_cells[!(norm_cells %in% c("lumhr"))]) #all nonlumhr are considered normal
  read_counts<-dat$nCount_RNA
  exp.rawdata <- LayerData(dat, assay = "RNA", layer = "counts")
  copykat.test <- copykat(rawmat=exp.rawdata, 
                          sam.name=sample_name[1], 
                          distance="euclidean", 
                          win.size=winsize,
                          norm.cell.names=norm_cells,
                          output.seg="FALSE", 
                          plot.genes="FALSE", 
                          genome="hg20",n.cores=cores) #hg20 is same as hg38
  CNA.test <- data.frame(copykat.test$CNAmat)
  #filter CNA.test[] to where there are predictions.
  dend <- t(CNA.test[,4:ncol(CNA.test)]) %>% 
  parallelDist::parDist(method = "euclidean", nbproc = cores) %>% 
    hclust(method="ward.D2") %>%  
    as.dendrogram

    k_optimal=find_k(dend, krange = 2:10)
    saveRDS(dend,file=paste0(output_dir,"/copykat/",sample_name[1],"/copykat.",sample_name[1],".",as.character(winsize),".dendrogram.rds"))

    print(paste("optimal k value for cutting hclust:", k_optimal$k))
    superclones=dendextend::cutree(dend,k=k_optimal$k+4)
    subclones=dendextend::cutree(dend,k=k_optimal$k+6)
    names(subclones)<-gsub("[.]","-",names(subclones))
    names(superclones)<-gsub("[.]","-",names(superclones))
    copykat.test$prediction$subclones<-subclones[copykat.test$prediction$cell.names]
    copykat.test$prediction$superclones<-superclones[copykat.test$prediction$cell.names]

    #pred.test is for plotting
    pred.test <- data.frame(copykat.test$prediction)
    cells_in=row.names(t(CNA.test)[4:ncol(CNA.test),])
    pred.test <- pred.test[cells_in,]
    

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
    'endo'='#f0b243',#'lymphatic'='#f0b243',
    'vascular'='#d0bd4a',
    'tcell'='#2e3fa3',
    'bcell'='#00adea',
    'myeloid'="#55157a",
    'plasma'='#006455',
    'basal'='#7200cc',
    'lumsec'='#af00af',
    'lumhr'='#d8007c')

    #plot heatmap
    ha = rowAnnotation(
        reads_annot=log10(read_counts[gsub("[.]","-",cells_in)]),
        celltype_annot=dat@meta.data[gsub("[.]","-",cells_in),]$coarse_celltype,
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

   #define colors based on data
    #log_col=colorRamp2(c(stats::quantile(unlist(t(CNA.test[,4:ncol(CNA.test)])),probs=c(0.01,0.1,0.25,0.5,0.75,0.9,0.99))),
    #                      c("#053061","#d9d9ff","white","white","white","#ffd9d9","#67001f")) 

    #define colors based on data
    #log_col=colorRamp2(c(stats::quantile(unlist(t(CNA.test[,4:ncol(CNA.test)])),probs=c(0.01,0.25,0.5,0.75,0.9))),
    #                      c("darkblue","blue","white","red","darkred")) 
    #copying color breaks from copykat
    #my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
    #copying how copykat does its color scaling
    my_palette <- colorRampPalette(col=c("darkblue","lightblue","white","pink","darkred"))(n = 999)

    col_breaks = unlist(c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=49)))
    log_col=colorRamp2(col=my_palette,breaks=col_breaks)

  plt<-Heatmap(
          t(CNA.test[,4:ncol(CNA.test)]),
          col=log_col,
          clustering_distance_rows = "manhattan",
          clustering_method_rows="ward.D2",
          cluster_columns=FALSE,
          left_annotation=ha,
          cluster_rows = dend,
          top_annotation=column_ha,
          column_split=CNA.test$chrom, 
          cluster_column_slices=FALSE,
          show_row_names=FALSE,
          #row_split=pred.test$superclones,
          show_column_names=FALSE)

  pdf(file=paste0(output_dir,"/copykat/",sample_name[1],"/copykat.",sample_name[1],".",as.character(winsize),".heatmap.pdf"),width=20)
  print(plt)
  dev.off()

  saveRDS(copykat.test,paste0(output_dir,"/copykat/",sample_name[1],"/copykat.",sample_name[1],".",as.character(winsize),".rds"))

}

winsize=25
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS05T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS07T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS102T-4h'),winsize=winsize) #need to correct name
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS124T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS22T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS28T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS32T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS35T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS41T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS49T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS52T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS65T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS66T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS70T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS74T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS80T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS79T_24hTis_DCIS','BCMDCIS79T_24hTis_IDC'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS80T_24hTis'),winsize=winsize)

opykat_per_sample(obj=obj,sample_name=c('BCMDCIS82T-24hTis'),winsize=winsize) #BCMDCIS82T-24hTis correct name
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS92T_24hTis'),winsize=winsize)

copykat_per_sample(obj=obj,sample_name=c('BCMDCIS94T-24hTis'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS97T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMDCIS99T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMHBCA03R'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMHBCA04R'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMHBCA09R-3h'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMHBCA12R-3h'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMHBCA16R-3h'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMHBCA17R-3h'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMHBCA22R-4h'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMHBCA26L-24hTis-4h'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMHBCA29L-2h'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMHBCA38L-3h'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMHBCA83L-3h'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('BCMHBCA85L-3h'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('ECIS25T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('ECIS26T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('ECIS36T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('ECIS48T'),winsize=winsize)
copykat_per_sample(obj=obj,sample_name=c('ECIS57T'),winsize=winsize)

```

Assign final clone names (and align to copykit output clone names).

Visually pairing RNA and methylation clones as first pass.

```R

output_dir="/data/rmulqueen/projects/scalebio_dcis/rna"
copykat_output<-list.files(path=paste0(output_dir,"/copykat/"),recursive=TRUE,full.names=TRUE,pattern=".100.rds")

#make clones a named list to collapse overclustering or low cell counts/cluster
winsize=100
sample_name<-"ECIS26T"
dat<-readRDS(file=paste0(output_dir,"/copykat/",sample_name[1],"/copykat.",sample_name[1],".",as.character(winsize),".rds"))
dend<-readRDS(file=paste0(output_dir,"/copykat/",sample_name[1],"/copykat.",sample_name[1],".",as.character(winsize),".dendrogram.rds"))

dend<-dat$hclustering


assign_copykit_aneuploid_clonename<-function(sample_name,cancer_clones,split_on="superclones",winsize=150){
    copykat.test<-readRDS(paste0(output_dir,"/copykat/",sample_name[1],"/copykat.",sample_name[1],".",as.character(winsize),".rds"))
    dend<-readRDS(file=paste0(output_dir,"/copykat/",sample_name[1],"/copykat.",sample_name[1],".",as.character(winsize),".dendrogram.rds"))
    CNA.test <- data.frame(copykat.test$CNAmat)

    copykat.test$prediction$subclones<-subclones[copykat.test$prediction$cell.names]
    copykat.test$prediction$superclones<-superclones[copykat.test$prediction$cell.names]

    #pred.test is for plotting
    pred.test <- data.frame(copykat.test$prediction)
    cells_in=row.names(t(CNA.test)[4:ncol(CNA.test),])
    pred.test <- pred.test[cells_in,]
    
    #define colors based on data
    my_palette <- colorRampPalette(col=c("darkblue","lightblue","white","pink","darkred"))(n = 999)

    col_breaks = unlist(c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=49)))
    log_col=colorRamp2(col=my_palette,breaks=col_breaks)
    reads_col=colorRamp2(c(min(log10(read_counts)),
                            max(log10(read_counts))),
                            c("white","black"))

    superclone_col=setNames(nm=unique(as.character(pred.test$superclones[!is.na(pred.test$superclones)])),
                            colorRampPalette(brewer.pal(9, "Pastel1"))(length(unique(as.character(pred.test$superclones[!is.na(pred.test$superclones)])))))
    subclone_col=setNames(nm=unique(as.character(pred.test$subclones[!is.na(pred.test$subclones)])),
                            colorRampPalette(brewer.pal(9, "Spectral"))(length(unique(as.character(pred.test$subclones[!is.na(pred.test$subclones)])))))

    celltype_col=c(
    'perivasc'='#c1d552',
    'fibro'='#7f1911',
    'endo'='#f0b243',#'lymphatic'='#f0b243',
    'vascular'='#d0bd4a',
    'tcell'='#2e3fa3',
    'bcell'='#00adea',
    'myeloid'="#55157a",
    'plasma'='#006455',
    'basal'='#7200cc',
    'lumsec'='#af00af',
    'lumhr'='#d8007c')

    #plot heatmap
    ha = rowAnnotation(
        reads_annot=log10(read_counts[gsub("[.]","-",cells_in)]),
        celltype_annot=dat@meta.data[gsub("[.]","-",cells_in),]$coarse_celltype,
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

  pdf(file=paste0(output_dir,"/copykat/",sample_name[1],"/copykat.",sample_name[1],".",as.character(winsize),".heatmap.pdf"),width=20)
  print(plt)
  dev.off()

    #plot heatmap
    ha = rowAnnotation(
        reads=log10(tmp@colData$unique_reads),
        cg_perc=tmp@colData$mcg_pct,
        celltype=tmp@colData$fine_celltype,
        superclones=as.character(tmp@colData$superclones),
        subclones=as.character(tmp@colData$subclones),
        cancerclone=as.character(tmp@colData$clonename),
        col= list(
            celltype=celltype_col,
            reads=reads_col,
            cg_perc=cg_perc_col,
            superclones=superclone_col,
            subclones=subclone_col,
            cancerclone_col
        ))

    cyto_overlap<-GenomicRanges::findOverlaps(tmp@rowRanges,
                                                makeGRangesFromDataFrame(cyto,keep=TRUE),
                                                select="first")
    tmp@rowRanges$stain <- cyto[cyto_overlap,]$stain

    arm_col=c("p"="grey","q"="darkgrey")
    band_col=c("acen"="#99746F","gneg"="white","gpos100"="black","gpos25"="lightgrey","gpos50"="grey","gpos75"="darkgrey","gvar"="#446879")
    dip_cov=colorRamp2(c(0.5,1,1.5), 
                            c("white","grey","black"))
    column_ha = HeatmapAnnotation(
        mappability=tmp@rowRanges$diploid_cov,
        arm = tmp@rowRanges$arm,
        band = tmp@rowRanges$stain,
        col=list(mappability=dip_cov,arm=arm_col,band=band_col))

    plt<-Heatmap(
        t(tmp@assays@data$logr),
        clustering_distance_rows = "manhattan", #CHANGING THIS
        left_annotation=ha,col=log_col,
        row_split=as.character(tmp@colData$clonename),
        show_column_names=FALSE,show_row_names=FALSE,
        top_annotation=column_ha,cluster_columns=FALSE,cluster_column_slices=FALSE,column_split=seqnames(tmp@rowRanges),
        name="logr")

    pdf(paste0("/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/copykit/",sample_name,"/copykit.",sample_name,".",resolution,".cancerclone.pdf"),width=20)
    print(plt)
    dev.off()

    saveRDS(tmp,file=paste0("/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/copykit/",sample_name,"/copykit.",sample_name,".",resolution,".rds"))
}


assign_copykit_aneuploid_clonename(sample_name="BCMDCIS05T",cancer_clones=c("c1"='2',"c2"='3',"c2"='4')) #done
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS07T',cancer_clones=c()) #done all diploid #x loss?
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS102T_24hTis',cancer_clones=c("c1"='3','c2'='4')) #done some x loss?
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS124T',cancer_clones=c("c1="='2',"c2"='4')) 
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS22T',cancer_clones=c("c1"='2','c2'='4')) 
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS28T',,cancer_clones=c("c1"='5',"c2"='3','c2'='6','c2'='7')) 
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS32T',cancer_clones=c()) #all diploid
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS35T',cancer_clones=c("c1"='1',"c1"='2',"c1"='3',"c1"='5',"c1"='6',"c1"='7',"c1"='8'))
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS41T',split_on='subclones',cancer_clones=c('c1'='4','c1'='6',
                                                                                        'c2'='1',
                                                                                        'c3'='3','c3'='7',
                                                                                        'c4'='2')) 
#assign_copykit_aneuploid_clonename(sample_name='BCMDCIS49T',cancer_clones=c()) #not run
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS52T',cancer_clones=c('c1'='3','c2'='2','c2'='4')) 
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS65T',cancer_clones=c('c1'='5')) 
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS66T',cancer_clones=c('c1'='2',
                                                                            'c2'='4','c2'='11',
                                                                            'c3'='14','c3'='12','c3'='7','c3'='15','c3'='13','c3'='8','c3'='5'))                                   
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS70T',cancer_clones=c('c1'='2','c1'='3','c1'='4')) 
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS74T',cancer_clones=c('c1'='11','c1'='8','c1'='9',
                                                                            'c2'='6',
                                                                            'c3'='10',
                                                                            'c4'='3',
                                                                            'c5'='2',
                                                                            'c6'='4',
                                                                            'c8'='12',
                                                                            'c9'='7')) #lotsa clear evolution in this one
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS79T_24hTis_DCIS',cancer_clones=c('c1'='2','c2'='4')) 
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS80T_24hTis',cancer_clones=c('c1'='3',
                                                                                    'c1'='5','c1'='4')) 
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS82T_24hTis',cancer_clones=c('c1'='3','c1'='4','c1'='2')) 
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS92T_24hTis',cancer_clones=c('c1'='3','c1'='4')) 
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS94T_24hTis',split_on='subclones',cancer_clones=c('c1'='5','c1'='4',
                                                                                                        'c2'='3')) 
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS97T',cancer_clones=c('c1'='11','c1'='5','c1'='10','c1'='12','c1'='17','c1'='7',
                                                                            'c2'='6','c2'='2','c2'='8',
                                                                            'c4'='15',
                                                                            'c3'='9')) 
assign_copykit_aneuploid_clonename(sample_name='BCMDCIS99T',cancer_clones=c('c1'='4')) 
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA03R',cancer_clones=c('c1'='6')) 
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA04R',split_on='subclones',cancer_clones=c('c1'='5','c2'='6','c2'='7')) 
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA09R-3h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA12R-3h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA16R-3h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA17R-3h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA19R-4h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA22R-4h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA26L-24hTis-4h',cancer_clones=c()) #i dont see any cancer clones previously split_on="subclones",cancer_clones=c('c1'='6')
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA29L-2h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA38L-3h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA83L-3h',cancer_clones=c('c1'='4')) 
assign_copykit_aneuploid_clonename(sample_name='BCMHBCA85L-3h',cancer_clones=c()) 
assign_copykit_aneuploid_clonename(sample_name='ECIS25T',split_on="subclones",cancer_clones=c('c1'='7','c1'='4','c1'='6',
                                                                                            'c2'='5','c2'='1','c2'='8',
                                                                                            'c3'='9','c3'='2')) 
assign_copykit_aneuploid_clonename(sample_name='ECIS26T',cancer_clones=c('c1'='3','c1'='17','c1'='15','c1'='9','c1'='10','c1'='5','c1'='11',
                                                                        'c1'='13','c1'='4','c1'='1','c1'='6','c1'='16','c1'='7','c1'='8',
                                                                        'c2'='14')) 

assign_copykit_aneuploid_clonename(sample_name='ECIS36T',,cancer_clones=c('c1'='1','c2'='4','c3'='2')) 
assign_copykit_aneuploid_clonename(sample_name='ECIS48T',cancer_clones=c('c1'='5')) #might have a cancer precursor in the diploid pop chr16 loss in some lumhr, 1q gain in lumsec?
assign_copykit_aneuploid_clonename(sample_name='ECIS57T',cancer_clones=c('c1'='2','c1'='3','c1'='4')) 
```
count telomere motifs using fuzzy matching on fastq per cell
perform normalization with nico
get telomere estimate sizes
```bash
singularity shell \
    --bind /data/rmulqueen/projects/scalebio_dcis \
    --bind ~/tools/ \
~/singularity/amethyst.sif

cd /data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/fastq
for plate in $(find -maxdepth 1 -type d -name "*plate*" | cut -c 3-);
do python /data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/telomere_fq_extract.py --plate $plate --cores 4 & done &

#NOTE in plates that were run multiple times, use ScaleMethylMerged fastqs, remove others.
#Just go through directories and check by hand for now.
cd /data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/telomere_fq

for plate in $(find -maxdepth 1 -type d -name "*plate*" | cut -c 3-); do
    if [ -f "./${plate}/ScaleMethylMerged_A01_S1_L001_R1_001.fastq.gz" ]; then
        cat ./${plate}/ScaleMethylMerged*R1*fastq.gz > ./${plate}/${plate}.merged.telo.R1.fq.gz
        cat ./${plate}/ScaleMethylMerged*R2*fastq.gz > ./${plate}/${plate}.merged.telo.R2.fq.gz
        cat ./${plate}/ScaleMethylMerged*readCounts.tsv > ./${plate}/${plate}_readCounts.telo.tsv
    else
        cat ./${plate}/ScaleMethyl*R1*fastq.gz > ./${plate}/${plate}.merged.telo.R1.fq.gz
        cat ./${plate}/ScaleMethyl*R2*fastq.gz > ./${plate}/${plate}.merged.telo.R2.fq.gz
        cat ./${plate}/ScaleMethyl*readCounts.tsv > ./${plate}/${plate}_readCounts.telo.tsv
    fi
done

#concatenate all plate-lane counts
cat /data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/telomere_fq/*/*_readCounts.telo.tsv > /data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/telomere_fq/telomere_readCounts.tsv

```

```R
library(dplyr)
library(ggplot2)
telo<-read.csv("/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/telomere_fq/telomere_readCounts.tsv",
                header=T,
                sep=",")
#clean up data and combine duplicate cells
telo<-telo[telo$cellID!="cellID",] #remove duplicate header files
telo$totalReads<-as.numeric(telo$totalReads)
telo$telomereReads<-as.numeric(telo$telomereReads)
telo<- telo %>% group_by(cellID) %>% summarise(across(c(totalReads, telomereReads), sum))
telo<-as.data.frame(telo %>% mutate(percent_telo=telomereReads/totalReads))
row.names(telo)<-telo$cellID

cellline_meta<-as.data.frame(read.csv(file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/merged_data/01_celllines.metadata.tsv",sep="\t"))
cellline_meta$fine_celltype<-cellline_meta$sample
cellline_meta$Group<-"cellline"
cellline_meta$cell_id<-row.names(cellline_meta)
dcis_meta<-as.data.frame(read.csv(file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/merged_data/05_scaledcis.fine_celltype.metadata.tsv",sep="\t"))
dcis_meta$cell_id<-row.names(dcis_meta)
library(plyr)
telomere_meta<-rbind.fill(cellline_meta, dcis_meta)
row.names(telomere_meta)<-telomere_meta$cell_id
telomere_meta$totalReads_fq<-telo[row.names(telomere_meta),]$totalReads
telomere_meta$telomereReads_fq<-telo[row.names(telomere_meta),]$telomereReads
telomere_meta$percent_telo<-telo[row.names(telomere_meta),]$percent_telo

write.table(telomere_meta,col.names=T,row.names=T,
            file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/telomere_fq/telomere_readCounts.summarized.tsv")
telomere_meta<-read.table(file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/telomere_fq/telomere_readCounts.summarized.tsv")

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
'lumhr'='#d8007c',
"MCF10A"="red",
"MCF7"="red",
"MDA-MB-231"="red")

group_col=c("HBCA"="#c6d8e8",
            "DCIS"="#8c86bc",
            "DCIS/LCIS"="#bfda9f",
            "IDC"="#e37f76")

telomere_meta$Group<-factor(telomere_meta$Group,levels=names(group_col))
telomere_meta$fine_celltype<-factor(telomere_meta$fine_celltype,levels=names(celltype_col))

plt<-ggplot(telomere_meta,aes(x=fine_celltype,y=percent_telo*100,color=fine_celltype,fill=Group))+
geom_boxplot(outlier.shape=NA)+scale_color_manual(values=celltype_col)+scale_fill_manual(values=group_col)+
theme_minimal()+coord_cartesian(ylim = c(0, 5))
ggsave(plt,width=10,file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/telomere_fq/telomere_percCount.summarized.pdf")

#plot celltypes by age 

plt<-ggplot(telomere_meta[!is.na(telomere_meta$Age) & telomere_meta$Group!="DCIS/LCIS",],aes(x=Age,y=percent_telo*100,color=fine_celltype,fill=Group))+
geom_point(size=0.5)+scale_color_manual(values=celltype_col)+scale_fill_manual(values=group_col)+
theme_minimal()+coord_cartesian(ylim = c(0, 5))+facet_wrap(fine_celltype~Group,ncol=3)+geom_smooth(method='lm', formula= y~x)

ggsave(plt,file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/telomere_fq/telomere_percCount.age.pdf",height=20,width=10)
```


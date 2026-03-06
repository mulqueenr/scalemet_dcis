transposable element methylation analysis

copy these to ref folder
- check cnv overlap with haploinsuffiency and triplosensitivity markers
- check methylation scores for repeats (expecting lower methylation in cancer for repeats)

## cgi
https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=regulation&hgta_track=cpgIslandExt&hgta_table=cpgIslandExt&hgta_doSchema=describe+table+schema
```bash
cd /data/rmulqueen/projects/scalebio_dcis/ref
wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz"
/data/rmulqueen/projects/scalebio_dcis/ref/cpgIslandExt.bed
```

## haploinsufficiency and triplosensitivity tracks
https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=3596528659_3AWPx4Dyqv9CtPn24flp4abaaA3K&db=hg38&c=chr6&g=dosageSensitivity
```bash
cd /data/rmulqueen/projects/scalebio_dcis/ref
wget "https://zenodo.org/records/6347673/files/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz?download=1"
mv Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz?download=1 Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz
/data/rmulqueen/projects/scalebio_dcis/ref/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz
#pHaplo scores ≥0.86 indicate that the average effect sizes of deletions are as strong as the loss-of-function of genes known to be constrained against protein truncating variants (average OR≥2.7) (Karczewski et al., 2020). pHaplo scores ≥0.55 indicate an odds ratio ≥2.

#pTriplo scores ≥0.94 indicate that the average effect sizes of deletions are as strong as the loss-of-function of genes known to be constrained against protein truncating variants (average OR≥2.7) (Karczewski et al., 2020). pHaplo scores ≥0.68 indicate an odds ratio ≥2.
#gene symbols, pHaplo, and pTriplo scores
```

## repeats, lines sines ltr etc
https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_doSchema=describe+table+schema

```bash
cd /data/rmulqueen/projects/scalebio_dcis/ref
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
/data/rmulqueen/projects/scalebio_dcis/ref/rmsk.txt.gz
```

# Reading in amethyst object and summarizing methylation over different annotations.

```R
library(copykit)
library(GenomicRanges)
library(amethyst)
library(rtracklayer)
library(ComplexHeatmap)
library(circlize)
#set environment and read in data
set.seed(111)
options(future.globals.maxSize= 200000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)

repeat_dir=paste(sep="/",project_data_directory,"repeat_analysis")
system(paste("mkdir -p",repeat_dir))

setwd(wd)
obj<-readRDS(file="06_scaledcis.celltype.amethyst.rds")

repeats<-read.table("/data/rmulqueen/projects/scalebio_dcis/ref/rmsk.txt.gz")

repeats<-repeats[,c("V6","V7","V8","V10","V11","V12")]
colnames(repeats)<-c("chr","start","end","strand","name","type")

repeats<-GRanges(repeats) %>% split(~type)

#make windows over copykit ranges
#actually i think i can just slice the smoothed 500bp windows I already have, per clone or per celltype
#get 500kb windows ranges
clone500bpwindows<-readRDS(file=paste0(project_data_directory,"/DMR_analysis/","cnv_clones_alllumhr/","dmr_analysis.cnv_clones_alllumhr.500bp_windows.rds"))
met <- GRanges(clone500bpwindows[["pct_matrix"]])


#calculate global methylation per clone
clone_met<-colMeans(as.data.frame(mcols(met)),na.rm=T)

#overlap with repeats
#mean percentage per track of overlap per type
#overlap with haplo and triplo tracks

#get percent methylation, per repeat region, per clone

repeat_methylation_perc<-lapply(1:length(repeats),function(x){
    overlaps<-findOverlaps(met,repeats[[x]])
    return(colMeans(as.data.frame(mcols(met[overlaps@from,])),na.rm=T))
})

repeat_methylation_perc<-do.call("cbind",repeat_methylation_perc)
colnames(repeat_methylation_perc)<-names(repeats)
repeat_methylation_perc<-as.data.frame(repeat_methylation_perc)

#add group metadata per sample
repeat_methylation_perc$Sample<-unlist(gsub(gsub(row.names(repeat_methylation_perc),pattern="_diploid",replacement=""),pattern="_c[0-9]",replacement=""))
repeat_methylation_perc$Sample<-gsub(repeat_methylation_perc$Sample,pattern="[.]",replacement="-")
group_meta<-unique(data.frame(Sample=obj@metadata$Sample,Group=obj@metadata$Group))
row.names(group_meta)<-group_meta$Sample
repeat_methylation_perc$Group<-group_meta[repeat_methylation_perc$Sample,]$Group

repeat_methylation_perc$Group<-factor(repeat_methylation_perc$Group,levels=c("NA","HBCA","DCIS","Synchronous","IDC"))

#https://emilhvitfeldt.github.io/r-color-palettes/discrete/NatParksPalettes/Acadia/
met_col=colorRamp2(breaks=seq(from = 50, to = 100, length.out = 9),colors=rev(c("#212E52FF", "#444E7EFF", "#8087AAFF", "#B7ABBCFF", "#F9ECE8FF", "#FCC893FF", "#FEB424FF", "#FD8700FF", "#D8511DFF")))
column_ha = HeatmapAnnotation(global_met = clone_met,col=list(global_met=met_col))


pdf(paste0(repeat_dir,"/","repeat_per_clone_methylation.pdf"),width=20)
Heatmap(t(repeat_methylation_perc[,c("LINE","SINE","Retroposon")]),
        column_split=repeat_methylation_perc$Group,
         bottom_annotation = column_ha, col=met_col,
        cluster_column_slices=FALSE)

dev.off()


#next up is to segment cnv elements so i can overlap with bed discretely for enrichment
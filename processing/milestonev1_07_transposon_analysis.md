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
library(GenomicRanges)
library(amethyst)
library(rtracklayer)
library(rhdf5)
library(parallel)
library(ComplexHeatmap)
library(circlize)
library(tidyr)

#set environment and read in data
set.seed(111)
options(future.globals.maxSize= 200000*1024^2) #80gb limit for parallelizing

project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"

#read in object from directory
task_cpus=300
processing_folder="04_repeat_analysis"
wd=paste(sep="/",project_data_directory,processing_folder)
system(paste0("mkdir -p ",wd))
setwd(wd)


obj<-readRDS(file=paste(project_data_directory,"03_fine_celltyping","03_scaledcis.final_celltypes.amethyst.rds",sep="/"))

repeats<-read.table("/data/rmulqueen/projects/scalebio_dcis/ref/rmsk.txt.gz")
repeats<-repeats[,c("V6","V7","V8","V10","V11","V12","V13")]
colnames(repeats)<-c("chr","start","end","strand","name","class","family")
repeats<-repeats[complete.cases(repeats),]
repeats$class_family<-paste(repeats$class,repeats$family)
repeats_fam<-GRanges(repeats) %>% split(~class_family,drop=TRUE)
length(repeats_fam)
#72 cases

#loop through all repeat element types (72 class_family) and summarize methylation over annotations
for(annot_index in 1:length(repeats_fam)){
    annot_name=mcols(repeats_fam[[annot_index]])$class_family[1]
    annot_name=gsub(annot_name,pattern=" ",repl="_")
    annot<-repeats_fam[[annot_index]]
    out<-mclapply(row.names(obj@metadata),function(cellid){
        print(paste("Calculating repeat type",annot_name,"for cellid",cellid))
        h5_path<-obj@metadata[cellid,]$h5_path
        cell_dat <- do.call("rbind",h5read(h5_path ,"/CG",cell_name))
        cell_dat$start<-cell_dat$pos
        cell_dat$end<-cell_dat$pos+1
        cell_dat<-GRanges(cell_dat)
        overlaps<-suppressWarnings(GenomicRanges::findOverlaps(query=annot,subject=cell_dat))
        out<-c(cellid=cellid,annotation=annot_name,colSums(as.matrix(mcols(cell_dat[overlaps@to,])[c("t","c")])))
        return(out)},
        mc.cores=200)
    out<-as.data.frame(do.call("rbind",out))
    saveRDS(out,
        file=paste0("04_repeat_analysis.",annot_name,".repeat_windows.rds"))
}

#after this runs, read in all rds files and summarize methylation percentage over repeat annotations for all cells
repeat_files<-list.files(pattern=".*.repeat_windows.rds")
repeat_elements<-do.call("rbind",lapply(repeat_files,function(x){
    in_dat<-readRDS(x)
}))

repeat_elements$c <- as.numeric(repeat_elements$c)
repeat_elements$t <- as.numeric(repeat_elements$t)

repeat_elements<-repeat_elements[complete.cases(repeat_elements),]
repeat_elements<-repeat_elements[which((repeat_elements$t+repeat_elements$c)>100),] #require 100 measurements
repeat_elements$methylation<-repeat_elements$c/(repeat_elements$c+repeat_elements$t)
repeat_mat <- repeat_elements %>% 
            pivot_wider(id_cols=cellid,names_from=annotation, values_from=methylation,values_fill=NA)

saveRDS(as.data.frame(repeat_mat),file="04_repetitive_elements_methylation.rds")

repeat_mat<-as.data.frame(repeat_mat)
row.names(repeat_mat)<-repeat_mat$cellid

library(dplyr)
retro_dat<-left_join(obj@metadata,repeat_mat)


repeat_elements_by_celltype_group<-retro_dat %>% 
    group_by(celltype_group) %>% 
    mutate(across(DNA_DNA:Unknown,as.numeric)) %>% 
    summarize(across(DNA_DNA:Unknown,mean,na.rm=T))

library(ComplexHeatmap)

repeat_elements_by_celltype_group <- as.data.frame(repeat_elements_by_celltype_group)
row.names(repeat_elements_by_celltype_group) <- repeat_elements_by_celltype_group$celltype_group
repeat_elements_by_celltype_group<-repeat_elements_by_celltype_group[,2:ncol(repeat_elements_by_celltype_group)]

repeat_elements_by_celltype_group<-repeat_elements_by_celltype_group %>% mutate(across(everything(), as.numeric))


#set colors
celltype_col=c(
    "pericyte"="#FF6600",
    "fibroblast"="#FF0066",
    "endothelial"="#FFCC00",
    "unknown"="#666666",

    "myeloid"="#00FFFF",
    "bcell"="#0099FF",
    "tcell"="#0033FF",

    "basal"="#6600FF",
    "lumsec"="#CC00FF",
    "lumhr"="#FF00CC",
    "cancer"="#00FF99")

cellgroups<-expand.grid(c("HBCA","DCIS","Synchronous","IDC"),names(celltype_col))
cellgroups<-paste(cellgroups$Var2,cellgroups$Var1,sep="_")
cellgroups<-cellgroups[which(cellgroups %in% row.names(repeat_elements_by_celltype_group))]
repeat_elements_by_celltype_group<-t(repeat_elements_by_celltype_group)
repeat_elements_by_celltype_group<-repeat_elements_by_celltype_group[,cellgroups]

plt<-Heatmap(t(scale(t(repeat_elements_by_celltype_group))),
    column_order=1:ncol(repeat_elements_by_celltype_group),
    cluster_column_slices=FALSE,cluster_columns=FALSE,
    column_split=factor(unlist(lapply(strsplit(cellgroups,"_"),"[",1)),levels=names(celltype_col)))

pdf("04_repeat_methylation.scaled.pdf",width=10,height=20)
print(plt)
dev.off()


library(circlize)
met_col=colorRamp2(c(min(unlist(repeat_elements_by_celltype_group),na.rm=T),
                    median(unlist(repeat_elements_by_celltype_group),na.rm=T),
                    max(unlist(repeat_elements_by_celltype_group),na.rm=T)),
    c("#ff70ff","#CCCCCC","#000000"))

plt<-Heatmap(repeat_elements_by_celltype_group,
    column_order=1:ncol(repeat_elements_by_celltype_group),
    col=met_col,
    cluster_column_slices=FALSE,cluster_columns=FALSE,
    column_split=factor(unlist(lapply(strsplit(cellgroups,"_"),"[",1)),levels=names(celltype_col)))

pdf("04_repeat_methylation.pdf",width=10,height=20)
print(plt)
dev.off()

#run some statistical tests
repeat_mat<-readRDS(file="04_repetitive_elements_methylation.rds")
retro_dat<-left_join(obj@metadata,repeat_mat,by="cellid")

repeat_elements_by_celltype_group <-retro_dat %>% 
    filter(celltype_group %in% c("lumhr_HBCA","cancer_DCIS","cancer_Synchronous","cancer_IDC")) %>% 
    mutate(across(DNA_DNA:Unknown,as.numeric))
repeat_elements_by_celltype_group$celltype_group<-factor(repeat_elements_by_celltype_group$celltype_group,
levels=c("lumhr_HBCA","cancer_DCIS","cancer_Synchronous","cancer_IDC"))

plt_list<-lapply(c("LINE","LINE_L1","Retroposon","SINE","Satellite_telo"),function(x){
    plt<-ggplot(repeat_elements_by_celltype_group,
                aes(x=celltype_group,y=repeat_elements_by_celltype_group[,x],color=Group))+
        geom_jitter(size=0.5)+
        geom_violin(fill=NA,color="black")+ylim(c(0,1))+ylab(x)+
        scale_color_manual(values=group_col)
    return(plt)
})

library(patchwork)
library(presto)
library(Matrix)
ggsave(wrap_plots(plt_list,ncol=1),file="04_repeat_methylation.multiple_elements.lumhrHBCA_v_cancerIDC.pdf",height=10)

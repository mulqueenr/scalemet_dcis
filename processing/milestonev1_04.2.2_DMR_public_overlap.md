
Overlapping public methylarray data on DCIS progression

Using illumina CG IDs to compare location to our DMRs

```R

set.seed(111)
options(future.globals.maxSize= 500000*1024^2) #80gb limit for parallelizing
library(amethyst)
library(data.table)
library(dplyr)
library(msigdbr)
library(fgsea)
library(ggplot2)
library(GenomicRanges)
library(Matrix)
library(parallel)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(GeneOverlap)
library(matrixStats)
library(sparseMatrixStats)
library(irlba)
library(rtracklayer)

project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"

#read in object from directory
task_cpus=300
processing_folder="04_dmr"
wd=paste(sep="/",project_data_directory,processing_folder)
system(paste0("mkdir -p ",wd))
setwd(wd)

```
Functions
```R


```

## Running on 500bp windows
## On celltype x group DMRs
```R

output_directory=paste0(wd,"/","dmr_celltype_group")
system(paste0("mkdir -p ",output_directory))

obj<-readRDS(file=paste(project_data_directory,"03_fine_celltyping","03_scaledcis.final_celltypes.amethyst.rds",sep="/"))

comparison_set<-"dmr_within_celltype_across_group"
input_dmr_folder<-paste(project_data_directory,"04_dmr",comparison_set,"dmr_out",sep="/")

output_plots_folder<-paste0(output_folder,"/plots")

group_order<-list()
for(i in names(celltype_col)){
  for(j in c("HBCA","DCIS","Synchronous","IDC")){
    group_order<-c(group_order,paste(i,j,"v","rest",sep="_"))
    for(k in c("HBCA","DCIS","Synchronous","IDC"))
      if(j!=k){
        group_order<-c(group_order,paste(i,j,"v",i,k,sep="_"))}
  }
}
group_order<-unlist(group_order)

#read in dmrs and set order for plotting
dmr=do.call("rbind",lapply(list.files(input_dmr_folder,full.names=T,pattern=".collapse.rds"),function(i){readRDS(i)}))
dmr$comparison_name=dmr$type


#read in public data
#from https://link.springer.com/article/10.1186/s13148-015-0094-0#Sec20
cgs_corr_to_idc_development<-readxl::read_excel("/data/rmulqueen/projects/scalebio_dcis/ref/13148_2015_94_MOESM4_ESM.xls",sheet=1,col_names=T,skip=1)
cgs_corr_to_idc_development$cg_name<-cgs_corr_to_idc_development$`Illumina CG ID`

#read in illumina array conversion
#https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=varRep&hgta_track=genotypeArrays&hgta_table=snpArrayIllumina450k&hgta_doSchema=describe+table+schema

illumina_450_array_data <- import("https://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/illumina/illumina450K.bb", format = "bigBed") %>% 
                            as.data.frame()
illumina_450_array_data$cg_name<-illumina_450_array_data$name


cgs_cor_to_idc_development<-left_join(cgs_corr_to_idc_development,illumina_450_array_data,by="cg_name")


cgs_cor_to_idc_development<-cgs_cor_to_idc_development %>% 
                select("cg_name","seqnames","start","end",
                "Gene", "Hazard Ratio", "P-value", 
                "Diff Meth Normal and DCIS") %>% as.data.frame()

johnson_etal<-cgs_cor_to_idc_development %>% filter(!is.na(cgs_cor_to_idc_development$seqnames)) %>% GRanges() #losing a lotta sites

dmr<-GRanges(dmr)
library(GenomicRanges)

overlaps_johnson_etal <- findOverlaps(query = dmr, subject = johnson_etal)
#total count of cgs that overlap at least one dmr site

length(unique(overlaps_johnson_etal@to))
#32
table(dmr[overlaps_johnson_etal@from,]$comparison_name,dmr[overlaps_johnson_etal@from,]$direction)

#interesting that its so low of an overlap

#fleischer et al
#https://pmc.ncbi.nlm.nih.gov/articles/PMC4514996/table/Tab2/


#read fleischer et al
fleischer<-read.csv("/data/rmulqueen/projects/scalebio_dcis/ref/PMC4514996.table2.tsv",skip=0,header=T,sep="\t")
fleischer$cg_name<-fleischer$`Illumina.cg.ID`

fleischer<-left_join(fleischer,illumina_450_array_data,by="cg_name")

fleischer<-fleischer %>% filter(!is.na(fleischer$seqnames)) %>% GRanges() #losing a lotta sites

overlaps_johnson_etal <- findOverlaps(query = dmr, subject = fleischer)

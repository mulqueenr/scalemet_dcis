
# Taking per sample copykit and copykat output to coembed and align clones

```R
library(Rsamtools)
library(GenomicRanges)
library(copykit)
library(copykat)
library(circlize)
library(dendextend)
library(RColorBrewer)
library(ComplexHeatmap)
library(parallel)
library(BiocParallel)
library(Matrix)
library(rliger)
library(ggplot2)
library(amethyst)
set.seed(111)

project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)

met_amethyst<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/merged_data/deprocated_analysis/07_scaledcis.integrated_celltyping.amethyst.rds")
rna_seurat<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")

#prepare metadata
rna_metadata<-rna_seurat@meta.data
met_metadata<-met_amethyst@metadata

rna_metadata$cnv_clonename<-NA
rna_metadata$Sample<-rna_metadata$sample
columns_to_keep<-intersect(colnames(rna_metadata),colnames(met_metadata))
row.names(rna_metadata)<-paste0("rna_",row.names(rna_metadata))
row.names(met_metadata)<-paste0("met_",row.names(met_metadata))

metadata<-rbind(rna_metadata[,columns_to_keep],
                met_metadata[,columns_to_keep])

#set metadata
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
sample_name="BCMDCIS41T"
resolution="500kb"

met<-readRDS(file=paste0("/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/copykit/",sample_name,"/copykit.",sample_name,".",resolution,".rds"))
rna<-readRDS(file=paste0("/data/rmulqueen/projects/scalebio_dcis/rna","/copykat/",sample_name,"/copykat.",sample_name,".100.rds"))

#copykit saves location as bed file (kinda) just get overlapping position value per cell for met
rna_bed<-rna$CNAmat[,1:3]
rna_cnv<-rna$CNAmat[,4:ncol(rna$CNAmat)]
colnames(rna_bed)<-c("chr","start","end")
rna_bed$end<-rna_bed$start+1
rna_bed$chr<-paste0("chr",rna_bed$chr)
rna_bed<-makeGRangesFromDataFrame(rna_bed)

#assign same regions as rna_bed for methylation cnvs
cnv_overlap<-findOverlaps(rna_bed,met@rowRanges,select="first")

#filter rna mat
rna_cnv<-rna_cnv[!is.na(cnv_overlap),]
cnv_overlap<-cnv_overlap[!is.na(cnv_overlap)]

#align met mat
met_cnv<-met@assays@data$logr[cnv_overlap,]

nrow(rna_cnv)==nrow(met_cnv)
row.names(met_cnv)<-row.names(rna_cnv)

sample_metadata<-metadata[metadata$Sample %in% sample_name,]

colnames(rna_cnv)<- gsub(colnames(rna_cnv),pattern=".1$",replacement="-1") #rename copykat output

sample_metadata <- rbind(sample_metadata[paste0("met_",colnames(met_cnv)),],
                        sample_metadata[paste0("rna_",colnames(rna_cnv)),])



#make liger object
cnvList <- list(
    "rna" = Matrix(as.matrix(rna_cnv),sparse=TRUE),
    "met" = Matrix(as.matrix(met_cnv),sparse=TRUE)
)

ligerObj <- createLiger(cnvList,removeMissing=TRUE)

#now follow liger for integration
ligerObj <- ligerObj %>% 
            rliger::normalize() %>%
            rliger::selectGenes(useDatasets = "met",thresh=0.01) %>%
            scaleNotCenter()

ligerObj <- rliger::runIntegration(ligerObj, k = 5)
ligerObj<- rliger::quantileNorm(ligerObj)
ligerObj <- rliger::runCluster(ligerObj, nNeighbors = 20)
ligerObj <- rliger::runUMAP(ligerObj,distance="cosine",minDist=0.5,n_neighbors=50)

#to do add cell metadata from met and rna obj
ligerObj@cellMeta$cnv_clonename<-NA
ligerObj@cellMeta[row.names(sample_metadata),]$cnv_clonename<-sample_metadata$cnv_clonename

plt1<-plotDatasetDimRed(ligerObj)
plt2<-plotClusterDimRed(ligerObj,"cnv_clonename")
#plt3<-plotClusterDimRed(ligerObj,"sample")
#plt4<-plotClusterDimRed(ligerObj)
ggsave(plt1|plt2,file=paste0(project_data_directory,"/test.cnv.integration2.pdf"),width=10,height=5)

#ggsave((plt1|plt2)/(plt3|plt4),file=paste0(project_data_directory,"/integration/","scaledcis.met_rna.integrated.umap.pdf"),width=25,height=25)

#coembed with liger for clones

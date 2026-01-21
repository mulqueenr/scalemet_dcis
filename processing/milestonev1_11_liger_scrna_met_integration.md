integrate with liger 
https://welch-lab.github.io/liger/articles/rna_methylation.html

https://pmc.ncbi.nlm.nih.gov/articles/PMC8132955/

BiocManager::install("HDF5Array")
install.packages('leidenAlg')
install.packages('rliger')

```R
#source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
library(amethyst)
library(Seurat)
library(rliger)
library(data.table)
library(dplyr)
library(msigdbr)
library(fgsea)
library(ggplot2)
library(GenomicRanges)
library(rliger)
library(Matrix)
library(parallel)
library(patchwork)
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)

obj<-readRDS(file="06_scaledcis.cnv_clones.amethyst.rds")

#################################
#Get DMR per Celltypes
#################################

#reading in pre-computed windows per cell type, to determine DMRs
#from scalemet_dcis/processing/milestonev1_08_DiffMet_analysis.md
dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")
dmr_celltype_outdir=paste(sep="/",dmr_outdir,"celltype")
collapsed_dmrs<-readRDS(file=paste0(dmr_celltype_outdir,"/","celltype_allcells",".dmr_filt_collapse.rds"))

#filter DMRS that are different by cell type to
# significant 
# gene overlap
# hypomethylated
# logFC atleast 1.5

table(collapsed_dmrs$type)
#   basal    bcell   cancer     endo    endo2   fibro1   fibro2    lumhr 
#   97038   174809   116983   128119   184160   130668   131439   134854 
#  lumsec myeloid1 myeloid2     peri    tcell 
#  126467    80008   130106   133283   129342 

#top 1000 per celltype
dmrs_for_integration <- collapsed_dmrs %>% 
                          filter(direction=="hypo") %>%
                          filter(dmr_padj < 0.05) %>%
                          filter(abs(dmr_logFC) > 1.5) %>% 
                          filter(dmr_length>5000) %>% 
                          filter(dmr_length<100000) %>% 
                          filter(gene_names != "NA") %>%
                          group_by(type) %>%
                          slice_min(dmr_logFC,n=5000)

table(dmrs_for_integration$type)
#   basal    bcell   cancer     endo    endo2   fibro1   fibro2    lumhr 
#    3948     4693     3110     4796     1062     5000     1777     4429 
#  lumsec myeloid1 myeloid2     peri    tcell 
#    5001     2455     1681     4397     1890 

genes_for_integration <- dmrs_for_integration %>% 
                          as.data.frame() %>%
                          select(gene_names) %>% 
                          unlist %>%
                          strsplit(split=",") %>% 
                          unlist() %>% 
                          stringr::str_replace_all(" ","") %>% 
                          unique() %>%
                          sort()


dmr_out<-dmrs_for_integration %>% 
        as.data.frame() %>% 
        select(chr,dmr_start,dmr_end) %>% 
        distinct(chr,dmr_start,dmr_end) %>% 
        makeGRangesFromDataFrame() 

dmr_bed<-data.frame(seqnames=seqnames(dmr_out),start=start(dmr_out),end=end(dmr_out))

#assign these windows to gene names
dmrs_for_integration$met_feature<-paste(dmrs_for_integration$chr,dmrs_for_integration$dmr_start,dmrs_for_integration$dmr_end,sep="_")

#key of meta features to gene names
dmrs_to_gene_names<-mclapply(1:nrow(dmrs_for_integration),function(i){
  met_feature<-dmrs_for_integration$met_feature[i]
  gene_names<-dmrs_for_integration$gene_names[i]
  gene_names_list<-gene_names %>% 
                  strsplit(split=",") %>% 
                  unlist() %>% 
                  stringr::str_replace_all(" ","")
  return(cbind(met_feature,gene_names_list))
  },mc.cores=50)

dmrs_to_gene_names<-as.data.frame(do.call("rbind",dmrs_to_gene_names))

#create new matrix from DMR sites for refined clustering
obj@genomeMatrices[["celltype_dmr_sites"]] <- makeWindows(obj, 
                                                     bed = dmr_bed,
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 50, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

#add filter to DMRs with atleast 10% coverage per cell
met_features<-which(rowSums(!is.na(obj@genomeMatrices[["celltype_dmr_sites"]]))>(ncol(obj@genomeMatrices[["celltype_dmr_sites"]])*0.1))

#set new matrix with average DMR score per gene using the key above
#add a weighted mean by size?????
gene_methylation<-mclapply(unique(dmrs_to_gene_names$gene_names_list), function(i){
  gene_name=i
  met_features<-dmrs_to_gene_names[dmrs_to_gene_names$gene_names_list==gene_name,1]
  return(colMeans(obj@genomeMatrices[["celltype_dmr_sites"]][row.names(obj@genomeMatrices[["celltype_dmr_sites"]]) %in% met_features,],na.rm=TRUE))
},mc.cores=100)

gene_met<-as.data.frame(do.call("rbind",gene_methylation))
row.names(gene_met)<-unique(dmrs_to_gene_names$gene_names_list)

# Create from raw score matrices MET
met_dat <- Matrix(as.matrix(gene_met), sparse = TRUE) # Convert the dense matrix to a sparse dgCMatrix

#set NA values to 0 (following amethyst example)
met_dat[which(is.na(met_dat),arr.ind=TRUE)]<-0

system(paste0("mkdir -p ",project_data_directory,"/integration/"))
saveRDS(met_dat,file=paste0(project_data_directory,"/integration/","scaledcis.gene_dmr_methylation.mat.rds"))

# Create from raw score matrices RNA (counts)
rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")
rna_dat<-JoinLayers(rna[["RNA"]])
rna_dat<-LayerData(rna_dat,features=row.names(met_dat),assay="RNA",layer="counts")
saveRDS(rna_dat,file=paste0(project_data_directory,"/integration/","scaledcis.gene_rna.mat.rds"))

#create shared metadata
meta<-rbind(obj@metadata[c("sample","fine_celltype")],rna@meta.data[c("sample","fine_celltype")])

#make liger object
met.liger <- createLiger(rawData=list(met=met_dat,rna=rna_dat), 
                        modal=c("meth","rna"),cellMeta=meta,
                        removeMissing=FALSE)

#now follow liger for integration
rna.met <- met.liger %>% 
            rliger::normalize() %>%
            rliger::selectGenes(useDatasets = "rna") %>%
            rliger::scaleNotCenter()

rna.met <- runIntegration(rna.met, k = 20)
rna.met <- quantileNorm(rna.met)
rna.met <- runCluster(rna.met, nNeighbors = 100)
rna.met <- runUMAP(rna.met)

plt1<-plotDatasetDimRed(rna.met)
plt2<-plotClusterDimRed(rna.met,"fine_celltype")
plt3<-plotClusterDimRed(rna.met,"sample")
plt4<-plotClusterDimRed(rna.met)

ggsave((plt1|plt2)/(plt3|plt4),file=paste0(project_data_directory,"/integration/","scaledcis.met_rna.integrated.umap.pdf"),width=25,height=25)

saveRDS(rna.met,file=paste0(project_data_directory,"/integration/","scaledcis.met_rna.integrated.liger.rds"))


```

Need to figure out best way to summarize methylation over genes. I think its 
cell type level DMRs > assign DMRs to genes by overlap > window per cell for those genes > filter RNA to just those genes


Important parameters of quantile_norm are as follows:

knn_k. This sets the number of nearest neighbors for the within-dataset KNN graph. The default is 20.

quantiles. This sets the number of quantiles to use for quantile normalization. The default is 50.

min_cells. This indicates the minimum number of cells to consider a cluster as shared across datasets. The default is 20.

dims.use. This sets the indices of factors to use for quantile normalization. The user can pass in a vector of indices indicating specific factors. This is helpful for excluding factors capturing biological signals such as the cell cycle or technical signals such as mitochondrial genes. The default is all k of the factors.

do.center. This indicates whether to center when scaling cell factor loadings. The default is FALSE. This option should be set to TRUE when cell factor loadings have a mean above zero, as with dense data such as DNA methylation.

max_sample. This sets the maximum number of cells used for quantile normalization of each cluster and factor. The default is 1000.

refine.knn. This indicates whether to increase robustness of cluster assignments using KNN graph. The default is TRUE. Disabling this option makes the function run faster but may increase the amount of spurious alignment in divergent datasets.

eps. This sets the error bound of the nearest neighbor search. The default is 0.9. Lower values give more accurate nearest neighbor graphs but take much longer to compute. We find that this parameter affects result quality very little.

ref_dataset. This indicates the name of the dataset to be used as a reference for quantile normalization. By default, the dataset with the largest number of cells is used.
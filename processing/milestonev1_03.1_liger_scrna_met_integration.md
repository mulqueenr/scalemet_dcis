Integration of seurat object with methylation object.

Could also potentially run DMRs on broad_celltypes assigned. But I'm going to start at the cluster level.

1. Using coarse cluster level DMRs (from scalemet_dcis/processing/milestonev1_01_amethyst_coarse_celltyping.md)
2. assign DMRs to genes by overlap 
3. window per cell for those genes 
4. filter RNA to just those genes 
5. integrate with liger following:
https://welch-lab.github.io/liger/articles/rna_methylation.html
https://pmc.ncbi.nlm.nih.gov/articles/PMC8132955/

Use integration to inform cell typing. Note some cell types were mislabelled based on small set of marker genes originally used. 
Using integration data, marker genes, clones and clustering for final calls. 

Add filter from RNA to limit genes that just overlap piggybacking?

```R
BiocManager::install("HDF5Array")
install.packages('leidenAlg')
install.packages('rliger')
```


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
library(ComplexHeatmap)
library(circlize)
set.seed(111)
options(future.globals.maxSize= 200000*1024^2) #80gb limit for parallelizing
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
dmr_celltype_outdir=paste(sep="/",dmr_outdir,"dmr_coarse_cluster")
collapsed_dmrs<-readRDS(file=paste0(dmr_celltype_outdir,"/","coarse_cluster",".dmr_filt_collapse.rds"))

#filter DMRS that are different by cell type to
# significant 
# gene overlap
# hypomethylated
# logFC atleast 1.5

table(collapsed_dmrs$type)
#   1     10     11     12     13     14     15     16     17     18     19 
#184447 185076 182009 177556 175150 173075 184027 167856 178465 157240 173844 
#     2     20      3      4      5      6      7      8      9 
#158315 161403 181188 176171 182814 171735 173040 137075 182240 

dmrs_for_integration <- collapsed_dmrs %>% 
                          filter(direction=="hypo") %>%
                          filter(dmr_padj < 0.05) %>%
                          filter(abs(dmr_logFC) > 1.5) %>% 
                          filter(dmr_length<50000) %>% 
                          filter(gene_names != "NA") %>%
                          group_by(type) 

table(dmrs_for_integration$type)
#  1   10   11   12   13   14   15   16   17   18   19    2   20    3    4    5 
#1385 3330 3714 4934 4394 3913 5510 6089 2707 2945  445 3208 4929 3931 3056 1426 
#   6    7    8    9 
#4879 1119 2147 3131 

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
        makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>%
        reduce()

dmr_out<-dmr_out[width(dmr_out)<50000,]
summary(width(dmr_out))

#assign genes to merged dmr_out features
reduced_overlap<-findOverlaps(dmr_out,
            makeGRangesFromDataFrame(dmrs_for_integration,keep.extra.columns=TRUE))


#add list of gene names and store original dmr sites in bed
dmr_out$gene_names<-NA
dmr_out$gene_names<-mclapply(unique(queryHits(reduced_overlap)), function(x){
    genes=dmrs_for_integration[subjectHits(reduced_overlap[queryHits(reduced_overlap)==x,]),]$gene_names
    genes <- genes %>% strsplit(split=",") %>% unlist() %>% stringr::str_replace_all(" ","") %>% unique() %>% paste(collapse=",")
    genes<-paste(unlist(genes),sep=",")
    dmr_out[x,]$gene_names<-genes
},mc.cores=100)
dmr_out$met_feature<-paste(seqnames(dmr_out),start(dmr_out),end(dmr_out),sep="_")


#create new matrix from DMR sites for refined clustering
obj@genomeMatrices[["coarse_cluster_dmr_sites"]] <- makeWindows(obj, 
                                                     bed = dmr_bed,
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 50, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

saveRDS(obj,file="07_scaledcis.rna_integrated.amethyst.rds")
obj<-readRDS(file="07_scaledcis.rna_integrated.amethyst.rds")

#add filter to DMRs with atleast 10% coverage per cell
met_features<-which(rowSums(!is.na(obj@genomeMatrices[["coarse_cluster_dmr_sites"]]))>(ncol(obj@genomeMatrices[["coarse_cluster_dmr_sites"]])*0.1))

#expand dmr_out so each row is its own gene name
#this will be used to summarize multiple DMRs to their genes
dmrs_to_gene_names<-lapply(1:length(dmr_out), function(x) {
  temp<-dmr_out[x,]
  genes<-unlist(temp$gene_names %>% unlist() %>% strsplit(split=",")) 
  return(cbind(genes=genes,met_feature=temp$met_feature))
})

dmrs_to_gene_key<-as.data.frame(do.call("rbind",dmrs_to_gene_names))

#set new matrix with average DMR score per gene using the key above
#add a weighted mean by size?????
gene_methylation<-mclapply(unique(dmrs_to_gene_key$genes), function(i){
  gene_name=i
  met_features<-dmrs_to_gene_key[dmrs_to_gene_key$genes==gene_name,2]
  return(colMeans(obj@genomeMatrices[["coarse_cluster_dmr_sites"]][row.names(obj@genomeMatrices[["coarse_cluster_dmr_sites"]]) %in% met_features,],na.rm=TRUE))
},mc.cores=100)

gene_met<-as.data.frame(do.call("rbind",gene_methylation))
row.names(gene_met)<-unique(dmrs_to_gene_key$genes)

# Create from raw score matrices MET
met_dat <- Matrix(as.matrix(gene_met), sparse = TRUE) # Convert the dense matrix to a sparse dgCMatrix

#set NA values to 0 (following amethyst example)
met_dat[which(is.na(met_dat),arr.ind=TRUE)]<-0

system(paste0("mkdir -p ",project_data_directory,"/integration/"))
saveRDS(met_dat,file=paste0(project_data_directory,"/integration/","scaledcis.gene_dmr_methylation.mat.rds"))

# Create from raw score matrices RNA (counts)
rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")
Idents(rna)<-rna$fine_celltype
rna_dat<-JoinLayers(rna[["RNA"]])
rna_dat<-subset(rna_dat,features=row.names(met_dat))
rna_dat <- NormalizeData(rna_dat, normalization.method = "LogNormalize", scale.factor = 10000)
rna_dat <- ScaleData(rna_dat, features = row.names(met_dat))
rna_markers<-FindAllMarkers(rna_dat,assay="RNA",only.pos=TRUE)

rna_markers_filt<-rna_markers %>% group_by(cluster) %>% filter(avg_log2FC > 3) %>% filter(p_val_adj <0.01) %>% as.data.frame()
rna_markers_filt<-rna_markers_filt[!duplicated(rna_markers_filt$gene),]

#filter met_dat to only genes in rna markers
met_dat<-met_dat[row.names(met_dat) %in% unique(rna_markers_filt$gene),]

#get raw RNA values on filtered methylation genes
rna_dat<-LayerData(rna_dat,features=row.names(met_dat),assay="RNA",layer="counts")
saveRDS(rna_dat,file=paste0(project_data_directory,"/integration/","scaledcis.gene_rna.mat.rds"))

rna$broad_celltype<-rna$fine_celltype
rna$coarse_cluster_id <- "NA"
rna$mcg_pct <- "NA"


#create shared metadata, add coarse_cluster_id to meta
meta<-rbind(obj@metadata[c("sample","broad_celltype","coarse_cluster_id","mcg_pct")],rna@meta.data[c("sample","broad_celltype","coarse_cluster_id","mcg_pct")])

#make liger object
met.liger <- createLiger(rawData=list(met=met_dat,rna=rna_dat), 
                        modal=c("meth","rna"),cellMeta=meta,
                        removeMissing=FALSE)

#now follow liger for integration
rna.met <- met.liger %>% 
            rliger::normalize() %>%
            rliger::selectGenes(useDatasets = "met") %>%
            rliger::scaleNotCenter()

rna.met <- rliger::runIntegration(rna.met, k = 20)
rna.met <- rliger::quantileNorm(rna.met)
rna.met <- runCluster(rna.met)
rna.met <- runUMAP(rna.met)

plt1<-plotDatasetDimRed(rna.met,splitBy="dataset")
ggsave(wrap_plots(plt1),file=paste0(project_data_directory,"/integration/","scaledcis.met_rna.","dataset",".integrated.umap.pdf"),width=20,height=10)

plt2<-plotClusterDimRed(rna.met,splitBy="dataset","broad_celltype")
ggsave(wrap_plots(plt2),file=paste0(project_data_directory,"/integration/","scaledcis.met_rna.","celltype",".integrated.umap.pdf"),width=20,height=10)

plt3<-plotClusterDimRed(rna.met,splitBy="dataset","sample")
ggsave(wrap_plots(plt3),file=paste0(project_data_directory,"/integration/","scaledcis.met_rna.","sample",".integrated.umap.pdf"),width=20,height=10)

plt4<-plotClusterDimRed(rna.met,splitBy="dataset","coarse_cluster_id")
ggsave(wrap_plots(plt4),file=paste0(project_data_directory,"R/integration/","scaledcis.met_rna.","coarse_cluster_id",".integrated.umap.pdf"),width=20,height=10)

plt5<-plotClusterDimRed(rna.met,slot="cellMeta","mcg_pct")
ggsave(wrap_plots(plt5),file=paste0(project_data_directory,"/integration/","scaledcis.met_rna.","mcg",".integrated.umap.pdf"),width=20,height=10)



plots <- plotGeneDimRed(rna.met, c("PTPRC", "RCOR3"), splitBy = "dataset",
                        titles = c(names(rna.met), names(rna.met)))

saveRDS(rna.met,file=paste0(project_data_directory,"/integration/","scaledcis.met_rna.integrated.liger.rds"))

rna.met<-readRDS(file=paste0(project_data_directory,"/integration/","scaledcis.met_rna.integrated.liger.rds"))


col_fun = colorRamp2(c(0, 0.5, 1), c("white", "red", "darkred"))

#plot confusion of clusters for met to rna
#plot 2 heatmaps next to each other
meta_rna<-rna.met@cellMeta[rna.met@cellMeta$dataset=="rna",]
meta_met<-rna.met@cellMeta[rna.met@cellMeta$dataset=="met",]

rna_confusion<-table(meta_rna$leiden_cluster,meta_rna$broad_celltype)
rna_confusion <- rna_confusion/rowSums(rna_confusion)
rna_confusion<-round(rna_confusion,2)

met_confusion<-table(meta_met$leiden_cluster,meta_met$broad_celltype)
met_confusion <- met_confusion/rowSums(met_confusion)
met_confusion<-round(met_confusion,2)

met_confusion_cluster<-table(meta_met$leiden_cluster,meta_met$coarse_cluster_id)
met_confusion_cluster <- met_confusion_cluster/rowSums(met_confusion_cluster)
met_confusion_cluster<-round(met_confusion_cluster,2)



#pivot to long format for plotting (or just use complexheatmap)\


plt1<-Heatmap(rna_confusion,col=col_fun,border = TRUE,rect_gp = gpar(col = "gray", lwd = 1),cluster_columns=FALSE)
plt2<-Heatmap(met_confusion,col=col_fun,border = TRUE,rect_gp = gpar(col = "gray", lwd = 1),cluster_columns=FALSE)
plt3<-Heatmap(met_confusion_cluster,col=col_fun,border = TRUE,rect_gp = gpar(col = "gray", lwd = 1),cluster_columns=FALSE)

pdf(paste0(project_data_directory,"/integration/","scaledcis.met_confusion.liger.pdf"))
print(plt1+plt2+plt3)
dev.off()


#correct met cell type info using confusion matrix
meta_met$corrected_celltype<-NA
meta_met[meta_met$leiden_cluster %in% c("13","7"),]$corrected_celltype<-"basal"
meta_met[meta_met$leiden_cluster %in% c("19","1","15","11","9"),]$corrected_celltype<-"lumhr"
meta_met[meta_met$leiden_cluster %in% c("16","21"),]$corrected_celltype<-"lumsec"
meta_met[meta_met$leiden_cluster %in% c("14"),]$corrected_celltype<-"endo_vein"
meta_met[meta_met$leiden_cluster %in% c("20"),]$corrected_celltype<-"plasma"
meta_met[meta_met$leiden_cluster %in% c("25","2"),]$corrected_celltype<-"tcell_cd4"
meta_met[meta_met$leiden_cluster %in% c("17"),]$corrected_celltype<-"tcell_nk"
meta_met[meta_met$leiden_cluster %in% c("0"),]$corrected_celltype<-"tcell_cd8"
meta_met[meta_met$leiden_cluster %in% c("12","6","4"),]$corrected_celltype<-"fibro"
meta_met[meta_met$leiden_cluster %in% c("5"),]$corrected_celltype<-"fibro_CAF"
meta_met[meta_met$leiden_cluster %in% c("10"),]$corrected_celltype<-"myeloid"
meta_met[meta_met$leiden_cluster %in% c("3"),]$corrected_celltype<-"endo"
meta_met[meta_met$leiden_cluster %in% c("3"),]$corrected_celltype<-"endo"
meta_met[meta_met$leiden_cluster %in% c("23"),]$corrected_celltype<-"epithelial"
meta_met[meta_met$leiden_cluster %in% c("18"),]$corrected_celltype<-"bcell"
meta_met[meta_met$leiden_cluster %in% c("24"),]$corrected_celltype<-"myeloid_mast"
meta_met[meta_met$leiden_cluster %in% c("22"),]$corrected_celltype<-"suspected_doublet"
meta_met[meta_met$leiden_cluster %in% c("8"),]$corrected_celltype<-"peri"


met_confusion_corrected<-table(meta_met$leiden_cluster,meta_met$corrected_celltype)
met_confusion_corrected <- met_confusion_corrected/rowSums(met_confusion_corrected)
met_confusion_corrected<-round(met_confusion_corrected,2)
plt3<-Heatmap(met_confusion_corrected,col=col_fun,border = TRUE,rect_gp = gpar(col = "gray", lwd = 1),cluster_columns=FALSE)


pdf(paste0(project_data_directory,"/integration/","scaledcis.met_confusion.corrected.liger.pdf"))
print(plt1+plt2+plt3)
dev.off()


#plot dim reduc
rna.met@cellMeta$corrected_met<-NA
rna.met@cellMeta[row.names(meta_met),]$corrected_met<-meta_met$corrected_celltype
plt1<-plotDatasetDimRed(rna.met)
plt2<-plotClusterDimRed(rna.met,"corrected_met")
plt3<-plotClusterDimRed(rna.met,"sample")
plt4<-plotClusterDimRed(rna.met)

ggsave((plt1|plt2)/(plt3|plt4),file=paste0(project_data_directory,"/integration/","scaledcis.met_rna.integrated.corrected_celltype.umap.pdf"),width=25,height=25)
saveRDS(rna.met,file=paste0(project_data_directory,"/integration/","scaledcis.met_rna.integrated.liger.rds"))


#correct celltypes in amethyst object
meta<-rna.met@cellMeta[rna.met@cellMeta$dataset=="met",]
row.names(meta)<-gsub(pattern="met_",replace="",row.names(meta))
obj@metadata$integrated_celltype<-NA
obj@metadata[row.names(meta),]$integrated_celltype<-meta$corrected_met
#and convert lumhr aneuploid back to cancer celltype
obj@metadata[obj@metadata$cnv_ploidy_500kb=="aneuploid",]$integrated_celltype<-"cancer"

saveRDS(obj,file="07_scaledcis.integrated_celltyping.amethyst.rds")


```


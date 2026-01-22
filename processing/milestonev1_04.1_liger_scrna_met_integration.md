Integration of seurat object with methylation object.

1. Using cell type level DMRs 
2. assign DMRs to genes by overlap 
3. window per cell for those genes 
4. filter RNA to just those genes 
5. integrate with liger following:
https://welch-lab.github.io/liger/articles/rna_methylation.html
https://pmc.ncbi.nlm.nih.gov/articles/PMC8132955/

Use integration to finalize cell typing. Note some cell types were mislabelled based on small set of marker genes originally used. 
Using integration data for final classifications.

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

rna.met<-readRDS(file=paste0(project_data_directory,"/integration/","scaledcis.met_rna.integrated.liger.rds"))


#plot confusion of clusters for met to rna
#plot 2 heatmaps next to each other
meta_rna<-rna.met@cellMeta[rna.met@cellMeta$dataset=="rna",]
meta_met<-rna.met@cellMeta[rna.met@cellMeta$dataset=="met",]

rna_confusion<-table(meta_rna$leiden_cluster,meta_rna$fine_celltype)
rna_confusion <- rna_confusion/rowSums(rna_confusion)
rna_confusion<-round(rna_confusion,2)

met_confusion<-table(meta_met$leiden_cluster,meta_met$fine_celltype)
met_confusion <- met_confusion/rowSums(met_confusion)
met_confusion<-round(met_confusion,2)
#pivot to long format for plotting (or just use complexheatmap)\

col_fun = colorRamp2(c(0, 0.5, 1), c("white", "red", "darkred"))

plt1<-Heatmap(rna_confusion,col=col_fun,border = TRUE,rect_gp = gpar(col = "gray", lwd = 1),cluster_columns=FALSE)
plt2<-Heatmap(met_confusion,col=col_fun,border = TRUE,rect_gp = gpar(col = "gray", lwd = 1),cluster_columns=FALSE)

pdf(paste0(project_data_directory,"/integration/","scaledcis.met_confusion.liger.pdf"))
print(plt1+plt2)
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


```R
library(Seurat)
library(hdf5r)
library(ggplot2)
library(patchwork)
setwd("/data/rmulqueen/projects/scalebio_dcis/rna")

#named vectors of sample name (as they are in the methylation, combined with 10x cellranger output)
rna_dat<-c(
'BCMDCIS05T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS05T',
'BCMDCIS07T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS07T',
'BCMDCIS102T_24hTis'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS102T_24hTis',
'BCMDCIS124T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS124T',
'BCMDCIS22T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS22T',
'BCMDCIS28T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS28T',
'BCMDCIS32T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS32T',
'BCMDCIS35T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS35T',
'BCMDCIS41T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS41T',
'BCMDCIS49T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS49T',
'BCMDCIS52T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS52T',
'BCMDCIS65T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS65T',
'BCMDCIS66T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS66T',
'BCMDCIS70T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS70T',
'BCMDCIS74T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS74T',
'BCMDCIS79T_24hTis_DCIS'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS79T_24hTis_DCIS',
'BCMDCIS79T_24hTis_IDC'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS79T_24hTis_IDC',
'BCMDCIS80T_24hTis'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS80T_24hTis',
'BCMDCIS82T_24hTis'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS82T_24hTis',
'BCMDCIS92T_24hTis'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS92T_24hTis',
'BCMDCIS94T_24hTis'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS94T_24hTis',
'BCMDCIS97T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS97T',
'BCMDCIS99T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMDCIS99T',
'ECIS25T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/ECIS25T',
'ECIS26T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/ECIS26T',
'ECIS36T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/ECIS36T',
'ECIS48T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/ECIS48T',
'ECIS57T'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/ECIS57T',
'BCMHBCA03R'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMHBCA03R',
'BCMHBCA04R'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMHBCA04R',
'BCMHBCA09R-3h'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMHBCA09R-3h',
'BCMHBCA12R-3h'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMHBCA12R-3h',
'BCMHBCA16R-3h'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMHBCA16R-3h',
'BCMHBCA16R-3h-nuc'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMHBCA16R-3h-nuc',
'BCMHBCA17R-3h'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMHBCA17R-3h',
'BCMHBCA17R-3h-nuc'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMHBCA17R-3h-nuc',
'BCMHBCA19R-4h'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMHBCA19R-4h',
'BCMHBCA19R-4h-nuc'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMHBCA19R-4h-nuc',
'BCMHBCA22R-4h'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMHBCA22R-4h',
'BCMHBCA22R-4h-nuc'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMHBCA22R-4h-nuc',
'BCMHBCA26L-24hTis-4h'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMHBCA26L-24hTis-4h',
'BCMHBCA26L-24hTis-4h-nuc'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMHBCA26L-24hTis-4h-nuc',
'BCMHBCA29L-2h'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMHBCA29L-2h',
'BCMHBCA38L-3h'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMHBCA38L-3h',
'BCMHBCA83L-3h'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMHBCA83L-3h',
'BCMHBCA85L-3h'='/data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output/BCMHBCA85L-3h')

make_seurat_object<-function(x){
    sample_name=names(rna_dat[x])
    sample_dir=paste0(unname(rna_dat[x]),"/","outs")
    print(paste("Processing...",sample_name))

    #read in data
    counts <- Read10X_h5(paste0(sample_dir,"/filtered_feature_bc_matrix.h5")) #count data
    # create a Seurat object containing the RNA data
    dat <- CreateSeuratObject(
        counts = counts,
        assay = "RNA")
    dat$sample<-sample_name
    dat<-RenameCells(dat,add.cell.id=sample_name)
  return(dat)
}

out<-lapply(1:length(rna_dat), make_seurat_object)

dat <- merge(out[[1]], y = as.list(out[2:length(out)]), project = "all_data")

#following standard seurat intial processing
dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^MT-")
dat <- subset(dat, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 5)
saveRDS(dat, file = "tenx_dcis.rds")
```


```R
library(Seurat)
library(hdf5r)
library(ggplot2)
library(patchwork)
setwd("/data/rmulqueen/projects/scalebio_dcis/rna")
obj<-readRDS(file = "tenx_dcis.rds")

#filtering to only cellular RNA for now
obj@meta.data$assay="cellular_RNA"
obj@meta.data[endsWith(obj@meta.data$sample,suffix="nuc"),]$assay="nuclear_RNA"
obj<-subset(obj,assay=="cellular_RNA")

prefix="coarse_celltype"
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = rownames(obj))
obj <- RunPCA(obj, features = VariableFeatures(object = obj))
plt<-ElbowPlot(obj)
ggsave(plt,file=paste("tenx_dcis","elbow",prefix,"sample.pdf",sep="."))
saveRDS(obj, file = "tenx_dcis.rds")

#function for umap clustering
cluster_object<-function(obj=obj,dims=1:15,prefix="coarse_celltype",res=0.2){
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    obj <- RunPCA(obj, features = VariableFeatures(object = obj))
    obj <- FindNeighbors(obj, dims = dims)
    obj <- FindClusters(obj, resolution = res)
    obj <- RunUMAP(obj, dims = dims)
    return(obj)
}
obj<-cluster_object(obj=obj,dims=1:15,prefix="coarse_celltype",res=0.2)
saveRDS(obj, file = "tenx_dcis.rds")

plt<-DimPlot(obj,group.by="sample",raster=T,label=TRUE)
ggsave(plt,file=paste("tenx_dcis","umap",prefix,"sample.pdf",sep="."),width=30,height=20,limitsize=FALSE)
plt<-DimPlot(obj,group.by="seurat_clusters",raster=T,label=TRUE)
ggsave(plt,file=paste("tenx_dcis","umap",prefix,"clusters.pdf",sep="."),width=30,height=20,limitsize=FALSE)

#plot hbca vs cancer
obj@meta.data$diag="cancer"
obj@meta.data[grepl(obj@meta.data$sample,pattern="HBCA"),]$diag="HBCA"
plt<-DimPlot(obj,group.by="diag",raster=T,label=TRUE)
ggsave(plt,file=paste("tenx_dcis","umap",prefix,"diag.pdf",sep="."),width=30,height=20,limitsize=FALSE)

#same cell markers list as used for methylation
cell_markers<-list()
cell_markers[["basal"]]<-c("CARMN","ACTA2","KRT17","KRT14","DST","KRT5")
cell_markers[["lumhr"]]<-c("AREG","AZGP1","KRT18","AGR2","PIP","ANKRD30A")
cell_markers[["lumsec"]]<-c("GABRP","ELF5","CCL28","KRT15","KIT","MMP7","LTF","SLPI")

cell_markers[["fibro"]]<-c("DCN","APOD","LUM","COL1A2","COL1A1","FAP")
cell_markers[["endo"]]=c("CCL21","TFF3","MMRN1","CLDN5","AL357507.1","PKHD1L1","KLHL4","LINC02147","RHOJ","ST6GALNAC3","MECOM","BTNL9","MCTP1","PTPRB","VWF","ADGRL4","LDB2")
cell_markers[["perivasc"]]=c("RGS6","KCNAB1","COL25A1","ADGRL3","PRKG1","MYL9","ADIRF","NR2F2-AS1","AC012409.2")

cell_markers[["myeloid"]]<-c("HLA-DRA","HLA-DPA1","CD74")
cell_markers[["tcell"]]<-c("PTPRC","IKZF1","IL7R","GNLY")
cell_markers[["bcell"]]=c("CD37","TCL1A","LTB","HLA-DPB1")
cell_markers[["plasma"]]=c("IGHA2","IGHA1","JCHAIN","IGHM","IGHG1","IGHG4","IGHG3","IGHG2")

cell_markers[["adipo"]]=c("PDE3B","ACACB","WDPCP","PCDH9","CLSTN2","ADIPOQ","TRHDE")

plt<-DotPlot(obj, features = cell_markers, group.by = "seurat_clusters",cluster.idents=TRUE) + RotatedAxis()
ggsave(plt,file="tenx_dcis.coarse_celltype.dotplot.pdf",width=20,height=20,limitsize=FALSE)

#pretty much everything over 25ish is suspected double
obj@meta.data$coarse_celltype<-"suspected_doublet" #32, 30, 21, 26, 31, 10, 29, 27, 33, 25, 23, 22, 11
obj@meta.data[obj@meta.data$seurat_clusters %in% c("5","3","14","15","17"),]$coarse_celltype<-"lumhr"
obj@meta.data[obj@meta.data$seurat_clusters %in% c("2"),]$coarse_celltype<-"basal"
obj@meta.data[obj@meta.data$seurat_clusters %in% c("9","19"),]$coarse_celltype<-"lumsec"

obj@meta.data[obj@meta.data$seurat_clusters %in% c("8","1","13"),]$coarse_celltype<-"fibro" #13 suspected CAFs
obj@meta.data[obj@meta.data$seurat_clusters %in% c("4","23","24"),]$coarse_celltype<-"endo"
obj@meta.data[obj@meta.data$seurat_clusters %in% c("6"),]$coarse_celltype<-"perivasc"

obj@meta.data[obj@meta.data$seurat_clusters %in% c("0"),]$coarse_celltype<-"tcell"
obj@meta.data[obj@meta.data$seurat_clusters %in% c("12"),]$coarse_celltype<-"bcell"
obj@meta.data[obj@meta.data$seurat_clusters %in% c("7","28","18"),]$coarse_celltype<-"myeloid"
obj@meta.data[obj@meta.data$seurat_clusters %in% c("16"),]$coarse_celltype<-"plasma"

plt<-DimPlot(obj,group.by="coarse_celltype",raster=F)
ggsave(plt,file="tenx_dcis.umap.coarse_celltype.pdf",width=20,height=20,limitsize=FALSE)
saveRDS(obj, file = "tenx_dcis.rds")

#some more umap groupings look like doublets to filter more during compartment fine celltyping

```

# Stromal celltyping

## Fibro
```R
dat_sub<-subset(obj,coarse_celltype %in% c("fibro"))
dat_sub<-cluster_object(obj=dat_sub,dims=1:15,prefix="stromal_fibro",res=0.2)

stromal_markers<-list()
stromal_markers[["fibro_major"]]<-c('KDM6B', 'HMOX1', 'GPRC5A', 'TNFAIP6', 'TNFRSF12A', 'HMGA1', 'CXCL1', 'NOP16', 'ATP13A3', 'CXCL2', 'GEM', 'MAT2A', 'DDX21', 'KLHL21', 'EFNB2', 'AMD1', 'PNO1', 'GCLM', 'WDR43', 'ATP1B3')
stromal_markers[["fibro_prematrix"]]<-c('PCOLCE2', 'WISP2', 'MFAP5', 'CHRDL1', 'GPX3', 'EFEMP1', 'ADH1B', 'TNXB', 'MGST1', 'FBLN2', 'PTGIS', 'GPC3', 'NFIB', 'C3', 'ABCA8', 'TGFBR2', 'PODN', 'ITIH5', 'PRELP', 'FHL1')
stromal_markers[["fibro_CAF"]]<-c('POSTN', 'COL10A1', 'INHBA', 'COL8A1', 'CTHRC1', 'COL5A1', 'MMP11', 'MFAP2', 'ITGBL1', 'BGN', 'SULF1', 'FAP', 'ASPN', 'COMP', 'SULF2', 'COL11A1', 'COL5A2', 'NREP', 'LOXL1', 'RAB31')
stromal_markers[["fibro_matrix"]]<-c('DLK1', 'IGFBP2', 'APOE', 'TAC1', 'ADAM12', 'GPM6B', 'TGFBI', 'COL15A1', 'EDNRB', 'DCX', 'CAPN6', 'ITM2A', 'VEGFD', 'IGF1', 'RELN', 'MMP16', 'OLFML2A', 'EFNB2', 'COL4A2', 'PLXDC1')
stromal_markers[["fibro_SFRP4"]]<-c('SFRP4', 'CD9', 'CLU', 'LTBP1', 'GREM1', 'IGFBP5', 'LEPR', 'OGN', 'KCNMA1', 'ABI3BP', 'CRYAB', 'ADIRF', 'RRAD', 'C12orf75', 'CAV1', 'PAPPA2', 'PTGIS', 'IGFBP6', 'PLP2', 'NFIB')
stromal_markers[["fibro_stress"]]<-c('TXNIP', 'FOS', 'HSPA1B', 'NR2F1', 'HSPA6', 'SPRY1', 'JUN', 'ABCA10', 'SOCS3', 'EGR1', 'ZFP36', 'HSPA1A', 'APOD', 'BTG2', 'DDIT3', 'PTGDS', 'DDIT4', 'MAFB', 'MRPL18', 'ABCA6')

dat_sub<-FindClusters(dat_sub,res=0.25)
dat_sub$stromal_fibro_subclusters<-dat_sub$seurat_clusters
plt_dim<-DimPlot(dat_sub,group.by="stromal_fibro_subclusters",label=TRUE)

plt_list<-lapply(1:length(stromal_markers),function(x) {
    plt<-DotPlot(dat_sub, 
        features = stromal_markers[x], 
        group.by = "stromal_fibro_subclusters",
        cluster.idents=FALSE) + 
    RotatedAxis() + ggtitle(names(stromal_markers)[x])
    return(plt)})

plt<-wrap_plots(plt_list,nrow=1,guides="collect")

layout <- "
A######
BBBBBBB"

ggsave(plt_dim+plt+plot_layout(design=layout),file="tenx_dcis.stromal.fibro.fine_celltype.dotplot.pdf",width=0.25*length(unlist(stromal_markers)),height=10,limitsize=FALSE)

#labelling fibro subtypes
dat_sub@meta.data$fine_celltype<-"suspected_doublet"
dat_sub@meta.data[dat_sub@meta.data$stromal_fibro_subclusters %in% c("0"),]$fine_celltype<-"fibro_major"
dat_sub@meta.data[dat_sub@meta.data$stromal_fibro_subclusters %in% c("2"),]$fine_celltype<-"fibro_prematrix"
dat_sub@meta.data[dat_sub@meta.data$stromal_fibro_subclusters %in% c("1"),]$fine_celltype<-"fibro_CAF"
dat_sub@meta.data[dat_sub@meta.data$stromal_fibro_subclusters %in% c("4"),]$fine_celltype<-"fibro_matrix"
dat_sub@meta.data[dat_sub@meta.data$stromal_fibro_subclusters %in% c("6"),]$fine_celltype<-"fibro_SFRP4"
dat_sub@meta.data[dat_sub@meta.data$stromal_fibro_subclusters %in% c("3"),]$fine_celltype<-"fibro_SFRP4"

plt_dim<-DimPlot(dat_sub,group.by=c("stromal_fibro_subclusters","fine_celltype"),label=TRUE)
ggsave(plt_dim,file="tenx_dcis.stromal.fibro.fine_celltype.umap.pdf",width=8)

saveRDS(dat_sub, file = "tenx_dcis.stromal.fibro.rds")
```


## Endo

```R
dat_sub<-subset(obj,coarse_celltype %in% c("endo"))
dat_sub
dat_sub<-cluster_object(obj=dat_sub,dims=1:15,prefix="stromal_endo",res=0.2)

stromal_markers<-list()
stromal_markers[["endo_artery"]]<-c('SEMA3G', 'CXCL12', 'IGFBP3', 'DEPP1', 'HEY1', 'LTBP4', 'PPP1R14A', 'FN1', 'PCSK5', 'ARL15', 'RHOB', 'FBLN5', 'TSPAN2', 'GUCY1A1', 'EGR1', 'SLC9A3R2', 'ICAM2', 'GJA4', 'NEBL', 'ANXA3')
stromal_markers[["endo_vein"]]<-c('ACKR1', 'SELE', 'SELP', 'IL1R1', 'ZNF385D', 'CLU', 'CNKSR3', 'VCAN', 'CCL14', 'NCOA7', 'EPB41L3', 'TLL1', 'OLFM1', 'PRCP', 'VCAM1', 'GNG12', 'NPC2', 'FAM84B', 'CYP1B1', 'ACTN1')
stromal_markers[["endo_capillary"]]<-c('BTNL9', 'CD300LG', 'MT1M', 'RGCC', 'MT1A', 'C11orf96', 'MT1E', 'ADGRF5', 'CA4', 'CD36', 'RBP7', 'MSX1', 'GPIHBP1', 'ITGA1', 'SGK1', 'MT1X', 'SCARB1', 'LITAF', 'KDR', 'MLEC')
stromal_markers[["endo_TEC"]]<-c('MMP2', 'INSR', 'HECW2', 'PODXL', 'ESM1', 'RGCC', 'NOTCH4', 'MAP1B', 'HTRA1', 'MLEC', 'COL4A2', 'PMEPA1', 'VWA1', 'PTPRG', 'COL6A2', 'NID1', 'TCIM', 'PDGFD', 'PLPP1', 'THY1')
stromal_markers[["endo_cycling"]]<-c('TYMS', 'NUSAP1', 'PCLAF', 'TOP2A', 'CENPF', 'MKI67', 'TK1', 'TPX2', 'GGH', 'CDK1', 'BIRC5', 'ASPM', 'CENPK', 'DHFR', 'DIAPH3', 'ZWINT', 'GTSE1', 'CENPW', 'UBE2C', 'KNL1')
stromal_markers[["endo_lymphatic"]]<-c('CCL21', 'TFF3', 'MMRN1', 'EFEMP1', 'PDPN', 'MRC1', 'NR2F1', 'LYVE1', 'PROX1', 'LAPTM5', 'COLEC12', 'IGF1', 'SEMA3D', 'MPP7', 'AKAP12', 'PKHD1L1', 'ABI3BP', 'CCND1', 'TFPI', 'FXYD6')

dat_sub<-FindClusters(dat_sub,res=0.25)
dat_sub$stromal_endo_subclusters<-dat_sub$seurat_clusters
plt_dim<-DimPlot(dat_sub,group.by="stromal_endo_subclusters",label=TRUE)

plt_list<-lapply(1:length(stromal_markers),function(x) {
    plt<-DotPlot(dat_sub, 
        features = stromal_markers[x], 
        group.by = "stromal_endo_subclusters",
        cluster.idents=FALSE) + 
    RotatedAxis() + ggtitle(names(stromal_markers)[x])
    return(plt)})

plt<-wrap_plots(plt_list,nrow=1,guides="collect")

layout <- "
A######
BBBBBBB"

ggsave(plt_dim+plt+plot_layout(design=layout),file="tenx_dcis.stromal.endo.fine_celltype.dotplot.pdf",width=0.25*length(unlist(stromal_markers)),height=10,limitsize=FALSE)

#labelling endo subtypes
dat_sub@meta.data$fine_celltype<-"suspected_doublet"
dat_sub@meta.data[dat_sub@meta.data$stromal_endo_subclusters %in% c("5"),]$fine_celltype<-"endo_unknown"
dat_sub@meta.data[dat_sub@meta.data$stromal_endo_subclusters %in% c("4"),]$fine_celltype<-"endo_artery"
dat_sub@meta.data[dat_sub@meta.data$stromal_endo_subclusters %in% c("7","2","0"),]$fine_celltype<-"endo_vein"
dat_sub@meta.data[dat_sub@meta.data$stromal_endo_subclusters %in% c("1"),]$fine_celltype<-"endo_capillary"
dat_sub@meta.data[dat_sub@meta.data$stromal_endo_subclusters %in% c("3"),]$fine_celltype<-"endo_TEC"
dat_sub@meta.data[dat_sub@meta.data$stromal_endo_subclusters %in% c("6"),]$fine_celltype<-"endo_lymphatic"
plt_dim<-DimPlot(dat_sub,group.by=c("stromal_endo_subclusters","fine_celltype"),label=TRUE)
ggsave(plt_dim,file="tenx_dcis.stromal.endo.fine_celltype.umap.pdf",width=8)

saveRDS(dat_sub, file = "tenx_dcis.stromal.endo.rds")
```


## Peri/VSMC

```R
dat_sub<-subset(obj,coarse_celltype %in% c("perivasc"))
dat_sub<-cluster_object(obj=dat_sub,dims=1:15,prefix="stromal_periVSMC")
saveRDS(dat_sub,file = "tenx_dcis.stromal.periVSMC.rds")

stromal_markers<-list()
stromal_markers[["peri_THY1"]]<-c('CD36', 'CYGB', 'THY1', 'HIGD1B', 'COL3A1', 'RGS5', 'IGFBP2', 'COL6A3', 'ADGRF5', 'ARHGDIB', 'FAM213A', 'COX4I2', 'FAM162B', 'LPL', 'GUCY1A2', 'MYO1B', 'TXNIP', 'CYBA', 'PLXDC1')
stromal_markers[["peri_CCL21"]]<-c('CTSC', 'STEAP4', 'C1S', 'GGT5', 'C1R', 'COL6A3', 'CFH', 'S100A10', 'FGF7', 'CCL21', 'EMP1', 'CFD', 'TMEM176B', 'CXCL12', 'CHRDL1', 'TMEM176A', 'CLSTN2', 'FHL2', 'CCL19', 'CD44')
stromal_markers[["peri_stress"]]<-c('HSPA6', 'CACYBP', 'HSPH1', 'CHORDC1', 'FKBP4', 'AHSA1', 'STIP1', 'MIR4435-2HG', 'ABL2', 'ZFAND2A', 'THY1', 'HSPA4', 'BAG3', 'TENT5A', 'MRPL18', 'TCP1', 'SERPINH1', 'CYTOR', 'RUNX1', 'DNAJA4')
stromal_markers[["VSMC_CREM"]]<-c('MYH11', 'ZNF331', 'CREM', 'NR4A2', 'CD9', 'CNN1', 'SYNM', 'PLEKHO1', 'SORBS2', 'CCDC107', 'FAM107B', 'ELL2', 'ATP1B3', 'NDUFA4', 'ZFHX3', 'C12orf75', 'BCAM', 'RERGL', 'EIF4A3')
stromal_markers[["VSMC_PLN"]]<-c('NET1', 'TXNIP', 'RERGL', 'CASQ2', 'BCAM', 'PLN', 'KLHL23', 'NDUFA4', 'MEF2C', 'RCAN2', 'NTRK2', 'LBH', 'PLAC9', 'HES1', 'WFDC1', 'TBX2', 'AC023157.3', 'TBX2-AS1', 'ACTA2', 'KCNAB1')
stromal_markers[["VSMC_stress"]]<-c('HSPA6', 'MYH11', 'HSPA1B', 'HSPH1', 'CRYAB', 'BAG3', 'CHORDC1', 'DNAJB1', 'SORBS2', 'GADD45G', 'HSPB8', 'NET1', 'MRPL18', 'PLN', 'SERPINH1', 'HSPA1A', 'BCAM', 'PHLDA2', 'ID2', 'DNAJB4')

dat_sub<-FindClusters(dat_sub,res=0.25)
dat_sub$stromal_peri_subclusters<-dat_sub$seurat_clusters
plt_dim<-DimPlot(dat_sub,group.by="stromal_peri_subclusters",label=TRUE)

plt_list<-lapply(1:length(stromal_markers),function(x) {
    plt<-DotPlot(dat_sub, 
        features = stromal_markers[x], 
        group.by = "stromal_peri_subclusters",
        cluster.idents=FALSE) + 
    RotatedAxis() + ggtitle(names(stromal_markers)[x])
    return(plt)})

plt<-wrap_plots(plt_list,nrow=1,guides="collect")

layout <- "
A######
BBBBBBB"

ggsave(plt_dim+plt+plot_layout(design=layout),file="tenx_dcis.stromal.peri.fine_celltype.dotplot.pdf",width=0.25*length(unlist(stromal_markers)),height=10,limitsize=FALSE)

#labelling endo subtypes
dat_sub@meta.data$fine_celltype<-"suspected_doublet"
dat_sub@meta.data[dat_sub@meta.data$stromal_peri_subclusters %in% c("1","2"),]$fine_celltype<-"peri"
dat_sub@meta.data[dat_sub@meta.data$stromal_peri_subclusters %in% c("0","4"),]$fine_celltype<-"VSMC"
dat_sub@meta.data[dat_sub@meta.data$stromal_peri_subclusters %in% c("3"),]$fine_celltype<-"periVSMC_unknown"

plt_dim<-DimPlot(dat_sub,group.by=c("stromal_peri_subclusters","fine_celltype"),label=TRUE)
ggsave(plt_dim,file="tenx_dcis.stromal.peri.fine_celltype.umap.pdf",width=8)

saveRDS(dat_sub, file = "tenx_dcis.stromal.peri.rds")
```

# Immune celltyping

## Myeloid

```R
obj<-readRDS(file = "tenx_dcis.rds")

dat_sub<-subset(obj,coarse_celltype %in% c("myeloid"))
dat_sub<-cluster_object(obj=dat_sub,dims=1:15,prefix="immune_myeloid")

immune_markers<-list()
immune_markers[["myeloid_mono_classical"]]<-c('SERPINB2', 'AQP9', 'FCN1', 'VCAN', 'AC245128.3', 'EREG', 'TNIP3', 'CD300E', 'S100A8', 'S100A9', 'THBS1', 'CYP1B1', 'SLC39A8', 'MAP3K20', 'EHD1', 'LUCAT1', 'FCAR', 'S100A12', 'PID1', 'ACOD1')
immune_markers[["myeloid_mono_nonclassical"]]<-c('CCDC68', 'CLEC12A', 'CD48', 'OLIG1', 'CDKN1C', 'S100A4', 'KRTAP2-3', 'CD52', 'WARS', 'LINC02432', 'CFP', 'UNC45B', 'APOBEC3A', 'TKT', 'RIPOR2', 'CDA', 'PLAC8', 'FFAR2', 'PAPSS2', 'GK5')
immune_markers[["myeloid_macro_C3"]]<-c('C3', 'GPR158', 'C1QC', 'C1QB', 'AC079760.2', 'FCGR3A', 'SDS', 'A2M', 'AL590705.1', 'MT3', 'TREM2', 'HLA-DQB2', 'OLR1', 'ZBED9', 'IL4I1', 'DOK5', 'IGLC7', 'PMEPA1', 'IBSP', 'CD9')
immune_markers[["myeloid_macro_LYVE1"]]<-c('PDK4', 'SELENOP', 'FOLR2', 'LYVE1', 'LILRB5', 'MAF', 'KLF4', 'RNASE1', 'SLC40A1', 'RHOB', 'F13A1', 'MRC1', 'TSC22D3', 'PLTP', 'STAB1', 'NR4A2', 'MS4A4A', 'DAB2', 'CD163', 'CTSC')
immune_markers[["myeloid_macro_CCL4"]]<-c('CCL4', 'SELENOP', 'FOLR2', 'RNASE1', 'MRC1', 'F13A1', 'LYVE1', 'CD163', 'STAB1', 'CXCL3', 'HMOX1', 'CCL3', 'CXCL2', 'AC010980.2', 'GCLM', 'CD93', 'CTSL', 'PLTP', 'MAFB', 'THBD')
immune_markers[["myeloid_macro_APOC1"]]<-c('ACP5', 'CYP27A1', 'LIPA', 'PLD3', 'TREM2', 'PLA2G7', 'CD36', 'FABP4', 'GPNMB', 'CCL18', 'NCEH1', 'RARRES1', 'FUCA1', 'APOC1', 'LPL', 'OTOA', 'BLVRB', 'HS3ST2', 'NR1H3', 'CHIT1')
immune_markers[["myeloid_macro_interferon"]]<-c( 'EPSTI1', 'CXCL10', 'GBP1', 'STAT1', 'CXCL11', 'CXCL9', 'GBP4', 'IFI44L', 'FCGR1A', 'SERPING1', 'IFI6', 'IGSF6', 'IFIH1', 'IDO1', 'RASSF4', 'XAF1', 'APOL6', 'ISG15', 'FYB1', 'PSME2')
immune_markers[["myeloid_macro_stress"]]<-c( 'HSPA1B', 'BAG3', 'DNAJB1', 'HSPA6', 'MRPL18', 'ZFAND2A', 'CHORDC1', 'SDS', 'CACYBP', 'FKBP4', 'DNAJA4', 'HSPE1', 'RGS1', 'HSPH1', 'HSPA1A', 'OLR1', 'HSPB1', 'CLEC2B', 'CKS2', 'PLEK')
immune_markers[["myeloid_macro_CTSK"]]<-c( 'AP000904.1', 'CTSK', 'ACP5', 'SLC9B2', 'MMP9', 'DGKI', 'SPP1', 'TCIRG1', 'ATP6V0D2', 'CKB', 'CST3', 'CD109', 'TIMP2', 'JDP2', 'SPARC', 'C9orf16', 'GRN', 'SIGLEC15', 'AK5', 'SNX10')
immune_markers[["myeloid_cycling"]]<-c( 'MKI67', 'PCLAF', 'CKS1B', 'DIAPH3', 'TYMS', 'TOP2A', 'HIST1H1B', 'ASPM', 'TK1', 'DHFR', 'NUSAP1', 'NCAPG', 'DLGAP5', 'CENPM', 'TPX2', 'RRM2', 'CEP55', 'MYBL2', 'BIRC5', 'CENPF')
immune_markers[["myeloid_neutrophil"]]<-c('CRISP3', 'ANKRD34B', 'ABCA13', 'CD177', 'ORM1', 'SERPINB4', 'DEFA3', 'DPYS', 'C7', 'AZU1', 'CLC', 'FOXQ1', 'RAB6B', 'SFRP1', 'TWIST1', 'MAP1B', 'SLPI', 'ADH1B', 'FOLR3', 'BTNL9')
immune_markers[["myeloid_mast"]]<-c('CPA3', 'TPSAB1', 'MS4A2', 'TPSB2', 'HDC', 'IL1RL1', 'GATA2', 'HPGDS', 'GCSAML', 'HPGD', 'ADCYAP1', 'KIT', 'KRT1', 'PRG2', 'CTSG', 'CLU', 'TPSD1', 'CALB2', 'RAB27B', 'SLC18A2')
immune_markers[["myeloid_cDC1"]]<-c('CLEC9A', 'XCR1', 'DNASE1L3', 'IDO1', 'LGALS2', 'CPNE3', 'WDFY4', 'C1orf54', 'DAPP1', 'RAB11FIP1', 'PPA1', 'CPVL', 'C9orf135', 'FOXB1', 'AC016717.2', 'CST3', 'PLEKHA5', 'GPR157', 'FOXN2', 'NCOA7')
immune_markers[["myeloid_cDC2"]]<-c('FCER1A', 'IL1R2', 'CLEC10A', 'CD1C', 'CST7', 'GPAT3', 'DAPP1', 'CFP', 'EREG', 'CCL22', 'LGALS2', 'JAML', 'PID1', 'AREG', 'IL7R', 'AC020656.1', 'IL1R1', 'MARCKSL1', 'ADAM19', 'SLC7A11')
immune_markers[["myeloid_mDC"]]<-c('LAMP3', 'BIRC3', 'LACRT', 'NUB1', 'CCR7', 'MARCKSL1', 'IDO1', 'DAPP1', 'POGLUT1', 'LINC01539', 'GPR157', 'IL12B', 'LAD1', 'KIF2A', 'FSCN1', 'IL7R', 'TXN', 'DUSP5', 'FOXD4L1', 'CD200', 'RAB9A')
immune_markers[["myeloid_pDC"]]<-c('CD5', 'LYPD8', 'PRL', 'AC136475.3', 'LINC01087', 'AC006058.1', 'BDKRB2', 'POTEI', 'DSP', 'AC026369.3', 'GZMB', 'TNFSF4', 'CD2', 'TIGIT', 'LTB', 'TMEM45A', 'PGR', 'CD3G', 'AC015936.1', 'ACOT7')
immune_markers[["myeloid_TAM"]]<-c('SDS', 'LAIR1', 'FCGR3A', 'SH3PXD2B', 'C1QB', 'FGL2', 'SGPL1', 'ADA2', 'AXL', 'TTYH3', 'TREM2', 'AOAH', 'ACP5', 'RAB20', 'SLC16A10', 'SATB1', 'FPR3', 'HLA-DOA', 'OLFML2B', 'CCDC107', 'MMP9', 'CALHM6', 'PLA2G7', 'GNA13', 'ARL4C', 'ZNF331', 'JMY', 'C2', 'A2M', 'STAT1')

#TAM list from Message Aatish Thennavan

dat_sub<-FindClusters(dat_sub,res=0.4)
dat_sub$immune_myeloid_subclusters<-dat_sub$seurat_clusters
plt_dim<-DimPlot(dat_sub,group.by="immune_myeloid_subclusters",label=TRUE)

library(Seurat)
plt_list<-lapply(1:length(immune_markers),function(x) {
    plt<-DotPlot(dat_sub, 
        features = immune_markers[x], 
        group.by = "immune_myeloid_subclusters",
        cluster.idents=FALSE) + 
    RotatedAxis() + ggtitle(names(immune_markers)[x])
    return(plt)})

plt<-wrap_plots(plt_list,nrow=1,guides="collect")

layout <- "
A######
BBBBBBB"

detach("package:Seurat",unload=TRUE) #seurat and patchwork have some incompatibility
library(patchwork)
ggsave(plt_dim+plt+patchwork::plot_layout(design=layout),file="tenx_dcis.immune.myeloid.fine_celltype.dotplot.pdf",width=0.25*length(unlist(immune_markers)),height=10,limitsize=FALSE) 


#labelling myeloid subtypes
dat_sub@meta.data$fine_celltype<-"suspected_doublet"
dat_sub@meta.data[dat_sub@meta.data$immune_myeloid_subclusters %in% c("4"),]$fine_celltype<-"myeloid_mono"

dat_sub@meta.data[dat_sub@meta.data$immune_myeloid_subclusters %in% c("2","1","10"),]$fine_celltype<-"myeloid_macro" #APOC INTERFERON
dat_sub@meta.data[dat_sub@meta.data$immune_myeloid_subclusters %in% c("7","6"),]$fine_celltype<-"myeloid_TAM" #LYVE1 CCL4

dat_sub@meta.data[dat_sub@meta.data$immune_myeloid_subclusters %in% c("9"),]$fine_celltype<-"myeloid_cycling"
dat_sub@meta.data[dat_sub@meta.data$immune_myeloid_subclusters %in% c("11"),]$fine_celltype<-"myeloid_neutrophil"
dat_sub@meta.data[dat_sub@meta.data$immune_myeloid_subclusters %in% c("8"),]$fine_celltype<-"myeloid_mast"
dat_sub@meta.data[dat_sub@meta.data$immune_myeloid_subclusters %in% c("0"),]$fine_celltype<-"myeloid_DC"

plt_dim<-DimPlot(dat_sub,group.by=c("immune_myeloid_subclusters","fine_celltype"),label=TRUE)
ggsave(plt_dim,file="tenx_dcis.immune.myeloid.fine_celltype.umap.pdf",width=15)


saveRDS(dat_sub,file = "tenx_dcis.immune.myeloid.rds")

```

```R
#nondiegetic code update, adding TAMs (based mostly on TREM2) into PF seurat object
rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")
rna@meta.data[row.names(dat_sub@meta.data),]$fine_celltype<-dat_sub$fine_celltype
rna$myeloid_clus<-NA
rna@meta.data[row.names(dat_sub@meta.data),]$myeloid_clus<-dat_sub$immune_myeloid_subclusters

dat_sub<-AddMetaData(dat_sub,rna@meta.data)
saveRDS(rna,"tenx_dcis.pf.rds")

plt_dim<-DimPlot(dat_sub,group.by=c("immune_myeloid_subclusters","fine_celltype","Group"),label=TRUE)
ggsave(plt_dim,file="tenx_dcis.immune.myeloid.fine_celltype.umap.pdf",width=15)


```

## T cells

```R
dat_sub<-subset(obj,coarse_celltype %in% c("tcell"))
dat_sub<-cluster_object(obj=dat_sub,dims=1:15,prefix="immune_tcell")

immune_markers<-list()
immune_markers[["tcell_cd4_memory"]]<-c('ADAM23', 'NEFL', 'LINC02273', 'ANTXR2', 'MFHAS1', 'AP3M2', 'GPR183', 'PASK', 'S1PR1', 'FTH1', 'GLIPR1', 'SESN3', 'IL7R', 'CCR7', 'SLC2A3', 'DDIT4', 'CD28', 'ANXA1', 'TRAT1', 'KLF2')
immune_markers[["tcell_cd4_stress"]]<-c('MYADM', 'LMNA', 'ANXA1', 'KLF6', 'HSP90AA1', 'RGCC', 'FAM107B', 'AHNAK', 'HSPH1', 'VIM', 'HSPA8', 'GPR183', 'HSP90AB1', 'S100A10', 'S100A11', 'EZR', 'HSPD1', 'IL7R', 'ANKRD12', 'TUBB4B')
immune_markers[["tcell_cd4_naive"]]<-c('CA6', 'ADTRP', 'GTSCR1', 'LINC00402', 'LINC01481', 'LEF1', 'MAL', 'SELL', 'AL391097.1', 'ANKRD55', 'NEXMIF', 'CCR7', 'LINC00861', 'TESPA1', 'ANK1', 'LTB', 'SESN3', 'TSHZ2', 'PASK', 'LDHB')
immune_markers[["tcell_treg"]]<-c('CD177', 'FANK1', 'LINC02099', 'IL1R2', 'CCR8', 'FOXP3', 'PTGIR', 'LAYN', 'IKZF4', 'F5', 'RTKN2', 'CARD16', 'MAGEH1', 'IL2RA', 'LINC01943', 'BATF', 'IKZF2', 'HTATIP2', 'GK', 'CTLA4')
immune_markers[["tcell_cd4_Tfh"]]<-c('IGFL2', 'CXCL13', 'CPM', 'CD200', 'KSR2', 'NMB', 'HMSD', 'TUBA3D', 'PTPN13', 'GNG4', 'SERPINE2', 'AC008011.2', 'HMHB1', 'LINC01229', 'TSHZ2', 'GK', 'PGM2L1', 'SESN3', 'MFHAS1', 'FKBP5')
immune_markers[["tcell_cd8_TEM"]]<-c('DKK3', 'ENC1', 'GZMK', 'EOMES', 'YBX3', 'GGA2', 'CMC1', 'CD8A', 'CST7', 'CD8B', 'SH2D1A', 'HLA-DPB1', 'CRTAM', 'DUSP2', 'LYST', 'COTL1', 'ZEB2', 'PIK3R1', 'RNF19A', 'TRAT1')
immune_markers[["tcell_cd8_TRMZNF683"]]<-c('CSN2', 'AL590434.1', 'CCL4', 'OR5AU1', 'FGF10-AS1', 'CCL4L2', 'TNFSF9', 'AL133163.2', 'DUSP6', 'FOSB', 'AC020916.1', 'IER5L', 'MZF1-AS1', 'NR4A2', 'IFNG', 'RASGEF1B', 'PLCH1', 'POLQ', 'GZMA', 'GZMK')
immune_markers[["tcell_cd8_TRMAUTS"]]<-c('AUTS2', 'PRSS3', 'GFPT2', 'DAPK2', 'PAGE5', 'ITGA1', 'KLRC1', 'IGFBP3', 'AC020571.1', 'HPSE2', 'CTH', 'COL25A1', 'LINC01871', 'FGF9', 'PARD6G', 'SPRY1', 'SNTG2', 'CD8A', 'ADAMTS6', 'TAC1')
immune_markers[["tcell_cd8_stress"]]<-c( 'HSPA1A', 'DNAJB1', 'HSPH1', 'HSP90AA1', 'HSPB1', 'HSPE1', 'DNAJA1', 'TUBA4A', 'HSP90AB1', 'HSPD1', 'CDKN1C', 'B3GNT5', 'CSF1', 'HSPA8', 'KLF6', 'EGR1', 'POU3F1', 'ANKRD37', 'ZNF80', 'MIR222HG')
immune_markers[["tcell_gdT"]]<-c('KIR2DL3', 'KLRC2', 'KIR2DL1', 'KIR3DL1', 'CHPT1', 'SMC4', 'TRGC2', 'KLRD1', 'XCL1', 'CEMIP2', 'CD7', 'METRNL', 'IKZF2', 'AHI1', 'CTSW', 'KIR2DL4', 'CAST', 'LINC01871', 'ID3', 'PRKX')
immune_markers[["tcell_cycling"]]<-c('DLGAP5', 'E2F8', 'HJURP', 'SPC25', 'RRM2', 'GTSE1', 'AC007240.1', 'BIRC5', 'MYBL2', 'MKI67', 'MCM10', 'KIF4A', 'CENPA', 'SPC24', 'CDC45', 'PCLAF', 'ESCO2', 'KIF20A', 'UBE2C', 'DIAPH3')
immune_markers[["tcell_cd8_interferon"]]<-c( 'IFIT3', 'IFIT1', 'IFI44L', 'HERC5', 'MX2', 'MX1', 'IFI6', 'XAF1', 'ISG15', 'STAT1', 'EIF2AK2', 'LY6E', 'EPSTI1', 'IFI44', 'PARP14', 'OASL', 'OAS2', 'RSAD2', 'SAMD9L', 'MT2A')
immune_markers[["tcell_cd8_exhausted"]]<-c( 'REG4', 'RYR2', 'HAVCR2', 'NPW', 'VCAM1', 'ETV1', 'PNOC', 'KRT86', 'TNS3', 'MYO7A', 'PHEX', 'LRRN3', 'AC243829.4', 'FNDC9', 'BUB1B', 'BEND4', 'AKAP5', 'ASB2', 'NHS', 'RGS13', 'ENTPD1' )
immune_markers[["tcell_nk_FCGR3A"]]<-c('MGAM', 'MYOM2', 'FCGR3A', 'KLRF1', 'PRSS57', 'MLC1', 'AKR1C3', 'FGFBP2', 'SPON2', 'RBPMS2', 'FCER1G', 'KIR3DL1', 'CX3CR1', 'KIR2DL1', 'SH2D1B', 'S1PR5', 'CHST2', 'CCL3', 'CFAP97D2', 'KIR2DL3')
immune_markers[["tcell_nk_cytotoxic"]]<-c('LILRB1', 'CES1', 'PROK2', 'FGFBP2', 'ASCL2', 'CX3CR1', 'FCRL6', 'LINC02384', 'GZMH', 'CFAP97D2', 'LINC01936', 'FCGR3A', 'S1PR5', 'AL590434.1', 'PLEK', 'NKG7', 'ADGRG1', 'AC011933.3', 'PRSS23', 'KIR2DL1')
immune_markers[["tcell_ilc1"]]<-c('ADGRG3', 'CCNJL', 'LINC00996', 'KIR2DL4', 'CSF2', 'B3GNT7', 'FAM166B', 'FES', 'SCT', 'SH2D1B', 'CFAP46', 'AL157895.1', 'KRT86', 'ATP8B4', 'FCER1G', 'GNLY', 'ZNF660', 'TYROBP', 'GOLIM4', 'HOXA6')
immune_markers[["tcell_ilc23"]]<-c('CD300LF', 'KIT', 'SOST', 'LINC00299', 'HOXA6', 'PLA2G4A', 'NPTX2', 'VSTM2L', 'FAM30A', 'SCT', 'XCL1', 'SOX4', 'XCL2', 'AREG', 'TNFRSF18', 'SCN1B', 'TNFSF4', 'MAP3K8', 'SSBP2', 'PLCH1')
immune_markers[["tcell_MAIT"]]<-c('KLRB1', 'RORA', 'KISS1', 'CEBPD', 'TRDV2', 'AC022392.1', 'AC026461.3', 'FEZ1', 'REL', 'GPR65', 'MYBL1', 'LTB', 'MT3', 'SLC4A10', 'NFKB1', 'IGHV3-23', 'CCL20', 'RORC', 'NCR3', 'DPP4')
immune_markers[["tcell_GZMK"]]<-c('ZMAT4', 'KRT13', 'COL24A1', 'CD69', 'COX6A2', 'AL157895.1', 'CCNJL', 'SOST', 'MLC1', 'ADGRG3', 'SMIM24', 'HOXA6', 'SH2D1B', 'TYROBP', 'XCL2', 'GNLY', 'TAL1', 'XCL1', 'AREG', 'FCER1G')


dat_sub<-FindClusters(dat_sub,res=0.25)
dat_sub$immune_tcell_subclusters<-dat_sub$seurat_clusters
plt_dim<-DimPlot(dat_sub,group.by="immune_tcell_subclusters",label=TRUE)

plt_list<-lapply(1:length(immune_markers),function(x) {
    plt<-DotPlot(dat_sub, 
        features = immune_markers[x], 
        group.by = "immune_tcell_subclusters",
        cluster.idents=FALSE) + 
    RotatedAxis() + ggtitle(names(immune_markers)[x])
    return(plt)})

plt<-wrap_plots(plt_list,nrow=1,guides="collect")

layout <- "
A######
BBBBBBB"

ggsave(plt_dim+plt+plot_layout(design=layout),file="tenx_dcis.immune.tcell.fine_celltype.dotplot.pdf",width=0.25*length(unlist(immune_markers)),height=10,limitsize=FALSE)

#labelling tcell subtypes
dat_sub@meta.data$fine_celltype<-"suspected_doublet"
dat_sub@meta.data[dat_sub@meta.data$immune_tcell_subclusters %in% c("9"),]$fine_celltype<-"tcell_nk"
dat_sub@meta.data[dat_sub@meta.data$immune_tcell_subclusters %in% c("5","4","3","1"),]$fine_celltype<-"tcell_cd8"
dat_sub@meta.data[dat_sub@meta.data$immune_tcell_subclusters %in% c("8"),]$fine_celltype<-"tcell_cycling"
dat_sub@meta.data[dat_sub@meta.data$immune_tcell_subclusters %in% c("4"),]$fine_celltype<-"tcell_interferon"
dat_sub@meta.data[dat_sub@meta.data$immune_tcell_subclusters %in% c("6"),]$fine_celltype<-"tcell_gdT"
dat_sub@meta.data[dat_sub@meta.data$immune_tcell_subclusters %in% c("2"),]$fine_celltype<-"tcell_treg"
dat_sub@meta.data[dat_sub@meta.data$immune_tcell_subclusters %in% c("0","8","7"),]$fine_celltype<-"tcell_cd4"

plt_dim<-DimPlot(dat_sub,group.by=c("immune_tcell_subclusters","fine_celltype"),label=TRUE)
ggsave(plt_dim,file="tenx_dcis.immune.tcell.fine_celltype.umap.pdf",width=15)

saveRDS(dat_sub,file = "tenx_dcis.immune.tcell.rds")

```

## B cell

```R
dat_sub<-subset(obj,coarse_celltype %in% c("bcell"))
#have to remove samples from layers list with no b cells
layersList <- lapply(dat_sub@assays$RNA@layers,function(x){dim(x)})
dat_sub@assays$RNA@layers[names(layersList[sapply(layersList, is.null)])] <- NULL

dat_sub<-cluster_object(obj=dat_sub,dims=1:15,prefix="immune_bcell")

immune_markers<-list()
immune_markers[["bcell_naive"]]<-c('IGHD', 'YBX3', 'TCL1A', 'SPRY1', 'CD83', 'BCL7A', 'TXNIP', 'FCER2', 'CLEC2D', 'HLA-DQA2', 'MYB', 'KLF2', 'ESCO2', 'STAG3', 'HLA-DRB5', 'GPR18', 'CMSS1', 'ACSM3', 'GBP4')
immune_markers[["bcell_memory"]]<-c('REL', 'LTB', 'CD69', 'GPR183', 'CCR7', 'CD48', 'ID3', 'BCL2A1', 'KYNU',  'PIKFYVE', 'SAMSN1', 'IER5', 'CD80', 'CHORDC1', 'SCIMP', 'DUSP2', 'KLF6', 'NFKB1')
immune_markers[["bcell_interferon"]]<-c('CCL17', 'STAT1', 'IFI44L',  'ISG15', 'SERPINA9', 'MX1', 'AICDA', 'DUSP4', 'XAF1', 'GBP1', 'ZBED2', 'GBP5', 'GBP2', 'AIM2', 'SAMD9L', 'IFNG', 'MX2', 'ITGAX', 'PLIN2')
immune_markers[["bcell_stress"]]<-c('HSPA1B', 'RRBP1', 'DNAJB1', 'HSPA1A', 'FOS', 'ANKRD28', 'SAMD7', 'H1FX', 'IGKV1-37', 'RGS1', 'HSPB1', 'AC026369.3', 'AMPD1', 'IGKV2-26', 'FAM92B', 'AC021074.3', 'IGHV3-41', 'LINC01405', 'CPEB4', 'IGKV2-4')

immune_markers[["myeloid_mono_classical"]]<-c('SERPINB2', 'AQP9', 'FCN1', 'VCAN', 'AC245128.3', 'EREG', 'TNIP3', 'CD300E', 'S100A8', 'S100A9', 'THBS1', 'CYP1B1', 'SLC39A8', 'MAP3K20', 'EHD1', 'LUCAT1', 'FCAR', 'S100A12', 'PID1', 'ACOD1')
immune_markers[["myeloid_mono_nonclassical"]]<-c('CCDC68', 'CLEC12A', 'CD48', 'OLIG1', 'CDKN1C', 'S100A4', 'KRTAP2-3', 'CD52', 'WARS', 'LINC02432', 'CFP', 'UNC45B', 'APOBEC3A', 'TKT', 'RIPOR2', 'CDA', 'PLAC8', 'FFAR2', 'PAPSS2', 'GK5')
immune_markers[["myeloid_macro_C3"]]<-c('C3', 'GPR158', 'C1QC', 'C1QB', 'AC079760.2', 'FCGR3A', 'SDS', 'A2M', 'AL590705.1', 'MT3', 'TREM2', 'HLA-DQB2', 'OLR1', 'ZBED9', 'IL4I1', 'DOK5', 'IGLC7', 'PMEPA1', 'IBSP', 'CD9')
immune_markers[["myeloid_macro_LYVE1"]]<-c('PDK4', 'SELENOP', 'FOLR2', 'LYVE1', 'LILRB5', 'MAF', 'KLF4', 'RNASE1', 'SLC40A1', 'RHOB', 'F13A1', 'MRC1', 'TSC22D3', 'PLTP', 'STAB1', 'NR4A2', 'MS4A4A', 'DAB2', 'CD163', 'CTSC')
immune_markers[["myeloid_macro_CCL4"]]<-c('CCL4', 'SELENOP', 'FOLR2', 'RNASE1', 'MRC1', 'F13A1', 'LYVE1', 'CD163', 'STAB1', 'CXCL3', 'HMOX1', 'CCL3', 'CXCL2', 'AC010980.2', 'GCLM', 'CD93', 'CTSL', 'PLTP', 'MAFB', 'THBD')
immune_markers[["myeloid_macro_APOC1"]]<-c('ACP5', 'CYP27A1', 'LIPA', 'PLD3', 'TREM2', 'PLA2G7', 'CD36', 'FABP4', 'GPNMB', 'CCL18', 'NCEH1', 'RARRES1', 'FUCA1', 'APOC1', 'LPL', 'OTOA', 'BLVRB', 'HS3ST2', 'NR1H3', 'CHIT1')
immune_markers[["myeloid_macro_interferon"]]<-c( 'EPSTI1', 'CXCL10', 'GBP1', 'STAT1', 'CXCL11', 'CXCL9', 'GBP4', 'IFI44L', 'FCGR1A', 'SERPING1', 'IFI6', 'IGSF6', 'IFIH1', 'IDO1', 'RASSF4', 'XAF1', 'APOL6', 'ISG15', 'FYB1', 'PSME2')
immune_markers[["myeloid_macro_stress"]]<-c( 'HSPA1B', 'BAG3', 'DNAJB1', 'HSPA6', 'MRPL18', 'ZFAND2A', 'CHORDC1', 'SDS', 'CACYBP', 'FKBP4', 'DNAJA4', 'HSPE1', 'RGS1', 'HSPH1', 'HSPA1A', 'OLR1', 'HSPB1', 'CLEC2B', 'CKS2', 'PLEK')
immune_markers[["myeloid_macro_CTSK"]]<-c( 'AP000904.1', 'CTSK', 'ACP5', 'SLC9B2', 'MMP9', 'DGKI', 'SPP1', 'TCIRG1', 'ATP6V0D2', 'CKB', 'CST3', 'CD109', 'TIMP2', 'JDP2', 'SPARC', 'C9orf16', 'GRN', 'SIGLEC15', 'AK5', 'SNX10')
immune_markers[["myeloid_cycling"]]<-c( 'MKI67', 'PCLAF', 'CKS1B', 'DIAPH3', 'TYMS', 'TOP2A', 'HIST1H1B', 'ASPM', 'TK1', 'DHFR', 'NUSAP1', 'NCAPG', 'DLGAP5', 'CENPM', 'TPX2', 'RRM2', 'CEP55', 'MYBL2', 'BIRC5', 'CENPF')
immune_markers[["myeloid_neutrophil"]]<-c('CRISP3', 'ANKRD34B', 'ABCA13', 'CD177', 'ORM1', 'SERPINB4', 'DEFA3', 'DPYS', 'C7', 'AZU1', 'CLC', 'FOXQ1', 'RAB6B', 'SFRP1', 'TWIST1', 'MAP1B', 'SLPI', 'ADH1B', 'FOLR3', 'BTNL9')
immune_markers[["myeloid_mast"]]<-c('CPA3', 'TPSAB1', 'MS4A2', 'TPSB2', 'HDC', 'IL1RL1', 'GATA2', 'HPGDS', 'GCSAML', 'HPGD', 'ADCYAP1', 'KIT', 'KRT1', 'PRG2', 'CTSG', 'CLU', 'TPSD1', 'CALB2', 'RAB27B', 'SLC18A2')
immune_markers[["myeloid_cDC1"]]<-c('CLEC9A', 'XCR1', 'DNASE1L3', 'IDO1', 'LGALS2', 'CPNE3', 'WDFY4', 'C1orf54', 'DAPP1', 'RAB11FIP1', 'PPA1', 'CPVL', 'C9orf135', 'FOXB1', 'AC016717.2', 'CST3', 'PLEKHA5', 'GPR157', 'FOXN2', 'NCOA7')
immune_markers[["myeloid_cDC2"]]<-c('FCER1A', 'IL1R2', 'CLEC10A', 'CD1C', 'CST7', 'GPAT3', 'DAPP1', 'CFP', 'EREG', 'CCL22', 'LGALS2', 'JAML', 'PID1', 'AREG', 'IL7R', 'AC020656.1', 'IL1R1', 'MARCKSL1', 'ADAM19', 'SLC7A11')
immune_markers[["myeloid_mDC"]]<-c('LAMP3', 'BIRC3', 'LACRT', 'NUB1', 'CCR7', 'MARCKSL1', 'IDO1', 'DAPP1', 'POGLUT1', 'LINC01539', 'GPR157', 'IL12B', 'LAD1', 'KIF2A', 'FSCN1', 'IL7R', 'TXN', 'DUSP5', 'FOXD4L1', 'CD200', 'RAB9A')
immune_markers[["myeloid_pDC"]]<-c('CD5', 'LYPD8', 'PRL', 'AC136475.3', 'LINC01087', 'AC006058.1', 'BDKRB2', 'POTEI', 'DSP', 'AC026369.3', 'GZMB', 'TNFSF4', 'CD2', 'TIGIT', 'LTB', 'TMEM45A', 'PGR', 'CD3G', 'AC015936.1', 'ACOT7')
#added myeloid because they cluster near eachother and there was an unaccounted for b cell cluster

dat_sub<-FindClusters(dat_sub,res=0.25)
dat_sub$immune_bcell_subclusters<-dat_sub$seurat_clusters
plt_dim<-DimPlot(dat_sub,group.by="immune_bcell_subclusters",label=TRUE)
table(dat_sub$immune_bcell_subclusters,useNA="ifany")
dat_sub<-subset(dat_sub,cells=Cells(dat_sub)[!is.na(dat_sub$immune_bcell_subclusters)])
plt_list<-lapply(1:length(immune_markers),function(x) {
    plt<-DotPlot(dat_sub, 
        features = immune_markers[x], 
        group.by = "immune_bcell_subclusters",
        cluster.idents=FALSE) + 
    RotatedAxis() + ggtitle(names(immune_markers)[x])
    return(plt)})

plt<-wrap_plots(plt_list,nrow=1,guides="collect")

layout <- "
A######
BBBBBBB"

ggsave(plt_dim+plt+plot_layout(design=layout),file="tenx_dcis.immune.bcell.fine_celltype.dotplot.pdf",width=0.25*length(unlist(immune_markers)),height=10,limitsize=FALSE)

#labelling bcell subtypes
dat_sub@meta.data$fine_celltype<-"bcell"

#labelling tcell subtypes
dat_sub@meta.data$fine_celltype<-"suspected_doublet"
dat_sub@meta.data[dat_sub@meta.data$immune_bcell_subclusters %in% c("0","1"),]$fine_celltype<-"bcell"
dat_sub@meta.data[dat_sub@meta.data$immune_bcell_subclusters %in% c("2","3"),]$fine_celltype<-"plasma"

plt_dim<-DimPlot(dat_sub,group.by=c("immune_bcell_subclusters","fine_celltype"),label=TRUE)
ggsave(plt_dim,file="tenx_dcis.immune.bcell.fine_celltype.umap.pdf",width=8)

saveRDS(dat_sub,file = "tenx_dcis.immune.bcell.rds")
```

# Plasma cells

```R
dat_sub<-subset(obj,coarse_celltype %in% c("plasma"))
layersList <- lapply(dat_sub@assays$RNA@layers,function(x){dim(x)})
dat_sub@assays$RNA@layers[names(layersList[sapply(layersList, is.null)])] <- NULL
dat_sub<-cluster_object(obj=dat_sub,dims=1:15,prefix="immune_plasma")

immune_markers<-list()
immune_markers[["plasma_IGA"]]<-c('JCHAIN', 'IGHA2', 'IGHA1', 'LGALSL', 'NCAM1', 'CPVL', 'CCR10', 'AC106028.4', 'ZNF215', 'IGHV7-81', 'IGKV2-28', 'CADPS2', 'LHX8', 'IGLV2-18', 'IGHV3OR16-8', 'OTP', 'CTSW', 'AC136428.1', 'IGHV3-72')
immune_markers[["plasma_IGG"]]<-c('MZB1', 'IGHG3', 'DERL3', 'FKBP11', 'IGHG4', 'ITM2C', 'IGHG1', 'PRDX4', 'SEC11C', 'XBP1', 'SSR4', 'SSR3', 'JSRP1', 'ERLEC1', 'CYTOR', 'IGKV3D-11', 'SELENOS', 'DPEP1', 'CD38', 'MYDGF')

dat_sub<-FindClusters(dat_sub,res=0.25)
dat_sub$immune_plasma_subclusters<-dat_sub$seurat_clusters
plt_dim<-DimPlot(dat_sub,group.by="immune_plasma_subclusters",label=TRUE)

table(dat_sub$immune_plasma_subclusters,useNA="ifany")
#dat_sub<-subset(dat_sub,cells=Cells(dat_sub)[!is.na(dat_sub$immune_plasma_subclusters)])


plt_list<-lapply(1:length(immune_markers),function(x) {
    plt<-DotPlot(dat_sub, 
        features = immune_markers[x], 
        group.by = "immune_plasma_subclusters",
        cluster.idents=FALSE) + 
    RotatedAxis() + ggtitle(names(immune_markers)[x])
    return(plt)})

plt<-wrap_plots(plt_list,nrow=1,guides="collect")

layout <- "
A######
BBBBBBB"

ggsave(plt_dim+plt+plot_layout(design=layout),file="tenx_dcis.immune.plasma.fine_celltype.dotplot.pdf",width=0.25*length(unlist(immune_markers)),height=10,limitsize=FALSE)

#labelling endo subtypes
dat_sub@meta.data$fine_celltype<-"plasma"
plt_dim<-DimPlot(dat_sub,group.by=c("immune_plasma_subclusters","fine_celltype"),label=TRUE)
ggsave(plt_dim,file="tenx_dcis.immune.plasma.fine_celltype.umap.pdf",width=8)

saveRDS(dat_sub,file = "tenx_dcis.immune.plasma.rds")

```


# Combine output of fine clustering back in main object

```R

dat<-readRDS("tenx_dcis.rds")

fine_celltype_list<-lapply(
    c("tenx_dcis.stromal.fibro.rds", "tenx_dcis.stromal.endo.rds","tenx_dcis.stromal.peri.rds","tenx_dcis.immune.myeloid.rds","tenx_dcis.immune.tcell.rds","tenx_dcis.immune.bcell.rds","tenx_dcis.immune.plasma.rds"),
    function(x){
        tmp<-readRDS(x)
        return(tmp$fine_celltype)
    }
)

fine_celltype<-unlist(fine_celltype_list)

dat<-AddMetaData(dat,col.name="fine_celltype",metadata=fine_celltype)
dat@meta.data[dat@meta.data$coarse_celltype=="lumhr",]$fine_celltype<-"lumhr"
dat@meta.data[dat@meta.data$coarse_celltype=="lumsec",]$fine_celltype<-"lumsec"
dat@meta.data[dat@meta.data$coarse_celltype=="basal",]$fine_celltype<-"basal"

saveRDS(dat,"tenx_dcis.rds")

#go through and label doublet clusters more
#obj<-cluster_object(obj=dat,dims=1:15,res=0.8,prefix="coarse_celltype")
#obj@meta.data[obj@meta.data$seurat_clusters %in% c("33","30","35","21","25","27","28","24","26"),]$broad_celltype<-"suspected_doublet"
#obj@meta.data[obj@meta.data$seurat_clusters %in% c("33","30","35","21","25","27","28","24","26"),]$fine_celltype<-"suspected_doublet"

plt<-DimPlot(dat,group.by=c("seurat_clusters","coarse_celltype","fine_celltype"),label=T)
ggsave(plt,file="tenx_dcis.celltypes.final.pdf",width=40,height=20,limitsize=F)
saveRDS(dat,"tenx_dcis.rds")


#clean umap with no doublets
obj<-subset(dat,coarse_celltype != "suspected_doublet")
obj<-cluster_object(obj=obj,dims=1:15,prefix="fine_celltype")

plt<-DimPlot(obj,group.by=c("seurat_clusters","coarse_celltype","fine_celltype"),label=T)
ggsave(plt,file="tenx_dcis.celltypes.final.pf.pdf",width=40,height=20,limitsize=F)
saveRDS(obj,"tenx_dcis.pf.rds")


```


Additional final plots for manuscript

```R
library(Seurat)
library(hdf5r)
library(ggplot2)
library(patchwork)
library(amethyst)
library(dplyr)
setwd("/data/rmulqueen/projects/scalebio_dcis/rna")

obj<-readRDS("tenx_dcis.pf.rds")

#this metadata is apriori, but im just loading from amethyst object because i'm lazy, will just replace with hard code for timeline consistency
#nondiegetic file loading (from 04 step)
met<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/merged_data/07_scaledcis.integrated_celltyping.amethyst.rds")
meta<-met@metadata %>% group_by(Sample) %>% slice_head(n=1) %>% as.data.frame() %>% select(sample,Sample,Dx,Age,Race_Ethnicity,Menopause,BRCA,ER,PR,HER2,Grade,Group)
#read in sample metadata from amethyst methylation data
#add diagnosis meta

unique(obj@meta.data$sample) %in% meta$Sample
row.names(meta)<-meta$Sample
obj@meta.data$Dx<-meta[obj@meta.data$sample,]$Dx
obj@meta.data$Age<-meta[obj@meta.data$sample,]$Age
obj@meta.data$Race_Ethnicity<-meta[obj@meta.data$sample,]$Race_Ethnicity
obj@meta.data$Menopause<-meta[obj@meta.data$sample,]$Menopause
obj@meta.data$BRCA<-meta[obj@meta.data$sample,]$BRCA
obj@meta.data$ER<-meta[obj@meta.data$sample,]$ER
obj@meta.data$PR<-meta[obj@meta.data$sample,]$PR
obj@meta.data$HER2<-meta[obj@meta.data$sample,]$HER2
obj@meta.data$Grade<-meta[obj@meta.data$sample,]$Grade
obj@meta.data$Group<-meta[obj@meta.data$sample,]$Group

res=0.2
dims=1:20

obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- RunPCA(obj, features = VariableFeatures(object = obj))
obj <- FindNeighbors(obj, dims = dims)
obj <- RunUMAP(obj, dims = dims)

for(i in c("coarse_celltype","fine_celltype","Group")){
    plt_dim<-DimPlot(obj,group.by=i,label=TRUE,raster=FALSE)
    ggsave(plt_dim,file=paste0("tenx_dcis.pf.final.",i,".umap.pdf"),width=10)
    }


saveRDS(obj,"tenx_dcis.pf.rds")

table(rna$fine_celltype,rna$Group)
#confirms strong enrichment for CAFs TECs and TAMs in IDC/DCIS
```

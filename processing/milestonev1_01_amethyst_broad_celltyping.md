```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
```

```R
source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing

task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)

system(paste("mkdir -p",wd))
setwd(wd)

#make merged object
plate_obj=paste(sep="/",project_data_directory,"scalemethyl_pipeline_out/amethyst_plate_obj")
amethyst_files=list.files(path=plate_obj,pattern="*.amethyst.rds",recursive=TRUE,full.names=TRUE)

########################################
## Read in all amethyst files
########################################

window_name="cg_100k_score"
dat_list<-mclapply(amethyst_files, function(x) {
    obj<-readRDS(x)
    return(obj)},mc.cores=20)

dat <- combineObject(objList = dat_list, genomeMatrices=window_name)
```

########################################
## Run nakshatri atac peaks on all samples to cluster
########################################

```R
nakshatri_peaks<-readRDS("/data/rmulqueen/projects/scalebio_dcis/ref/celltype_peaks.500bp.rds")
nakshatri_bed<-data.frame(
        seqnames=seqnames(nakshatri_peaks),
        starts=start(nakshatri_peaks),
        ends=end(nakshatri_peaks))

nakshatri_bed<-nakshatri_bed[nakshatri_bed$seqnames %in% unique(dat@ref$seqid),]
dat@genomeMatrices[["nakshatri_peaks"]] <- makeWindows(dat, 
                                                     bed = nakshatri_bed,
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 300, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

saveRDS(dat,file="00_merged.amethyst.rds")

dat<-readRDS(dat,file="00_merged.amethyst.rds")
#subset to cell lines
dat_celllines<-subsetObject(dat,cells=row.names(dat@metadata)[dat@metadata$sample %in% c("MCF10A","MCF7","MDA-MB-231")])
saveRDS(dat_celllines,file="01_celllines.amethyst.rds")

#subset to just patient samples
dat<-subsetObject(dat,cells=row.names(dat@metadata)[!(dat@metadata$sample %in% c("MCF10A","MCF7","MDA-MB-231"))])
patient_metadata<-read.csv("/data/rmulqueen/projects/scalebio_dcis/sample_selection/simplified_patient_metadata.csv")
merge_meta<-merge(dat@metadata,patient_metadata,by.x="sample",by.y="Processing_Name",all.x=TRUE)
table(merge_meta$Sample,useNA="ifany")

dat@metadata<-merge_meta
row.names(dat@metadata)<-paste(dat@metadata$cell_id,dat@metadata$plate_info,sep="+")
saveRDS(dat,file="01_scaledcis.amethyst.rds")

```

########################################
## QC AND FILTER
########################################

```R

dat@metadata %>% group_by(Sample) %>% dplyr::summarize(cell_count=n(),mean_reads=mean(unique_reads),mean_mcg=mean(cg_cov)) 

#simple filters
dat<-subsetObject(dat,cells=row.names(dat@metadata[dat@metadata$cg_cov>=25000,]))
dat<-subsetObject(dat,cells=row.names(dat@metadata[dat@metadata$mcg_pct>=50,]))

#replace cov with cg specific cov 
dat@metadata$ch_cov<-dat@metadata$cov
dat@metadata$cov<-dat@metadata$cg_cov

#final hyperparamters for initial grouping
c=6;m=200;n=0.01;o=200
dat<-cluster_by_windows(obj=dat,
                        window_name="nakshatri_peaks",
                        outname=paste("nakshatri_init",l,m,n,o,sep="."),
                        threads=task_cpus,est_dim=l,neighbors=m,dist=n,k_pheno=o)


###### 1kb window generation, and DMR Calculation ####
dat<-dmr_and_1kb_window_gen(dat,prefix="02_scaledcis.nakshatri_peaks")

#saves DMRs as separate file, returns updated object with smoothed windows
saveRDS(dat,file="02_scaledcis.nakshatri_cluster.amethyst.rds")

```

########################################
## From Nakshatri ATAC peaks, generate list of DMRs and recluster on those.
### Nakshatri is only normal breast tissue, so clustering is hampered for DCIS/IDC epithelia
########################################

```R
#read from [prefix].dmr.[group_by].collapsed.rds
dat<-readRDS(file="02_scaledcis.nakshatri_cluster.amethyst.rds")
dmr_out<-readRDS("02_scaledcis.nakshatri_peaks.dmr.cluster_id.collapsed.rds")

#filter by significance <0.01 and logFC > 1, merge regions that overlap
dmr_out<-dmr_out |> dplyr::filter(dmr_padj<0.01) |> dplyr::filter(dmr_logFC < -1 | dmr_logFC > 1) |> select(chr,dmr_start,dmr_end) |> distinct(chr,dmr_start,dmr_end) |> makeGRangesFromDataFrame() |> reduce()
dmr_bed<-data.frame(seqnames=seqnames(dmr_out),start=start(dmr_out),end=end(dmr_out))

#create new matrix from DMR sites for refined clustering
dat@genomeMatrices[["nakshatri_dmr_sites"]] <- makeWindows(dat, 
                                                     bed = dmr_bed,
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 20, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

# Using Nakshatri initial clusters, define DMRs to recluster.

d=9 #9 dimensions
neigh=15 #15 neighbors
q=0.001
task_cpus=100

dat<-cluster_by_windows(
    obj=dat,
    window_name="nakshatri_dmr_sites",
    outname=paste("dmr_sites_initial_clustering",neigh,sep="."),
    threads=task_cpus,est_dim=d,neighbors=neigh,dist=q,
    k_pheno=neigh)

 dat@metadata$nakshatri_clusters<-dat@metadata$cluster_id

 saveRDS(dat,file="03_scaledcis.dmr_cluster.amethyst.rds")

#summarize clusters and find DMRs per cluster to help with celltyping
 out<-dmr_and_1kb_window_gen(obj=dat,prefix="03_scaledcis.dmr_cluster",groupBy = "nakshatri_clusters",threads=40)
 #saves DMRs as separate file, returns updated object with smoothed windows
 saveRDS(out,file="03_scaledcis.dmr_cluster.amethyst.rds")

# #read from [prefix].dmr.[group_by].collapsed.rdsa
# dmr_out<-readRDS("dmr_clusters.dmr.cluster_id.collapsed.rds")
# dat<-readRDS(file="03_scaledcis.dmr_cluster.amethyst.rds")
# ```

# ########################################
# ## Broad cell typing through marker genes
# ########################################

# ```R

# cell_markers<-list()
# cell_markers[["basal"]]<-c("CARMN","ACTA2","KRT17","KRT14","DST","KRT5")
# cell_markers[["lumhr"]]<-c("AREG","AZGP1","KRT18","AGR2","PIP","ANKRD30A")
# cell_markers[["lumsec"]]<-c("GABRP","ELF5","CCL28","KRT15","KIT","MMP7","LTF","SLPI")

# cell_markers[["fibro"]]<-c("DCN","APOD","LUM","COL1A2","COL1A1","FAP")
# cell_markers[["endo"]]=c("CCL21","TFF3","MMRN1","CLDN5","AL357507.1","PKHD1L1","KLHL4","LINC02147","RHOJ","ST6GALNAC3","MMRN1","MECOM","BTNL9","MCTP1","PTPRB","VWF","ADGRL4","LDB2")
# cell_markers[["perivasc"]]=c("RGS6","KCNAB1","COL25A1","ADGRL3","PRKG1","MYL9","ADIRF","NR2F2-AS1","AC012409.2")

# cell_markers[["myeloid"]]<-c("HLA-DRA","HLA-DPA1","CD74")
# cell_markers[["tcell"]]<-c("PTPRC","IKZF1","IL7R","GNLY")
# cell_markers[["mast"]]<-c("NTM","SYTL3","SLC24A3","TPSB2","HDC")
# cell_markers[["bcell"]]=c("CD37","TCL1A","LTB","HLA-DPB1","HLA-DRA","HLA-DPA1")
# cell_markers[["plasma"]]=c("IGHA2","IGHA1","JCHAIN","IGHM","IGHG1","IGHG4","IGHG3","IGHG2")

# cell_markers[["adipo"]]=c("PDE3B","ACACB","WDPCP","PCDH9","CLSTN2","ADIPOQ","TRHDE")

# cell_markers[["dcis_diff"]]=c("HOBX13","EN1","TBX15","DLX4")

# cell_colors=c(
# "basal"="#844c9d",
# "lumhr"="#e23e96",
# "lumsec"="#ff6498",
# "fibro"="#f58e90",
# "lymphatic"="#b5d564",
# "vascular"="#5bbb5a",
# "perivasc"="#eaba67",
# "myeloid"="#8088c2",
# "tcell"="#1d87c8",
# "mast"="#dcd0ff",
# "bcell"="#65cbe4",
# "plasma"="#7ecdc2",
# "adipo"="#b48454")

# nakshatri<-list()
# nakshatri[["basal"]]<-c("CELSR2","TP73","UQCRC1")
# nakshatri[["lumhr"]]<-c("ANKRD30A","CELSR1","STC2","FOXA1")
# nakshatri[["lumsec"]]<-c("TFAP2C","EPB41L1","ZMYND8","SLC52A3","WFDC3")
# nakshatri[["fibro"]]<-c("GAS7","SH3PXD2B","COL5A1","TWIST2","ISLR","GABRB3")
# nakshatri[["adipo"]]=c("PPARG","TMEM100","C14orf180","PCK1","ACACB")
# nakshatri[["endothelial"]]=c("CLEC14A")
# nakshatri[["macrophage"]]=c("HHEX","LYL1","CD300LB","C1QB","TYROBP","FPR3")
# nakshatri[["tcell"]]=c("SCML4","FYN","BCL11B","RUNX3")


# nakshatri_cell_colors=c(
# "basal"="#844c9d",
# "lumhr"="#e23e96",
# "lumsec"="#ff6498",
# "fibro"="#f58e90",
# "endothelial"="#b5d564",
# "macrophage"="#8088c2",
# "tcell"="#1d87c8",
# "adipo"="#b48454")


# print("Plotting...")
# p1 <- dimFeature(dat, colorBy = cluster_id, reduction = "umap") + ggtitle("Clusters")
# p2 <- dimFeature(dat, colorBy = sample, reduction = "umap") + ggtitle("Samples")
# p3 <- dimFeature(dat, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red"),guide="colourbar") + ggtitle("Coverage distribution")
# p4 <- dimFeature(dat, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
# p5 <- dimFeature(dat, colorBy = Group, reduction = "umap") + ggtitle("Group")
# p6 <- dimFeature(dat, colorBy = batch, reduction = "umap") + ggtitle("Batch")

# plt<-plot_grid(p1, p2,p3, p4,p5,p6,ncol=2)
# ggsave(plt,file="03_scaledcis.dmr_cluster.umap.pdf",width=20,height=30)  


# #prepare cgi
# cgisland="/data/rmulqueen/projects/scalebio_dcis/ref/cpgIslandExt.bed"
# cgi<-rtracklayer::import(cgisland)
# cgi<-as.data.frame(cgi)
# colnames(cgi)<-c("chr","start","end","strand")

# plot_histogram_page<-function(celltype){
#     plt<-histograModified(dat, 
#         baseline="mean",
#         genes = unlist(cell_markers[celltype]),
#         colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
#         matrix = "cg_cluster_id_tracks", arrowScale = .03, trackScale = .5,
#         legend = F, cgisland=cgi)
#     return(plt)
# }

# #rough ordering by what clustered together
# order<-c(
# "29","14","15","9","1",
# "24","8","35","2","16","30","17","7",
# "21","27","32","3","34","13","18",
# "19","11","4","22","5","25",
# "33","23","26","20","31",
# "6","12","28","10")
    
# plt_list<-lapply(names(cell_colors),
# function(celltype){
#     genes<-unlist(cell_markers[celltype])
#     genes<-genes[genes %in% dat@ref$gene_name]
#     print(paste("Plotting:",celltype))
#     print(paste("Genes to plot:",genes))
#     plt<-histograModified(dat, 
#         baseline="mean",
#         genes = unlist(cell_markers[celltype]),
#         colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
#         matrix = "cg_cluster_id_tracks", arrowScale = .03, trackScale = .5,
#         legend = F, cgisland=cgi,order=order) + ggtitle(celltype)
#     return(plt)
# })

# plt_out<-wrap_plots(plt_list,ncol=1) 
# ggsave(plt_out ,file="test.marker.pdf",width=30,height=length(cell_colors)*30,limitsize=F)

# plt_list<-lapply(names(nakshatri_cell_colors),
# function(celltype){
#     genes<-unlist(nakshatri[celltype])
#     genes<-genes[genes %in% dat@ref$gene_name]
#     print(paste("Plotting:",celltype))
#     print(paste("Genes to plot:",genes))
#     plt<-histograModified(dat, 
#         baseline="mean",
#         genes = unlist(nakshatri[celltype]),
#         colors= c(nakshatri_cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
#         matrix = "cg_cluster_id_tracks", arrowScale = .03, trackScale = .5,
#         legend = F, cgisland=cgi,order=order) + ggtitle(celltype)
#     return(plt)
# })

# plt_out<-wrap_plots(plt_list,ncol=1) 
# ggsave(plt_out ,file="nakshatri_test.marker.pdf",width=30,height=length(cell_colors)*30,limitsize=F)

# #set celltype plotting order
# dat@metadata$celltype<-factor(dat@metadata$celltype,levels=c(
# "basal",
# "lumhr",
# "lumsec",
# "fibro",
# "endothelial",
# "myeloid",
# "tcell"))

# cell_colors=c(
# "basal"="#844c9d",
# "lumhr"="#e23e96",
# "lumsec"="#f68d3e",
# "fibro"="#68b1af",
# "lymphatic"="#b5d564",
# "endothelial"="#b5d564",
# "vascular"="#5bbb5a",
# "perivasc"="#eaba67",
# "myeloid"="#8088c2",
# "tcell"="#1d87c8",
# "mast"="#dcd0ff",
# "bcell"="#65cbe4",
# "plasma"="#7ecdc2",
# "adipo"="#b48454",
# "dcis_diff"="black"
# )

# "basal"="#844c9d",
# "lumhr"="#e23e96",
# "lumsec"="#ff6498",
# "fibro"="#f58e90",
# "endothelial"="#b5d564",
# "macrophage"="#8088c2",
# "tcell"="#1d87c8",
# "adipo"="#b48454"

# dat@metadata$broad_celltype<-"lumhr_cancer" #20, 15, 8, 17, 11, 3, 9, 19, 12
# dat@metadata[dat@metadata$cluster_id %in% c("10"),]$broad_celltype<-"lumhr"
# dat@metadata[dat@metadata$cluster_id %in% c("4"),]$broad_celltype<-"basal"
# dat@metadata[dat@metadata$cluster_id %in% c("5"),]$broad_celltype<-"lumsec"
# dat@metadata[dat@metadata$cluster_id %in% c("6","1"),]$broad_celltype<-"stromal"
# dat@metadata[dat@metadata$cluster_id %in% c("2","18","7","16"),]$broad_celltype<-"immune"

# #run DMR and 1kb windows on reclustered data
# #running with 500 step so I can use the same chunks for methyltree
# out<-dmr_and_1kb_window_gen(dat,
#     prefix="broad_celltype",
#     groupBy="broad_celltype",
#     threads=10,step=500)

# dat<-out
# saveRDS(dat,file="04_scaledcis.amethyst.broad_celltype.rds")


# #run 3d plotting just for fun


# ```

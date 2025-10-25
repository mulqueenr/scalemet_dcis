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
                                                     threads = 50, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

dim(dat@genomeMatrices[["nakshatri_dmr_sites"]])
saveRDS(dat,file="03_scaledcis.dmr_cluster.amethyst.rds")

# Using Nakshatri initial clusters, define DMRs to recluster.

d=9 #9 dimensions
neigh=15 #15 neighbors
q=0.001
task_cpus=100

dat<-cluster_by_windows(
    obj=dat,
    window_name="nakshatri_dmr_sites",
    outname=paste("dmr_sites_initial_clustering",neigh,sep="."),
    threads=task_cpus,est_dim=d,neighbors=neigh,dist=q,k_pheno=60)

table(dat@metadata$cluster_id)

 dat@metadata$nakshatri_clusters<-dat@metadata$cluster_id

 saveRDS(dat,file="03_scaledcis.dmr_cluster.amethyst.rds")

dat<-readRDS(file="03_scaledcis.dmr_cluster.amethyst.rds")

#down to 21 clusters
out<-dmr_and_1kb_window_gen(obj=dat,prefix="03_scaledcis.dmr_cluster",groupBy = "nakshatri_clusters",threads=40) 
#saves DMRs as separate file, returns updated object with smoothed windows
saveRDS(out,file="03_scaledcis.dmr_cluster.1kb.amethyst.rds")

#read from [prefix].dmr.[group_by].collapsed.rdsa
dmr_out<-readRDS("03_scaledcis.dmr_cluster.dmr.nakshatri_clusters.collapsed.rds")
dat<-readRDS(file="03_scaledcis.dmr_cluster.1kb.amethyst.rds")
```

########################################
## Broad cell typing through marker genes
########################################

```R


print("Plotting...")
p1 <- dimFeature(dat, colorBy = cluster_id, reduction = "umap") + ggtitle("Clusters")
p2 <- dimFeature(dat, colorBy = sample, reduction = "umap") + ggtitle("Samples")
p3 <- dimFeature(dat, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red"),guide="colourbar") + ggtitle("Coverage distribution")
p4 <- dimFeature(dat, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
plt<-plot_grid(p1, p2,p3, p4,ncol=2)
ggsave(plt,file="03_scaledcis.dmr_cluster.broad_final.umap.pdf",width=20,height=20)  

cell_markers<-list()
cell_markers[["basal"]]<-c("CARMN","ACTA2","KRT17","KRT14","DST","KRT5")
cell_markers[["lumhr"]]<-c("AREG","AZGP1","KRT18","AGR2","PIP","ANKRD30A")
cell_markers[["lumsec"]]<-c("GABRP","ELF5","CCL28","KRT15","KIT","MMP7","LTF","SLPI")
cell_markers[["fibro"]]<-c("DCN","APOD","LUM","COL1A2","COL1A1","FAP")
cell_markers[["endo"]]=c("CCL21","TFF3","MMRN1","CLDN5","AL357507.1","PKHD1L1","KLHL4","LINC02147","RHOJ","ST6GALNAC3","MMRN1","MECOM","BTNL9","MCTP1","PTPRB","VWF","ADGRL4","LDB2")
cell_markers[["perivasc"]]=c("RGS6","KCNAB1","COL25A1","ADGRL3","PRKG1","MYL9","ADIRF","NR2F2-AS1","AC012409.2")
cell_markers[["myeloid"]]<-c("HLA-DRA","HLA-DPA1","CD74")
cell_markers[["tcell"]]<-c("PTPRC","IKZF1","IL7R","GNLY")
cell_markers[["mast"]]<-c("NTM","SYTL3","SLC24A3","TPSB2","HDC")
cell_markers[["bcell"]]=c("CD37","TCL1A","LTB","HLA-DPB1","HLA-DRA","HLA-DPA1")
cell_markers[["plasma"]]=c("IGHA2","IGHA1","JCHAIN","IGHM","IGHG1","IGHG4","IGHG3","IGHG2")
cell_markers[["adipo"]]=c("PDE3B","ACACB","WDPCP","PCDH9","CLSTN2","ADIPOQ","TRHDE")

cell_colors=c(
"basal"="#844c9d",
"lumhr"="#e23e96",
"lumsec"="#ff6498",
"fibro"="#f58e90",
"endo"="#5bbb5a",
"perivasc"="#eaba67",
"myeloid"="#8088c2",
"tcell"="#1d87c8",
"mast"="#dcd0ff",
"bcell"="#65cbe4",
"plasma"="#7ecdc2",
"adipo"="#b48454")


#prepare cgi
cgisland="/data/rmulqueen/projects/scalebio_dcis/ref/cpgIslandExt.bed"
cgi<-rtracklayer::import(cgisland)
cgi<-as.data.frame(cgi)
colnames(cgi)<-c("chr","start","end","strand")

plot_histogram_page<-function(celltype){
    plt<-histograModified(dat, 
        baseline="mean",
        genes = unlist(cell_markers[celltype]),
        colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
        matrix = "cg_cluster_id_tracks", arrowScale = .03, trackScale = .5,
        legend = F, cgisland=cgi)
    return(plt)
}

#rough ordering by what clustered together
order<-c(
"10","4","5","18",
"16","21","19","14",
"6",
"7","20",
"11",
"8","2","17",
"9","1","12",
"3","13","15")
    
plt_list<-lapply(names(cell_colors),
function(celltype){
    genes<-unlist(cell_markers[celltype])
    genes<-genes[genes %in% dat@ref$gene_name]
    print(paste("Plotting:",celltype))
    print(paste("Genes to plot:",genes))
    plt<-histograModified(dat, 
        baseline="mean",
        genes = unlist(cell_markers[celltype]),
        colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
        matrix = "cg_nakshatri_clusters_tracks", arrowScale = .03, trackScale = .5,
        legend = F, cgisland=cgi,order=order) + ggtitle(celltype)
    return(plt)
})

plt_out<-wrap_plots(plt_list,ncol=1) 
ggsave(plt_out ,file="03_scaledcis.dmr_cluster.test.marker.pdf",width=30,height=length(cell_colors)*30,limitsize=F)


dat@metadata$broad_celltype<-"lumhr" #10, 4, 5, 18, 16, 21, 19, 14, 6
dat@metadata[dat@metadata$cluster_id %in% c("7","20","11"),]$broad_celltype<-"basal"
dat@metadata[dat@metadata$cluster_id %in% c("17","8","2"),]$broad_celltype<-"lumsec"
dat@metadata[dat@metadata$cluster_id %in% c("12","1","9"),]$broad_celltype<-"fibro"
dat@metadata[dat@metadata$cluster_id %in% c("15","3","13"),]$broad_celltype<-"immune"

saveRDS(dat,file="04_scaledcis.broad_celltype.amethyst.rds")


print("Plotting...")

cell_colors=c(
"basal"="#844c9d",
"lumhr"="#e23e96",
"lumsec"="#ff6498",
"fibro"="#f58e90",
"endo"="#5bbb5a",
"myeloid"="#8088c2",
"lymphoid"="#1d87c8")

p1 <- dimFeature(dat, colorBy = broad_celltype, colors=cell_colors, reduction = "umap",) + ggtitle("Broad cell types")
ggsave(p1,file="03_scaledcis.broad_celltypes.umap.pdf",width=20,height=20)  
```

```R

#run 3d plotting just for fun
dim3d<-uwot::umap(
  X=dat@reductions$nakshatri_dmr_sites_irlba_regressed,
  n_neighbors = 15,
  n_components = 3,
  metric = "euclidean",
  seed = 123,
  n_threads=50,
)

dim3d<-as.data.frame(dim3d)
colnames(dim3d)<-c("X","Y","Z")
dim3d$cellname<-row.names(dim3d)
dim3d$celltype<-as.data.frame(dat@metadata[row.names(dim3d),])$broad_celltype
dim3d$hex_color<-cell_colors[dim3d$celltype]
dim3d$r_col<-unlist(col2rgb(dim3d$hex_color)["red",])
dim3d$g_col<-unlist(col2rgb(dim3d$hex_color)["green",])
dim3d$b_col<-unlist(col2rgb(dim3d$hex_color)["blue",])

write.table(dim3d,col.names=T,file="03_scaledcis.broad_celltypes.3dumap.csv",sep=",",row.names=F)

```

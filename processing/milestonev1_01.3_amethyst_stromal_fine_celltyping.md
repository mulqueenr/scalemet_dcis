```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
```


```R
source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
library(Seurat)
library(ComplexHeatmap)
library(GeneOverlap)
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)

system(paste0("mkdir -p ",project_data_directory,"/fine_celltyping"))
obj<-readRDS(file="04_scaledcis.broad_celltype.amethyst.rds")

#prepare cgi for plotting
cgisland="/data/rmulqueen/projects/scalebio_dcis/ref/cpgIslandExt.bed"
cgi<-rtracklayer::import(cgisland)
cgi<-as.data.frame(cgi)
colnames(cgi)<-c("chr","start","end","strand")

```

# Subcluster to assess more specific cell types
Clustering endothelial and fibroblast broad cell types for finer celltyping resolution.

```R

broad_celltype="stromal"

dat<-obj
output_directory=paste0(project_data_directory,"/fine_celltyping/",broad_celltype)
system(paste0("mkdir -p ",output_directory))

print(paste("Subsetting object based on broad celltype:",broad_celltype))
dat<-subsetObject(dat,cells=row.names(dat@metadata[dat@metadata$broad_celltype %in% c("endo","fibro"),]))

#overcluster by existing DMR
print(paste("Reclustering on only",broad_celltype,"cells"))
dat<-cluster_by_windows(
    dat, 
    window_name="nakshatri_dmr_sites",
    outname=paste(output_directory,
            paste(broad_celltype,"broad_dmr_sites_clustering",sep="."),
            sep="/"),
    threads=10, 
    est_dim=15, #estimated dims
    neighbors=50, 
    dist=0.01,
    k_pheno=80)

dat@metadata[[paste(broad_celltype,"clusters",sep=".")]]<-dat@metadata$cluster_id  #set celltype clusters
table(dat@metadata[[paste(broad_celltype,"clusters",sep=".")]])
#cell counts per cluster
#  1   2   3   4   5   6   7   8 
#803 566 698 468 315 288 841 262 

#clusters 
#find broad celltype specific DMRs
print(paste("Identifying DMRs on",broad_celltype,"cells"))
dat<-dmr_and_1kb_window_gen(dat,
    prefix=paste(output_directory,broad_celltype,sep="/"),
    groupBy=paste(broad_celltype,"clusters",sep="."),
    threads=10,step=1000)

collapsed_dmrs<-readRDS(file=paste0(output_directory,"/",
    broad_celltype,".dmr.",paste(broad_celltype,"clusters",sep="."),".collapsed.rds"))

table(collapsed_dmrs$celltype)
#dmr per cluster
#   1     2     3     4     5     6     7     8
#29347   623  1504  5957  4555  3571  3222  3206

#recluster on broad_celltype DMRs
print(paste("Recluster",broad_celltype,"cells on cell type DMRs"))
dmr_bed<-as.data.frame(cbind(collapsed_dmrs$chr,
                        collapsed_dmrs$dmr_start,
                        collapsed_dmrs$dmr_end))
colnames(dmr_bed)<-c("chr","start","end")
dmr_bed$start<-as.numeric(dmr_bed$start)
dmr_bed$end<-as.numeric(dmr_bed$end)

dat@genomeMatrices[["fine_dmr_sites"]] <- makeWindows(dat, 
                                                    bed = dmr_bed,
                                                    type = "CG", 
                                                    metric = "score", 
                                                    threads = 20, 
                                                    index = "chr_cg", 
                                                    nmin = 2) 

#cluster on defined DMRs within broad celltype
dat<-cluster_by_windows(
    dat, 
    window_name="fine_dmr_sites",
    outname=paste0(output_directory,"/",paste(broad_celltype,"broad_dmr_sites_clustering",sep=".")),
    threads=10, 
    est_dim=18, 
    neighbors=12, 
    dist=1E-8,
    k_pheno=65)

table(dat@metadata$cluster_id )

#overwrite temporary clusters to celltype-dmr defined clusters
dat@metadata[[paste(broad_celltype,"clusters",sep=".")]]<-dat@metadata$cluster_id  #set celltype clusters

print(paste("Make 500bp windows on",broad_celltype,"clusters"))

#summarize celltype-dmr clusters over windows
dat<-dmr_and_1kb_window_gen(dat,
    prefix=paste(output_directory,broad_celltype,sep="/"),
    groupBy=paste(broad_celltype,"clusters",sep="."),
    threads=10,step=500)

print(paste("Saving final amethyst object of",broad_celltype))
saveRDS(dat,paste0(output_directory,"/",paste("04_scaledcis",broad_celltype,"fine_celltype.amethyst.rds",sep=".")))

#read in subset amethyst object and dmrs
dat_fine<-readRDS(
    paste(project_data_directory,"fine_celltyping",broad_celltype,
    paste("04_scaledcis",broad_celltype,"fine_celltype.amethyst.rds",sep="."),
    sep="/"))

#UMAP PLOTTING
p1 <- dimFeature(dat_fine, colorBy = cluster_id, reduction = "umap") + ggtitle(broad_celltype)
p2 <- dimFeature(dat_fine, colorBy = mcg_pct, reduction = "umap") + ggtitle("%metCG")+scale_color_gradientn(colors = c("black", "turquoise", "gold", "red"),guide="colourbar")
p3 <- dimFeature(dat_fine, colorBy = log10(cov), reduction = "umap") + ggtitle("cov")+ scale_color_gradientn(colors = c("black", "turquoise", "gold", "red"),guide="colourbar")
p4 <- dimFeature(dat_fine, colorBy = Group, reduction = "umap") + ggtitle("Group")
p5 <- dimFeature(dat_fine, colorBy = Race_Ethnicity, reduction = "umap") + ggtitle("Race/Ethnicity")
p6 <- dimFeature(dat_fine, colorBy = Age, reduction = "umap") + ggtitle("Age")
ggsave((p1|p2)/(p3|p4)/(p5|p6),file=paste0(output_directory,"/",paste("04_scaledcis",broad_celltype,"fine_celltype.umap.pdf",sep=".")),width=20,height=20)  


#set order by umap groupings
order<-c("1","3","9","7","2","8","10","6","4","5")


#STACKED BARPLOT PER CLUSTER PLOTTING
p1<-ggplot(dat_fine@metadata, aes(fill=Group,x=cluster_id)) + 
    geom_bar(position="fill", stat="count")+
    scale_x_discrete(limits=factor(order))
p2<-ggplot(dat_fine@metadata, aes(fill=batch,x=cluster_id)) + 
    geom_bar(position="fill", stat="count")+
    scale_x_discrete(limits=factor(order))
p3<-ggplot(dat_fine@metadata, aes(fill=method,x=cluster_id)) + 
    geom_bar(position="fill", stat="count")+
    scale_x_discrete(limits=factor(order))
p4<-ggplot(dat_fine@metadata, aes(fill=Sample,x=cluster_id)) + 
    geom_bar(position="fill", stat="count")+
    scale_x_discrete(limits=factor(order))
p5<-ggplot(dat_fine@metadata, aes(fill=Menopause,x=cluster_id)) + 
    geom_bar(position="fill", stat="count")+
    scale_x_discrete(limits=factor(order))
p6<-ggplot(dat_fine@metadata, aes(color=cluster_id,y=mcg_pct,x=cluster_id)) + 
    geom_jitter()+geom_violin()+ylim(c(0,100))+
    scale_x_discrete(limits=factor(order))
p7<-ggplot(dat_fine@metadata, aes(color=cluster_id,y=cov,x=cluster_id)) + 
    geom_jitter()+geom_violin()+
    scale_x_discrete(limits=factor(order))
ggsave(p1/p2/p3/p4/p5/p6/p7,file=paste0(output_directory,"/",paste("04_scaledcis",broad_celltype,"fine_celltype.stacked_barplots.pdf",sep=".")),width=10,height=40)  

#more broad cell markers
cell_markers<-list()
cell_markers[["fibro"]]<-c("DCN","APOD","LUM","COL1A2","COL1A1","FAP")
cell_markers[["lymphatic"]]=c("PKHD1L1","KLHL4","LINC02147","RHOJ","ST6GALNAC3","MMRN1") #"AL357507.1",
cell_markers[["vascular"]]=c("MECOM","BTNL9","MCTP1","PTPRB","VWF","ADGRL4","LDB2")
cell_markers[["perivasc"]]=c("RGS6","KCNAB1","COL25A1","ADGRL3","PRKG1","MYL9","ADIRF","NR2F2-AS1","AC012409.2")
cell_markers[["adipo"]]=c("PDE3B","ACACB","WDPCP","PCDH9","CLSTN2","ADIPOQ","TRHDE")

cell_markers<-cell_markers
plt_list<-lapply(names(cell_markers),
    function(celltype){
            genes<-unlist(cell_markers[celltype])
            print(paste("Plotting:",celltype))
            print(paste("Genes to plot:",genes))
            plt<-histograModified(obj=dat_fine, 
                baseline="mean",
                genes = unlist(cell_markers[celltype]),
                order=order, trackOverhang=10000, #promoter_focus=TRUE,
                matrix = paste0("cg_",broad_celltype,".clusters_tracks"), arrowScale = .03, trackScale = .5,
                legend = F, cgisland=cgi) + ggtitle(celltype)
            ggsave(plt,
            file=paste(project_data_directory,"fine_celltyping",broad_celltype,
            paste("04_scaledcis",broad_celltype,celltype,"broad_celltype.markers.pdf",sep="."),
            sep="/"),
            width=length(genes)*2,
            height=length(order)*2,
            limitsize=F)
        })


#fine stromal markers
stromal_markers<-list()
stromal_markers[["fibro_major"]]<-c('KDM6B', 'HMOX1', 'GPRC5A', 'TNFAIP6', 'TNFRSF12A', 'HMGA1', 'CXCL1', 'NOP16', 'ATP13A3', 'CXCL2', 'GEM', 'MAT2A', 'DDX21', 'KLHL21', 'EFNB2', 'AMD1', 'PNO1', 'GCLM', 'WDR43', 'ATP1B3')
stromal_markers[["fibro_prematrix"]]<-c('PCOLCE2', 'WISP2', 'MFAP5', 'CHRDL1', 'GPX3', 'EFEMP1', 'ADH1B', 'TNXB', 'MGST1', 'FBLN2', 'PTGIS', 'GPC3', 'NFIB', 'C3', 'ABCA8', 'TGFBR2', 'PODN', 'ITIH5', 'PRELP', 'FHL1')
stromal_markers[["fibro_CAF"]]<-c('POSTN', 'COL10A1', 'INHBA', 'COL8A1', 'CTHRC1', 'COL5A1', 'MMP11', 'MFAP2', 'ITGBL1', 'BGN', 'SULF1', 'FAP', 'ASPN', 'COMP', 'SULF2', 'COL11A1', 'COL5A2', 'NREP', 'LOXL1', 'RAB31')
stromal_markers[["fibro_matrix"]]<-c('DLK1', 'IGFBP2', 'APOE', 'TAC1', 'ADAM12', 'GPM6B', 'TGFBI', 'COL15A1', 'EDNRB', 'DCX', 'CAPN6', 'ITM2A', 'VEGFD', 'IGF1', 'RELN', 'MMP16', 'OLFML2A', 'EFNB2', 'COL4A2', 'PLXDC1')
stromal_markers[["fibro_SFRP4"]]<-c('SFRP4', 'CD9', 'CLU', 'LTBP1', 'GREM1', 'IGFBP5', 'LEPR', 'OGN', 'KCNMA1', 'ABI3BP', 'CRYAB', 'ADIRF', 'RRAD', 'C12orf75', 'CAV1', 'PAPPA2', 'PTGIS', 'IGFBP6', 'PLP2', 'NFIB')
stromal_markers[["fibro_stress"]]<-c('TXNIP', 'FOS', 'HSPA1B', 'NR2F1', 'HSPA6', 'SPRY1', 'JUN', 'ABCA10', 'SOCS3', 'EGR1', 'ZFP36', 'HSPA1A', 'APOD', 'BTG2', 'DDIT3', 'PTGDS', 'DDIT4', 'MAFB', 'MRPL18', 'ABCA6')

stromal_markers[["endo_artery"]]<-c('SEMA3G', 'CXCL12', 'IGFBP3', 'DEPP1', 'HEY1', 'LTBP4', 'PPP1R14A', 'FN1', 'PCSK5', 'ARL15', 'RHOB', 'FBLN5', 'TSPAN2', 'GUCY1A1', 'EGR1', 'SLC9A3R2', 'ICAM2', 'GJA4', 'NEBL', 'ANXA3')
stromal_markers[["endo_vein"]]<-c('ACKR1', 'SELE', 'SELP', 'IL1R1', 'ZNF385D', 'CLU', 'CNKSR3', 'VCAN', 'CCL14', 'NCOA7', 'EPB41L3', 'TLL1', 'OLFM1', 'PRCP', 'VCAM1', 'GNG12', 'NPC2', 'FAM84B', 'CYP1B1', 'ACTN1')
stromal_markers[["endo_capillary"]]<-c('BTNL9', 'CD300LG', 'MT1M', 'RGCC', 'MT1A', 'C11orf96', 'MT1E', 'ADGRF5', 'CA4', 'CD36', 'RBP7', 'MSX1', 'GPIHBP1', 'ITGA1', 'SGK1', 'MT1X', 'SCARB1', 'LITAF', 'KDR', 'MLEC')
stromal_markers[["endo_TEC"]]<-c('MMP2', 'INSR', 'HECW2', 'PODXL', 'ESM1', 'RGCC', 'NOTCH4', 'MAP1B', 'HTRA1', 'MLEC', 'COL4A2', 'PMEPA1', 'VWA1', 'PTPRG', 'COL6A2', 'NID1', 'TCIM', 'PDGFD', 'PLPP1', 'THY1')
stromal_markers[["endo_cycling"]]<-c('TYMS', 'NUSAP1', 'PCLAF', 'TOP2A', 'CENPF', 'MKI67', 'TK1', 'TPX2', 'GGH', 'CDK1', 'BIRC5', 'ASPM', 'CENPK', 'DHFR', 'DIAPH3', 'ZWINT', 'GTSE1', 'CENPW', 'UBE2C', 'KNL1')
stromal_markers[["endo_lymphatic"]]<-c('CCL21', 'TFF3', 'MMRN1', 'EFEMP1', 'PDPN', 'MRC1', 'NR2F1', 'LYVE1', 'PROX1', 'LAPTM5', 'COLEC12', 'IGF1', 'SEMA3D', 'MPP7', 'AKAP12', 'PKHD1L1', 'ABI3BP', 'CCND1', 'TFPI', 'FXYD6')

stromal_markers[["peri_THY1"]]<-c('CD36', 'CYGB', 'THY1', 'HIGD1B', 'COL3A1', 'RGS5', 'IGFBP2', 'COL6A3', 'ADGRF5', 'ARHGDIB', 'FAM213A', 'COX4I2', 'FAM162B', 'LPL', 'GUCY1A2', 'MYO1B', 'TXNIP', 'CYBA', 'PLXDC1')
stromal_markers[["peri_CCL21"]]<-c('CTSC', 'STEAP4', 'C1S', 'GGT5', 'C1R', 'COL6A3', 'CFH', 'S100A10', 'FGF7', 'CCL21', 'EMP1', 'CFD', 'TMEM176B', 'CXCL12', 'CHRDL1', 'TMEM176A', 'CLSTN2', 'FHL2', 'CCL19', 'CD44')
stromal_markers[["peri_stress"]]<-c('HSPA6', 'CACYBP', 'HSPH1', 'CHORDC1', 'FKBP4', 'AHSA1', 'STIP1', 'MIR4435-2HG', 'ABL2', 'ZFAND2A', 'THY1', 'HSPA4', 'BAG3', 'TENT5A', 'MRPL18', 'TCP1', 'SERPINH1', 'CYTOR', 'RUNX1', 'DNAJA4')

stromal_markers[["VSMC_CREM"]]<-c('MYH11', 'ZNF331', 'CREM', 'NR4A2', 'CD9', 'CNN1', 'SYNM', 'PLEKHO1', 'SORBS2', 'CCDC107', 'FAM107B', 'ELL2', 'ATP1B3', 'NDUFA4', 'ZFHX3', 'C12orf75', 'BCAM', 'RERGL', 'EIF4A3')
stromal_markers[["VSMC_PLN"]]<-c('NET1', 'TXNIP', 'RERGL', 'CASQ2', 'BCAM', 'PLN', 'KLHL23', 'NDUFA4', 'MEF2C', 'RCAN2', 'NTRK2', 'LBH', 'PLAC9', 'HES1', 'WFDC1', 'TBX2', 'AC023157.3', 'TBX2-AS1', 'ACTA2', 'KCNAB1')
stromal_markers[["VSMC_stress"]]<-c('HSPA6', 'MYH11', 'HSPA1B', 'HSPH1', 'CRYAB', 'BAG3', 'CHORDC1', 'DNAJB1', 'SORBS2', 'GADD45G', 'HSPB8', 'NET1', 'MRPL18', 'PLN', 'SERPINH1', 'HSPA1A', 'BCAM', 'PHLDA2', 'ID2', 'DNAJB4')


#set order by umap groupings
cell_markers<-stromal_markers
plt_list<-lapply(names(cell_markers),
    function(celltype){
            genes<-unlist(cell_markers[celltype])
            print(paste("Plotting:",celltype))
            print(paste("Genes to plot:",genes))
            plt<-histograModified(obj=dat_fine, 
                baseline="mean",
                genes = unlist(cell_markers[celltype]),
                order=order, promoter_focus=TRUE,trackOverhang=10000,
                matrix = paste0("cg_",broad_celltype,".clusters_tracks"), arrowScale = .03, trackScale = .5,
                legend = F, cgisland=cgi) + ggtitle(celltype)
            ggsave(plt,
            file=paste(project_data_directory,"fine_celltyping",broad_celltype,
            paste("04_scaledcis",broad_celltype,celltype,"fine_celltype.markers.pdf",sep="."),
            sep="/"),
            width=length(genes)*2,
            height=length(order)*2,
            limitsize=F)
        })


#dmr set
dmr<-readRDS(
        paste(project_data_directory,"fine_celltyping",broad_celltype,
        paste(broad_celltype,"dmr",broad_celltype,"clusters.collapsed.rds",sep="."),
        sep="/"))

#rna markers defined by paired 10x RNA analysis
rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.rds")
table(rna$broad_celltype)
rna<-subset(rna,broad_celltype %in% c("fibro","endo","perivasc","vascular","lymphatic")) #using all stromal types
rna<-subset(rna,fine_celltype %in% c("suspected_doublet"),invert=TRUE)
table(rna$fine_celltype)
#rna<-subset(rna,fine_celltype %in% c("fibro_SFRP4","fibro_major","fibro_matrix","fibro_prematrix","fibro_CAF"))

rna<-JoinLayers(rna)
markers<-FindAllMarkers(rna,group.by="fine_celltype")

#plot overlap of RNA markers for fine cell types and our cluster hypo DMRs
plot_dmr_RNAmarker_overlap(dmr=dmr,
                            rna_markers=markers,
                            celltype=broad_celltype,
                            clus_order=order,
                            plot_value="pval")

plot_dmr_RNAMarker_hypoM(dat=dat_fine,
                            dmr=dmr,
                            order=order,
                            rna_markers=markers,
                            broad_celltype,
                            gtf_file="/container_ref/gencode.v43.annotation.gtf.gz")

#grab cluster 7 dmr hypo markers
clu7_markers<-dmr %>% filter(celltype=="7") %>% filter(direction=="hypo") %>% slice_min(dmr_logFC,n=100)
clu7_markers<-unlist(
        lapply(
            unlist(
                strsplit(
                    unlist(clu7_markers$gene_names),",")),
                function(x) gsub("[[:space:]]", "", x)
            )
        )
#overlap with RNA markers
markers %>% filter(gene %in% clu7_markers) %>% filter(avg_log2FC>2) %>% select(cluster) %>% table()

######FINAL STROMAL CELL TYPING###########
dat_fine@metadata$fine_celltypes<-"fibro1"
dat_fine@metadata[dat_fine@metadata[[paste(broad_celltype,"clusters",sep=".")]] %in% c("3","1"),]$fine_celltypes<-"endo"
dat_fine@metadata[dat_fine@metadata[[paste(broad_celltype,"clusters",sep=".")]] %in% c("9"),]$fine_celltypes<-"endo2"
dat_fine@metadata[dat_fine@metadata[[paste(broad_celltype,"clusters",sep=".")]] %in% c("5","6","10"),]$fine_celltypes<-"fibro2"
dat_fine@metadata[dat_fine@metadata[[paste(broad_celltype,"clusters",sep=".")]] %in% c("7"),]$fine_celltypes<-"peri"

#summarize celltype-dmr clusters over windows
dat_fine<-dmr_and_1kb_window_gen(dat_fine,
    prefix=paste(output_directory,broad_celltype,sep="/"),
    groupBy="fine_celltypes",
    threads=10,step=500)

print(paste("Saving final amethyst object of",broad_celltype))
saveRDS(dat_fine,paste0(output_directory,"/",paste("04_scaledcis",broad_celltype,"fine_celltype.amethyst.rds",sep=".")))



# #plot with atac peak overlaps
# atac_markers<-read.table("/data/rmulqueen/projects/scalebio_dcis/ref/CEDAR_multiome.celltypes.markers.peaks.df.tsv",header=T)
# atac_markers$chr<-unlist(lapply(strsplit(atac_markers$feature,"-"),"[",1))
# atac_markers$start<-unlist(lapply(strsplit(atac_markers$feature,"-"),"[",2))
# atac_markers$end<-unlist(lapply(strsplit(atac_markers$feature,"-"),"[",3))
# clus_order<-order

# plot_dmr_ATACmarker_overlap<-function(dmr,atac_markers,celltype,clus_order,plot_value="Jaccard"){
#     #plot_value can be one of pval, odds.ratio, intersection, union, Jaccard
#     #flow through:
#     # - filter DMRs to HYPO
#     # - filter RNA markers to overexpressed
#     # - table dmr to marker count, and plot odds ratio
#     # - plot gene body methylation hypoM plots
#     output_directory=paste0(project_data_directory,"/fine_celltyping/",broad_celltype)

#     #PLOT HYPO
#     #filter to hypo DMRs (more likely to increase expression)
#     dmr_hypo<- dmr %>% group_by(test) %>% filter(direction=="hypo") 
#     dmr_granges<-makeGRangesFromDataFrame(dmr_hypo,keep.extra.columns=TRUE)
#     dmr_granges<-split(dmr_granges,dmr_granges$test)

#     atac_markers <- atac_markers %>% filter(logFC>0) #look only for increased expression markers
#     atac_markers_granges<-makeGRangesFromDataFrame(atac_markers,keep.extra.columns=TRUE)
#     atac_markers_granges<-split(atac_markers_granges,atac_markers_granges$group)

#     overlap<-as.data.frame(do.call("rbind",lapply(1:length(dmr_granges),calc_overlap_stats)))
#     colnames(overlap)<-c("cluster","celltype","odds.ratio","pval","jaccard")

#     dmr_count<-dmr_hypo %>% group_by(test,.drop=FALSE) %>% tally() 
#     dmr_count<-setNames(dmr_count$n,nm=dmr_count$test)

#     atac_count<-atac_markers %>% group_by(group,.drop=FALSE) %>% tally() 
#     atac_count<-setNames(atac_count$n,nm=atac_count$group)

#     or_mat<-as.data.frame(data.table::dcast(data=as.data.table(overlap),formula=cluster~celltype,value.var=plot_value))
#     row.names(or_mat)<-or_mat$cluster
#     or_mat<-or_mat[,2:ncol(or_mat)]
#     or_mat <- or_mat %>%
#                 mutate_all(~as.numeric(as.character(.)))    
    
#     print("Plotting hypo overlap of DMRs...")
#     pdf(paste0(output_directory,"/","04_scaledcis.",celltype,".dmr_ATAC_overlap.heatmap.pdf"))
#     print(
#         Heatmap(or_mat,
#             bottom_annotation=HeatmapAnnotation(bar=anno_barplot(atac_count)),
#             row_order=clus_order,cluster_rows=FALSE,
#             show_row_names = TRUE, show_column_names = TRUE, row_names_side="right",column_names_side="bottom",
#             column_title = "Hypomethylated and ATAC Peaks",
#             left_annotation=rowAnnotation(bar=anno_barplot(dmr_count)),
#             name=plot_value,
#             cell_fun = function(j, i, x, y, width, height, fill) {
#                                 grid.text(format(round(or_mat,digits=2),nsmall=2)[i, j], x, y, gp = gpar(fontsize = 10))
#                                 }
#             ))
#     dev.off()
# }


# plot_dmr_ATACmarker_overlap(dmr=dmr,
#                             clus_order=order,
#                             atac_markers=atac_markers,
#                             celltype=broad_celltype,
#                             plot_value="pval")
   
# #plot percentage cells in fine celltype per patient(rna), and percentage cells in cluster per patient (met)

# # #plot RNA defined markers that overlap with large DMRs
#  plot_dmr_RNAMarker_hypoM(dat=dat_fine,
#                              dmr=dmr,
#                              rna_markers=markers,
#                              broad_celltype=broad_celltype,
#                              order=c("1","3","5","9","4","14","2","12","11","7","6","13","10","8")) #loose ordering by umap grouping

# # #ASSIGN CELL TYPES BY PLOT
# # dat_fine@metadata$fine_celltypes<-"unknown"
# # dat_fine@metadata[dat_fine@metadata[[paste(celltype,"clusters",sep=".")]] %in% c("0"),]$fine_celltypes<-"unknown"

# # #SAVE FINE CELL TYPES
# # saveRDS(dat_fine,paste("04_scaledcis",celltype,"fine_celltype.amethyst.rds",sep="."))

# ```


# Integrating RNA and Methylation based on promoter methylation

# Proper way to do this is correlate RNA expression with different methylation metrics per gene (whole body, promoter, promoter extended etc) but for now just doing promoter +5kb into gene body
# ```R

# library(Seurat)
# library(Signac)

# #rna markers defined by paired 10x RNA analysis
# rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.rds")
# table(rna$broad_celltype)
# rna<-subset(rna,broad_celltype %in% c("fibro","endo","perivasc","vascular","lymphatic")) #using all stromal types
# rna<-subset(rna,fine_celltype %in% c("suspected_doublet"),invert=TRUE)
# table(rna$fine_celltype)
# #rna<-subset(rna,fine_celltype %in% c("fibro_SFRP4","fibro_major","fibro_matrix","fibro_prematrix","fibro_CAF"))

# rna<-JoinLayers(rna)
# markers<-FindAllMarkers(rna,group.by="fine_celltype")
# rna_markers<-markers %>% filter(p_val_adj<0.05 & avg_log2FC>2) %>% slice_max(order_by=avg_log2FC,by=cluster,n=1000)

# #read in subset amethyst object and dmrs
# dat_fine<-readRDS(
#     paste(project_data_directory,"fine_celltyping",broad_celltype,
#     paste("04_scaledcis",broad_celltype,"fine_celltype.amethyst.rds",sep="."),
#     sep="/"))

# rna_markers_bed <- dat_fine@ref %>% 
#                     filter(gene_name %in% rna_markers$gene) %>% 
#                     filter(type=="gene") %>% 
#                     filter(!duplicated(gene_name)) %>% 
#                     mutate(
#                         chr=seqid,
#                         promoter_start=start-5000,
#                         promoter_end=start+5000,
#                         window_name=paste(chr,promoter_start,promoter_end,sep="_")) %>% 
#                     select(chr,promoter_start,promoter_end,window_name,gene_name) 
# windows_name<-setNames(rna_markers_bed$gene_name,nm=rna_markers_bed$window_name)
# #make bed of promoters
# window_name="rna_markers"
# dat@genomeMatrices[[window_name]] <- makeWindows(dat_fine, 
#                                                     bed = rna_markers_bed[,1:3],
#                                                     type = "CG", 
#                                                     metric = "score", 
#                                                     threads = 20, 
#                                                     index = "chr_cg", 
#                                                     nmin = 2) 
# row.names(dat@genomeMatrices[[window_name]])<-windows_name[row.names(dat@genomeMatrices[[window_name]])]
# dat@genomeMatrices[[window_name]] <- dat@genomeMatrices[[window_name]][rowSums(!is.na(dat@genomeMatrices[[window_name]])) >= 45, ]

# dat@reductions[[paste(window_name,"irlba",sep="_")]] <- runIrlba(dat, genomeMatrices = c(window_name), dims = 8, replaceNA = c(0))

# dat@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- regressCovBias(dat, reduction = paste(window_name,"irlba",sep="_")) 


# met_obj<-CreateSeuratObject(counts=as.matrix(dat@genomeMatrices[["rna_markers"]]),
#                 data = as.matrix(dat@genomeMatrices[["rna_markers"]]), 
#                             assay="met", 
#                             project = "met",
#                             min.cells = 0, 
#                             min.features = 0)

# met_obj<-AddMetaData(met_obj,dat@metadata)

# transfer.anchors <- FindTransferAnchors(
#     reference = rna, 
#     query = met_obj, 
#     features = rna_markers$gene,
#     reference.assay = "RNA", query.assay = "met", reduction = "cca")
 
# #caused a memory core fault

# ```
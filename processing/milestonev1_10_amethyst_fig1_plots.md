
```R
library(amethyst)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(ggplot2)
library(patchwork)
#source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)
obj<-readRDS(file="07_scaledcis.cnv_clones.amethyst.rds")
rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")


# meta<-read.csv("/data/rmulqueen/projects/scalebio_dcis/sample_selection/simplified_patient_metadata.csv")
# meta<-meta[!duplicated(meta$Sample),]
# row.names(meta)<-meta$Sample

# obj@metadata$ER<-meta[obj@metadata$Sample,]$ER
# obj@metadata$PR<-meta[obj@metadata$Sample,]$PR
# obj@metadata$HER2<-meta[obj@metadata$Sample,]$HER2
# obj@metadata$Group<-meta[obj@metadata$Sample,]$Group
# obj@metadata$DCIS_Grade<-meta[obj@metadata$Sample,]$DCIS_Grade
# obj@metadata$Race_Ethnicity<-meta[obj@metadata$Sample,]$Race_Ethnicity
# obj@metadata$IDC_Differentiation<-meta[obj@metadata$Sample,]$IDC_Differentiation

# saveRDS(obj,file="07_scaledcis.cnv_clones.amethyst.rds")

# saveRDS(obj,file="06_scaledcis.celltype.amethyst.rds")
# write.table(obj@metadata,file="06_scaledcis.celltype.metadata.csv.",sep=",",col.names=T,row.names=T,quote=F)

####################################################
#           Fig 1 Sample Heatmap                  #
###################################################

met<-obj@metadata
#summarize per sample (heatmap categorical)
met<-met[!duplicated(met$Sample),]

meta_cat<-met %>% group_by(Sample) %>% select(Sample,Group,DCIS_Grade,IDC_Differentiation,Age,Race_Ethnicity,Menopause,ER,PR,HER2) %>% as.data.frame()
row.names(meta_cat)<-meta_cat$Sample
meta_cat[which(is.na(meta_cat),arr.ind=T)]<-"N/A"
meta_cat[which(meta_cat=="",arr.ind=T)]<-"N/A"

meta_cat$Group<-factor(meta_cat$Group,levels=c("HBCA","DCIS","Synchronous","IDC"))
row_order<-rev(meta_cat %>% arrange(Group,ER,PR,HER2,Sample) %>% pull(Sample))

race_ethnicity_col=c(
    "African_American"="#0F6FC6",
    "Asian"="#10CF9B",           
    "White"="#A5C249",  
    "Hispanic_or_Latino"="#DBEFF9")

age_col=colorRamp2(breaks=c(min(meta_cat$Age,na.rm=T),max(meta_cat$Age,na.rm=T)),c("#d789d7","#2a3d66"))

class_col=c("+"="black",
            "-"="grey",
            "N/A"="white")

grade_col=c("N/A"="white",
            "G1"="#CCDF92",
            "G2"="#8A9A5B",
            "G3"="#45503B")
differentiation_col=c("N/A"="white",
            "well"="#B7957C",
            "moderate"="#734939",
            "poor"="#A6432D")

menopause_col=c("Hysterectomy (perimenopausal)"="#E2BD6B",
                "Perimenopausal"="#E2BD6B",
                "Post-menopausal(Hysterectomy)"="#4D067B",
                "Post-menopausal"="#4D067B",
                "Pre-menopausal"="#B984DB",
                "Unknown"="white",
                "s/p hysterectomy"="#E2BD6B")

group_col=c("DCIS"="#278192",
            "HBCA"="#20223E",
            "Synchronous"="#00B089",
            "IDC"="#8FF7BD")

#plot metadata
ha = rowAnnotation(
  Group=meta_cat$Group,
  DCIS_Grade=meta_cat$DCIS_Grade,
  IDC_Differentiation=meta_cat$IDC_Differentiation,
  Age=meta_cat$Age,
  Race_Ethnicity=meta_cat$Race_Ethnicity,
  Menopause=meta_cat$Menopause,
  col = list(Group=group_col,
            DCIS_Grade=grade_col,
            IDC_Differentiation=differentiation_col,
            Age=age_col,
            Race_Ethnicity=race_ethnicity_col,
            Menopause=menopause_col))

plt<-Heatmap(meta_cat[,c("ER","PR","HER2")],
 row_split=meta_cat$group,
 cluster_row_slices=FALSE,
 left_annotation=ha,row_order=row_order)

 pdf("scaledcis.metadata.heatmap.pdf",width=10)
 print(plt)
 dev.off()

####################################################
#           Fig 1 Stacked Celltype ID             #
###################################################
#Make stacked barplot on identities per cluster


#color assignment is fluor as cancer associated cell type, rest muted versions
celltype_col=c(
"cancer"="#ff00ff",
"basal"="#92278F",
"lumsec"="#AD8CFF",
"lumhr"="#FF6AD5",

"TEC"="#ff8000",
"endothelial"="#ffab5f",
"endo"="#ffab5f",

"CAF"="#ff2222",
"pericyte_VSMC"="#FA7876",
"perivasc"="#FA7876",

"fibroblast"="#9b1c31",
"fibro"="#9b1c31",

"TAM"="#ccff00",
"TAM_2"="#00ff80",
"monocyte"="#98d3b9",
"macrophage"="#00af5f",
"DC"="#008080",
"myeloid"="#00af5f",


"nk_tnk"="#00ffff",
"tcell_cd4"="#00bae5",
"tcell_cd8"="#1800ff",
"tcell_cd8_2"="#0016b7",
"bcell"="#87ceeb",
"plasma"="#00ffff",
"tcell"="#1800ff"
)

group<-setNames(nm=row.names(meta_cat),meta_cat[,"Group"])
row_order<-rev(row_order)
#met cell counts
met_count<-obj@metadata %>% as.data.frame() %>% dplyr::count(Sample, .drop=FALSE)
met_count$Group<-group[met_count$Sample]
plt1<-ggplot(met_count,aes(x=Sample,fill=Group,y=n))+geom_bar(stat="identity")+theme_minimal()+scale_x_discrete(limits=c(row_order))+scale_fill_manual(values=group_col)

#met cell types
met_cell<-obj@metadata %>% as.data.frame() %>% dplyr::count(Sample,celltype, .drop=FALSE)
met_cell$Group<-group[met_cell$Sample]
met_cell$celltype<-factor(met_cell$celltype,levels=names(celltype_col))
plt2<-ggplot(met_cell,aes(x=Sample,fill=celltype,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()+scale_x_discrete(limits=c(row_order))+scale_fill_manual(values=celltype_col)


#met cell types percent
met_cell<-obj@metadata %>% as.data.frame() %>% dplyr::count(Group,celltype, .drop=FALSE)
met_cell$Group<-group[met_cell$Group]
met_cell$celltype<-factor(met_cell$celltype,levels=names(celltype_col))
plt5<-ggplot(met_cell,aes(x=celltype,fill=Group,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()+scale_x_discrete(limits=names(celltype_col))+scale_fill_manual(values=group_col)

ggsave(plt5,file="scaledcis.assigned_celltype_by_group.barplots.pdf",width=50,limitsize=F)

#rna cell counts
rna_counts<-rna@meta.data %>% as.data.frame() %>% dplyr::count(sample, .drop=FALSE)
rna_counts$Group<-group[rna_counts$sample]
plt3<-ggplot(rna_counts,aes(x=sample,fill=Group,y=n))+geom_bar(stat="identity")+theme_minimal()+scale_x_discrete(limits=c(row_order))+scale_fill_manual(values=group_col)

#rna cell types
rna_cell<-rna@meta.data %>% as.data.frame() %>% dplyr::count(sample,coarse_celltype, .drop=FALSE)
rna_cell$Group<-group[rna_cell$sample]
rna_cell$coarse_celltype<-factor(rna_cell$coarse_celltype,levels=names(celltype_col))
plt4<-ggplot(rna_cell,aes(x=sample,fill=coarse_celltype,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()+scale_x_discrete(limits=c(row_order))+scale_fill_manual(values=celltype_col)

plt<-patchwork::wrap_plots(plt1,plt2,plt3,plt4,plt5,ncol=1,guides='collect')
ggsave(plt,file="scaledcis.assigned_celltype_barplots.pdf",width=50,limitsize=F)


clone_cov<-obj@metadata %>% filter(!endsWith(cnv_ploidy_500kb,suffix="diploid")) %>% filter(!is.na(cnv_clonename) & cnv_clonename!="NA") %>% group_by(Sample,cnv_clonename) %>% summarize(count=n(),cov=sum(cg_cov)/28.3e6)

clone_cov %>% filter(cov>0.5)
select(Sample,cnv_clonename,cov) 


markers<-list()
markers[["basal"]]<-c("KRT14","KRT17","SAA1","SFN")
markers[["lumsec"]]<-c("LTF","KRT15","MMP7")
markers[["lumhr"]]<-c("AREG","AZGP1","KRT18","ANKRD30A")

markers[["endothelial"]]<-c("PECAM1","CDH5","CD93")
markers[["TEC"]]<-c("APLNR","VWA1","COL6A2","IGFBP7")
markers[["fibroblast"]]<-c("COL1A1","DCN","COL1A2")
markers[["CAF"]]<-c("FAP","CXCL12","ACTA2","COL10A1")
markers[["perivasc"]]<-c("RGS5","ABCC9","KCNJ8","MYH11","TAGLN")

markers[["immune"]]<-c("PTPRC","IKZF1")
markers[["tcell"]]<-c("IL7R","CD3D","CD3G","CD4","CD3E")
markers[["tcell_cd4"]]<-c("ICOS","CD28")
markers[["tcell_cd8"]]<-c("CCL5","NKG7","KLRK1","CXCR3")
markers[["tcell_nk"]]<-c("FGFBP2","CCL4","GZMB","XCL2")
markers[["bcell/plasma"]]<-c("CD79A","CD19","TCL1A","MS4A1","CD74","JCHAIN","SDC1","TNFRSF17")

markers[["DC"]]<-c("BCL11A","GRLH2","RORA","IRF7","PACSIN")

markers[["myeloid"]]<-c("SPI1","LYZ","CD68")
markers[["macrophage"]]<-c("AIF1","CSF1R","S100A8")
markers[["monocyte"]]<-c("S100A9","FCGR3A","CD36","FCGR3B","MPO","MPEG1")
markers[["TAM"]]<-c("TREM2","C1QA","PLA2G4A","HACD1")


Idents(rna)<-factor(rna$coarse_celltype,levels=c("basal","lumsec","lumhr","cancer","perivasc","fibro","endo","myeloid","tcell","bcell","plasma"))
plt<-DotPlot(rna, features = markers, cluster.idents=FALSE,cols=c("lightgrey","darkred")) + RotatedAxis()
ggsave(plt,file="tenx_dcis.coarse_celltype.dotplot.pdf",width=30,height=10,limitsize=FALSE)

		

					
									
####################################################
#           Fig 1 Geom Jitter/Box plot             #
###################################################
#Make stacked barplot on identities per cluster

#uniq reads
plt1<-ggplot(obj@metadata,aes(x=Sample,y=log10(unique_reads*2),color=Group))+#geom_jitter(alpha=0.2)
geom_boxplot(fill=NA,outlier.shape=NA)+scale_color_manual(values=group_col)+scale_x_discrete(limits=c(row_order))+theme_minimal()

#%Mito
plt2<-ggplot(obj@metadata,aes(x=Sample,y=mito_reads,color=Group))+#geom_jitter(alpha=0.2)+
geom_boxplot(fill=NA,outlier.shape=NA)+scale_color_manual(values=group_col)+scale_x_discrete(limits=c(row_order))+theme_minimal()

#%telomere

#mCG perc
plt3<-ggplot(obj@metadata,aes(x=Sample,y=mcg_pct,color=Group))+#geom_jitter(alpha=0.2)+
geom_boxplot(fill=NA,outlier.shape=NA)+scale_color_manual(values=group_col)+scale_x_discrete(limits=c(row_order))+theme_minimal()

#CG Covered
plt4<-ggplot(obj@metadata,aes(x=Sample,y=log10(cg_cov),color=Group))+#geom_jitter(alpha=0.2)+
geom_boxplot(fill=NA,outlier.shape=NA)+scale_color_manual(values=group_col)+scale_x_discrete(limits=c(row_order))+theme_minimal()

plt<-wrap_plots(plt1,plt2,plt3,plt4,ncol=1,guides='collect')
ggsave(plt,file="scaledcis.assigned_qc_boxplots.pdf",width=50,limitsize=F)

####################################################
#           Fig 1 UMAP            RNA               #
###################################################
#UMAP by celltype for RNA
library(Seurat)
rna@meta.data$Group<-group[rna@meta.data$sample]
rna_celltype_plt<-DimPlot(rna,group.by="coarse_celltype")+scale_color_manual(values=celltype_col)+theme_void()
rna_group_plt<-DimPlot(rna,group.by="Group")+scale_color_manual(values=group_col)+theme_void()

ggsave(rna_celltype_plt,file="scaledcis.rna.celltype.umap.pdf")
ggsave(rna_group_plt,file="scaledcis.rna.group.umap.pdf")


####################################################
#           Fig 1 UMAP            MET               #
###################################################
library(parallel)
library(dplyr)

options(future.globals.maxSize= 2000000*1024^2) #200gb limit for parallelizing

library(ggplot2)
library(patchwork)
p1 <- ggplot(data=obj@metadata,aes(x=umap_x,y=umap_y,color=celltype))+geom_point(size=0.25)+theme_void()+ggtitle(paste0("UMAP All cells"))+scale_color_manual(values=celltype_col)
ggsave(p1,file=paste0("07_scaledcis.umap.celltype.pdf"),width=10,height=10,limitsize=FALSE)

p2 <- ggplot(data=obj@metadata,aes(x=umap_x,y=umap_y,color=Group))+geom_point(size=0.25)+theme_void()+ggtitle(paste0("UMAP All cells"))+scale_color_manual(values=group_col)
ggsave(p2,file=paste0("07_scaledcis.umap.group.pdf"),width=10,height=10,limitsize=FALSE)

met_cg_plt<- dimFeature(obj, colorBy = mcg_pct, reduction = "umap")+scale_color_gradientn(colors = c("black", "turquoise", "gold", "red"),guide="colourbar")
ggsave(met_cg_plt,file="07_scaledcis.umap.mcg.pdf")

met_cgcov_plt <- dimFeature(obj, colorBy = log10(cg_cov), reduction = "umap") + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red"),guide="colourbar")
ggsave(met_cgcov_plt,file="07_scaledcis.umap.mcg.pdf")

```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
```


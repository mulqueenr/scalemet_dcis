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
setwd(wd)

obj<-readRDS(file="05_scaledcis.fine_celltype.amethyst.rds")
rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")
table(rna@meta.data$sample) %in% table(obj@metadata$Sample)
rna$original_rna_sample_name<-rna$sample
rna$Sample<-rna$sample

#set names to be the same in RNA data
rna@meta.data[rna@meta.data$Sample=="BCMDCIS102T-4h",]$Sample<-"BCMDCIS102T_24hTis"
rna@meta.data[rna@meta.data$Sample=="BCMDCIS35T-3h",]$Sample<-"BCMDCIS35T"
rna@meta.data[rna@meta.data$Sample=="BCMDCIS49T-24hTis",]$Sample<-"BCMDCIS49T"
rna@meta.data[rna@meta.data$Sample=="DCIS-79T",]$Sample<-"BCMDCIS79T_24hTis_DCIS"
rna@meta.data[rna@meta.data$Sample=="IDC-79T",]$Sample<-"BCMDCIS79T_24hTis_IDC"
rna@meta.data[rna@meta.data$Sample=="BCMDCIS80T",]$Sample<-"BCMDCIS80T_24hTis"
rna@meta.data[rna@meta.data$Sample=="BCMDCIS82T-24hTis",]$Sample<-"BCMDCIS82T_24hTis"
rna@meta.data[rna@meta.data$Sample=="BCMDCIS94T-24hTis",]$Sample<-"BCMDCIS94T_24hTis"   
rna@meta.data[rna@meta.data$Sample=="BCMDCIS99T",]$Sample<-"BCMDCIS99T"  
rna@meta.data[rna@meta.data$Sample=="DCIS-92T",]$Sample<-"BCMDCIS92T_24hTis"  
rna@meta.data[rna@meta.data$Sample=="HBCA-16R",]$Sample<-"BCMHBCA16R-3h"
rna@meta.data[rna@meta.data$Sample=="HBCA-19T",]$Sample<-"BCMHBCA19R-4h"
rna@meta.data[rna@meta.data$Sample=="HBCA-83L",]$Sample<-"BCMHBCA83L-3h"

saveRDS(rna,"/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")

####################################################
#           Fig 1 Sample Heatmap                  #
###################################################
met<-obj@metadata
#summarize per sample (heatmap categorical)
met<-met[!duplicated(met$Sample),]
age=setNames(nm=met$Sample,met$Age)
race_ethnicity=setNames(nm=met$Sample,met$Race_Ethnicity)
menopause=setNames(nm=met$Sample,met$Menopause)
er=setNames(nm=met$Sample,met$ER)
pr=setNames(nm=met$Sample,met$PR)
her2=setNames(nm=met$Sample,met$HER)
grade=setNames(nm=met$Sample,met$Grade)
group=setNames(nm=met$Sample,met$Group)
meta_cat<-cbind(group,age,race_ethnicity,menopause,er,pr,her2,grade)
meta_cat[which(is.na(meta_cat),arr.ind=T)]<-"N/A"
meta_cat<-as.data.frame(meta_cat)
meta_cat$age<-as.numeric(meta_cat$age)
row_order<-met %>% arrange(Group,ER,PR,HER2,Sample) %>% pull(Sample)

age_col=colorRamp2(breaks=c(min(meta_cat$age,na.rm=T),max(meta_cat$age,na.rm=T)),c("#d789d7","#2a3d66"))
class_col=c("+"="black",
            "-"="grey",
            "N/A"="white")
grade_col=c("N/A"="white",
            "G1"="#CCDF92",
            "G2"="#8A9A5B",
            "G2-int"="#6C7C59",
            "G3"="45503B")
menopause_col=c("Hysterectomy (perimenopausal)"="#E2BD6B",
                "Perimenopausal"="#E2BD6B",
                "Post-menopausal(Hysterectomy)"="#4D067B",
                "Post-menopausal"="#4D067B",
                "Pre-menopausal"="#B984DB",
                "Unknown"="white",
                "s/p hysterectomy"="#E2BD6B")
group_col=c("DCIS"="#8c86bc",
            "DCIS/LCIS"="#bfda9f",
            "HBCA"="#c6d8e8",
            "IDC"="#e37f76")
#plot metadataq
ha = rowAnnotation(
  group=meta_cat$group,
  age=meta_cat$age,
  race_ethnicity=meta_cat$race_ethnicity,
  menopause=meta_cat$menopause,
  col = list(group=group_col,
              age=age_col,
              menopause=menopause_col))
plt<-Heatmap(meta_cat[,c("er","pr","her2")],
 cluster_columns=F,cluster_rows=F,
 left_annotation=ha,row_order=row_order)
 #col=class_col)

 pdf("scaledcis.metadata.heatmap.pdf")
 print(plt)
 dev.off()

####################################################
#           Fig 1 Stacked Celltype ID             #
###################################################
#Make stacked barplot on identities per cluster

celltype_col=c(
'peri'='#c1d552',
'perivasc'='#c1d552',
'fibro1'='#7f1911',
'fibro'='#7f1911',
'fibro2'='#e791f9',
'endo'='#f0b243',
'endo2'='#d0bd4a',

'tcell'='#2e3fa3',
'bcell'='#00adea',
'myeloid1'='#00a487',
'myeloid'='#00a487',
'myeloid2'='#006455',
'plasma'='#006455',

'basal'='#7200cc',
'lumsec'='#af00af',
'lumhr'='#d8007c')

#met cell counts
met_count<-obj@metadata %>% as.data.frame() %>% dplyr::count(Sample, .drop=FALSE)
met_count$Group<-group[met_count$Sample]
plt1<-ggplot(met_cell,aes(x=Sample,fill=Group,y=n))+geom_bar(stat="identity")+theme_minimal()+scale_x_discrete(limits=c(row_order))+scale_fill_manual(values=group_col)

#met cell types
met_cell<-obj@metadata %>% as.data.frame() %>% dplyr::count(Sample,fine_celltype, .drop=FALSE)
met_cell$Group<-group[met_cell$Sample]
met_cell$fine_celltype<-factor(met_cell$fine_celltype,levels=names(celltype_col))
plt2<-ggplot(met_cell,aes(x=Sample,fill=fine_celltype,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()+scale_x_discrete(limits=c(row_order))+scale_fill_manual(values=celltype_col)

#rna cell counts
rna_counts<-rna@meta.data %>% as.data.frame() %>% dplyr::count(Sample, .drop=FALSE)
rna_counts$Group<-group[rna_counts$Sample]
plt3<-ggplot(rna_counts,aes(x=Sample,fill=Group,y=n))+geom_bar(stat="identity")+theme_minimal()+scale_x_discrete(limits=c(row_order))+scale_fill_manual(values=group_col)

#rna cell types
rna_cell<-rna@meta.data %>% as.data.frame() %>% dplyr::count(Sample,coarse_celltype, .drop=FALSE)
rna_cell$Group<-group[rna_cell$Sample]
rna_cell$coarse_celltype<-factor(rna_cell$coarse_celltype,levels=names(celltype_col))
plt4<-ggplot(rna_cell,aes(x=Sample,fill=coarse_celltype,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()+scale_x_discrete(limits=c(row_order))+scale_fill_manual(values=celltype_col)

plt<-wrap_plots(plt1,plt2,plt3,plt4,ncol=1,guides='collect')
ggsave(,file="scaledcis.assigned_celltype_barplots.pdf",width=50,limitsize=F)


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
ggsave(,file="scaledcis.assigned_qc_boxplots.pdf",width=50,limitsize=F)

####################################################
#           Fig 1 UMAP            RNA               #
###################################################
#UMAP by celltype for RNA
library(Seurat)
rna@meta.data$Group<-group[rna@meta.data$Sample]
rna_celltype_plt<-DimPlot(rna,group.by="coarse_celltype")+scale_color_manual(values=celltype_col)
rna_group_plt<-DimPlot(rna,group.by="Group")+scale_color_manual(values=group_col)

ggsave(rna_celltype_plt,file="scaledcis.rna.celltype.umap.pdf")
ggsave(rna_group_plt,file="scaledcis.rna.group.umap.pdf")


####################################################
#           Fig 1 UMAP            MET               #
###################################################
#UMAP by celltype for MET
obj<-cluster_by_windows(
    obj, 
    window_name="nakshatri_dmr_sites",
    outname="scaledcis.met.umap",
    threads=10, 
    est_dim=12, 
    neighbors=80, 
    dist=1E-8,
    k_pheno=65)

#UMAP PLOTTING
met_celltype_plt <- dimFeature(obj, colorBy = fine_celltype, reduction = "umap")+scale_color_manual(values=celltype_col)
ggsave(met_celltype_plt,file="scaledcis.met.celltype.umap.pdf")

met_group<- dimFeature(obj, colorBy = Group, reduction = "umap") +scale_color_manual(values=group_col)
ggsave(met_group,file="scaledcis.met.group.umap.pdf")

met_cg_plt<- dimFeature(obj, colorBy = mcg_pct, reduction = "umap")+scale_color_gradientn(colors = c("black", "turquoise", "gold", "red"),guide="colourbar")
ggsave(met_cg_plt,file="scaledcis.met.cgperc.umap.pdf")

met_cgcov_plt <- dimFeature(obj, colorBy = log10(cg_cov), reduction = "umap") + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red"),guide="colourbar")
ggsave(met_cgcov_plt,file="scaledcis.met.cgcov.umap.pdf")

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
setwd(wd)

obj<-readRDS(file="04_scaledcis.broad_celltype.amethyst.rds")
obj@metadata$fine_celltype<-obj@metadata$broad_celltype

#stromal
stromal<-"/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/fine_celltyping/stromal/04_scaledcis.stromal.fine_celltype.amethyst.rds"
stromal<-readRDS(stromal)
stromal_celltypes<-setNames(stromal@metadata[which(row.names(stromal@metadata) %in% row.names(obj@metadata)),]$fine_celltypes,nm=row.names(stromal@metadata))
obj@metadata[names(stromal_celltypes),]$fine_celltype<-stromal_celltypes

#immune
immune<-"/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/fine_celltyping/immune/04_scaledcis.immune.fine_celltype.amethyst.rds"
immune<-readRDS(immune)
immune_celltypes<-setNames(immune@metadata[which(row.names(immune@metadata) %in% row.names(obj@metadata)),]$fine_celltypes,nm=row.names(immune@metadata))
obj@metadata[names(immune_celltypes),]$fine_celltype<-immune_celltypes

#summarize celltype-dmr clusters over windows
obj<-dmr_and_1kb_window_gen(obj,
    prefix="fine_celltypes",
    groupBy="fine_celltype",
    threads=10,step=500)

saveRDS(obj,file="05_scaledcis.fine_celltype.amethyst.rds")
write.table(obj@metadata,file="05_scaledcis.fine_celltype.metadata.tsv",quote=F,col.names=T,row.names=T,sep="\t")

obj<-readRDS(file="01_celllines.amethyst.rds")
write.table(obj@metadata,file="01_celllines.metadata.tsv",quote=F,col.names=T,row.names=T,sep="\t")
```


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
obj<-readRDS(file="09_scaledcis.final_ploidy.amethyst.rds")
rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")

####################################################
#           Fig 1 Sample Heatmap                  #
###################################################
outdir=paste0(project_data_directory,"/","figure_1")
system(paste0("mkdir -p ", outdir))

met<-obj@metadata
#summarize per sample (heatmap categorical)
met<-met[!duplicated(met$Sample),]
met<-met[!is.na(met$Sample),]
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

age_col=colorRamp2(breaks=c(min(meta_cat$Age,na.rm=T),median(meta_cat$Age,na.rm=T),max(meta_cat$Age,na.rm=T)),c("#53FFFF","#90A2FF","#FF7BFF"))

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


#color assignment is fluor as cancer associated cell type, rest muted versions
celltype_col=c(
"perivascular"="#FF9900",
"fibroblast"="#FF0000",
"endothelial"="#FFFF66",
"unknown"="#FF6699",

"monocyte"="#99FFFF",
"macrophage"="#0066FF",
"bcell"="#0099CC",
"tcell_treg"="#99FF99",
"tcell_cd4"="#009966",
"tcell_cd8"="#66FF00",

"basal"="#990099",
"lumsec"="#CC0066",
"lumhr"="#FF00CC",
"cancer"="#00FF99")

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

 pdf(paste0(outdir,"/","scaledcis.metadata.heatmap.pdf"),width=10)
 print(plt)
 dev.off()

####################################################
#           Fig 1 Stacked Celltype ID             #
###################################################
#Make stacked barplot on identities per cluster

group<-setNames(nm=row.names(meta_cat),meta_cat[,"Group"])
row_order<-rev(row_order)

#met cell counts
met_count<-obj@metadata %>% as.data.frame() %>% dplyr::count(Sample, .drop=TRUE)
#met_count<-met_count[!is.na(met_count$Sample),]
met_count$Group<-group[met_count$Sample]
plt1<-ggplot(met_count,aes(x=Sample,fill=Group,y=n))+geom_bar(stat="identity")+theme_minimal()+scale_x_discrete(limits=c(row_order))+scale_fill_manual(values=group_col)
ggsave(plt1,file=paste0(outdir,"/","scaledcis.counts_per_sample.barplots.pdf"),width=50,limitsize=F)

#met cell types percent
met_cell<-obj@metadata %>% as.data.frame() %>% dplyr::count(Sample,celltype, .drop=TRUE)
met_cell$celltype<-factor(met_cell$celltype,levels=names(celltype_col))
plt5<-ggplot(met_cell,aes(x=Sample,fill=celltype,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()+scale_x_discrete(limits=c(row_order))+scale_fill_manual(values=celltype_col)
ggsave(plt2,file=paste0(outdir,"/","scaledcis.assigned_celltype_by_sample.barplots.pdf"),width=50,limitsize=F)

#rna cell counts
rna_counts<-rna@meta.data %>% as.data.frame() %>% dplyr::count(sample, .drop=TRUE)
rna_counts$Group<-group[rna_counts$sample]
plt3<-ggplot(rna_counts,aes(x=sample,fill=Group,y=n))+geom_bar(stat="identity")+theme_minimal()+scale_x_discrete(limits=c(row_order))+scale_fill_manual(values=group_col)
ggsave(plt3,file=paste0(outdir,"/","rna.assigned_celltype_by_sample.barplots.pdf"),width=50,limitsize=F)

#rna cell types
#matching celltype assignment to met level
#perivascular   fibroblast  endothelial      unknown     monocyte   macrophage 
#          34           40           38           29           28           41 
#       bcell   tcell_treg    tcell_cd4    tcell_cd8        basal       lumsec 
#          32           35           41           31           40           39 
#       lumhr       cancer 
#          41           27 

#note that this is based on cluster proximity, so where we think lower resolution methylation data will group cells together
rna@meta.data[rna@meta.data$fine_celltype %in% c("bcell","bcell_stress","plasma"),]$coarse_celltype<-"bcell"
rna@meta.data[rna@meta.data$fine_celltype %in% c("endo_artery","endo_capillary","endo_lymphatic","endo_TEC","endo_unknown","endo_vein"),]$coarse_celltype<-"endothelial"
rna@meta.data[rna@meta.data$fine_celltype %in% c("fibro_CAF","fibro_major","fibro_matrix","fibro_prematrix","fibro_SFRP4"),]$coarse_celltype<-"fibroblast"
rna@meta.data[rna@meta.data$fine_celltype %in% c("fibro_CAF","fibro_major","fibro_matrix","fibro_prematrix","fibro_SFRP4"),]$coarse_celltype<-"fibroblast"

rna@meta.data[rna@meta.data$fine_celltype %in% c("myeloid_3","myeloid_cycling","myeloid_DC","myeloid_macro","myeloid_mast","myeloid_neutrophil","myeloid_TAM"),]$coarse_celltype<-"macrophage"
rna@meta.data[rna@meta.data$fine_celltype %in% c("myeloid_mono"),]$coarse_celltype<-"monocyte"
rna@meta.data[rna@meta.data$fine_celltype %in% c("peri","periVSMC_unknown","VSMC"),]$coarse_celltype<-"perivascular"

rna@meta.data[rna@meta.data$fine_celltype %in% c("tcell_cd4"),]$coarse_celltype<-"tcell_cd4"
rna@meta.data[rna@meta.data$fine_celltype %in% c("tcell_cd8","tcell_gdT","tcell_nk","tcell_interferon"),]$coarse_celltype<-"tcell_cd8"
rna@meta.data[rna@meta.data$fine_celltype %in% c("tcell_treg"),]$coarse_celltype<-"tcell_treg"

rna_cell<-rna@meta.data %>% as.data.frame() %>% filter(coarse_celltype %in% names(celltype_col)) %>% dplyr::count(sample,coarse_celltype, .drop=FALSE)
rna_cell$Group<-group[rna_cell$sample]
rna_cell$coarse_celltype<-factor(rna_cell$coarse_celltype,levels=names(celltype_col))
plt4<-ggplot(rna_cell,aes(x=sample,fill=coarse_celltype,y=n))+geom_bar(position="fill",stat="identity")+theme_minimal()+scale_x_discrete(limits=c(row_order))+scale_fill_manual(values=celltype_col)

plt<-patchwork::wrap_plots(plt1,plt2,plt3,plt4,ncol=1,guides='collect')
ggsave(plt,file=paste0(outdir,"/","scaledcis_and_rna.assigned_celltype_barplots.pdf"),width=50,limitsize=F)

rna<-subset(rna,coarse_celltype %in% names(celltype_col))
Idents(rna)<-factor(rna$coarse_celltype,levels=c("perivascular","fibroblast","endothelial","tcell_treg","tcell_cd4","tcell_cd8","bcell","macrophage","monocyte","basal","lumsec","lumhr","cancer"))
genes<-c("RGS5","COL1A1","PECAM1","PTPRC","CD8A","CD8B","CD19","AIF1","KRT17","KRT15","ANKRD30A")
plt<-DotPlot(rna,features=factor(genes,levels=rev(genes)),cols=c("lightgrey","#e57528"))
ggsave(plt,file=paste0(outdir,"/","rna.gene_dotplot.pdf"),width=10,limitsize=F)

####################################################
#           Fig 1 UMAPS                          #
###################################################
#UMAP by celltype for RNA
library(Seurat)
rna@meta.data$Group<-group[rna@meta.data$sample]
rna_meta<-rna@meta.data
rna_meta<-rna_meta[row.names(rna_meta[rna_meta$coarse_celltype %in% names(celltype_col),]),]
umap_dat<-as.data.frame(cbind(rna@reductions$umap@cell.embeddings[,1],rna@reductions$umap@cell.embeddings[,2]))
colnames(umap_dat)<-c("umap_x","umap_y")
umap_dat<-umap_dat[row.names(rna_meta),]
rna_meta$umap_x<-umap_dat$umap_x
rna_meta$umap_y<-umap_dat$umap_y

rna_celltype_plt<-ggplot(rna_meta,aes(x=umap_x,y=umap_y,color=coarse_celltype))+geom_point(color="black",size=2)+geom_point(size=1,alpha=0.1)+scale_color_manual(values=celltype_col)+theme_void()
ggsave(rna_celltype_plt,file=paste0(outdir,"/","scaledcis.rna.celltype.umap.pdf"))

rna_group_plt<-ggplot(rna_meta,aes(x=umap_x,y=umap_y,color=Group))+geom_point(color="black",size=2)+geom_point(size=1,alpha=0.1)+scale_color_manual(values=group_col)+theme_void()
ggsave(rna_group_plt,file=paste0(outdir,"/","scaledcis.rna.group.umap.pdf"))



####################################################
# Box plot of %met per clone #
####################################################
met_meta_celltype<-obj@metadata %>% filter(Group=="HBCA")
met_meta_celltype$celltype<-factor(met_meta_celltype$celltype,levels=names(celltype_col))

met_meta<-obj@metadata %>% filter(!endsWith(cnv_clonename,suffix="diploid")) %>% filter(!is.na(cnv_clonename)) %>% filter(!is.na(scquantum_ploidy))  %>% filter(Group!="HBCA")
met_meta_clone<- met_meta[!duplicated(met_meta$cnv_clonename),] 
clone_order<-met_meta_clone %>% arrange(scquantum_ploidy) %>% select(cnv_clonename) %>% unlist()

met_meta$cnv_clonename<-factor(met_meta$cnv_clonename,levels=clone_order)
met_meta_clone$cnv_clonename<-factor(met_meta_clone$cnv_clonename,levels=clone_order)

plt1<-ggplot(met_meta, aes(x=cnv_clonename,y=mcg_pct,fill=Group)) +geom_boxplot(outlier.shape = NA) +ylim(c(50,100))
plt2<-ggplot(met_meta_clone, aes(x=cnv_clonename,y=scquantum_ploidy,fill=Group)) + geom_bar(stat="identity")
plt3<-ggplot(met_meta_celltype, aes(x=celltype,y=mcg_pct,fill=celltype)) +geom_boxplot(outlier.shape = NA) +ylim(c(50,100)) +scale_fill_manual(values=celltype_col)
ggsave(plt1/plt2/plt3,file=paste0(outdir,"/","scaledcis.cnv_clonename.percmet.pdf"))



####################################################
#           Fig 1 UMAP            MET               #
###################################################

library(amethyst)
library(ggplot2)
library(patchwork)

#plot all methylation cells
met_celltype_plt<-ggplot(data=obj@metadata,aes(x=umap_x,y=umap_y,color=celltype_color))+geom_point(size=2,color="black")+geom_point(size=1,alpha=0.5)+scale_color_identity()+theme_void()
ggsave(plt1,file=paste0(outdir,"/","scaledcis.methylation.celltype.umap.pdf"))

met_group_plt<-ggplot(data=obj@metadata,aes(x=umap_x,y=umap_y,color=Group))+geom_point(size=2,color="black")+geom_point(size=1,alpha=0.5)+scale_color_manual(values=group_col)+theme_void()
ggsave(plt2,file=paste0(outdir,"/","scaledcis.methylation.group.umap.pdf"))

percent_met_col<-colorRampPalette(c(min(obj@metadata$mcg_pct,na.rm=T),max(obj@metadata$mcg_pct,na.rm=T)),c("grey","#FF00FF"))

met_percmet_plot<-ggplot(data=obj@metadata,aes(x=umap_x,y=umap_y,color=mcg_pct))+geom_point(size=2,color="black")+geom_point(size=1,alpha=0.5)+scale_colour_gradient2(high="black",mid="white",low="#FF00FF",midpoint=75)+theme_void()

ggsave(plt3,file=paste0(outdir,"/","scaledcis.methylation.percmet.umap.pdf"))

plt<-patchwork::wrap_plots(met_celltype_plt,met_group_plt,met_percmet_plot,rna_celltype_plt,rna_group_plt,ggplot(),ncol=3,guides='collect')
ggsave(plt,file=paste0(outdir,"/","scaledcis.umap.pdf"),width=14.5,height=9)
ggsave(plt,file=paste0(outdir,"/","scaledcis.umap.png"),width=14.5,height=9,dpi = 600)


############################################################
# RNA and methylation gene dosage
############################################################

#subset to clones
rna_cnv<-subset(rna, cells=row.names(rna@meta.data)[!endsWith(rna@meta.data$rna_clonename, suffix="diploid")])
rna_cnv@meta.data[rna_cnv@meta.data$rna_clonename=="BCMDCIS124T_c2=",]$rna_clonename<-"BCMDCIS124T_c2"
met_cnv<-subsetObject(obj,cells=row.names(obj@metadata)[!endsWith(obj@metadata$cnv_clonename, suffix="diploid")])
met_cnv<-subsetObject(met_cnv,cells=row.names(met_cnv@metadata)[!is.na(met_cnv@metadata$cnv_clonename)])
 
#correlate rna expression to integer copy number per region
#genomic ranges to overlap with cnv sites and ref genes, assign integer copy per cell
rna_clones<-setNames(nm=row.names(rna_cnv@meta.data),rna_cnv@meta.data$rna_clonename)
met_clones<-setNames(nm=row.names(met_cnv@metadata),met_cnv@metadata$cnv_clonename)

unique(rna_clones)[!(unique(rna_clones) %in% unique(met_clones))]
#clones in RNA not in met "BCMDCIS07T_c1" "BCMDCIS74T_c5"   "ECIS57T_c1"
#BCMDCIS07T_c1 not in met
#BCMDCIS74T_c5 looks different than met clones
#ECIS57T_c1 low cell count in met

#shared clone list
common_clones<-intersect(unique(rna_clones),unique(met_clones))
met_cnv<-subsetObject(met_cnv,cells=row.names(met_cnv@metadata)[met_cnv@metadata$cnv_clonename %in% common_clones])
rna_cnv<-subset(rna_cnv, cells=row.names(rna_cnv@meta.data)[rna_cnv@meta.data$rna_clonename %in% common_clones])

#make pseudocells per clone????
rna_clone<-AggregateExpression(rna_cnv,group.by="rna_clonename",return.seurat=TRUE)
rna_clone<-JoinLayers(rna_clone)
rna_clone<-LayerData(rna_clone,assay="RNA",layer="data") #get raw counts for summarizing over feature windows
colnames(rna_clone)<-gsub(colnames(rna_clone),pattern="-",replace="_")

#get cnv windows
cnv_win <- row.names(obj@genomeMatrices$cg_cnv_segments)
cnv_win<- GRanges(data.frame(
    "chr"=unlist(lapply(strsplit(cnv_win,"_"),"[",1)),
    "start"=unlist(lapply(strsplit(cnv_win,"_"),"[",2)),
    "end"=unlist(lapply(strsplit(cnv_win,"_"),"[",3))))

#get list of genes per window
rna_win <- obj@ref %>% 
  filter(type=="gene") %>% 
  filter(gene_type %in% c("protein_coding")) %>% #,"lncRNA")) %>% 
  filter(gene_name %in% row.names(rna_clone))
rna_win <- rna_win[!duplicated(rna_win$gene_name),]
row.names(rna_win) <- rna_win$gene_name

rna_feat<-intersect(row.names(rna_win),row.names(rna_clone))
rna_win<-rna_win[rna_feat,]
rna_clone<-rna_clone[rna_feat,]
rna_clone <- as.data.frame(rna_clone)
rna_clone$chr<-unlist(lapply(rna_win$seqid,as.character))
rna_clone$start<-as.character(rna_win$start)
rna_clone$end<-as.character(rna_win$end)

#run a calcSmooth function to merge over clones
cluster1kbwindows <- calcSmoothedWindows(met_cnv, 
                                         type = "CG", 
                                         threads = 50,
                                         step = 1000, # change to 500 for real data unless you have really low coverage
                                         smooth = 1,
                                         genome = "hg38",
                                         index = "chr_cg",
                                         groupBy = "cnv_clonename",
                                         returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                         returnPctMatrix = TRUE)
saveRDS(cluster1kbwindows,file=paste0(outdir,"/","scaledcis.methylation.cnv_clone.1kb.windows.rds"))

met_clone<-cluster1kbwindows[["pct_matrix"]]
met_clone_keep<-c("chr","start","end",colnames(met_clone)[which(colnames(met_clone) %in% common_clones)])
met_clone<-met_clone %>% select(all_of(met_clone_keep))

#split RNA by CNV window overlaps
rna_cnv_overlaps <- findOverlaps(cnv_win, GRanges(rna_clone %>% select(all_of(c("chr","start","end")))), select="all")
met_cnv_overlaps <- findOverlaps(cnv_win, GRanges(met_clone %>% select(all_of(c("chr","start","end")))), select="all")

met_clone <- met_clone[met_cnv_overlaps@to,]
met_clone$win <- met_cnv_overlaps@from

rna_clone<-rna_clone[rna_cnv_overlaps@to,]
rna_clone$win <- rna_cnv_overlaps@from

#add integer copy per met/rna feature
cnv<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/copykit/all_samples.aneuploid_only.copykit.Rds")
cnv_clone<-cnv@consensus[,common_clones]
row.names(cnv_clone)<-1:nrow(cnv_clone) #matching windownames

#methylation features to cnv integer number
met_win_cor<-mclapply(1:nrow(met_clone),
  function(i){
    met<-as.numeric(met_clone[i] %>% select(all_of(common_clones)))
    cnv<-as.numeric(cnv_clone[met_clone[i]$win,common_clones])
    if(sum(is.na(met))<length(met)){
      met_cor<-cor(cnv,met,use="complete",method="pearson")
      cnv_var<-var(cnv,na.rm=T)
      cnv_mean<-mean(cnv,na.rm=T)
      met_mean<-mean(met,na.rm=T)
    return(c(i,cnv_var,cnv_mean,met_mean,met_cor))
    }
},mc.cores=300)


cor_dat<-as.data.frame(do.call("rbind",met_win_cor))
colnames(cor_dat)<-c("win","cnv_var","cnv_mean","met_mean","met_cor")
cor_dat<-cor_dat[cor_dat$cnv_var>0.25,]

plt<-ggplot(cor_dat,aes(x=as.numeric(met_mean),y=as.numeric(met_cor),color=as.numeric(cnv_mean),alpha=abs(met_cor)))+geom_point(size=0.2)+theme_minimal()+scale_color_gradient(low="blue",high="red")
ggsave(plt,file=paste0(outdir,"/","test.met_cor.png"),width=10,height=10)


#methylation features to cnv integer number
rna_win_cor<-mclapply(1:nrow(rna_clone),
  function(i){
    met<-as.numeric(rna_clone[i,common_clones])
    cnv<-as.numeric(cnv[rna_clone[i,]$win,common_clones])
    met_cor<-cor(cnv,met,use="complete",method="spearman")
    cnv_var<-var(cnv)
    cnv_mean<-mean(cnv)
    met_mean<-mean(met)
  return(c(i,
    cnv_var,
    cnv_mean,
    met_mean,
    rna_mean,
    met_cor,
    rna_cor))
},mc.cores=50)


met_clone<-split(met_clone,f=met_clone$win)
rna_clone<-split(rna_clone,f=rna_clone$win)

#win
win<-intersect(names(met_clone),names(rna_clone))

win_cor<-mclapply(win,function(i){
  cnv<-as.numeric(cnv[i,common_clones])
  met<-as.numeric(colMeans(met_clone[[i]]%>% select(all_of(common_clones)),na.rm=T))
  rna<-as.numeric(colMeans(rna_clone[[i]]%>% select(all_of(common_clones)),na.rm=T))
  met_cor<-cor(cnv,met,use="complete",method="spearman")
  rna_cor<-cor(cnv,rna,use="complete",method="spearman")
  cnv_var<-var(cnv)
  cnv_mean<-mean(cnv)
  met_mean<-mean(met)
  rna_mean<-mean(met)
  return(c(i,
    cnv_var,
    cnv_mean,
    met_mean,
    rna_mean,
    met_cor,
    rna_cor))
},mc.cores=50)

cor_dat<-as.data.frame(do.call("rbind",win_cor))
colnames(cor_dat)<-c("win","cnv_var","cnv_mean","met_mean","rna_mean","met_cor","rna_cor")
cor_dat<-cor_dat[cor_dat$cnv_var>0.5,]

plt<-ggplot(cor_dat,aes(x=as.numeric(met_cor),y=as.numeric(rna_cor),color=as.numeric(cnv_var)))+geom_point()
ggsave(plt,file=paste0(outdir,"/","test.dotplot.png"))


cnv_clone_tmp <- obj@genomeMatrices$aneuploid_cnv_clonelevel[names(met_clones[names(met_clones==clone)][1])]
rna_clone$cnv<-as.numeric(cnv_clone_tmp[rna_clone_tmp$win,])
met_clone$cnv<-as.numeric(cnv_clone_tmp[met_clone_tmp$win,])

#plot methylation-cnv correlation by mean cnv
met_clone
ref <- ref[rna_cnv_overlaps@to,]
ref$window <- gene_overlaps@from

#correlation between methylation clone and integer copy number
#per clone doesnt make much sense, should do it per feature instead

cnv_cor<-lapply(common_clones,function(clone){
  if(clone %in% colnames(met_clone) & clone %in% colnames(rna_clone)){
  #assign window per met clone methylation
  met_clone_tmp <- met_clone %>% select(all_of(c("chr","start","end",clone))) %>% GRanges()
  met_clone_tmp <- met_clone_tmp[met_cnv_overlaps@to,]
  met_clone_tmp$win <- met_cnv_overlaps@from

  #assign window per rna feature
  rna_clone_tmp<-rna_clone %>% select(all_of(c("chr","start","end",clone))) %>% GRanges()
  rna_clone_tmp<-rna_clone_tmp[rna_cnv_overlaps@to,]
  rna_clone_tmp$win <- rna_cnv_overlaps@from
  

  met_cor = cor(unlist(mcols(met_clone_tmp)[clone]),as.numeric(met_clone_tmp$cnv),method="spearman",use="complete")
  rna_cor = cor(unlist(mcols(rna_clone_tmp)[clone]),as.numeric(rna_clone_tmp$cnv),method="spearman",use="complete")
  return(c(clone,met_cor,rna_cor))}
  #add overlap between rna and methylation?

})

as.data.frame(do.call("rbind",cnv_cor))

############################################################
#################### Methylation Tracks ####################
############################################################

#read in bigwig files,

#modified from histograM
histograModified <- function(obj,
    genes = NULL,
    matrix,colors = NULL,
    ncol = length(genes),
    trackOverhang = 5000,arrowOverhang = 3000,trackScale = 1.5,arrowScale = 0.025,colorMax = 100,
    legend = TRUE,removeNA = TRUE,
    order = NULL,baseline = "mean",
    promoter_focus=FALSE,
    cgisland=NULL) {

    if (!is.null(colors)) {pal <- colors} else {pal <- c("#FF0082", "#dbdbdb", "#cccccc", "#999999")}
    
    genes<-genes[genes %in% obj@ref$gene_name] #ensure genes are in ref
    p <- vector("list", length(genes)) # empty plot list

    for (i in 1:length(genes)) {
            group_count<-ncol(obj@genomeMatrices[[matrix]])-3 #-3 for chr start end
            ngroup<-vector("list",group_count+1)  #+1 for gene track
            if (!is.null(order)) {
                names(ngroup)<-c("gene_track",order)
            }else{
                names(ngroup)<-c("gene_track",colnames(obj@genomeMatrices[[matrix]])[4:ncol(obj@genomeMatrices[[matrix]])])
            }

            ref <- obj@ref |> dplyr::filter(gene_name == genes[i])
            aggregated <- obj@genomeMatrices[[matrix]]

            toplot <- aggregated[c((aggregated$chr == ref$seqid[ref$type == "gene"] &
                                aggregated$start > (ref$start[ref$type == "gene"] - trackOverhang) &
                                aggregated$end < (ref$end[ref$type == "gene"] + trackOverhang))), ]

            promoter <- ref |> dplyr::filter(type == "gene") |>
                  dplyr::mutate(promoter_start = ifelse(strand == "+", (start - 2000), (end + 2000)),
                  promoter_end = ifelse(strand == "+", (promoter_start+2000), (promoter_start-2000)))
            
            if(promoter_focus){
              if(ref[ref$type == "gene",]$strand[1] == "+"){
                plot_left<-promoter$promoter_start-trackOverhang
                plot_right<-promoter$promoter_end+trackOverhang
              }else{
                plot_right<-promoter$promoter_start+trackOverhang
                plot_left<-promoter$promoter_end-trackOverhang
              }
            }
            exon <- ref |> dplyr::filter(type == "exon")
            
            trackHeight <- group_count * trackScale
            toplot <- tidyr::pivot_longer(toplot, cols = c(4:ncol(toplot)), names_to = "group", values_to = "pct_m") |> dplyr::rowwise() |> dplyr::mutate(middle = mean(c(start, end), na.rm = TRUE))
        
        if (removeNA) {toplot <- toplot |> dplyr::filter(!is.na(pct_m))}

        if (baseline == "mean") {
        glob_m <- data.frame(group = colnames(aggregated[, 4:ncol(aggregated)]),glob_m = colMeans(aggregated[, 4:ncol(aggregated)], na.rm = T))
        toplot <- dplyr::left_join(toplot, glob_m, by = "group")}

        if (!is.null(order)) {toplot$group <- factor(toplot$group, levels = order)}
        
        if(!is.null(cgisland)) {cgi <- cgisland[cgisland$chr == as.character(toplot$chr[1]) &
                            cgisland$start > min(toplot$start) &
                            cgisland$end < max(toplot$end),]}

        for (j in names(ngroup)[2:length(ngroup)]){
            toplot_sub<-toplot[toplot$group==j,]
            # plotting histogram
            p[[i]][[j]] <- ggplot2::ggplot() +
            {if (baseline == 0) ggplot2::geom_col(data = toplot_sub, ggplot2::aes(x = middle, y = pct_m, fill = pct_m), width = mean(toplot_sub$end - toplot_sub$start))} +
            {if (baseline == "mean")ggplot2::geom_rect(data = toplot_sub, ggplot2::aes(xmin = start, xmax = end, ymin = glob_m, ymax = pct_m, fill = pct_m))} + 
            {if(promoter_focus) ggplot2::coord_cartesian(xlim=c(plot_left,plot_right),expand=FALSE)}+
            ggplot2::scale_fill_gradientn(colors = pal, limits = c(0,colorMax), oob = scales::squish) + theme_minimal()+ ylab(j) + ggplot2::theme(legend.position="none") +
            ylim(c(0,100))+
            ggplot2::theme(panel.background = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank(), panel.grid.major.y = element_line(color = "#dbdbdb", linetype = "dashed"))
            }
        #last track is annot
        p[[i]][["gene_track"]] <-ggplot()+
        ggplot2::geom_rect(fill = "pink",alpha=0.8, data = promoter, ggplot2::aes(xmin = promoter_start, xmax = promoter_end, ymin = 0, ymax = trackHeight/4)) +
        ggplot2::geom_rect(alpha=0.8,fill = "gray", data = exon,ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = trackHeight/8)) +
        ggplot2::geom_rect(alpha=0.8,fill = "green", data = cgi, ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = trackHeight/4)) +
        ggplot2::geom_segment(data = ref, color=ifelse(ref$strand == "+","gray","black"),
                            aes(x = ifelse(strand == "+", (min(start) - arrowOverhang), (max(end)) + arrowOverhang),
                            y = trackHeight/8,
                            xend = ifelse(strand == "+", (max(end) + arrowOverhang), (min(start)) - arrowOverhang),
                            yend = trackHeight/8), arrow = arrow(length = unit(trackHeight/40, "cm"))) + 
                            labs(title=paste(genes[i]),subtitle=paste(toplot$chr[1],toplot$start[1],"-",toplot$end[nrow(toplot)]),size=3) +
                            theme(axis.title.x = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank(),
                            axis.text.y = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.border = element_blank(),
                            panel.background = element_blank(),
                            title=element_text(size=8))+
                            theme(legend.position="none") #+xlim(c(xmin_pos,xmax_pos))

    if(promoter_focus){ #focus to plot region to promoter + track overhang for either direction
      p[[i]][["gene_track"]]<-p[[i]][["gene_track"]]+ggplot2::coord_cartesian(xlim=c(plot_left,plot_right),expand=FALSE)
    }
    p[[i]]<-wrap_plots(p[[i]],guides='collect')+plot_layout(ncol=1)
    }
    all_genes_plot<-wrap_plots(p)+plot_layout(nrow=1)
    return(all_genes_plot)
}

#prepare cgi
cgisland="/data/rmulqueen/projects/scalebio_dcis/ref/cpgIslandExt.bed"
cgi<-rtracklayer::import(cgisland)
cgi<-as.data.frame(cgi)
colnames(cgi)<-c("chr","start","end","strand")


plt<-histograModified(obj=obj,matrix="cg_celltype_tracks",genes=c("PTPRC"),colors=celltype_col,cgisland=cgi)
ggsave(plt,file=paste0(outdir,"/","scaledcis.markertest.pdf"),width=14.5,height=9)


#plot stromal
stromal_sub<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/fine_celltyping/stromal_50kb/06_scaledcis.stromal_celltype_finecelltyping.amethyst.rds")
stromal_sub@metadata$celltype_color<-celltype_col[stromal_sub@metadata$celltype]
plt1<-ggplot(data=stromal_sub@metadata,aes(x=umap_x,y=umap_y,color=celltype_color))+geom_point(size=2,color="black")+geom_point(size=1,alpha=0.8)+scale_color_identity()+theme_minimal()
ggsave(plt1,file="08_scaledcis.final_celltype.stromal.umap.pdf")

#plot immune
immune_sub<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/fine_celltyping/immune_50kb/06_scaledcis.immune_50kb_finecelltyping.amethyst.rds")
immune_sub@metadata$celltype<-unlist(lapply(immune_sub@metadata$celltype,function(i){ifelse(i=="macro","macrophage",i)}))
saveRDS(immune_sub,file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/fine_celltyping/immune_50kb/06_scaledcis.immune_50kb_finecelltyping.amethyst.rds")

immune_sub@metadata$celltype_color<-celltype_col[immune_sub@metadata$celltype]
plt1<-ggplot(data=immune_sub@metadata,aes(x=umap_x,y=umap_y,color=celltype_color))+geom_point(size=2,color="black")+geom_point(size=1,alpha=0.8)+scale_color_identity()+theme_minimal()
ggsave(plt1,file="08_scaledcis.final_celltype.immune.umap.pdf")


saveRDS(dat,file="08_scaledcis.final_celltype.amethyst.rds") #just saving because of the macrophage name correction

##########################################################################################

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


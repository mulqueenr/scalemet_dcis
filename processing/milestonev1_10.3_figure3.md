#Plotting copy number calling and methylation clone comparisons

```R
library(amethyst)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(GenomicRanges)
#source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"

#read in object from directory
processing_folder="10_fig_plots"
wd=paste(sep="/",project_data_directory,processing_folder)
outdir=paste0(wd,"/","figure_3")

system(paste0("mkdir -p ",wd))
setwd(wd)

obj<-readRDS(file=paste(project_data_directory,"03_fine_celltyping","03_scaledcis.final_celltypes.amethyst.rds",sep="/"))
rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")

```

Set up color palettes

```R
class_col2=c("+"="black",
            "-"="grey",
            "N/A"="white",
            "Equivocal FISH pending"="#ebe8c7")

class_col=c("+"="black",
            "-"="grey",
            "N/A"="white")            
group_col=c("DCIS"="#278192",
            "HBCA"="#20223E",
            "Synchronous"="#00B089",
            "IDC"="#8FF7BD")


#color assignment is fluor as cancer associated cell type, rest muted versions
celltype_col=c(
#orange and reds
"pericyte"="#F1D302",
"fibroblast"="#780000",
"endothelial"="#F86624",

#blues
"myeloid"="#0B5563",
"bcell"="#98C1D9",
"tcell"="#43BCCD",

#purple and high contract green
"basal"="#f72585",
"lumsec"="#C490D1",
"lumhr"="#412854",
"cancer"="#99ffd3")


arm_col=c("p"="grey","q"="darkgrey")
band_col=c("acen"="#99746F","gneg"="white","gpos100"="black","gpos25"="lightgrey","gpos50"="grey","gpos75"="darkgrey","gvar"="#446879")
dip_cov=colorRamp2(c(0.5,1,1.5), c("white","grey","black"))
ploidy_col=colorRamp2(c(1,2,4), c("blue","white","darkred"))
                
percent_met_col<-colorRamp2(c(100,75,60),c("black","white","#FF00FF"))
           
reads_col=colorRamp2(c(min(log10(meta$unique_reads)),
                        max(log10(meta$unique_reads))),
                        c("white","black"))

int_col=c("0"="#053061","1"="#4393c3","2"="#f7f7f7","3"="#f4a582","4"="#b2182b","5"="#67001f","6"="#3d0229")
       

#######################################
############# Plotting #############
#######################################

cnv<-obj@genomeMatrices$scquantum_cnv[,!startsWith(prefix="metadata_",colnames(obj@genomeMatrices$scquantum_cnv))]
all_cnv_cells<-row.names(obj@metadata[!is.na(obj@metadata$cnv_clonename_500kb) & !endsWith(obj@metadata$cnv_clonename_500kb,suffix="_diploid"),])

cnv<-cnv[,colnames(cnv) %in% all_cnv_cells]
meta<-obj@metadata[colnames(cnv),]

cnv_ranges<-GRanges(data.frame(
    chr=unlist(lapply(strsplit(row.names(obj@genomeMatrices$scquantum_cnv),"_"),"[",1)),
    start=unlist(lapply(strsplit(row.names(obj@genomeMatrices$scquantum_cnv),"_"),"[",2)),
    end=unlist(lapply(strsplit(row.names(obj@genomeMatrices$scquantum_cnv),"_"),"[",3))))

cnv_ranges<-GRanges(cbind(as.data.frame(cnv_ranges),
                        obj@genomeMatrices$scquantum_cnv[,grep(pattern="metadata_",colnames(obj@genomeMatrices$scquantum_cnv))]))
                 
percent_met_col=colorRamp2(c(min(meta$mcg_pct),mean(meta$mcg_pct),
                        max(meta$mcg_pct)),
                        c("#FF00FF","white","black"))

clone_col=setNames(nm=unique(as.character(meta$cnv_clonename)),
                        colorRampPalette(brewer.pal(9, "Spectral"))(length(unique(as.character(meta$cnv_clonename)))))

#make column annotations
ha = rowAnnotation(
    group=meta$Group,
    er=meta$ER,
    pr=meta$PR,
    her2=meta$HER2,
    reads=log10(meta$unique_reads),
    cg_perc=meta$mcg_pct,
    ploidy=meta$cnv_scquantum_ploidy,
    clone=as.character(meta$cnv_clonename_500kb),
    col= list(
        group=group_col,
        er=class_col,
        pr=class_col,
        her2=class_col2,
        reads=reads_col,
        cg_perc=percent_met_col,
        ploidy=ploidy_col,
        clone=clone_col))

column_ha = HeatmapAnnotation(
    mappability=cnv_ranges$metadata_diploid_cov,
    arm = cnv_ranges$metadata_arm,
    band = cnv_ranges$metadata_stain,
    col=list(mappability=dip_cov,arm=arm_col,band=band_col))

print("Plotting heatmap...")
plt1<-Heatmap(t(cnv),
            col=int_col,
            cluster_columns=FALSE,
            cluster_rows=TRUE,
            show_row_names = FALSE, row_title_rot = 0,
            show_column_names = FALSE,
            cluster_row_slices = TRUE,
            top_annotation = column_ha, left_annotation = ha,
            name="Single cell copy number",
            row_split=meta[row.names(t(cnv)),]$cnv_clonename_500kb,
            column_split=cnv_ranges@seqnames,
            border = TRUE)

pdf_outname=paste0(outdir,"aneuploid.integer.heatmap.pdf")
pdf(pdf_outname,width=40,height=20)
print(plt1)
dev.off()

#ploidy by met (order clones by ploidy)
#wgd vs not : running (use DMR groups)
#met window by cnv : running


#```bash
#singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
#source activate
#update data.table
#```

# To Do: Also try applying Bartlett's test for differential variability instead of just diff mean
see DNA methylation outliers in normal breast tissue  identify field defects that are enriched in cancer 2015


```R
#install.packages("data.table")
#ensure 0.9 version of amethyst
#system('wget https://github.com/lrylaarsdam/amethyst/releases/download/v0.0.0.9000/amethyst_0.0.0.9000.tar.gz')
#system('tar -xf amethyst_0.0.0.9000.tar.gz')
#devtools::install("./amethyst")
#source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
#running locally

set.seed(111)
library(amethyst)
library(data.table)
library(dplyr)
library(msigdbr)
library(fgsea)
library(ggplot2)
library(data.table)
library(GenomicRanges)
library(GeneNMF)
library(ComplexHeatmap)
library(parallel)

options(future.globals.maxSize= 200000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)
obj<-readRDS(file="07_scaledcis.cnv_clones.amethyst.rds")


```

Running through all the 500bp window generation split by:
- Celltype (fine_celltype)
- CNV based clones (cnv_csampllonenames)
- CNV based clones against all lumhr (cnv_clonename_alllumhr)
- Celltype + Diagnosis (Group + fine_celltype)

Saving to folders for further processing.
For CNV based clones:
1. loop through cnv based clones samples to do DMR per sample per clone.
2. group together clones that have shared features (e.g. 1q amp) and compare to others
3. Plot DMR regions across genome per clone (see if clustered around major CNV changes or in trans)

For celltype by diagnosis
1. DMRs for TME cells per diagnosis (especially SYNCH vs non-SYNCH DCIS, and normal vs DCIS)

```R
##################FUNCTIONS##################


testDMR <- function(sumMatrix, eachVsAll = TRUE, comparisons = NULL, nminTotal = 3,nminGroup = 3) {

  if (!eachVsAll && is.null(comparisons)) {
    stop("Please either specify eachVsAll = TRUE or provide a data frame of comparisons to make.")
  }

  # filter counts table
  data.table::setDT(sumMatrix)
  counts <- data.table::copy(sumMatrix)
  counts <- counts[rowSums(counts[, .SD, .SDcols = patterns("_c$|_t$")], na.rm = TRUE) >= nminTotal]

  # fast fisher's exact test developed by @zellerivo; see https://github.com/al2na/methylKit/issues/96
  fast.fisher <- function (
    cntg_table) {
    q <- cntg_table[1, 1]
    m <- cntg_table[1, 1] + cntg_table[2, 1]
    n <- cntg_table[1, 2] + cntg_table[2, 2]
    k <- cntg_table[1, 1] + cntg_table[1, 2]
    pval_right <- phyper(q = q, m = m, n = n, k = k, lower.tail = FALSE) +
      (0.5 * dhyper(q, m, n, k))
    pval_left <- phyper(q = q - 1, m = m, n = n, k = k, lower.tail = TRUE) +
      (0.5 * dhyper(q, m, n, k))
    return(ifelse(test = pval_right > pval_left, yes = pval_left *
                    2, no = pval_right * 2))
  }

  if (is.null(comparisons)) {
    # get unique groups
    groups <- as.list(sub("_c$", "", colnames(sumMatrix)[grep("_c$", colnames(sumMatrix))]))

    for (gr in groups) {
      m_c <- paste0(gr, "_c") # m = member
      m_t <- paste0(gr, "_t")

      nm_c <- setdiff(grep("_c$", colnames(counts), value = TRUE), m_c) # nm = nonmember
      nm_t <- setdiff(grep("_t$", colnames(counts), value = TRUE), m_t)

      counts <- counts[, `:=`(
        member_c = get(paste0(gr, "_c")),
        member_t = get(paste0(gr, "_t")),
        nonmember_c = rowSums(.SD[, mget(nm_c)]),
        nonmember_t = rowSums(.SD[, mget(nm_t)])
      )]

      # don't test where the minimum observations per group is not met
      counts <- counts[member_c + member_t <= nminGroup | nonmember_c + nonmember_t <= nminGroup, c("member_c", "member_t", "nonmember_c", "nonmember_t") := .(NA, NA, NA, NA)]

      # apply fast fishers exact test
      counts <- counts[, paste0(gr, "_all_pval") := apply(.SD, 1, function(x) fast.fisher(matrix(x, nrow = 2, byrow = TRUE))), .SDcols = c("member_c", "member_t", "nonmember_c", "nonmember_t")]
      counts <- counts[, paste0(gr, "_all_logFC") := round(log2((member_c / (member_c + member_t)) / (nonmember_c / (nonmember_c + nonmember_t))), 4)]
      counts <- counts[, c("member_c", "member_t", "nonmember_c", "nonmember_t") := NULL] # this line used to be outside the loop in < v1.0.2, causing nonmember variable buildup :(
      cat(paste0("Finished group ", gr, "\n"))
    }

  } else if (!is.null(comparisons)) {
    for (i in 1:nrow(comparisons)) {
      m <- unlist(strsplit(comparisons[i, "A"], ','))
      nm <- unlist(strsplit(comparisons[i, "B"], ',', fixed = FALSE))
      name <- comparisons[i, "name"]

      m_c <- paste0(m, "_c") # m = member
      m_t <- paste0(m, "_t")

      nm_c <- paste0(nm, "_c") # n = nonmember
      nm_t <- paste0(nm, "_t")

      counts <- counts[, `:=`(
        member_c = rowSums(.SD[, mget(m_c)]),
        member_t = rowSums(.SD[, mget(m_t)]),
        nonmember_c = rowSums(.SD[, mget(nm_c)]),
        nonmember_t = rowSums(.SD[, mget(nm_t)])
      )]

      # don't test where the minimum observations per group is not met
      counts <- counts[member_c + member_t <= nminGroup | nonmember_c + nonmember_t <= nminGroup, c("member_c", "member_t", "nonmember_c", "nonmember_t") := .(NA, NA, NA, NA)]

      # apply fast fishers exact test
      counts <- counts[, paste0(name, "_pval") := apply(.SD, 1, function(x) fast.fisher(matrix(x, nrow = 2, byrow = TRUE))), .SDcols = c("member_c", "member_t", "nonmember_c", "nonmember_t")]
      counts <- counts[, paste0(name, "_logFC") := round(log2((member_c / (member_c + member_t)) / (nonmember_c / (nonmember_c + nonmember_t))), 4)]
      counts <- counts[, c("member_c", "member_t", "nonmember_c", "nonmember_t") := NULL]
      cat(paste0("\nFinished testing ", name, ": ", paste0(m, collapse = ", "), " vs. ", paste0(nm, collapse = ", ")))
    }
  }
  return(counts)
}

#gsea of top DMR genes
gsea_enrichment_againstlumhr<-function(dmrs,species="human",
                          category="C3",
                          subcategory="TFT:GTRD",
                          out_setname="TFT",sample_name=sample_name,
                          outdir=outdir,
                          obj=obj){

  pathwaysDF <- msigdbr(species=species, 
                        collection=category, 
                        subcollection = subcategory)

  #limit pathways to genes in our data
  pathwaysDF<-pathwaysDF[pathwaysDF$ensembl_gene %in% unlist(lapply(strsplit(unique(obj@ref$gene_id),"[.]"),"[",1)),]
  
  pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

  #run plotting per group order by logFC (higher logFC is hyper, lower logFC is hypo)
  plt_list<-lapply(unique(dmrs$type),function(group1){
    #treat multiple gene overlaps as same logFC
    #set -Inf to -3 and Inf to 3
    group_features<-dmrs %>%
      dplyr::filter(type == group1) %>%
      dplyr::filter(dmr_padj<0.05) %>% 
      dplyr::filter(gene_names!="NA") %>% 
      dplyr::arrange(dmr_logFC) %>%
      dplyr::select(gene_names, dmr_logFC) %>%
      tidyr::separate_rows(gene_names) %>%
      dplyr::mutate(across(where(is.numeric), ~ replace(., .==-Inf, -3))) %>%
      dplyr::mutate(across(where(is.numeric), ~ replace(., .==Inf, 3)))

  ranks<-setNames(nm=group_features$gene_names,group_features$dmr_logFC)
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = ranks,
                    minSize  = 5,
                    nproc = 1)
  saveRDS(fgseaRes,file=paste0(outdir,"/",paste0(sample_name,".",group1,".",out_setname,".GSEA_enrichment.clone_lumhrall.rds")))

  topPathwaysUp <- fgseaRes %>% filter(ES > 0) %>% slice_max(NES,n=10) %>% dplyr::select(pathway)
  topPathwaysDown <- fgseaRes %>% filter(ES < 0) %>% slice_max(abs(NES),n=10) %>% dplyr::select(pathway)
  topPathways <- unlist(c(topPathwaysUp, rev(topPathwaysDown)))

  plt1<-plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam=0.5)+
        theme(axis.text.y = element_text( size = rel(0.2)),
        axis.text.x = element_text( size = rel(0.2)))

  # only plot the top 20 pathways NES scores
  nes_plt_dat<-rbind(
    fgseaRes  %>% slice_max(NES,n= 10),
    fgseaRes  %>% slice_min(NES,n= 10))

  plt2<-ggplot(nes_plt_dat, aes(reorder(pathway, NES), NES)) +
    geom_point(aes(color= padj,size=size)) + scale_color_gradient2(low="darkred",high="grey",mid="red",midpoint=0.05)+
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score \n Compared to normal LumHR \n (Hyper<->Hypo)",
        title="Hallmark pathways NES from GSEA") + 
    theme_minimal()+scale_fill_identity()+ggtitle(paste(group1,out_setname))+ylim(c(5,-5))
  return(plt2)
  })

  patchwork::wrap_plots(plt_list,ncol=1)

}

plot_gsea<-function(obj, collapsed_dmrs,
                    sample_name,
                    outdir){
  print(paste("Loading DMRS for sample:",sample_name))

  #run gsea enrichment on different sets
  print(paste("Calculating TF Binding Enrichment Compared to LumHR"))
  tft_plt<-gsea_enrichment_againstlumhr(species="human",
              category="C3",
              subcategory="TFT:GTRD",
              out_setname="TFT",sample_name=sample_name,
              dmrs=collapsed_dmrs,obj=obj,outdir=outdir)

  print(paste("Calculating Position Enrichment Compared to LumHR"))
  position_plt<-gsea_enrichment_againstlumhr(species="human",
              category="C1",
              subcategory=NULL,
              out_setname="position",sample_name=sample_name,
              dmrs=collapsed_dmrs,obj=obj,outdir=outdir)

  print(paste("Calculating Hallmark Enrichment Compared to LumHR"))
  hallmark_plt<-gsea_enrichment_againstlumhr(species="human",
              category="H",
              subcategory=NULL,
              out_setname="hallmark",sample_name=sample_name,
              dmrs=collapsed_dmrs,obj=obj,outdir=outdir)

  print(paste("Calculating Cancer Cell Atlas Enrichment Compared to LumHR"))
  cancercellatlas_plt<-gsea_enrichment_againstlumhr(species="human",
              category="C4",
              subcategory="3CA",
              out_setname="3CA",sample_name=sample_name,
              dmrs=collapsed_dmrs,obj=obj,outdir=outdir)
  plt<-patchwork::wrap_plots(list((tft_plt),(position_plt),(hallmark_plt),(cancercellatlas_plt)),ncol=4,axes="collect")

  ggsave(plt,
          file=paste0(sample_name,".clones_lumhrall.GSEA_enrichment.pdf"),
          path=outdir,
          width=40,height=length(unique(collapsed_dmrs$type))*5,limitsize = FALSE)
  print(paste("Finished sample:",sample_name))
}

dmr_clone_by_lumhr<-function(obj,sample_name,dmr_clones_outdir,clone500bpwindows){
    #set dir
    clone_outdir=paste0(dmr_clones_outdir,"/",sample_name)
    system(paste("mkdir -p",clone_outdir))

    print(paste("Subset object for sample:",sample_name))
    
    #subset object to all lumhr OR sample specific clones
    cells=row.names(obj@metadata)[obj@metadata$Sample %in% c(sample_name) & obj@metadata$cnv_clonename!="NA" & !endsWith(obj@metadata$cnv_clonename,suffix="_diploid")]
    lumhr=row.names(obj@metadata)[obj@metadata$celltype %in% c("lumhr") & obj@metadata$Group=="HBCA"]
    cells<-unique(c(cells,lumhr))
    obj_sample<-subsetObject(obj,cells=cells)
    
    #subset to just sample clones
    print(paste("Subset 500bp windows for sample:",sample_name))
    pct_columns_to_keep=c(1,2,3,
        which(
          grepl(pattern=sample_name,
          unlist(lapply(strsplit(colnames(clone500bpwindows[["pct_matrix"]]),"_"),"[",1)))))
    sum_columns_to_keep=c(1,2,3,
        which(
          grepl(pattern=sample_name,
          unlist(lapply(strsplit(colnames(clone500bpwindows[["sum_matrix"]]),"_"),"[",1)))))
    #add lumhr columns to the keeprs
    pct_columns_to_keep<-c(pct_columns_to_keep,
          which(grepl(pattern="lumhr",(colnames(clone500bpwindows[["pct_matrix"]])))))

    sum_columns_to_keep<-c(sum_columns_to_keep,
          which(grepl(pattern="lumhr",(colnames(clone500bpwindows[["sum_matrix"]])))))

    pct_mat <-clone500bpwindows[["pct_matrix"]] %>% select(all_of(pct_columns_to_keep))
    sum_mat <-clone500bpwindows[["sum_matrix"]] %>% select(all_of(sum_columns_to_keep))

    print(paste("Saving amethyst object for clones:",sample_name))

    #save clone object for future genome track plotting
    obj_sample@genomeMatrices[["cg_clone_tracks"]] <- pct_mat #load it into amethyst object for plotting
    saveRDS(obj_sample,file=paste0(clone_outdir,"/",sample_name,".clones_lumhrall_amethyst.rds"))

    #run tests of clones against lumhr, and all clones v lumhr
    comparisons=colnames(pct_mat)[4:ncol(pct_mat)]
    comparisons<-comparisons[!grepl(comparisons,pattern="lumhr")]
    comparisons<-comparisons[!is.na(comparisons)]

    comparisons=data.frame(
      name=c(paste0(comparisons,"_lumhr"),paste0(sample_name,"_allclones_lumhr")),
      A=c(comparisons,paste(comparisons,collapse=",")),
      B=rep("lumhr",length(comparisons)+1))
    
    row.names(comparisons)<-comparisons$name

  lapply(row.names(comparisons), function(i) {
  print(paste("Running DMRs across clones for:",i))
  dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
          comparisons=comparisons[i,],
          nminTotal = 3, # Min number observations across all groups to include the region in calculations
          nminGroup = 3) # Min number observations across either members or nonmembers to include the region
  dmrs$type <- i
      saveRDS(dmrs,file=paste0(clone_outdir,"/",sample_name,".",i,".clones_lumhrall_dmr.unfilt.rds"))

  print(paste("Filtering DMRs across clones for:",i))
  dmrs <- filterDMR(dmrs, 
              method = "bonferroni", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
              filter = TRUE, # If TRUE, removes insignificant results
              pThreshold = 0.05 # Maxmimum adjusted p value to allow if filter = TRUE
              ) # Minimum absolute value of the log2FC to allow if filter = TRUE
  dmrs$type <- i
      saveRDS(dmrs,file=paste0(clone_outdir,"/",sample_name,".",i,".clones_lumhrall_dmr.rds"))
  collapsed_dmrs <- collapseDMR(obj, 
                          dmrs, 
                          maxDist = 2000, # Max allowable overlap between DMRs to be considered adjacent
                          minLength = 2000, # Min length of collapsed DMR window to include in the output
                          reduce = T, # Reduce results to unique observations (recommended)
                          annotate = T) # Add column with overlapping gene names

  collapsed_dmrs$type <- i
      saveRDS(collapsed_dmrs,file=paste0(clone_outdir,"/",sample_name,".",i,".clones_lumhrall_dmr_filt_collapse.rds"))
  })


}

```

# Calculate 500bp windows per clone and lumhr

```R
#############################################
#Make 500bp windows once and slice for comparisons
#############################################

#### USING ONLY CELLS WITH CNV CLONE NAME ASSIGNED
#### AND CELL TYPE LABELLED AS LUMHR OR CANCER

dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")
system(paste("mkdir -p", dmr_outdir))

#clones (all lumhr treated as same group)
dmr_clones_outdir=paste(sep="/",dmr_outdir,"cnv_clones_alllumhr")
system(paste("mkdir -p", dmr_clones_outdir))

obj<-readRDS(file="07_scaledcis.cnv_clones.amethyst.rds")

#select cells used in analysis
clone_cells <-row.names(obj@metadata)[(obj@metadata$celltype %in% c("lumhr","cancer")) & (obj@metadata$cnv_clonename != "NA")]
lumhr_hbca <-row.names(obj@metadata)[(obj@metadata$celltype %in% c("lumhr")) & obj@metadata$Group %in% c("HBCA")]

obj_lumhr_all<-subsetObject(obj,cells=unique(c(clone_cells,lumhr_hbca)))

obj_lumhr_all@metadata$cnv_clonename_alllumhr<-obj_lumhr_all@metadata$cnv_clonename
obj_lumhr_all@metadata[lumhr_hbca,]$cnv_clonename_alllumhr<-"lumhr"
table(obj_lumhr_all@metadata$cnv_clonename_alllumhr)

obj_lumhr_all<-subsetObject(obj_lumhr_all,cells=row.names(obj_lumhr_all@metadata[obj_lumhr_all@metadata$cnv_clonename_alllumhr!="NA",]))


#RERUN, USE ONLY LUMHR FROM HBCA
clone500bpwindows <- amethyst::calcSmoothedWindows(obj_lumhr_all, 
                                        type = "CG", 
                                        threads = 20,
                                        step = 500,
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "cnv_clonename_alllumhr",
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)

saveRDS(clone500bpwindows,file=paste0(dmr_clones_outdir,"/","dmr_analysis.cnv_clones_alllumhr.500bp_windows.rds"))

```

# Calculating clones compared to lumhr

```R
#add plotting functions for DMR 
#https://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/pbmc_vignette/pbmc_vignette.html
#add dmr locations per clone over heatmap cnv plot

obj<-readRDS(file="07_scaledcis.cnv_clones.amethyst.rds")

#output to DMR analysis folder
#clones (all lumhr treated as same group)
dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")
dmr_clones_outdir=paste(sep="/",dmr_outdir,"cnv_clones_alllumhr")
clone500bpwindows<-readRDS(file=paste0(dmr_clones_outdir,"/","dmr_analysis.cnv_clones_alllumhr.500bp_windows.rds"))

#run dmr analysis clone with all lumhr
#use dmr_analysis.cnv_clones_alllumhr.500bp_windows.rds matrix
#note some of these also have the coverage to do specific clone vs lumhr analysis, can merge matrices for this case
clone_v_lumhr_analysis<-c(
'BCMHBCA83L-3h',
'BCMDCIS05T',
'BCMDCIS82T',
'BCMDCIS99T',
'BCMHBCA26L',
'BCMDCIS22T',
'BCMDCIS79T',
'BCMDCIS28T',
'BCMDCIS52T',
'BCMDCIS80T',
'ECIS48T',
'BCMDCIS102T',
'BCMDCIS92T',
'BCMDCIS70T',
'BCMDCIS124T',
'ECIS36T',
'BCMDCIS74T',
'BCMDCIS97T',
'BCMDCIS94T',
'BCMDCIS41T',
'ECIS25T',
'ECIS26T',
'BCMDCIS66T',
'BCMDCIS35T'
)

mclapply(clone_v_lumhr_analysis,function(x) {
  dmr_clone_by_lumhr(obj,sample_name=x,dmr_clones_outdir,clone500bpwindows)}
  ,mc.cores=50)
```

Plot GO terms per clone as heatmaps

```R
dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")
dmr_clones_outdir=paste(sep="/",dmr_outdir,"cnv_clones_alllumhr")


gsea_set<-c("TFT","3CA","hallmark","position")


lapply(c("TFT","3CA","hallmark","position"),function(i){
gsea_list<-list.files(path=dmr_clones_outdir,full.names=T,recursive=T,pattern=paste0(i,".GSEA_enrichment.clone_lumhrall.rds"))

gsea_clones<-lapply(gsea_list,function(x){
  sample_name<-strsplit(basename(x),"[.]")[[1]][1]
  clone_name<-strsplit(basename(x),"[.]")[[1]][2]
  gsea<-readRDS(x)
  gsea$sample_name<-sample_name
  gsea$clone_name<-clone_name
  return(gsea)
})

gsea_all_clones<-do.call("rbind",gsea_clones)

nes<-gsea_all_clones %>% tidyr::pivot_wider(names_from=clone_name,values_from=NES,id_cols=pathway,id_expand=TRUE) %>% as.data.frame()
pval<-gsea_all_clones %>% tidyr::pivot_wider(names_from=clone_name,values_from=padj,id_cols=pathway,id_expand=TRUE) %>% as.data.frame()

row.names(nes)<-nes$pathway
nes<-nes[,2:ncol(nes)]

row.names(pval)<-pval$pathway
pval<-pval[,2:ncol(pval)]
pval[which(is.na(pval),arr.ind=T)]<-1

column_sample<-unlist(lapply(strsplit(colnames(nes),"_"),"[[",1))
sample_info<-obj@metadata[!duplicated(obj@metadata$Sample),]

sample_info$matched_names<-unlist(lapply(strsplit(sample_info$Sample,"_"),"[[",1))
sample_info<-sample_info[!duplicated(sample_info$matched_names),]
row.names(sample_info)<-sample_info$matched_names

arm_col=c("p"="grey","q"="darkgrey")
band_col=c("acen"="#99746F","gneg"="white","gpos100"="black","gpos25"="lightgrey","gpos50"="grey","gpos75"="darkgrey","gvar"="#446879")

col_fun=colorRamp2(c(1,0.05,0),c("white","red","darkred"))

column_ha = HeatmapAnnotation(
                            Group = sample_info[column_sample,]$Group,
                            ER = sample_info[column_sample,]$ER,
                            PR = sample_info[column_sample,]$PR,
                            HER2 = sample_info[column_sample,]$HER2)

plt1<-Heatmap(pval,
  column_split=factor(sample_info[column_sample,]$Group,levels=c("HBCA","DCIS","Synchronous","IDC")),
  show_row_names = TRUE, row_title_rot = 0,
  col=col_fun,
  row_names_gp = gpar(fontsize = 5),
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10),
  top_annotation=column_ha,
  border = TRUE,
  width = 10)
  
pdf(paste0(dmr_clones_outdir,"/","clones_v_lumhr.",i,".GSEA.pdf"),height=10,width=10)
print(plt1)
dev.off()
})




```

# Metaprogram per clone
Using same subsetting of cells
This isnt working all that well.

```R

obj<-readRDS(file="07_scaledcis.cnv_clones.amethyst.rds")

#select cells used in analysis
clone_cells <-row.names(obj@metadata)[(obj@metadata$celltype %in% c("lumhr","cancer")) & (obj@metadata$cnv_clonename != "NA")]
lumhr_hbca <-row.names(obj@metadata)[(obj@metadata$celltype %in% c("lumhr")) & obj@metadata$Group %in% c("HBCA")]

obj_lumhr_all<-subsetObject(obj,cells=unique(c(clone_cells,lumhr_hbca)))

obj_lumhr_all@metadata$cnv_clonename_alllumhr<-obj_lumhr_all@metadata$cnv_clonename
obj_lumhr_all@metadata[lumhr_hbca,]$cnv_clonename_alllumhr<-"lumhr"
table(obj_lumhr_all@metadata$cnv_clonename_alllumhr)

obj_lumhr_all<-subsetObject(obj_lumhr_all,cells=row.names(obj_lumhr_all@metadata[obj_lumhr_all@metadata$cnv_clonename_alllumhr!="NA",]))

#Adding metaprogram (running pca on dmrs per clone, or per sample?)
clone500bpwindows<-readRDS(file=paste0(dmr_clones_outdir,"/","dmr_analysis.cnv_clones_alllumhr.500bp_windows.rds"))

#calc dmr for all
#calc over window per cell
sample_name="all_clones"
sum_mat<-clone500bpwindows[["sum_matrix"]]

print(paste("Running DMRs across clones for:",sample_name))
dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
        eachVsAll=TRUE,
        nminTotal = 3, # Min number observations across all groups to include the region in calculations
        nminGroup = 3) # Min number observations across either members or nonmembers to include the region

print(paste("Filtering DMRs across clones for:",sample_name))
dmrs <- filterDMR(dmrs, 
            method = "bonferroni", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
            filter = TRUE, # If TRUE, removes insignificant results
            pThreshold = 0.05 # Maxmimum adjusted p value to allow if filter = TRUE
            ) 

rename_dmr_output<-setNames(nm=1:length(unique(dmrs$test)),gsub(colnames(sum_mat)[grepl(colnames(sum_mat),pattern="_t")],pattern="_t",replacement=""))
dmrs$type <- rename_dmr_output[dmrs$test]
saveRDS(dmrs,file=paste0(dmr_clones_outdir,"/",sample_name,".clones_lumhrall_dmr.rds"))

collapsed_dmrs <- collapseDMR(obj_lumhr_all, 
                        dmrs, 
                        maxDist = 2000, # Max allowable overlap between DMRs to be considered adjacent
                        minLength = 2000, # Min length of collapsed DMR window to include in the output
                        reduce = T, # Reduce results to unique observations (recommended)
                        annotate = T) # Add column with overlapping gene names

rename_dmr_output<-setNames(nm=1:length(unique(collapsed_dmrs$test)),gsub(colnames(sum_mat)[grepl(colnames(sum_mat),pattern="_t")],pattern="_t",replacement=""))
collapsed_dmrs$type <- rename_dmr_output[collapsed_dmrs$test]
saveRDS(collapsed_dmrs,file=paste0(dmr_clones_outdir,"/",sample_name,".clones_lumhrall_dmr_filt_collapse.rds"))


collapsed_dmrs<-readRDS(file=paste0(dmr_clones_outdir,"/",sample_name,".clones_lumhrall_dmr_filt_collapse.rds"))
bed<-collapsed_dmrs %>% filter(dmr_padj<0.05) %>% GRanges() %>% reduce()
dmr_bed<-bed %>% as.data.frame() %>% mutate(chr=seqnames) %>% select(chr,start,end)



#create new matrix from DMR sites for refined metaprograms
window_name="clone_dmr_sites"
obj_lumhr_all@genomeMatrices[[window_name]] <- makeWindows(obj_lumhr_all, 
                                                     bed = dmr_bed,
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 50, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 
saveRDS(obj_lumhr_all,file=paste0(dmr_clones_outdir,"/",sample_name,".clones_lumhrall_dmr.amethyst.rds"))

mat<-obj_lumhr_all@genomeMatrices[[window_name]]
mat[which(is.na(mat),arr.ind=T)]<-0
mat<-mat+1 #scaled and centered around 1 (regular is -1 to 1,now is 0 to 2)
mat<- as.data.frame(t(mat))  %>% split(obj_lumhr_all@metadata[colnames(mat),]$cnv_clonename_alllumhr)

#similar to https://github.com/carmonalab/GeneNMF/blob/master/R/main.R


k=c(10,15)
#run nmf by sample
nmf.res <- lapply(names(mat), function(this) {
    if(nrow(as.data.frame(mat[[this]]))>50){
    res.k <- lapply(k, function(k.this) {
      model <- RcppML::nmf(t(mat[[this]]), k = k.this, L1 = c(0,0), verbose=FALSE, seed=123)      
      rownames(model$h) <- paste0("pattern",1:nrow(model$h))
      colnames(model$h) <- colnames(t(mat[[this]]))
      rownames(model$w) <- rownames(t(mat[[this]]))
      colnames(model$w) <- paste0("pattern",1:ncol(model$w))
      model
    })
    names(res.k) <- paste0("k",k)
    res.k
    }
  })

nmf.res <- unlist(nmf.res, recursive = FALSE)
nmf.genes<-getNMFgenes(nmf.res)
metaprograms<-getMetaPrograms(nmf.res)

plt<-plotMetaPrograms(metaprograms,downsample=1,)
library(pheatmap)
library(viridis)


mp.res<-metaprograms
similarity.cutoff=c(0.2,1)
scale = "none"
showtree = TRUE
colfun = colorRamp2(breaks=c(0,1),c("white","darkred"))
col<-colfun(seq(0,1,length=100))

main = "Clustered Heatmap"
show_rownames = FALSE
show_colnames = FALSE


  J <- mp.res[["programs.similarity"]]
  tree <- mp.res[["programs.tree"]]
  cl_members <- mp.res[["programs.clusters"]]
  labs.order <- labels(as.dendrogram(tree))

  
  cl_names <- names(cl_members)
  cl_members[!is.na(cl_members)] <- paste0("MP",cl_members[!is.na(cl_members)])
  names(cl_members) <- cl_names

  #Recover order of MP clusters
  cluster.order <- unique(cl_members[labs.order])
  nMP <- length(cluster.order)
  
  #Gaps for heatmap
  diffs <- diff(as.integer(as.factor(cl_members))) 
  gaps <- which(diffs != 0)
  
  #Annotation column
  annotation_col <- as.data.frame(cl_members)
  colnames(annotation_col) <- "Metaprogram"
  annotation_col[["Metaprogram"]] <- factor(cl_members, levels=cluster.order)

  #Apply trimming to similarity for plotting  
  J[J<similarity.cutoff[1]] <- similarity.cutoff[1]
  J[J>similarity.cutoff[2]] <- similarity.cutoff[2]
  
  
  ph <- pheatmap(J,
                 scale = scale,
                 color = col,
                 main = main,
                 cluster_rows = tree,
                 cluster_cols = tree,
                 cutree_rows = nMP,
                 cutree_cols = nMP,
                 annotation_names_col = FALSE,
                 annotation_names_row = FALSE,
                 show_rownames = show_rownames,
                 show_colnames = show_colnames)

ggsave(ph,file=paste0(dmr_clones_outdir,"/",sample_name,".clones_lumhrall_dmr.metaprogram.pdf"))


  

#to add:


    #ADD TF ENRICHMENT HERE LOLA
    #universe
    #use motifmatchr to make list of motifs per 500bp windows?
#    universe<-sum_mat %>% GRanges()
#    mcols(universe)<-NULL
    #split dmrs
#    dmr_list<-dmrs %>% filter(padj<0.05)  %>% GRanges()  %>% split(~type+direction)
    #get regionDB
#    jaspar<-rtracklayer::import.bb("/data/rmulqueen/projects/scalebio_dcis/ref/JASPAR2026_hg38.bb") #WAITING ON DOWNLOAD
#    locResults = runLOLA(dmr_list, universe, regionDB, cores=50)





  chromvar_methylation<-function(obj,counts,prefix="allcells",threads){
    if(dim(counts)[2]>200){
      print("Treating counts matrix as single cell input.")
    }else {
      print("Treating counts matrix as summarized cluster input.")
    }
    #prepare summarized experiment for chromvar
    peaks<-GenomicRanges::makeGRangesFromDataFrame(data.frame(
      seqnames=unlist(lapply(strsplit(row.names(counts),"_"),"[",1)),
      start=unlist(lapply(strsplit(row.names(counts),"_"),"[",2)),
      end=unlist(lapply(strsplit(row.names(counts),"_"),"[",3))))

    #prepare motifs
    opts <- list()
    opts[["species"]] <- "Homo sapiens"
    opts[["collection"]] <- "CORE"
    opts[["all_versions"]] <- FALSE
    motifs <- TFBSTools::getMatrixSet(JASPAR2020,opts)

  #split peaks evenly into chunks so we can multicore the motif scanning
  motif_matches<-mclapply(split(peaks,  cut(seq_along(peaks), threads, labels = FALSE)),
                          function(x){
                          matchMotifs(motifs, x, genome = BSgenome.Hsapiens.UCSC.hg38, p.cutoff=0.01)},
                          mc.cores=threads)
  motif_ix<-do.call("rbind",motif_matches)

  #create summarized experiment
  if(dim(counts)[2]<200){ colnames(counts)<-paste("cluster",colnames(counts),sep="_")}

  rse <- SummarizedExperiment::SummarizedExperiment(
                                  assays=list(counts=as(counts, "sparseMatrix")),
                                  rowRanges=peaks)
  colData(rse)<-as(obj@metadata[colnames(counts),],"DataFrame")
  rse <- addGCBias(rse, genome = BSgenome.Hsapiens.UCSC.hg38)
  dev <- computeDeviations(object = rse, annotations = motif_ix)
  saveRDS(dev,file=paste0(prefix,".chromvar.rds"))

  #calculate variability
  variability <- computeVariability(dev)
  ggsave(plotVariability(variability, use_plotly = FALSE),file=paste0(prefix,".chromvar_variability.pdf"))

  #Differential motif analysis (for single cell)
  if(dim(dev)[2]>200){
  diff_acc <- differentialDeviations(dev, "cluster_id")
  diff_var <- differentialVariability(dev, "cluster_id")
  }

  #differential tfbs by highest variability
  var_cut<-ifelse(dim(dev)[2]>200,1.5,0.3)

  diff_tfbs<-row.names(variability[variability$variability>var_cut,])
  devs<-deviationScores(dev)
  devs[is.na(devs)]<-0 #fill in NA for dev scores
  #dim_out<-irlba::irlba(devs[diff_tfbs,], 30)
  dim_out<-t(devs[diff_tfbs,])
  dim<-uwot::umap(dim_out,n_neighbors=2)
  dim<-as.data.frame(dim)
  colnames(dim)<-c("chromvar_umap_x","chromvar_umap_y")
  row.names(dim)<-colnames(devs)
  if(dim(dev)[2]>200){
  dim$cluster_id<-obj@metadata[row.names(dim),]$cluster_id
  rowannot<-as.data.frame(colData(dev)[c("cluster_id","sample")])
  } else {
    dim$cluster_id<-colnames(counts)
    rowannot<-as.data.frame(colnames(counts))
    row.names(rowannot)<-row.names(dim)
  }

  plt<-ggplot(dim,aes(x=chromvar_umap_x,y=chromvar_umap_y,color=cluster_id))+geom_point()+theme_minimal()
  ggsave(plt,file=paste0(prefix,".chromvar_umap.pdf"))

  sample_cor <- getSampleCorrelation(dev,threshold=var_cut)
  sample_cor[is.na(sample_cor)]<-0 #fill in na as 0 for sites with no overlap
  plt<-pheatmap(as.dist(sample_cor), 
          annotation_row = rowannot,
          clustering_distance_rows = as.dist(1-sample_cor), 
          clustering_distance_cols = as.dist(1-sample_cor))
  ggsave(plt,file=paste0(prefix,".chromvar_motifs.correlation.heatmap.pdf"))

    devs<-deviationScores(dev)
    colnames(devs)<-colnames(counts)
    row.names(devs)<-rowData(dev)$name
    diff_tfbs_names<-rowData(dev)[diff_tfbs,]$name
    devs<-devs[diff_tfbs_names,]
    devs<-scale(devs)
    devs_row_order<-hclust(dist(devs))
    plt<-pheatmap(devs,cluster_rows=devs_row_order,fontsize=4,angle_col="90",treeheight_col=0,color = colorRampPalette(c("black", "white", "magenta"))(100))
    ggsave(plt,file=paste0(prefix,".chromvar_motifs.heatmap.pdf"))

###Add motifs to bottom of chromvar plots, so we can see some CG enrichment!
library(ggplotify)

    motifs_to_plot <-motifs[diff_tfbs] 
    motifs_to_plot<-motifs_to_plot[devs_row_order$order]

    plt_motifs<-ggseqlogo::ggseqlogo(data=Matrix(motifs_to_plot),seq_type="dna",method="bits",ncol=1)
    #plt_out<-patchwork::wrap_plots(as.ggplot(plt),plt_motif)+patchwork::plot_layout(ncol = 2, heights = c(10, 10), widths = c(10,2))
    #ggsave(plt_out,file="test_motif.pdf",width=10,height=10)
   #plt_out<-plot_grid(as.ggplot(plt), plt_motif, ncol = 2, align = "h",axis="rt",rel_widths = c(1, 0.1))
    ggsave(plt_motifs,file="test_motif.pdf",width=5,height=100,limitsize=F)



    plt<-MotifPlot(object = x,assay="peaks",motifs = motif_list[get_order(o_rows,1)],ncol=1)+theme_void()+theme(strip.text = element_blank())
    ggsave(plt,file=paste0(prefix,".tf.heatmap.motif.pdf"),height=100,width=2,limitsize=F)
  saveRDS(dev,file=paste0(prefix,".chromvar.rds"))
}


```

# DMR Per CNV Clone Per Sample    

```R
#output to DMR analysis folder
dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")
dmr_clones_outdir=paste(sep="/",dmr_outdir,"cnv_clones")

#run 500bp smoothed windows by clone split once, save and subset
#run above, note that amethyst also allows for group by group specified comparisons now too
clone500bpwindows<-readRDS(file=paste0(dmr_clones_outdir,"/","dmr_analysis.cnv_clones.500bp_windows.rds"))

#samples with >= 2 clones and >=5% coverage per clone
#use dmr_analysis.cnv_clones.500bp_windows.rds matrix
clone_v_clone_analysis<-c('BCMDCIS05T',
'BCMDCIS124T',
'BCMDCIS28T',
'BCMDCIS35T',
'BCMDCIS41T',
'BCMDCIS52T',
'BCMDCIS66T',
'BCMDCIS74T',
'BCMDCIS79T',
'BCMDCIS82T',
'BCMDCIS97T',
'ECIS25T',
'ECIS26T',
'ECIS36T',
'ECIS48T')

lapply(clone_v_clone_analysis,function(x) dmr_clone_by_clone(obj,sample_name=x,dmr_clones_outdir,clone500bpwindows))

#add bedfile of dmrs, bigWig of CNVs per clone, stacked bigwig of methylation per clone
#split bw into 4 files per cell type by methylation/average methylation
for(i in colnames(clone500bpwindows[["pct_matrix"]])[4:ncol(clone500bpwindows[["pct_matrix"]])]){
    sample<-strsplit(i,"_")[[1]]
    sample<-paste(sample[1:length(sample)-1],collapse="_")
    sample_dir=paste0(dmr_clones_outdir,"/",sample,"/bigwig")
    system(paste("mkdir -p",sample_dir))

    out_dat<-clone500bpwindows[["pct_matrix"]] %>% select(chr,start,end,i) 

    hg38_seq_info<-Seqinfo(genome="hg38")
    out_dat<-GRanges(out_dat[complete.cases(out_dat),]) #filter NA
    out_dat<-out_dat[out_dat@seqnames %in% hg38_seq_info@seqnames,] #filter chr
    out_dat<-out_dat[out_dat@seqnames %in% paste0("chr",c(1:22,"X")),] #filter chr
    seqlevels(out_dat)<-paste0("chr",c(1:22,"X"))

    out_dat<-resize(out_dat,width=500)
    names(out_dat@elementMetadata)<-"score"
    mean_score<-mean(mcols(out_dat)$score)

    #bin to 100-mean, mean-50, 50-20, 20-0  
    #color black, grey, lightgrey, celltypecol
    #333333, #444444, #bcbcbc, celltypecol
    #subtract score-meanscorevalue

    out_dat_hypermet <- out_dat %>% 
                        as.data.frame() %>% 
                        filter(mcols(out_dat)$score > mean_score) %>% 
                        mutate(score=score-mean_score) %>% 
                        GRanges()
    names(out_dat_hypermet@elementMetadata)<-"score"
    genome(out_dat_hypermet)<-"hg38"
    seqlengths(out_dat_hypermet)<-seqlengths(hg38_seq_info)[seqlevels(out_dat_hypermet)] #filter by seqlengths

    out_dat_met_mid <- out_dat %>% 
                      as.data.frame() %>% 
                      filter(mcols(out_dat)$score <= mean_score & mcols(out_dat)$score > 50) %>% 
                      mutate(score=score-mean_score) %>% 
                      GRanges()
    names(out_dat_met_mid@elementMetadata)<-"score"
    genome(out_dat_met_mid)<-"hg38"
    seqlengths(out_dat_met_mid)<-seqlengths(hg38_seq_info)[seqlevels(out_dat_met_mid)] #filter by seqlengths

    out_dat_met_low <- out_dat %>% 
                        as.data.frame() %>% 
                        filter(mcols(out_dat)$score <= 50 & mcols(out_dat)$score > 20) %>% 
                        mutate(score=score-mean_score) %>% 
                        GRanges()
    names(out_dat_met_low@elementMetadata)<-"score"
    genome(out_dat_met_low)<-"hg38"
    seqlengths(out_dat_met_low)<-seqlengths(hg38_seq_info)[seqlevels(out_dat_met_low)] #filter by seqlengths

    out_dat_met_hypomet <- out_dat %>% 
                            as.data.frame() %>% filter(mcols(out_dat)$score <= 20) %>%
                            mutate(score=score-mean_score) %>% 
                            GRanges()
    names(out_dat_met_hypomet@elementMetadata)<-"score"
    genome(out_dat_met_hypomet)<-"hg38"
    seqlengths(out_dat_met_hypomet)<-seqlengths(hg38_seq_info)[seqlevels(out_dat_met_hypomet)] #filter by seqlengths

    print(paste("Saving bedgraphs for...",i))
    rtracklayer::export(out_dat_hypermet,con=paste0(sample_dir,"/",paste(i,"hypermet","bw",sep=".")))
    rtracklayer::export(out_dat_met_mid,con=paste0(sample_dir,"/",paste(i,"midmet","bw",sep=".")))
    rtracklayer::export(out_dat_met_low,con=paste0(sample_dir,"/",paste(i,"lowmet","bw",sep=".")))
    rtracklayer::export(out_dat_met_hypomet,con=paste0(sample_dir,"/",paste(i,"hypomet","bw",sep=".")))
    print(paste("Completed ",i, " for sample ",sample))
    print(paste("Saved in ",sample_dir))
    }

#move some files around to make igv tracks per sample (including the bigwig from the cnv bedgraph)
```



Run DMR results through gene ontology for interpretation 

- For 500bp windows, look for TF binding site enrichment
- For collapsed DMRS, look for ontology and positional enrichments

```R
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")
dmr_clones_outdir=paste(sep="/",dmr_outdir,"cnv_clones_alllumhr")

dmr_list <- list.files(path=dmr_outdir,full.names=TRUE,recursive=TRUE,pattern=".clones_lumhrall_dmr.rds")

#USING JASPAR AND LOLA PACKAGE FOR ENRICHMENT


#gsea of top DMR genes
gsea_enrichment_againstlumhr<-function(dmr_indx,species="human",
                          category="C3",
                          subcategory="TFT:GTRD",
                          out_setname="TFT",sample_name=sample_name,
                          outdir=outdir,
                          obj=obj){

  dmr<-readRDS(dmr_list[dmr_indx]) #read in 500bp dmr
  #run collapse, but without join 
  #split by hypo and hyper and which is significant
  #run LOLA with the universe as all 500bp sets on JASPAR

  dmrfilt <- filterDMR(dmr, 
              method = "bonferroni", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
              filter = FALSE) # Minimum absolute value of the log2FC to allow if filter = TRUE
  dmr$cluster<-
  dmr<-dmr %>% split(direction,test)select(chr,start,end,ends_with("_pval"))  %>% filter(if_any(where(is.numeric), ~ .x < 0.05)) #filter to pval column and only those significant
     
  collapsed_dmrs <- collapseDMR(obj, 
                          dmr, 
                          maxDist = 0, # Max allowable overlap between DMRs to be considered adjacent
                          minLength = 1, # Min length of collapsed DMR window to include in the output
                          reduce = F, # Reduce results to unique observations (recommended)
                          annotate = T) # Add column with overlapping gene names

  pathwaysDF <- msigdbr(species=species, 
                        collection=category, 
                        subcollection = subcategory)

  #limit pathways to genes in our data
  pathwaysDF<-pathwaysDF[pathwaysDF$ensembl_gene %in% unlist(lapply(strsplit(unique(obj@ref$gene_id),"[.]"),"[",1)),]
  
  pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

  #run plotting per group and per hyper and hypo? i think just using logFC is sufficient

  plt_list<-lapply(unique(dmrs$type),function(group1){
  #treat multiple gene overlaps as same logFC
  #set -Inf to -3 and Inf to 3
  group_features<-dmrs %>%
    dplyr::filter(type == group1) %>%
    dplyr::filter(dmr_padj<0.05) %>% 
    dplyr::filter(gene_names!="NA") %>% 
    dplyr::arrange(dmr_logFC) %>%
    dplyr::select(gene_names, dmr_logFC) %>%
    tidyr::separate_rows(gene_names) %>%
    dplyr::mutate(across(where(is.numeric), ~ replace(., .==-Inf, -3))) %>%
    dplyr::mutate(across(where(is.numeric), ~ replace(., .==Inf, 3)))

  ranks<-setNames(nm=group_features$gene_names,group_features$dmr_logFC)
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = ranks,
                    minSize  = 5,
                    nproc = 1)

  topPathwaysUp <- fgseaRes %>% filter(ES > 0) %>% slice_max(NES,n=10) %>% dplyr::select(pathway)
  topPathwaysDown <- fgseaRes %>% filter(ES < 0) %>% slice_max(abs(NES),n=10) %>% dplyr::select(pathway)
  topPathways <- unlist(c(topPathwaysUp, rev(topPathwaysDown)))
  saveRDS(fgseaRes,file=paste0(outdir,sample_name,".",out_setname,".GSEA_enrichment.clone_lumhrall.rds"))

  #not returned
  plt1<-plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam=0.5)+
        theme(axis.text.y = element_text( size = rel(0.2)),
        axis.text.x = element_text( size = rel(0.2)))

  # only plot the top 20 pathways NES scores
  nes_plt_dat<-rbind(
    fgseaRes  %>% slice_max(NES,n= 10),
    fgseaRes  %>% slice_min(NES,n= 10))

  plt2<-ggplot(nes_plt_dat, aes(reorder(pathway, NES), NES)) +
    geom_point(aes(color= padj,size=size)) + scale_color_gradient2(low="darkred",high="grey",mid="red",midpoint=0.05)+
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score \n Compared to normal LumHR \n (Hyper<->Hypo)",
        title="Hallmark pathways NES from GSEA") + 
    theme_minimal()+scale_fill_identity()+ggtitle(paste(group1,out_setname))+ylim(c(5,-5))
  return(plt2)
  })
  patchwork::wrap_plots(plt_list,ncol=1)

}

```R
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")
dmr_clones_outdir=paste(sep="/",dmr_outdir,"cnv_clones_alllumhr")

collapsed_dmr_list <- list.files(path=dmr_outdir,full.names=TRUE,recursive=TRUE,pattern=".clones_lumhrall_dmr_filt_collapse.rds")



```

#coregulation by gsea https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/geseca-tutorial.html



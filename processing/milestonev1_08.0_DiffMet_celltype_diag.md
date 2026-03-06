#```bash
#singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
#source activate
#update data.table
#```

# To Do: Also try applying Bartlett's test for differential variability instead of just diff mean
see DNA methylation outliers in normal breast tissue  identify field defects that are enriched in cancer 2015

# Applying bartlett.test()
result = bartlett.test(len ~ interaction(supp, dose)
, 
                                  data = ToothGrowth)

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
library(LOLA)
library(parallel)

options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)
obj<-readRDS(file="07_scaledcis.cnv_clones.amethyst.rds")
```

Running through all the 500bp window generation split by:
- Celltype (fine_celltype)
- CNV based clones (cnv_clonenames)
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
  saveRDS(fgseaRes,file=paste0(outdir,"/",paste0(sample_name,".",out_setname,".GSEA_enrichment.clone_lumhrall.rds")))

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

dmr_celltype_by_all<-function(obj){
    #read 500bp windows
    celltype500bpwindows<-readRDS(paste0(dmr_celltype_outdir,"/","dmr_analysis.celltype.500bp_windows.rds"))
    
    pct_mat<-celltype500bpwindows[["pct_matrix"]] 
    sum_mat<-celltype500bpwindows[["sum_matrix"]] 

    #save clone object for future genome track plotting
    obj@genomeMatrices[["cg_celltype_tracks"]] <- pct_mat #load it into amethyst object for plotting

    dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
            eachVsAll = TRUE, # If TRUE, each group found in the sumMatrix will be tested against all others
            nminTotal = 3, # Min number observations across all groups to include the region in calculations
            nminGroup = 3) # Min number observations across either members or nonmembers to include the region
    saveRDS(dmrs,file=paste0(dmr_celltype_outdir,"/","celltype_allcells",".dmr.rds"))
   
    dmrs <- filterDMR(dmrs, 
                method = "bonferroni", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
                filter = FALSE, # If TRUE, removes insignificant results
                pThreshold = 0.05, # Maxmimum adjusted p value to allow if filter = TRUE
                logThreshold = 1.5) # Minimum absolute value of the log2FC to allow if filter = TRUE
    collapsed_dmrs <- collapseDMR(obj, 
                            dmrs, 
                            maxDist = 2000, # Max allowable overlap between DMRs to be considered adjacent
                            minLength = 2000, # Min length of collapsed DMR window to include in the output
                            reduce = T, # Reduce results to unique observations (recommended)
                            annotate = T) # Add column with overlapping gene names
    saveRDS(collapsed_dmrs,file=paste0(dmr_celltype_outdir,"/","celltype_allcells",".dmr_filt_collapse.rds"))

    rename_dmr_output<-setNames(nm=1:length(unique(collapsed_dmrs$test)),gsub(colnames(sum_mat)[grepl(colnames(sum_mat),pattern="_t$")],pattern="_t",replacement=""))
    collapsed_dmrs$type <- rename_dmr_output[collapsed_dmrs$test]
    saveRDS(collapsed_dmrs,file=paste0(dmr_celltype_outdir,"/","celltype_allcells",".dmr_filt_collapse.rds"))
    
    #plot dmr counts per clone
    pal <- colorRampPalette(colors = c("#F9AB60", "#E7576E", "#630661", "#B5DCA5"))
    COLS <- pal(length(unique(collapsed_dmrs$type)))

    plt<-ggplot(collapsed_dmrs |> dplyr::group_by(type, direction) |> dplyr::summarise(n = n()), 
        aes(y = type, x = n, fill = type)) + 
        geom_col() + 
        facet_grid(vars(direction), scales = "free_y") + 
        scale_fill_manual(values = COLS) + 
        theme_classic()
    ggsave(plt,file=paste0(dmr_celltype_outdir,"/","celltype_allcellsm",".dmr_counts.pdf"))
}

```

Calculate met per 500bp window for celltype by diagnosis (can add back together if needed for DMR)

```R
#############################################
#Make 500bp windows once and slice for comparisons
#############################################
dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")
system(paste("mkdir -p", dmr_outdir))

#celltype by diagnosis
dmr_celltype_by_diag_outdir=paste(sep="/",dmr_outdir,"celltype_by_diagnosis")
system(paste("mkdir -p", dmr_celltype_by_diag_outdir))
obj@metadata$celltype_diag<-paste(sep="_",obj@metadata$Group,obj@metadata$celltype)
celltypediag500bpwindows <- calcSmoothedWindows(obj, 
                                        type = "CG", 
                                        threads = 100,
                                        step = 500, 
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "celltype_diag",
                                        returnSumMatrix = TRUE, 
                                        returnPctMatrix = TRUE)
saveRDS(celltypediag500bpwindows,file=paste0(dmr_celltype_by_diag_outdir,"/","dmr_analysis.celltype_by_diag.500bp_windows.rds"))

```


# Run specific DMR comparisons

"chr"                       "start"
 [3] "end"                       "HBCA_TEC"
 [5] "HBCA_basal"                "HBCA_tcell_cd8"
 [7] "DCIS_cancer"               "HBCA_lumsec"
 [9] "Synchronous_cancer"        "IDC_cancer"
[11] "IDC_fibroblast"            "IDC_DC"
[13] "IDC_endothelial"           "DCIS_tcell_cd8"
[15] "HBCA_fibroblast"           "IDC_basal"
[17] "HBCA_CAF"                  "HBCA_endothelial"
[19] "DCIS_endothelial"          "DCIS_fibroblast"
[21] "HBCA_pericyte_VSMC"        "IDC_tcell_cd4"
[23] "HBCA_lumhr"                "DCIS_monocyte"
[25] "IDC_lumhr"                 "DCIS_basal"
[27] "Synchronous_fibroblast"    "HBCA_cancer"
[29] "DCIS_lumsec"               "DCIS_TEC"
[31] "Synchronous_lumhr"         "DCIS_CAF"
[33] "IDC_CAF"                   "HBCA_nk_tnk"
[35] "IDC_lumsec"                "DCIS_lumhr"
[37] "DCIS_DC"                   "IDC_monocyte"
[39] "HBCA_macrophage"           "IDC_bcell"
[41] "IDC_tcell_cd8"             "Synchronous_lumsec"
[43] "Synchronous_basal"         "IDC_nk_tnk"
[45] "DCIS_tcell_cd8_2"          "IDC_TAM_2"
[47] "DCIS_tcell_cd4"            "Synchronous_TEC"
[49] "DCIS_bcell"                "HBCA_monocyte"
[51] "HBCA_DC"                   "IDC_macrophage"
[53] "DCIS_macrophage"           "Synchronous_tcell_cd8"
[55] "HBCA_tcell_cd4"            "Synchronous_tcell_cd8_2"
[57] "Synchronous_CAF"           "DCIS_pericyte_VSMC"
[59] "Synchronous_TAM_2"         "DCIS_TAM"
[61] "Synchronous_endothelial"   "Synchronous_bcell"
[63] "IDC_TAM"                   "Synchronous_tcell_cd4"
[65] "IDC_tcell_cd8_2"           "IDC_pericyte_VSMC"
[67] "Synchronous_macrophage"    "HBCA_tcell_cd8_2"
[69] "Synchronous_pericyte_VSMC" "HBCA_TAM"
[71] "Synchronous_monocyte"      "Synchronous_TAM"
[73] "HBCA_bcell"                "DCIS_nk_tnk"
[75] "Synchronous_nk_tnk"        "DCIS_TAM_2"
[77] "Synchronous_DC"            "IDC_TEC"
[79] "HBCA_TAM_2"


```R

celltypediag500bpwindows<-readRDS(file=paste0(dmr_celltype_by_diag_outdir,"/","dmr_analysis.celltype_by_diag.500bp_windows.rds"))

sum_mat<-celltypediag500bpwindows[["sum_matrix"]]

sample_name="diagnosis_celltype"


#TEC vs Endo
#CAF vs Fibro
#basal DCIS vs HBCA
#basal DCIS vs Synchronous
#basal DCIS vs IDC

comparisons<-as.data.frame(rbind(
    caf_v_fibro=c(name="fibro_v_caf",A=paste(c("DCIS_CAF","Synchronous_CAF","IDC_CAF"),collapse=","),B="HBCA_fibroblast"),
    tec_v_endo=c(name="endo_v_tec",A=paste(c("DCIS_TEC","Synchronous_TEC","IDC_TEC"),collapse=","),B="HBCA_endothelial"),
    cancerBasal_v_basal=c(name="lumhr_v_dcisBasal",A=paste(c("DCIS_basal","Synchronous_basal","IDC_basal"),collapse=","),B="HBCA_basal"),
    cBasal_v_basal=c(name="lumhr_v_dcisBasal",A="Synchronous_basal",B="HBCA_basal"),
    dcisCancer_v_lumhr=c(name="lumhr_v_dcisCancer",A="DCIS_cancer",B="HBCA_lumhr"),
    syncCancer_v_dcisCancer=c(name="dcisCancer_v_syncCancer",A="Synchronous_cancer",B="DCIS_cancer"),
    idcCancer_v_dcisCancer=c(name="dcisCancer_v_syncCancer",A="IDC_cancer",B="DCIS_cancer"),
    idcCancer_v_lumhr=c(name="dcisCancer_v_syncCancer",A="IDC_cancer",B="HBCA_lumhr")

))

mclapply(row.names(comparisons), function(i) {
    print(paste("Running DMRs across clones for:",i))
    system(paste0("mkdir -p ",dmr_celltype_by_diag_outdir,"/",i))
    dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
            comparisons=comparisons[i,],
            nminTotal = 3, # Min number observations across all groups to include the region in calculations
            nminGroup = 3) # Min number observations across either members or nonmembers to include the region
    dmrs$type <- i
    saveRDS(dmrs,file=paste0(dmr_celltype_by_diag_outdir,"/",i,"/",i,".",sample_name,".celltype_diagnosis.dmrs.unfilt.rds"))

    print(paste("Filtering DMRs across clones for:",i))
    dmrs <- filterDMR(dmrs, 
                method = "bonferroni", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
                filter = TRUE, # If TRUE, removes insignificant results
                pThreshold = 0.05 # Maxmimum adjusted p value to allow if filter = TRUE
                ) # Minimum absolute value of the log2FC to allow if filter = TRUE
    dmrs$type <- i
    saveRDS(dmrs,file=paste0(dmr_celltype_by_diag_outdir,"/",i,"/",i,".",sample_name,".celltype_diagnosis.dmrs.filt.rds"))
    collapsed_dmrs <- collapseDMR(obj, 
                            dmrs, 
                            maxDist = 2000, # Max allowable overlap between DMRs to be considered adjacent
                            minLength = 2000, # Min length of collapsed DMR window to include in the output
                            reduce = T, # Reduce results to unique observations (recommended)
                            annotate = T) # Add column with overlapping gene names

    collapsed_dmrs$type <- i
    saveRDS(collapsed_dmrs,file=paste0(dmr_celltype_by_diag_outdir,"/",i,"/",i,".",sample_name,".celltype_diagnosis.dmrs.collapse.rds"))
},mc.cores=20)


celltype_dmr<-list.files(path=dmr_celltype_by_diag_outdir,full.name=T,recursive=T,pattern=".dmrs.collapse.rds")

celltype_dmr<-do.call("rbind",lapply(celltype_dmr,function(x) readRDS(x)))


#plot dmr counts per clone
pal <- colorRampPalette(colors = c("#F9AB60", "#E7576E", "#630661", "#B5DCA5"))
COLS <- pal(length(unique(celltype_dmr$type)))

plt<-ggplot(celltype_dmr |> dplyr::group_by(type, direction) |> dplyr::summarise(n = n()), 
    aes(y = type, x = n, fill = type)) + 
    geom_col() + 
    facet_grid(vars(direction,type), scales = "free") + 
    scale_fill_manual(values = COLS) + 
    theme_classic()
ggsave(plt,file=paste0(dmr_celltype_by_diag_outdir,"/","celltype_comparisons",".dmr_counts.pdf"))

#plot width of DMRs
celltype_dmr$width<-width(GRanges(celltype_dmr))

plt<-ggplot(celltype_dmr |> dplyr::group_by(type, direction) |> dplyr::summarise(n = n())), 
    aes(y = type, x = n, fill = type)) + 
    geom_col() + 
    facet_grid(vars(direction,type), scales = "free") + 
    scale_fill_manual(values = COLS) + 
    theme_classic()
ggsave(plt,file=paste0(dmr_celltype_by_diag_outdir,"/","celltype_comparisons",".dmr_counts.pdf"))


#ADD HYPERGEO ENRICHMENT HERE GSEA
plot_gsea(obj_sample,
        sample_name=sample_name,collapsed_dmrs=collapsed_dmrs,
        outdir=clone_outdir)


```



```R
#add plotting functions for DMR 
#https://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/pbmc_vignette/pbmc_vignette.html
#add dmr locations per clone over heatmap cnv plot

obj<-readRDS(file="07_scaledcis.cnv_clones.amethyst.rds")

#output to DMR analysis folder
#clones (all lumhr treated as same group)
dmr_clones_outdir=paste(sep="/",dmr_outdir,"cnv_clones_alllumhr")
system(paste("mkdir -p", dmr_clones_outdir))

obj_lumhr_all<-subsetObject(obj,cells=row.names(obj@metadata)[(obj@metadata$celltype %in% c("lumhr","cancer"))])
obj_lumhr_all@metadata$cnv_clonename_alllumhr<-obj_lumhr_all@metadata$cnv_clonename
obj_lumhr_all@metadata[obj_lumhr_all@metadata$celltype=="lumhr",]$cnv_clonename_alllumhr<-"lumhr"
clone500bpwindows<-readRDS(file=paste0(dmr_clones_outdir,"/","dmr_analysis.cnv_clones_alllumhr.500bp_windows.rds"))

#run dmr analysis clone with all lumhr
#use dmr_analysis.cnv_clones_alllumhr.500bp_windows.rds matrix
#note some of these also have the coverage to do specific clone vs lumhr analysis, can merge matrices for this case
clone_v_lumhr_analysis<-c(
'BCMHBCA04R',
'BCMHBCA03R',
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

lapply(clone_v_lumhr_analysis,function(x) dmr_clone_by_lumhr(obj,sample_name=x,dmr_clones_outdir,clone500bpwindows))


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

plot_gsea<-function(obj,
                    sample_name,
                    outdir){
  print(paste("Loading DMRS for sample:",sample_name))

  collapsed_dmrs<-readRDS(file=paste0(dmr_outdir,"/cnv_clones_alllumhr/",sample_name,"/",sample_name,".clones_lumhrall_dmr_filt_collapse.rds"))

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
  hallmark_plt<-gsea_enrichment_againstlumhr(species="human",
              category="C4",
              subcategory="3CA",
              out_setname="3CA",sample_name=sample_name,
              dmrs=collapsed_dmrs,obj=obj,outdir=outdir)
  plt<-patchwork::wrap_plots(list((tft_plt),(position_plt),(hallmark_plt)),ncol=4,axes="collect")

  ggsave(plt,
          file=paste0(sample_name,".clones_lumhrall.GSEA_enrichment.pdf"),
          path=outdir,
          width=40,height=length(unique(collapsed_dmrs$type))*5,limitsize = FALSE)
  print(paste("Finished sample:",sample_name))
}

lapply(clone_v_lumhr_analysis, function(sample_name){
  plot_gsea(obj,
          sample_name=sample_name,
          outdir=paste0(dmr_outdir,"/cnv_clones_alllumhr/",sample_name,"/"))
})


```

#coregulation by gsea https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/geseca-tutorial.html



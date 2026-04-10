Running differential methylation across cell types (one v all strategy)

Running DMR comparisons across different sets.

Celltype comparisons
1. Cell type vs rest (all cells)
2. Cell type vs rest (HBCA only)

Clonal comparisons
4. Clone vs lumhr (per clone, for clones > 30 cells)
5. Cancer vs lumhr (per sample, sum of clones > 30 cells)

Group comparisons
6. Endothelial (HBCA) vs TEC (dcis, synch, idc)
7. Fibroblast (HBCA) vs CAF (dcis, synch, idc)
8. Macrophage (HBCA) vs TAM (dcis, synch, idc)



```R

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
system(paste0("mkdir -p ",project_data_directory,"/DMR_analysis"))

#set colors
celltype_col=c(
"perivascular"="#FF9900",
"fibroblast"="#FF0000",
"endothelial"="#FFFF66",
"unknown"="#FF6699",

"monocyte"="#99FFFF",
"macrophage"="#0066FF",
"macro"="#0066FF",
"bcell"="#0099CC",
"tcell_treg"="#006666",
"tcell_cd4"="#009966",
"tcell_cd8"="#66FF00",

"basal"="#990099",
"lumsec"="#CC0066",
"lumhr"="#FF00CC",
"cancer"="#00FF99")

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

find_cluster_markers<-function(dat,celltype500bp_windows,comp,prefix){
  #comparisons: If eachVsAll is not desired, provide a data frame
   #       describing which tests to run. The data.frame should have
   #       three columns with rows describing conditions of each test.
   #       "name" determines the name of the test in the output; "A"
   #       lists group members, and "B" lists group nonmembers.

  pct_mat<-celltype500bpwindows[["pct_matrix"]] 
  sum_mat<-celltype500bpwindows[["sum_matrix"]] 

  #i dont like how it doesnt keep name in tact, running per row
  dmrs<-testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
                        comparisons = comp, # If TRUE, each group found in the sumMatrix will be tested against all others
                        nminTotal = 3, # Min number observations across all groups to include the region in calculations
                        nminGroup = 3) # Min number observations across either members or nonmembers to include the region
  saveRDS(dmrs,file=paste0(prefix,".unfiltered.500bp.dmrs.rds"))

  dmrs <- filterDMR(dmrs, 
              method = "bonferroni", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
              filter = FALSE, # If TRUE, removes insignificant results
              pThreshold = 0.05, # Maxmimum adjusted p value to allow if filter = TRUE
              logThreshold = 1) # Minimum absolute value of the log2FC to allow if filter = TRUE
  test_name<-setNames(nm=as.data.frame(comp)$name,1:nrow(comp))
  dmrs$test<-names(test_name[dmrs$test])
  saveRDS(dmrs,file=paste0(prefix,".filtered.500bp.dmrs.rds"))

  collapsed_dmrs <- collapseDMR(dat, 
                        dmrs, 
                          maxDist = 1000, # Max allowable overlap between DMRs to be considered adjacent
                          minLength = 500, # Min length of collapsed DMR window to include in the output
                          reduce = T, # Reduce results to unique observations (recommended)
                          annotate = T) # Add column with overlapping gene names
  saveRDS(collapsed_dmrs,file=paste0(prefix,".collapsed.dmrs.rds"))
  return(collapsed_dmrs)
}

#gsea of top DMR genes (hypomet)
gsea_enrichment<-function(dmrs,species="human",
                          category="C3",
                          subcategory="TFT:GTRD",
                          out_setname="TFT",
                          prefix=prefix,
                          obj=obj,
                          sample_name=sample_name){

  print(paste("Calculating DMR overlap with:",out_setname))
  pathwaysDF <- msigdbr(species=species, 
                        collection=category, 
                        subcollection = subcategory)

  #limit pathways to genes in our data
  pathwaysDF<-pathwaysDF[pathwaysDF$ensembl_gene %in% unlist(lapply(strsplit(unique(obj@ref$gene_id),"[.]"),"[",1)),]
  
  pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

  #run plotting per group order by logFC (higher logFC is hyper, lower logFC is hypo)

  fgsea_list<-lapply(unique(dmrs$test),function(group1){
    print(paste("Running:",group1,"..."))
    #treat multiple gene overlaps as same logFC
    #set -Inf to -3 and Inf to 3
    group_features<-dmrs %>%
      dplyr::filter(test == group1) %>%
      dplyr::filter(dmr_padj<0.05) %>% 
      dplyr::filter(gene_names!="NA") %>% 
      dplyr::arrange(dmr_logFC) %>%
      dplyr::select(gene_names, dmr_logFC) %>%
      tidyr::separate_rows(gene_names) %>%
      dplyr::mutate(across(where(is.numeric), ~ replace(., .==-Inf, -3))) %>%
      dplyr::mutate(across(where(is.numeric), ~ replace(., .==Inf, 3))) %>% 
      group_by(gene_names) %>%
      slice_min(order_by=dmr_logFC,n = 1) %>%
      ungroup()

  ranks<-setNames(nm=group_features$gene_names,group_features$dmr_logFC)
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = ranks,
                    minSize  = 5,
                    nproc = 1)

  fgseaRes$test <- group1
  return(fgseaRes)

  })
  fgseaRes<-do.call("rbind",fgsea_list)
  fgseaRes$set_name <- out_setname
  saveRDS(fgseaRes,file=paste0(prefix,".GSEA_enrichment.",sample_name,".",out_setname,".rds"))
}

gsea_across_sets<-function(obj, 
                    collapsed_dmrs,
                    sample_name, 
                    prefix){
  print(paste("Loading DMRS for sample:",sample_name))

  #run gsea enrichment on different sets
  print(paste("Calculating TF Binding Enrichment"))
  tft_gsea<-gsea_enrichment(species="human",
                          category="C3",
                          subcategory="TFT:GTRD",
                          out_setname="TFT",
                          prefix=prefix,
                          sample_name=sample_name,
                          dmrs=collapsed_dmrs,obj=obj)


  print(paste("Calculating Position Enrichment"))
  position_gsea<-gsea_enrichment(species="human",
              category="C1",
              subcategory=NULL,
              out_setname="position",
              prefix=prefix,
              sample_name=sample_name,
              dmrs=collapsed_dmrs,obj=obj)

  print(paste("Calculating Hallmark Enrichment"))
  hallmark_gsea<-gsea_enrichment(species="human",
              category="H",
              subcategory=NULL,
              out_setname="hallmark",              
              prefix=prefix,
              sample_name=sample_name,
              dmrs=collapsed_dmrs,obj=obj)

  print(paste("Calculating Cancer Cell Atlas Enrichment"))
  cancercellatlas_gsea<-gsea_enrichment(species="human",
              category="C4",
              subcategory="3CA",
              out_setname="3CA",
              prefix=prefix,
              sample_name=sample_name,
              dmrs=collapsed_dmrs,obj=obj)

  print(paste("Finished sample:",sample_name))
}

plot_gsea<-function(gsea=hallmark_dmr,out_setname="hallmark",prefix=prefix,dmr_hypo_count,dmr_hyper_count){
  gsea_nes <- gsea %>% tidyr::pivot_wider(names_from=test, id_cols=pathway, values_from=NES)  %>% as.data.frame()
  row.names(gsea_nes) <- gsea_nes$pathway
  gsea_nes<-gsea_nes[,2:ncol(gsea_nes)]

  gsea_pval <- gsea %>% tidyr::pivot_wider(names_from=test, id_cols=pathway, values_from=padj) %>% as.data.frame()
  row.names(gsea_pval) <- gsea_pval$pathway
  gsea_pval<-gsea_pval[,2:ncol(gsea_pval)]
  feature_to_keep<- gsea_pval %>% filter(if_any(everything(), ~ .x < 0.05, na.rm=T)) %>% row.names() #filter hallmark to just columns with signficance
  if(length(feature_to_keep)>1){

  gsea_pval <- -log10(gsea_pval)
  col_fun = circlize::colorRamp2(c(-5, 0, 5), c("#b84d9c","white","#2c50a3"))

  gsea_pval<-gsea_pval[feature_to_keep,]
  gsea_nes<-gsea_nes[feature_to_keep,]
  gsea_nes[which(is.na(gsea_nes),arr.ind=T)]<-0

  #cap pval size for visualization
  max_size=quantile(unlist(gsea_pval),na.rm=T,probs=0.75)
  gsea_pval[which(gsea_pval>max_size,arr.ind=T)]<-max_size

  column_ha = HeatmapAnnotation(
    hyper_n = anno_barplot(dmr_hyper_count[colnames(gsea_nes)],gp = gpar(fill="black",col="black")),
    hypo_n = anno_barplot(dmr_hypo_count[colnames(gsea_nes)],gp = gpar(fill="#FF00FF",col="#FF00FF")))

  if(out_setname=="position"){
  #for plotting position in order
  row.names(gsea_nes)<-gsub(row.names(gsea_nes),pattern="chr",replacement="")
  row_order_chr<-unlist(lapply(strsplit(row.names(gsea_nes),split="p|q"),"[",1))
  row_order_arm<-gsub("\\d", "", row.names(gsea_nes))
  row_order_band<-unlist(lapply(strsplit(row.names(gsea_nes),split="p|q"),"[",2))
  row_order<-data.frame(chr=row_order_chr,arm=row_order_arm,band=row_order_band)
  row_order$chr<-factor(row_order$chr,c(1:22,"X"))
  row_order$arm<-factor(row_order$arm,c("p","q","Xp","Xq"))
  row.names(row_order)<-row.names(gsea_nes)
  row_order<-row_order %>% arrange(chr,arm,band) 
  gsea_nes<-gsea_nes[row.names(row_order),]

  plt<-Heatmap(gsea_nes,
      bottom_annotation=column_ha,
      col = col_fun,rect_gp = gpar(type = "none"),
      cell_fun = function(j, i, x, y, width, height, fill) {
          #draw a rectangle at all sites
          #grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
          if(!is.na(gsea_nes[i, j]) & !is.na(gsea_pval[i, j])){
            #draw a circle sized by NES and colored by pval if significant
            grid.circle(x = x, y = y, 
                  r = (abs(gsea_pval[i, j])/max_size)/2 * min(unit.c(width, height)), 
                  gp = gpar(fill = col_fun(gsea_nes[i, j]), col = NA))}},
    row_order=1:nrow(gsea_nes),
    cluster_rows=FALSE,cluster_columns=TRUE)
  } else {
  plt<-Heatmap(gsea_nes,
      bottom_annotation=column_ha,
      col = col_fun,rect_gp = gpar(type = "none"),
      cell_fun = function(j, i, x, y, width, height, fill) {
          #draw a rectangle at all sites
          #grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
          if(!is.na(gsea_nes[i, j]) & !is.na(gsea_pval[i, j])){
            #draw a circle sized by NES and colored by pval if significant
            grid.circle(x = x, y = y, 
                  r = (abs(gsea_pval[i, j])/max_size)/2 * min(unit.c(width, height)), 
                  gp = gpar(fill = col_fun(gsea_nes[i, j]), col = NA))}},
    cluster_rows=TRUE,cluster_columns=TRUE)
  }

  pdf(paste0(prefix,".",out_setname,".NES.heatmap.pdf"),width=10,height=10)
  print(plt)
  dev.off()
    } else { print(paste("No significant findings in",out_setname))}

}


```


1. DMR characterization per cell type for all cells

```R
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
dmr_output_dir<-paste0(project_data_directory,"/DMR_analysis/","/DMR_celltype")
system(paste0("mkdir -p ",dmr_output_dir))

prefix=paste0(dmr_output_dir,"/","08_scaledcis.final_celltype.onevrest.celltype")

dat<-readRDS(file="08_scaledcis.final_celltype.amethyst.rds") #just saving because of the macrophage name correction
celltype500bpwindows<-readRDS(file="08_scaledcis.final_celltype.500bp_windows.rds") #pre-generated during fine cell typing

comparisons_set=colnames(celltype500bpwindows[["pct_matrix"]])[4:ncol(celltype500bpwindows[["pct_matrix"]])]

#compare each celltype to all other cell types
comp<-lapply(comparisons_set,function(i){
  comparison_A=i
  comparison_B=comparisons_set[!(comparisons_set==i)]
  comparison_name=paste0(i)
  comp<-cbind("name"=comparison_name,
            "A"=i,
            "B"=paste(comparison_B,collapse=","))
  return(comp)}
)
comp<-as.data.frame(do.call("rbind",comp))
row.names(comp)<-comp$name

# Fine DMRs across comparisons
find_cluster_markers(dat=dat,
                  prefix=prefix,
                  celltype500bp_windows=celltype500bp_windows,
                  comp=comp)

#Analyze DMRs
collapsed_dmrs <- readRDS(file=paste0(prefix,".collapsed.dmrs.rds"))

#plot dmr counts per test
plt<-ggplot(collapsed_dmrs |> dplyr::group_by(test, direction) |> dplyr::summarise(n = n()), 
    aes(y = test, x = n, fill = test)) + 
    geom_col() + 
    facet_grid(vars(direction), scales = "free_y") + 
    scale_fill_manual(values = celltype_col) + 
    theme_classic() +  theme(legend.position = "none")
ggsave(plt,file=paste0(prefix,".collapsed.barplot.pdf"))

#gsea for dmr uses rank with hypomet down and hypermet up
gsea_across_sets(obj=dat, 
                  collapsed_dmrs=collapsed_dmrs,
                  sample_name="allcells", 
                  prefix=prefix)


#plot gsea
collapsed_dmrs <- readRDS(file=paste0(prefix,".collapsed.dmrs.rds"))
dmr_hypo_count <- collapsed_dmrs %>% dplyr::filter(dmr_padj < 0.05) %>% dplyr::filter(direction %in% c("hypo")) %>% dplyr::group_by(test, direction) %>% dplyr::summarise(n = n()) %>% as.data.frame()
dmr_hyper_count <- collapsed_dmrs %>% dplyr::filter(dmr_padj < 0.05) %>% dplyr::filter(direction %in% c("hyper")) %>% dplyr::group_by(test, direction) %>% dplyr::summarise(n = n()) %>% as.data.frame()
dmr_hyper_count <- setNames(nm=dmr_hyper_count$test,dmr_hyper_count$n)
dmr_hypo_count <- setNames(nm=dmr_hypo_count$test,dmr_hypo_count$n)


hallmark_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.allcells.","hallmark",".rds"))
plot_gsea(gsea=hallmark_dmr,out_setname="hallmark",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

tft_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.allcells.","TFT",".rds"))
plot_gsea(gsea=tft_dmr,out_setname="TFT",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

pos_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.allcells.","position",".rds"))
plot_gsea(gsea=pos_dmr,out_setname="position",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

cancermeta_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.allcells.","3CA",".rds"))
plot_gsea(gsea=cancermeta_dmr,out_setname="3CA",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)



```

2. HBCA only cells (to find cell type differences without confounding cancer signatures)

```R
#set up directories
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
dmr_output_dir<-paste0(project_data_directory,"/DMR_analysis/","/DMR_celltype_HBCAonly")
system(paste0("mkdir -p ",dmr_output_dir))
prefix=paste0(dmr_output_dir,"/","08_scaledcis.final_celltype.onevrest.HBCAonly_celltype")

#read in files
dat<-readRDS(file="08_scaledcis.final_celltype.amethyst.rds") #just saving because of the macrophage name correction
celltype500bpwindows<-readRDS(file=paste0(project_data_directory,"/merged_data/","08_scaledcis.final_celltype_by_group.500bp_windows.rds")) #pre-generated during fine cell typing

comparisons_set=colnames(celltype500bpwindows[["pct_matrix"]])[4:ncol(celltype500bpwindows[["pct_matrix"]])]
comparisons_set<-comparisons_set[endsWith(comparisons_set,suffix=".HBCA")]
comparisons_set<-comparisons_set[!(comparisons_set=="cancer.HBCA")]

#compare each celltype to all other cell types
comp<-lapply(comparisons_set,function(i){
  comparison_A=i
  comparison_B=comparisons_set[!(comparisons_set==i)]
  comparison_name=paste0(i)
  comp<-cbind("name"=comparison_name,
            "A"=i,
            "B"=paste(comparison_B,collapse=","))
  return(comp)}
)
comp<-as.data.frame(do.call("rbind",comp))
row.names(comp)<-gsub(comp$name,pattern=".HBCA",replace="") #replacing .hbca just for coloring ease

collapsed_dmrs <- find_cluster_markers(dat=dat,
                  prefix=prefix,
                  celltype500bp_windows=celltype500bp_windows,
                  comp=comp)

#plot dmr counts per test
plt<-ggplot(collapsed_dmrs |> dplyr::group_by(test, direction) |> dplyr::summarise(n = n()), 
    aes(y = test, x = n, fill = test)) + 
    geom_col() + 
    facet_grid(vars(direction), scales = "free_y") + 
    scale_fill_manual(values = celltype_col) + 
    theme_classic() +  theme(legend.position = "none")

ggsave(plt,file=paste0(prefix,".collapsed.barplot.pdf"))
  

```

3. Clone vs LumHR comparisons (per clone)

```R
#set up directories
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
dmr_output_dir<-paste0(project_data_directory,"/DMR_analysis/","/DMR_clone_v_lumhr")
system(paste0("mkdir -p ",dmr_output_dir))
prefix=paste0(dmr_output_dir,"/","08_scaledcis.final_celltype.clone_v_lumhr")

#read in files
dat<-readRDS(file="09_scaledcis.final_ploidy.amethyst.rds")

lumhr_hbca_cells <- row.names(dat@metadata[dat@metadata$Group=="HBCA" & dat@metadata$celltype=="lumhr",])

clones_passing_filter<-names(which(table(dat@metadata$cnv_clonename_500kb)>1)) #apply clone filters after the fact
clones_passing_filter <- clones_passing_filter[!(endsWith(clones_passing_filter,suffix="_diploid"))]
clones_passing_filter <- clones_passing_filter[clones_passing_filter!="NA"]

clone_cells <- row.names(dat@metadata[dat@metadata$Group!="HBCA" & 
                                                dat@metadata$celltype %in% c("lumhr","cancer") &
                                                dat@metadata$cnv_clonename_500kb %in% clones_passing_filter,])

dat<-subsetObject(dat,cells=c(lumhr_hbca_cells,clone_cells))
dat@metadata[lumhr_hbca_cells,]$cnv_clonename_500kb<-"HBCA_lumhr"
table(dat@metadata$cnv_clonename_500kb)

celltype500bpwindows <- calcSmoothedWindows(dat, 
                                        type = "CG", 
                                        threads = 200,
                                        step = 500, 
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "cnv_clonename_500kb",
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)
saveRDS(celltype500bpwindows,file=paste0(prefix,".500bp_windows.rds"))


celltype500bpwindows<-readRDS(file=paste0(prefix,".500bp_windows.rds"))
comparisons_set=colnames(celltype500bpwindows[["pct_matrix"]])[4:ncol(celltype500bpwindows[["pct_matrix"]])]

#compare each celltype to HBCA_lumhr
comp<-lapply(comparisons_set,function(i){
  if(i!="HBCA_lumhr"){
    comparison_A=i
    comparison_B="HBCA_lumhr"
    comparison_name=paste0(i)
    comp<-cbind("name"=comparison_name,
              "A"=i,
              "B"=comparison_B)
  return(comp)}}
)
comp<-as.data.frame(do.call("rbind",comp))
row.names(comp)<-comp$name #replacing .hbca just for coloring ease

collapsed_dmrs <- find_cluster_markers(dat=dat,
                  prefix=prefix,
                  celltype500bp_windows=celltype500bp_windows,
                  comp=comp)
                  #running


#Analyze DMRs
collapsed_dmrs <- readRDS(file=paste0(prefix,".collapsed.dmrs.rds"))

#plot dmr counts per test
plt<-ggplot(collapsed_dmrs |> dplyr::group_by(test, direction) |> dplyr::summarise(n = n()), 
    aes(y = test, x = n, fill = test)) + 
    geom_col() + 
    facet_grid(vars(direction), scales = "free_y") + 
    theme_classic() +  theme(legend.position = "none")
ggsave(plt,file=paste0(prefix,".collapsed.barplot.pdf"))

#gsea for dmr uses rank with hypomet down and hypermet up
gsea_across_sets(obj=dat, 
                  collapsed_dmrs=collapsed_dmrs,
                  sample_name="clone_v_lumhr", 
                  prefix=prefix)

#plot gsea
collapsed_dmrs <- readRDS(file=paste0(prefix,".collapsed.dmrs.rds"))
dmr_hypo_count <- collapsed_dmrs %>% dplyr::filter(dmr_padj < 0.05) %>% dplyr::filter(direction %in% c("hypo")) %>% dplyr::group_by(test, direction) %>% dplyr::summarise(n = n()) %>% as.data.frame()
dmr_hyper_count <- collapsed_dmrs %>% dplyr::filter(dmr_padj < 0.05) %>% dplyr::filter(direction %in% c("hyper")) %>% dplyr::group_by(test, direction) %>% dplyr::summarise(n = n()) %>% as.data.frame()
dmr_hyper_count <- setNames(nm=dmr_hyper_count$test,dmr_hyper_count$n)
dmr_hypo_count <- setNames(nm=dmr_hypo_count$test,dmr_hypo_count$n)


hallmark_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.clone_v_lumhr.","hallmark",".rds"))
hallmark_dmr <- hallmark_dmr %>% filter(test %in% c("BCMDCIS41T_c1","BCMDCIS41T_c2","BCMDCIS41T_c3"))
plot_gsea(gsea=hallmark_dmr,out_setname="hallmark",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

tft_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.clone_v_lumhr.","TFT",".rds"))
plot_gsea(gsea=tft_dmr,out_setname="TFT",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

pos_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.clone_v_lumhr.","position",".rds"))
pos_dmr <- pos_dmr %>% filter(test %in% c("BCMDCIS41T_c1","BCMDCIS41T_c2","BCMDCIS41T_c3"))
plot_gsea(gsea=pos_dmr,out_setname="position",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

cancermeta_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.clone_v_lumhr.","3CA",".rds"))
plot_gsea(gsea=cancermeta_dmr,out_setname="3CA",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

```


```R
#set up directories
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
dmr_output_dir<-paste0(project_data_directory,"/DMR_analysis/","/DMR_BCMDCIS41T")
system(paste0("mkdir -p ",dmr_output_dir))
prefix=paste0(dmr_output_dir,"/","08_scaledcis.final_celltype.BCMDCIS41T")

#using data split by clone
dat<-readRDS(file="09_scaledcis.final_ploidy.amethyst.rds")
celltype500bpwindows<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/DMR_analysis//DMR_clone_v_lumhr/08_scaledcis.final_celltype.clone_v_lumhr.500bp_windows.rds")


comp_c1<-as.data.frame(cbind("name"="BCMDCIS41T_c1",
          "A"="BCMDCIS41T_c1",
          "B"=paste(c("BCMDCIS41T_c2","BCMDCIS41T_c3"),collapse=",")))
          
comp_c2<-as.data.frame(cbind("name"="BCMDCIS41T_c2",
          "A"="BCMDCIS41T_c2",
          "B"=paste(c("BCMDCIS41T_c1","BCMDCIS41T_c3"),collapse=",")))

comp_c3<-as.data.frame(cbind("name"="BCMDCIS41T_c3",
          "A"="BCMDCIS41T_c3",
          "B"=paste(c("BCMDCIS41T_c1","BCMDCIS41T_c2"),collapse=",")))

comp<-rbind(comp_c1,comp_c2,comp_c3)

collapsed_dmrs <- find_cluster_markers(dat=dat,
                  prefix=prefix,
                  celltype500bp_windows=celltype500bpwindows,
                  comp=comp)
                  #running

#Analyze DMRs
collapsed_dmrs <- readRDS(file=paste0(prefix,".collapsed.dmrs.rds"))

#plot dmr counts per test
plt<-ggplot(collapsed_dmrs |> dplyr::group_by(test, direction) |> dplyr::filter(dmr_padj<0.05)|> dplyr::summarise(n = n()), 
    aes(y = test, x = n, fill = test)) + 
    geom_col() + 
    facet_grid(vars(direction), scales = "free_y") + 
    theme_classic() +  theme(legend.position = "none")
ggsave(plt,file=paste0(prefix,".collapsed.barplot.pdf"))

#gsea for dmr uses rank with hypomet down and hypermet up
sample_name="BCMDCIS41T"
gsea_across_sets(obj=dat, 
                  collapsed_dmrs=collapsed_dmrs,
                  sample_name=sample_name, 
                  prefix=prefix)

#plot gsea requires multiple comparisons
collapsed_dmrs <- readRDS(file=paste0(prefix,".collapsed.dmrs.rds"))
dmr_hypo_count <- collapsed_dmrs %>% dplyr::filter(dmr_padj < 0.05) %>% dplyr::filter(direction %in% c("hypo")) %>% dplyr::group_by(test, direction) %>% dplyr::summarise(n = n()) %>% as.data.frame()
dmr_hyper_count <- collapsed_dmrs %>% dplyr::filter(dmr_padj < 0.05) %>% dplyr::filter(direction %in% c("hyper")) %>% dplyr::group_by(test, direction) %>% dplyr::summarise(n = n()) %>% as.data.frame()
dmr_hyper_count <- setNames(nm=dmr_hyper_count$test,dmr_hyper_count$n)
dmr_hypo_count <- setNames(nm=dmr_hypo_count$test,dmr_hypo_count$n)


hallmark_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.",sample_name,".","hallmark",".rds"))
lot_gsea(gsea=hallmark_dmr,out_setname="hallmark",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

tft_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.",sample_name,".","TFT",".rds"))
plot_gsea(gsea=tft_dmr,out_setname="TFT",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

pos_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.",sample_name,".","position",".rds"))
plot_gsea(gsea=pos_dmr,out_setname="position",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

cancermeta_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.",sample_name,".","3CA",".rds"))
plot_gsea(gsea=cancermeta_dmr,out_setname="3CA",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)


```

Additional comparison of clones based on large scale shared events

```R
#set up directories
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
dmr_output_dir<-paste0(project_data_directory,"/DMR_analysis/","/DMR_WGDclone_v_nonWGDclone")
system(paste0("mkdir -p ",dmr_output_dir))
prefix=paste0(dmr_output_dir,"/","08_scaledcis.final_celltype.WGDclone_v_nonWGDclone")

#using data split by clone
dat<-readRDS(file="09_scaledcis.final_ploidy.amethyst.rds")
celltype500bpwindows<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/DMR_analysis//DMR_clone_v_lumhr/08_scaledcis.final_celltype.clone_v_lumhr.500bp_windows.rds")

lumhr_hbca_cells <- row.names(dat@metadata[dat@metadata$Group=="HBCA" & dat@metadata$celltype=="lumhr",])
clones_passing_filter<-names(which(table(dat@metadata$cnv_clonename_500kb)>1))
clones_passing_filter <- clones_passing_filter[!(endsWith(clones_passing_filter,suffix="_diploid"))]
clones_passing_filter <- clones_passing_filter[clones_passing_filter!="NA"]
clone_cells <- row.names(dat@metadata[dat@metadata$Group!="HBCA" & 
                                                dat@metadata$celltype %in% c("lumhr","cancer") &
                                                dat@metadata$cnv_clonename_500kb %in% clones_passing_filter,])
dat<-subsetObject(dat,cells=c(lumhr_hbca_cells,clone_cells))
dat@metadata[lumhr_hbca_cells,]$cnv_clonename_500kb<-"HBCA_lumhr"
table(dat@metadata$cnv_clonename_500kb)

clones_passing_filter <- colnames(celltype500bpwindows[["pct_matrix"]])[4:ncol(celltype500bpwindows[["pct_matrix"]])]


#wgd
wgd_clones=c("BCMDCIS124T_c1","BCMDCIS124T_c2","BCMDCIS124T_c3","BCMDCIS80T_24hTis_c1","ECIS25T_c1","ECIS25T_c2","ECIS25T_c3")
wgd_clones <- wgd_clones[wgd_clones %in% clones_passing_filter]
nonwgd_clones<-clones_passing_filter[which(!(clones_passing_filter %in% wgd_clones))]
nonwgd_clones<-nonwgd_clones[nonwgd_clones %in% clones_passing_filter]
nonwgd_clones<-nonwgd_clones[!grepl(pattern="HBCA",nonwgd_clones)]

#compare each celltype to HBCA_lumhr
comparison_A=paste(wgd_clones,collapse=",")
comparison_B=paste(nonwgd_clones,collapse=",")
comp_wgd<-as.data.frame(cbind("name"="WGD_v_nonWGD",
          "A"=comparison_A,
          "B"=comparison_B))

#1q gain
clones_1qgain<-c("BCMDCIS28T_c2", "BCMDCIS22T_c1", "BCMDCIS52T_c1", "ECIS36T_c3", "ECIS36T_c1", "ECIS36T_c2", "BCMDCIS102T_24hTis_c2", "BCMDCIS102T_24hTis_c3", "BCMDCIS74T_c3", "BCMDCIS28T_c1", "BCMDCIS70T_c1", "BCMDCIS35T_c2", "BCMDCIS35T_c1", "BCMDCIS74T_c2", "BCMDCIS94T_24hTis_c2", "BCMDCIS94T_24hTis_c1", "ECIS26T_c1", "ECIS26T_c2", "BCMDCIS102T_24hTis_c1", "BCMHBCA03R_c1", "BCMDCIS79T_24hTis_DCIS_c1", "BCMDCIS99T_c1", "BCMDCIS41T_c1", "BCMDCIS41T_c2", "BCMDCIS97T_c5", "BCMDCIS97T_c3", "BCMDCIS97T_c4", "BCMDCIS70T_c2", "BCMDCIS79T_24hTis_DCIS_c2", "BCMDCIS92T_24hTis_c1", "BCMDCIS97T_c1", "BCMDCIS97T_c2", "BCMDCIS80T_24hTis_c2")
clones_1qgain<-clones_1qgain[which(clones_1qgain %in% clones_passing_filter)]

nonclones_1qgain<-clones_passing_filter[which(!(clones_passing_filter %in% clones_1qgain))]
#compare each celltype to HBCA_lumhr
comparison_A=paste(clones_1qgain,collapse=",")
comparison_B=paste(nonclones_1qgain,collapse=",")
comp_1q<-as.data.frame(cbind("name"="1qgain_v_non1qgain",
          "A"=comparison_A,
          "B"=comparison_B))

#16q loss
clones_16qloss<-c("ECIS36T_c1", "ECIS36T_c2", "BCMDCIS74T_c3", "BCMDCIS28T_c1", "BCMDCIS70T_c1", "BCMDCIS35T_c2", "BCMDCIS35T_c1", "BCMDCIS74T_c2", "BCMDCIS94T_24hTis_c2", "BCMDCIS94T_24hTis_c1", "ECIS26T_c1", "ECIS26T_c2", "BCMDCIS74T_c1", "BCMDCIS41T_c3", "BCMDCIS66T_c1", "BCMDCIS41T_c1", "BCMDCIS41T_c2", "BCMDCIS97T_c5", "BCMDCIS97T_c3", "BCMDCIS97T_c4", "BCMDCIS70T_c2", "BCMDCIS79T_24hTis_DCIS_c2", "BCMDCIS92T_24hTis_c1", "BCMDCIS97T_c1")
clones_16qloss<-clones_16qloss[which(clones_16qloss %in% clones_passing_filter)]

nonclones_16qloss<-clones_passing_filter[which(!(clones_passing_filter %in% clones_16qloss))]
comparison_A=paste(clones_16qloss,collapse=",")
comparison_B=paste(nonclones_16qloss,collapse=",")
comp_16q<-as.data.frame(cbind("name"="16qloss_v_non16qloss",
          "A"=comparison_A,
          "B"=comparison_B))

comp<-as.data.frame(rbind(comp_wgd,comp_1q,comp_16q))

collapsed_dmrs <- find_cluster_markers(dat=dat,
                  prefix=prefix,
                  celltype500bp_windows=celltype500bpwindows,
                  comp=comp)
                  #running

#Analyze DMRs
collapsed_dmrs <- readRDS(file=paste0(prefix,".collapsed.dmrs.rds"))

#plot dmr counts per test
plt<-ggplot(collapsed_dmrs |> dplyr::group_by(test, direction) |> dplyr::filter(dmr_padj<0.05)|> dplyr::summarise(n = n()), 
    aes(y = test, x = n, fill = test)) + 
    geom_col() + 
    facet_grid(vars(direction), scales = "free_y") + 
    theme_classic() +  theme(legend.position = "none")
ggsave(plt,file=paste0(prefix,".collapsed.barplot.pdf"))

#gsea for dmr uses rank with hypomet down and hypermet up
sample_name="cnv_groups"
gsea_across_sets(obj=dat, 
                  collapsed_dmrs=collapsed_dmrs,
                  sample_name=sample_name, 
                  prefix=prefix)

#plot gsea requires multiple comparisons
collapsed_dmrs <- readRDS(file=paste0(prefix,".collapsed.dmrs.rds"))
dmr_hypo_count <- collapsed_dmrs %>% dplyr::filter(dmr_padj < 0.05) %>% dplyr::filter(direction %in% c("hypo")) %>% dplyr::group_by(test, direction) %>% dplyr::summarise(n = n()) %>% as.data.frame()
dmr_hyper_count <- collapsed_dmrs %>% dplyr::filter(dmr_padj < 0.05) %>% dplyr::filter(direction %in% c("hyper")) %>% dplyr::group_by(test, direction) %>% dplyr::summarise(n = n()) %>% as.data.frame()
dmr_hyper_count <- setNames(nm=dmr_hyper_count$test,dmr_hyper_count$n)
dmr_hypo_count <- setNames(nm=dmr_hypo_count$test,dmr_hypo_count$n)


hallmark_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.",sample_name,".","hallmark",".rds"))
plot_gsea(gsea=hallmark_dmr,out_setname="hallmark",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

tft_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.",sample_name,".","TFT",".rds"))
plot_gsea(gsea=tft_dmr,out_setname="TFT",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

pos_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.",sample_name,".","position",".rds"))
plot_gsea(gsea=pos_dmr,out_setname="position",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

cancermeta_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.",sample_name,".","3CA",".rds"))
plot_gsea(gsea=cancermeta_dmr,out_setname="3CA",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)


```


3. Cancer vs LumHR comparisons (per sample)

TO MODIFY AND RUN
RUN BASED ON SAMPLE+PLOIDY

Of 41 samples
-13 HBCA
-9 DCIS
-6 SYNCH
-12 IDC

Samples with enough cells for cancer v lumhr (min 30 cancer cells by CNV)
-5 DCIS
-3 SYNCH
-6 IDC

```R
#set up directories
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
dmr_output_dir<-paste0(project_data_directory,"/DMR_analysis/","/DMR_cancer_v_lumhr")
system(paste0("mkdir -p ",dmr_output_dir))
prefix=paste0(dmr_output_dir,"/","08_scaledcis.final_celltype.cancer_v_lumhr")

#read in files
dat<-readRDS(file="09_scaledcis.final_ploidy.amethyst.rds")

lumhr_hbca_cells <- row.names(dat@metadata[dat@metadata$Group=="HBCA" & dat@metadata$celltype=="lumhr",])

dat@metadata$cnv_sampleploidy_500kb <- paste0(dat@metadata$Sample,"_",dat@metadata$ploidy)

cancer_passing_filter <- names(which(table(dat@metadata$cnv_sampleploidy_500kb)>30))
cancer_passing_filter <- cancer_passing_filter[!(endsWith(cancer_passing_filter,suffix="_diploid"))]
cancer_passing_filter <- cancer_passing_filter[!(endsWith(cancer_passing_filter,suffix="NA"))]

cancer_cells <- row.names(dat@metadata[dat@metadata$celltype %in% c("lumhr","cancer") &
                                      dat@metadata$cnv_sampleploidy_500kb %in% cancer_passing_filter,])

dat<-subsetObject(dat,cells=c(lumhr_hbca_cells,cancer_cells))
dat@metadata[lumhr_hbca_cells,]$cnv_sampleploidy_500kb<-"HBCA_lumhr"

celltype500bpwindows <- calcSmoothedWindows(dat, 
                                        type = "CG", 
                                        threads = 200,
                                        step = 500, 
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "cnv_sampleploidy_500kb",
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)
saveRDS(celltype500bpwindows,file=paste0(prefix,".500bp_windows.rds"))

comparisons_set=colnames(celltype500bpwindows[["pct_matrix"]])[4:ncol(celltype500bpwindows[["pct_matrix"]])]

#compare each celltype to HBCA_lumhr
comp<-lapply(comparisons_set,function(i){
  if(i!="HBCA_lumhr"){
    comparison_A=i
    comparison_B="HBCA_lumhr"
    comparison_name=paste0(i)
    comp<-cbind("name"=comparison_name,
              "A"=i,
              "B"=comparison_B)
  return(comp)}}
)
comp<-as.data.frame(do.call("rbind",comp))
row.names(comp)<-comp$name #replacing .hbca just for coloring ease

collapsed_dmrs <- find_cluster_markers(dat=dat,
                  prefix=prefix,
                  celltype500bp_windows=celltype500bp_windows,
                  comp=comp)
                  #running

#Analyze DMRs
collapsed_dmrs <- readRDS(file=paste0(prefix,".collapsed.dmrs.rds"))

#plot dmr counts per test
plt<-ggplot(collapsed_dmrs |> dplyr::group_by(test, direction) |> dplyr::summarise(n = n()), 
    aes(y = test, x = n, fill = test)) + 
    geom_col() + 
    facet_grid(vars(direction), scales = "free_y") + 
    theme_classic() +  theme(legend.position = "none")
ggsave(plt,file=paste0(prefix,".collapsed.barplot.pdf"))

#gsea for dmr uses rank with hypomet down and hypermet up
gsea_across_sets(obj=dat, 
                  collapsed_dmrs=collapsed_dmrs,
                  sample_name="cancer_v_lumhr", 
                  prefix=prefix)

#plot gsea
collapsed_dmrs <- readRDS(file=paste0(prefix,".collapsed.dmrs.rds"))
dmr_hypo_count <- collapsed_dmrs %>% dplyr::filter(dmr_padj < 0.05) %>% dplyr::filter(direction %in% c("hypo")) %>% dplyr::group_by(test, direction) %>% dplyr::summarise(n = n()) %>% as.data.frame()
dmr_hyper_count <- collapsed_dmrs %>% dplyr::filter(dmr_padj < 0.05) %>% dplyr::filter(direction %in% c("hyper")) %>% dplyr::group_by(test, direction) %>% dplyr::summarise(n = n()) %>% as.data.frame()
dmr_hyper_count <- setNames(nm=dmr_hyper_count$test,dmr_hyper_count$n)
dmr_hypo_count <- setNames(nm=dmr_hypo_count$test,dmr_hypo_count$n)


hallmark_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.cancer_v_lumhr.","hallmark",".rds"))
plot_gsea(gsea=hallmark_dmr,out_setname="hallmark",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

tft_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.cancer_v_lumhr.","TFT",".rds"))
plot_gsea(gsea=tft_dmr,out_setname="TFT",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

pos_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.cancer_v_lumhr.","position",".rds"))
plot_gsea(gsea=pos_dmr,out_setname="position",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

cancermeta_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.cancer_v_lumhr.","3CA",".rds"))
plot_gsea(gsea=cancermeta_dmr,out_setname="3CA",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

```

6-8. Group comparisons
Group comparisons
6. Endothelial (HBCA) vs TEC (dcis, synch, idc)
7. Fibroblast (HBCA) vs CAF (dcis, synch, idc)
8. Macrophage (HBCA) vs TAM (dcis, synch, idc)

```R
#set up directories
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
dmr_output_dir<-paste0(project_data_directory,"/DMR_analysis/","/DMR_celltype_bygroup")
system(paste0("mkdir -p ",dmr_output_dir))
prefix=paste0(dmr_output_dir,"/","08_scaledcis.final_celltype.onevrest.celltype_bygroup")

#read in files
dat<-readRDS(file="09_scaledcis.final_ploidy.amethyst.rds")
celltype500bpwindows<-readRDS(file=paste0(project_data_directory,"/merged_data/","08_scaledcis.final_celltype_by_group.500bp_windows.rds"))

comparisons_set=colnames(celltype500bpwindows[["pct_matrix"]])[4:ncol(celltype500bpwindows[["pct_matrix"]])]

endothelial_comp<-cbind("name"="endothelial_HBCA_v_cancer",
            "A"="endothelial.HBCA",
            "B"=paste(c("endothelial.IDC","endothelial.DCIS","endothelial.Synchronous"),collapse=","))
fibroblast_comp<-cbind("name"="fibroblast_HBCA_v_cancer",
            "A"="fibroblast.HBCA",
            "B"=paste(c("fibroblast.IDC","fibroblast.DCIS","fibroblast.Synchronous"),collapse=","))
macro_comp<-cbind("name"="macro_HBCA_v_cancer",
            "A"="macro.HBCA",
            "B"=paste(c("macro.IDC","macro.DCIS","macro.Synchronous"),collapse=","))
basal_comp<-cbind("name"="basal_HBCA_v_DCIS",
            "A"="basal.HBCA",
            "B"="basal.DCIS")
lumhr_to_dcis_comp<-cbind("name"="lumhr_HBCA_v_cancer_DCIS",
            "A"="lumhr.HBCA",
            "B"="cancer.DCIS")
dcis_to_synch_comp<-cbind("name"="cancer_DCIS_v_cancer_SYNCH",
            "A"="cancer.DCIS",
            "B"="cancer.Synchronous")
synch_to_idc_comp<-cbind("name"="cancer_SYNCH_v_cancer_IDC",
            "A"="cancer.Synchronous",
            "B"="cancer.IDC")
dcis_to_idc_comp<-cbind("name"="cancer_DCIS_v_cancer_IDC",
            "A"="cancer.DCIS",
            "B"="cancer.IDC")
lumhr_to_idc_comp<-cbind("name"="lumhr_HBCA_v_cancer_IDC",
            "A"="lumhr.HBCA",
            "B"="cancer.IDC")
comp<-rbind(endothelial_comp,fibroblast_comp,macro_comp,basal_comp,
lumhr_to_dcis_comp,dcis_to_synch_comp,synch_to_idc_comp,dcis_to_idc_comp,lumhr_to_idc_comp)

collapsed_dmrs <- find_cluster_markers(dat=dat,
                  prefix=prefix,
                  celltype500bp_windows=celltype500bp_windows,
                  comp=comp)

#plot dmr counts per test
collapsed_dmrs <- collapsed_dmrs %>% dplyr::filter(test %in% c("basal_HBCA_v_DCIS","endothelial_HBCA_v_cancer","fibroblast_HBCA_v_cancer","macro_HBCA_v_cancer")) 

collapsed_dmrs %>% dplyr::filter(test %in% c("basal_HBCA_v_DCIS","endothelial_HBCA_v_cancer","fibroblast_HBCA_v_cancer","macro_HBCA_v_cancer")) %>% dplyr::filter(dmr_padj<0.05) %>% group_by(test,direction) %>% summarise(sum(dmr_length))
plt<-ggplot(collapsed_dmrs |> dplyr::filter(dmr_padj<0.05) |> dplyr::group_by(test, direction) |> dplyr::summarise(n = n()), 
    aes(y = test, x = n, fill = test)) + 
    geom_col() + 
    facet_grid(vars(direction), scales = "free_y") + 
    theme_classic() +  theme(legend.position = "none")
ggsave(plt,file=paste0(prefix,".collapsed.barplot.pdf"))

sample_name="celltype_by_group"
#gsea for dmr uses rank with hypomet down and hypermet up
gsea_across_sets(obj=dat, 
                  collapsed_dmrs=collapsed_dmrs,
                  sample_name=sample_name, 
                  prefix=prefix)

#plot gsea
collapsed_dmrs <- readRDS(file=paste0(prefix,".collapsed.dmrs.rds"))
collapsed_dmrs <- collapsed_dmrs %>% dplyr::filter(test %in% c("basal_HBCA_v_DCIS","endothelial_HBCA_v_cancer","fibroblast_HBCA_v_cancer","macro_HBCA_v_cancer")) 
dmr_hypo_count <- collapsed_dmrs %>% dplyr::filter(dmr_padj < 0.05) %>% dplyr::filter(direction %in% c("hypo")) %>% dplyr::group_by(test, direction) %>% dplyr::summarise(n = n()) %>% as.data.frame()
dmr_hyper_count <- collapsed_dmrs %>% dplyr::filter(dmr_padj < 0.05) %>% dplyr::filter(direction %in% c("hyper")) %>% dplyr::group_by(test, direction) %>% dplyr::summarise(n = n()) %>% as.data.frame()
dmr_hyper_count <- setNames(nm=dmr_hyper_count$test,dmr_hyper_count$n)
dmr_hypo_count <- setNames(nm=dmr_hypo_count$test,dmr_hypo_count$n)


hallmark_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.",sample_name,".","hallmark",".rds"))
plot_gsea(gsea=hallmark_dmr,out_setname="hallmark",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

tft_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.",sample_name,".","TFT",".rds"))
tft_dmr <- tft_dmr %>% dplyr::filter(test %in% c("basal_HBCA_v_DCIS","endothelial_HBCA_v_cancer","fibroblast_HBCA_v_cancer","macro_HBCA_v_cancer")) 
plot_gsea(gsea=tft_dmr,out_setname="TFT",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

pos_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.",sample_name,".","position",".rds"))
plot_gsea(gsea=pos_dmr,out_setname="position",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)

cancermeta_dmr <- readRDS(paste0(prefix,".GSEA_enrichment.",sample_name,".","3CA",".rds"))
plot_gsea(gsea=cancermeta_dmr,out_setname="3CA",prefix=prefix,dmr_hypo_count=dmr_hypo_count,dmr_hyper_count=dmr_hyper_count)


```




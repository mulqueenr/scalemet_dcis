
Running differential methylation across clones.

Running DMR comparisons across different sets.

1. For samples with multiple clones (clones > 50 cells):
- Run one v one and one v rest comparisons within sample

2. For all samples with >=1 clone (clones > 50 cells):
- Run clone vs lumhr_HBCA set

3. For classes of 1q gain, 16q gain, 16p loss:
- Run aggregate clones groups against lumhr_HBCA

4. For classes of WGD:
- Run WGD groups against lumhr_HBCA

```R

set.seed(111)
options(future.globals.maxSize= 500000*1024^2) #80gb limit for parallelizing
library(amethyst)
library(data.table)
library(dplyr)
library(msigdbr)
library(fgsea)
library(ggplot2)
library(GenomicRanges)
library(Matrix)
library(parallel)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(GeneOverlap)
library(matrixStats)
library(sparseMatrixStats)
library(irlba)
library(rtracklayer)

project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"

#read in object from directory
task_cpus=300
processing_folder="04_dmr"
wd=paste(sep="/",project_data_directory,processing_folder)
system(paste0("mkdir -p ",wd))
setwd(wd)


```
Functions
```R
#run DMR
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


comparison_set_dmr<-function(i){
  #calculate DMRs for group-wise comparisons across clones
  #run DMR analysis per row in comparisons
  print(paste("Running DMRs across cell types for:",i))
  dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
          comparisons=comparisons[i,],
          nminTotal = 3, # Min number observations across all groups to include the region in calculations
          nminGroup = 3) # Min number observations across either members or nonmembers to include the region
  dmrs$type <- i
  saveRDS(dmrs,paste0(output_folder,"/","04.1.3.dmr.",suffix,".",i,".dmr.rds"))
  print(paste("Filtering DMRs across clones for:",i))
  dmrs <- filterDMR(dmrs, 
              method = "bonferroni", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
              filter = TRUE, # If TRUE, removes insignificant results
              pThreshold = 0.05 # Maxmimum adjusted p value to allow if filter = TRUE
              ) # Minimum absolute value of the log2FC to allow if filter = TRUE
  dmrs$type <- i
  saveRDS(dmrs,paste0(output_folder,"/","04.1.3.dmr.",suffix,".",i,".dmr.filt.rds"))
  collapsed_dmrs <- collapseDMR(obj, 
                          dmrs, 
                          maxDist = 2000, # Max allowable overlap between DMRs to be considered adjacent
                          minLength = 500, # Min length of collapsed DMR window to include in the output
                          reduce = T, # Reduce results to unique observations (recommended)
                          annotate = T) # Add column with overlapping gene names
  collapsed_dmrs$type <- i
  saveRDS(collapsed_dmrs,paste0(output_folder,"/","04.1.3.dmr.",suffix,".",i,".dmr.collapse.rds"))
  return(collapsed_dmrs)
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
                          minLength = 1, # Min length of collapsed DMR window to include in the output
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

  fgsea_list<-lapply(unique(dmrs$comparison_name),function(group1){
    print(paste("Running:",group1,"..."))
    #treat multiple gene overlaps as same logFC
    #set -Inf to -3 and Inf to 3
    group_features<-dmrs %>%
      dplyr::filter(comparison_name == group1) %>%
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

plot_gsea<-function(gsea=hallmark_dmr,
                    out_setname="hallmark",
                    prefix=prefix){

  gsea_nes <- gsea %>% tidyr::pivot_wider(names_from=test, id_cols=pathway, values_from=NES)  %>% as.data.frame()
  row.names(gsea_nes) <- gsea_nes$pathway
  gsea_nes<-gsea_nes[,2:ncol(gsea_nes)]

  gsea_pval <- gsea %>% tidyr::pivot_wider(names_from=test, id_cols=pathway, values_from=padj) %>% as.data.frame()
  row.names(gsea_pval) <- gsea_pval$pathway
  gsea_pval<-gsea_pval[,2:ncol(gsea_pval)]
  feature_to_keep<- gsea_pval %>% filter(if_any(everything(), ~ .x < 0.05, na.rm=T)) %>% row.names() #filter hallmark to just columns with signficance
  
  if(length(feature_to_keep)>0){
  gsea_pval <- -log10(gsea_pval)
  col_fun = circlize::colorRamp2(c(-5, 0, 5), c("#b84d9c","white","#2c50a3"))

  gsea_pval<-gsea_pval[feature_to_keep,]
  gsea_nes<-gsea_nes[feature_to_keep,]
  gsea_nes[which(is.na(gsea_nes),arr.ind=T)]<-0

  #cap pval size for visualization
  max_size=quantile(unlist(gsea_pval),na.rm=T,probs=0.75)
  gsea_pval[which(gsea_pval>max_size,arr.ind=T)]<-max_size

  column_order=factor(colnames(gsea_nes), levels=comparison_order_for_plotting)
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
  plt<-Heatmap(t(gsea_nes), #adding transpose for plotting, so column and row swapped
      col = col_fun,rect_gp = gpar(type = "none"),
      cell_fun = function(j, i, x, y, width, height, fill) {
          #draw a rectangle at all sites
          #grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
          if(!is.na(gsea_nes[i, j]) & !is.na(gsea_pval[i, j])){
            #draw a circle sized by NES and colored by pval if significant
            grid.circle(x = x, y = y, 
                  r = (abs(gsea_pval[i, j])/max_size)/2 * min(unit.c(width, height)), 
                  gp = gpar(fill = col_fun(gsea_nes[i, j]), col = NA))}},
    column_order=1:nrow(gsea_nes),
    row_order=column_order)
  
  } else {

  plt<-Heatmap(t(gsea_nes), #adding transpose for plotting, so column and row swapped
      col = col_fun,rect_gp = gpar(type = "none"),
      cell_fun = function(j, i, x, y, width, height, fill) {
          #draw a rectangle at all sites
          #grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
          if(!is.na(gsea_nes[i, j]) & !is.na(gsea_pval[i, j])){
            #draw a circle sized by NES and colored by pval if significant
            grid.circle(x = x, y = y, 
                  r = (abs(gsea_pval[i, j])/max_size)/2 * min(unit.c(width, height)), 
                  gp = gpar(fill = col_fun(gsea_nes[i, j]), col = NA))}},
    cluster_columns=TRUE,
    row_order=column_order)
  }

  pdf(paste0(prefix,".",out_setname,".NES.heatmap.pdf"),width=10,height=10)
  print(plt)
  dev.off()
    } else { print(paste("No significant findings in",out_setname))}

}

```
Get amethyst object situated.

```R

output_directory=paste0(wd,"/","dmr_celltype_group")
system(paste0("mkdir -p ",output_directory))

obj<-readRDS(file=paste(project_data_directory,"03_fine_celltyping","03_scaledcis.final_celltypes.amethyst.rds",sep="/"))

dat<-subsetObject(obj,cells=row.names(obj@metadata[obj@metadata$celltype %in% c("cancer","lumhr"),]))
dat<-subsetObject(dat,cells=row.names(dat@metadata[!(dat@metadata$cnv_clonename_500kb %in% c("NA")),]))

#add in lumhr as a class to also act as baseline
dat@metadata[(dat@metadata$Group=="HBCA") & (dat@metadata$celltype=="lumhr"),]$cnv_clonename_500kb<-"HBCA_lumhr"

#require at least 50 cells for clones passing filter
clones_passing_filter<-table(dat@metadata$cnv_clonename_500kb)[table(dat@metadata$cnv_clonename_500kb)>50]

dat<-subsetObject(dat,cells=row.names(dat@metadata[dat@metadata$cnv_clonename_500kb %in% names(clones_passing_filter),]))
```

Read in precalculated matrix 
See scalemet_dcis/processing/milestonev1_04.1.1_DMR_generate_bin_summaries.md

```R
clone500kbwindows<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/04_dmr/bigwig_output_dmr_clones_500bp/03.1.VMR_umap.dmr_clones_500bp.fine_cluster.500bp_windows.rds")

pct_mat<-clone500kbwindows[["pct_matrix"]] 
sum_mat<-clone500kbwindows[["sum_matrix"]] 

```

1. For samples with multiple clones (clones > 50 cells):
- Run one v one and one v rest comparisons within sample

```R
clones_passing_filter

#get samples with 2 or more clones with >=50 cells
samples_with_multiple_clones<-table(dat@metadata$cnv_clonename_500kb,dat@metadata$Sample)
samples_with_multiple_clones<-colnames(samples_with_multiple_clones)[which(colSums(samples_with_multiple_clones>=50)>1)]

# 1. across cell types
comparison_generator_within_sample<-function(i){
    sample_comparisons <-colnames(pct_mat)[!(colnames(pct_mat) %in% c("chr","start","end"))]
    sample_comparisons <- sample_comparisons[grep(sample_comparisons,pattern=i)]

    #1. One vs rest (if there are more than 2)
    if(length(sample_comparisons)>2){
    comparisons_one_v_rest<-do.call("rbind",
      lapply(sample_comparisons,function(j){
        comparisons_one_v_rest <-  data.frame(
          name=paste0(j,"_v_rest"),
          A=j,
          B=paste(sample_comparisons[grep(sample_comparisons,invert=T,pattern=j)],collapse=","))
    }))
    } #otherwise its the same as one v one

    #2. One vs one 
    comparisons_one_v_one <-do.call("rbind",
      lapply(sample_comparisons,function(j){
          do.call("rbind",lapply(sample_comparisons,function(k){
            if(j!=k){
            data.frame(
              name=paste0(j,"_v_",k),
              A=j,
              B=k)
            }else{data.frame(name=NA,A=NA,B=NA)}}))}))

    if(length(sample_comparisons)>2){
    comparisons_sample <- do.call("rbind",list(
                                            comparisons_one_v_rest,
                                            comparisons_one_v_one))
    } else {comparisons_sample<-comparisons_one_v_one}
    return(comparisons_sample)
}

comparisons_within <- do.call("rbind",lapply(samples_with_multiple_clones,comparison_generator_within_sample))
comparisons_within <- comparisons_within[complete.cases(comparisons_within),]
row.names(comparisons_within)<-comparisons_within$name

comparisons<-comparisons_within
output_folder<-paste(project_data_directory,"04_dmr","dmr_within_sample_across_clones",sep="/")
dir.create(output_folder)
suffix="dmr_within_sample_across_clones"
dmr_out<-lapply(row.names(comparisons),comparison_set_dmr)
dmr_out<-do.call("rbind",dmr_out)
```

#2. for each clone (>50 cells) compare to lumhr_HBCA also 
compare all clones (sum > 50 cells) within the same sample merged to lumhr_HBCA
```R

#get samples with 1 or more clones with >=50 cells
samples_with_any_clones<-table(dat@metadata$cnv_clonename_500kb,dat@metadata$Sample)
samples_with_any_clones<-colnames(samples_with_any_clones)[which(colSums(samples_with_any_clones)>=50)]

# 2. across clones to compare clones or merged cancer per sample to lumhr_HBCA
comparison_generator_to_lumhr_hbca<-function(i){
    sample_comparisons <-colnames(pct_mat)[!(colnames(pct_mat) %in% c("chr","start","end"))]
    sample_comparisons <- sample_comparisons[grep(sample_comparisons,pattern=i)]

    #1. All vs lumhr (so sample vs lumhr essentially)
    if(length(sample_comparisons)>2){
      comparisons_all_vs_lumhr <-  data.frame(
          name=paste0(i,"_allclones_v_lumhr_HBCA"),
          A=paste(sample_comparisons,collapse=","),
          B="HBCA_lumhr")
    }  #otherwise its the same as one v lumhr

    #2. One vs one 
    comparisons_one_v_lumhr <-do.call("rbind",
      lapply(sample_comparisons,function(j){
        if(j %in% names(clones_passing_filter)){
          data.frame(
              name=paste0(j,"_v_","lumhr_HBCA"),
              A=j,
              B="HBCA_lumhr")
              }else{data.frame(name=NA,A=NA,B=NA)}
            }))

    if(length(sample_comparisons)>2){
    comparisons_sample <- do.call("rbind",list(
                                            comparisons_all_vs_lumhr,
                                            comparisons_one_v_lumhr))
    } else {comparisons_sample<-comparisons_one_v_lumhr}
    return(comparisons_sample)
}

comparisons_without <- do.call("rbind",lapply(samples_with_any_clones,comparison_generator_to_lumhr_hbca))
comparisons_without <- comparisons_without[complete.cases(comparisons_without),]
row.names(comparisons_without)<-comparisons_without$name

#run dmr
comparisons<-comparisons_without
output_folder<-paste(project_data_directory,"04_dmr","dmr_clones_to_lumhrHBCA",sep="/")
dir.create(output_folder)
suffix="dmr_clones_to_lumhrHBCA"
dmr_out<-lapply(row.names(comparisons),comparison_set_dmr)
dmr_out<-do.call("rbind",dmr_out)

```

4. DMR of clone groups
4.1. WGD vs not, WGD vs lumhr, notwgd vs lumhr

```R
clone500kbwindows<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/04_dmr/bigwig_output_dmr_clones_500bp/03.1.VMR_umap.dmr_clones_500bp.fine_cluster.500bp_windows.rds")

pct_mat<-clone500kbwindows[["pct_matrix"]] 
sum_mat<-clone500kbwindows[["sum_matrix"]] 

#1. WGD vs not
wgd_clones<-obj@metadata %>% group_by(cnv_clonename_500kb) %>% summarize(mean_ploidy=median(cnv_scquantum_ploidy)) %>% as.data.frame()
wgd_clones<-wgd_clones %>% filter(mean_ploidy>2.5) #cut off per van loo and hanghui ye paper
#    cnv_clonename_500kb mean_ploidy
#1        BCMDCIS124T_c2    4.636497
#2         BCMDCIS52T_c1    2.819645
#3  BCMDCIS80T_24hTis_c2    3.498913
#4  BCMDCIS80T_24hTis_c3    3.464776
#5            ECIS25T_c1    3.080883
#6            ECIS25T_c2    3.164176
#7            ECIS25T_c3    3.068273
#8            ECIS25T_c4    2.970785
#9            ECIS57T_c1    3.443908
#10           ECIS57T_c2    3.413622
nonwgd_clones<-unique(dat@metadata$cnv_clonename_500kb)[!(unique(dat@metadata$cnv_clonename_500kb) %in% wgd_clones)]
nonwgd_clones<-nonwgd_clones[!endsWith(nonwgd_clones,suffix="_diploid")]
nonwgd_clones<-nonwgd_clones[nonwgd_clones!="HBCA_lumr"]
wgd_clones<-wgd_clones$cnv_clonename_500kb
comparisons_wgd_vs_nonwgd <- data.frame(
  name="wgd_vs_nonwgd",
  A=paste(wgd_clones,collapse=","),
  B=paste(nonwgd_clones,collapse=",")
)

comparisons_wgd_vs_lumhr <- data.frame(
  name="wgd_vs_lumhr",
  A=paste(wgd_clones,collapse=","),
  B="HBCA_lumhr"
)

comparisons_nonwgd_vs_lumhr <- data.frame(
  name="nonwgd_vs_lumhr",
  A=paste(nonwgd_clones,collapse=","),
  B="HBCA_lumhr"
)

comparisons <- do.call("rbind",list(comparisons_wgd_vs_nonwgd,comparisons_wgd_vs_lumhr,comparisons_nonwgd_vs_lumhr))
comparisons <- comparisons[complete.cases(comparisons),]
row.names(comparisons)<-comparisons$name


output_folder<-paste(project_data_directory,"04_dmr","dmr_wgd_clones",sep="/")
dir.create(output_folder)
suffix="dmr_wgd_clones"
dmr_out<-lapply(row.names(comparisons),comparison_set_dmr)
dmr_out<-do.call("rbind",dmr_out)

```

#3. der(1,16) clones
```R

#diploid cells with clones that lead to the t(1,16) genotype
chr1q_2n_16p_2n_16q_2n=c(
  "BCMDCIS102T_24hTis_diploid",
  "BCMDCIS28T_diploid",
  "BCMDCIS35T_diploid",
  "BCMDCIS41T_diploid",
  "BCMDCIS66T_diploid",
  "BCMDCIS70T_diploid",
  "BCMDCIS74T_diploid",
  "BCMDCIS74T_diploid",
  "BCMDCIS79T_24hTis_DCIS_diploid",
  "BCMDCIS92T_24hTis_diploid",
  "BCMDCIS94T_24hTis_diploid",
  "BCMDCIS97T_diploid",
  "BCMHBCA03R_diploid",
  "BCMHBCA83L−3h_diploid",
  "ECIS26T_diploid",
  "ECIS36T_diploid")
  
chr1q_2n_16p_2n_16q_1n=c("BCMDCIS74T_c2",
                        "BCMDCIS41T_c1",
                        "BCMDCIS28T_c1",
                        "BCMHBCA83L−3h_c1")

chr1q_2n_16p_3n_16q_1n=c("BCMDCIS66T_c2")

chr1q_3n_16p_2n_16q_2n=c("BCMDCIS102T_24hTis_c1",
                        "BCMDCIS102T_24hTis_c5",
                        "BCMDCIS102T_24hTis_c2",
                        "BCMDCIS102T_24hTis_c3",
                        "BCMDCIS92T_24hTis_c2",
                        "BCMHBCA03R_c1")

chr1q_3n_16p_3n_16q_2n=c("BCMDCIS74T_c5")

chr1q_3n_16p_2n_16q_1n=c("BCMDCIS70T_c2",
                        "BCMDCIS65T_c1",
                        "BCMDCIS102T_24hTis_c4",
                        "BCMDCIS70T_c1",
                        "BCMDCIS35T_c1",
                        "ECIS26T_c1",
                        "BCMDCIS74T_c1",
                        "BCMDCIS94T_24hTis_c1",
                        "BCMDCIS79T_24hTis_DCIS_c1",
                        "BCMDCIS94T_24hTis_c2",
                        "ECIS26T_c2ECIS26T_c3",
                        "BCMDCIS97T_c2",
                        "ECIS36T_c2",
                        "ECIS36T_c3")

chr1q_3n_16p_3n_16q_1n=c("BCMDCIS74T_c3","BCMDCIS74T_c4")

chr1q_4n_16p_2n_16q_1n=c("BCMDCIS41T_c4","BCMDCIS41T_c3")

chr1q_4n_16p_3n_16q_1n=c("BCMDCIS41T_c5",
                          "BCMDCIS41T_c6",
                          "BCMDCIS70T_c3",
                          "BCMDCIS28T_c2",
                          "BCMDCIS79T_24hTis_DCIS_c2",
                          "BCMDCIS92T_24hTis_c3",
                          "BCMDCIS97T_c3",
                          "BCMDCIS97T_c4",
                          "BCMDCIS79T_24hTis_DCIS_c5",
                          "BCMDCIS92T_24hTis_c4")

dat@metadata$cnv_comparison_set<-NA
dat@metadata[(dat@metadata$Group=="HBCA") & (dat@metadata$celltype=="lumhr"),]$cnv_comparison_set<-"HBCA_lumhr"
dat@metadata[dat@metadata$cnv_clonename_500kb %in% chr1q_2n_16p_2n_16q_2n,]$cnv_comparison_set<-"chr1q_2n_16p_2n_16q_2n"
dat@metadata[dat@metadata$cnv_clonename_500kb %in% chr1q_2n_16p_2n_16q_1n,]$cnv_comparison_set<-"chr1q_2n_16p_2n_16q_1n"
dat@metadata[dat@metadata$cnv_clonename_500kb %in% chr1q_2n_16p_3n_16q_1n,]$cnv_comparison_set<-"chr1q_2n_16p_3n_16q_1n"
dat@metadata[dat@metadata$cnv_clonename_500kb %in% chr1q_3n_16p_2n_16q_2n,]$cnv_comparison_set<-"chr1q_3n_16p_2n_16q_2n"
dat@metadata[dat@metadata$cnv_clonename_500kb %in% chr1q_3n_16p_3n_16q_2n,]$cnv_comparison_set<-"chr1q_3n_16p_3n_16q_2n"
dat@metadata[dat@metadata$cnv_clonename_500kb %in% chr1q_3n_16p_2n_16q_1n,]$cnv_comparison_set<-"chr1q_3n_16p_2n_16q_1n"
dat@metadata[dat@metadata$cnv_clonename_500kb %in% chr1q_3n_16p_3n_16q_1n,]$cnv_comparison_set<-"chr1q_3n_16p_3n_16q_1n"
dat@metadata[dat@metadata$cnv_clonename_500kb %in% chr1q_4n_16p_2n_16q_1n,]$cnv_comparison_set<-"chr1q_4n_16p_2n_16q_1n"
dat@metadata[dat@metadata$cnv_clonename_500kb %in% chr1q_4n_16p_3n_16q_1n,]$cnv_comparison_set<-"chr1q_4n_16p_3n_16q_1n"
dat<-subsetObject(dat,cells=row.names(dat@metadata[!is.na(dat@metadata$cnv_comparison_set),]))

table(dat@metadata$cnv_comparison_set)


#major steps
comparisons_der1_16 <-do.call("rbind",list(
#t(1;16) maybe?
data.frame(
  name="chr1q_2n_16p_2n_16q_2n_v_hbca_lumhr",
  A=paste(chr1q_2n_16p_2n_16q_2n,collapse=","),
  B="HBCA_lumhr"),

#del(16q)
data.frame(
  name="chr1q_2n_16p_2n_16q_1n_v_chr1q_2n_16p_2n_16q_2n",
  A=paste(chr1q_2n_16p_2n_16q_1n,collapse=","),
  B=paste(chr1q_2n_16p_2n_16q_2n,collapse=",")),

#der(1;16)(q10;p10)
data.frame(
  name="chr1q_3n_16p_2n_16q_1n_v_chr1q_2n_16p_2n_16q_2n",
  A=paste(chr1q_3n_16p_2n_16q_1n,collapse=","),
  B=paste(chr1q_2n_16p_2n_16q_2n,collapse=",")),
  
#der(1;16)(q10;p10)x2
data.frame(
  name="chr1q_4n_16p_3n_16q_1n_v_chr1q_3n_16p_2n_16q_1n",
  A=paste(chr1q_4n_16p_3n_16q_1n,collapse=","),
  B=paste(chr1q_3n_16p_2n_16q_1n,collapse=","))
  ))

row.names(comparisons_der1_16)<-comparisons_der1_16$name

comparisons<-comparisons_der1_16

comparisons <- comparisons[complete.cases(comparisons),]
row.names(comparisons)<-comparisons$name
output_folder<-paste(project_data_directory,"04_dmr","dmr_der1_16_clones",sep="/")
dir.create(output_folder)
suffix="dmr_der1_16_clones"
dmr_out<-lapply(row.names(comparisons),comparison_set_dmr)
dmr_out<-do.call("rbind",dmr_out)

```
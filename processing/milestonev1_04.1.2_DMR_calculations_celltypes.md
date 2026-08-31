
Running differential methylation across cell types 

Running DMR comparisons across different sets.

Celltype comparisons
1. Cell type vs rest (all cells)
2. Cell type vs rest (HBCA only)

By Group comparisons
3. Celltype comparison one v one within group (i.e. DCIS fibroblast vs IDC fibroblast)
4. Celltype comparison one v rest within group 


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
  saveRDS(dmrs,paste0(output_folder,"/","04.1.2.dmr.",suffix,".",i,".dmr.rds"))
  print(paste("Filtering DMRs across clones for:",i))
  dmrs <- filterDMR(dmrs, 
              method = "bonferroni", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
              filter = TRUE, # If TRUE, removes insignificant results
              pThreshold = 0.05 # Maxmimum adjusted p value to allow if filter = TRUE
              ) # Minimum absolute value of the log2FC to allow if filter = TRUE
  dmrs$type <- i
  saveRDS(dmrs,paste0(output_folder,"/","04.1.2.dmr.",suffix,".",i,".dmr.filt.rds"))
  collapsed_dmrs <- collapseDMR(obj, 
                          dmrs, 
                          maxDist = 2000, # Max allowable overlap between DMRs to be considered adjacent
                          minLength = 500, # Min length of collapsed DMR window to include in the output
                          reduce = T, # Reduce results to unique observations (recommended)
                          annotate = T) # Add column with overlapping gene names
  collapsed_dmrs$type <- i
  saveRDS(collapsed_dmrs,paste0(output_folder,"/","04.1.2.dmr.",suffix,".",i,".dmr.collapse.rds"))
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

## Run 500 windows per celltype/group and per clone and lumhr

### Run per celltype
Running at 500bp and 250bp window sizes
Only run on cancer and lumhr

```R

output_directory=paste0(wd,"/","dmr_celltype_group")
system(paste0("mkdir -p ",output_directory))

obj<-readRDS(file=paste(project_data_directory,"03_fine_celltyping","03_scaledcis.final_celltypes.amethyst.rds",sep="/"))
celltype500kbwindows<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/04_dmr/bigwig_output_dmr_celltype_group_500bp/03.1.VMR_umap.dmr_celltype_group_500bp.fine_cluster.500bp_windows.rds")
pct_mat<-celltype500kbwindows[["pct_matrix"]] 
sum_mat<-celltype500kbwindows[["sum_matrix"]] 

comparison_columns<-colnames(pct_mat)
celltypes<-unique(unlist(lapply(strsplit(comparison_columns[4:length(comparison_columns)],"_"),"[",1)))

# 1. across cell types
comparison_generator_across_celltypes<-function(i){
      sample_comparisons <-colnames(pct_mat)[!(colnames(pct_mat) %in% c("chr","start","end"))]

    #1. One vs rest
    comparisons_one_v_rest <-  data.frame(
        name=paste0(i,"_v_rest"),
        A=paste(sample_comparisons[grep(sample_comparisons,pattern=i)],collapse=","),
        B=paste(sample_comparisons[grep(sample_comparisons,invert=T,pattern=i)],collapse=","))
    
    #2. One vs one (at cell type level)
    comparisons_one_v_one <-do.call("rbind",lapply(celltypes,function(k){
          if(i!=k){
            data.frame(name=paste0(i,"_v_",k),
              A=paste(sample_comparisons[grep(sample_comparisons,pattern=i)],collapse=","),
              B=paste(sample_comparisons[grep(sample_comparisons,pattern=k)],collapse=","))
              }else{
              data.frame(name=NA,A=NA,B=NA)}}))

    comparisons_sample <- do.call("rbind",list(
                                            comparisons_one_v_rest,
                                            comparisons_one_v_one))
    return(comparisons_sample)
}

comparisons_without <- do.call("rbind",lapply(celltypes,comparison_generator_across_celltypes))
comparisons_without <- comparisons_without[complete.cases(comparisons_without),]
row.names(comparisons_without)<-comparisons_without$name

#run dmr
comparisons<-comparisons_without
output_folder<-paste(project_data_directory,"04_dmr","dmr_across_celltype","dmr_out",sep="/")
dir.create(output_folder)
suffix="dmr_across_celltype"
dmr_out<-parallel::mclapply(row.names(comparisons),comparison_set_dmr,mc.cores=100)
dmr_out<-do.call("rbind",dmr_out)

# 2. across cell types HBCA only (for normal set)
comparison_generator_within_HBCA<-function(i){
      sample_comparisons <-colnames(pct_mat)[endsWith(colnames(pct_mat),suffix="_HBCA")]

    #1. One vs rest
    comparisons_one_v_rest <-  data.frame(
        name=paste0(i,"_v_rest_HBCAonly"),
        A=paste(sample_comparisons[grep(sample_comparisons,pattern=i)],collapse=","),
        B=paste(sample_comparisons[grep(sample_comparisons,invert=T,pattern=i)],collapse=","))
    
    #2. One vs one (at cell type level)
    comparisons_one_v_one <-do.call("rbind",lapply(celltypes,function(k){
          if(i!=k){
            data.frame(name=paste0(i,"_v_",k,"_HBCAonly"),
              A=paste(sample_comparisons[grep(sample_comparisons,pattern=i)],collapse=","),
              B=paste(sample_comparisons[grep(sample_comparisons,pattern=k)],collapse=","))
              }else{
              data.frame(name=NA,A=NA,B=NA)}}))

    comparisons_sample <- do.call("rbind",list(
                                            comparisons_one_v_rest,
                                            comparisons_one_v_one))
    return(comparisons_sample)
}

comparisons_HBCA_only <- do.call("rbind",lapply(celltypes,comparison_generator_within_HBCA))
comparisons_HBCA_only <- comparisons_HBCA_only[complete.cases(comparisons_HBCA_only),]
row.names(comparisons_HBCA_only)<-comparisons_HBCA_only$name


comparisons<-comparisons_HBCA_only
output_folder<-paste(project_data_directory,"04_dmr","dmr_across_celltype_HBCAonly","dmr_out",sep="/")
dir.create(output_folder)
suffix="dmr_across_celltype_HBCAonly"
dmr_out<-parallel::mclapply(row.names(comparisons),comparison_set_dmr,mc.cores=100)


#3-4. do all one v one and one v rest comparisons per cell type
comparison_generator_within_celltypes<-function(i){
    sample_comparisons <-colnames(pct_mat)[grep(colnames(pct_mat),pattern=i)]

    #1. One vs rest
    comparisons_one_v_rest <-do.call("rbind",lapply(sample_comparisons,function(j){
      data.frame(
        name=paste0(j,"_v_rest"),
        A=j,
        B=paste(sample_comparisons[grep(sample_comparisons,invert=T,pattern=j)],collapse=","))}))
    #2. One vs one
    comparisons_one_v_one <-do.call("rbind",
      lapply(sample_comparisons,function(j){
        do.call("rbind",lapply(sample_comparisons,function(k){
          if(j!=k){data.frame(name=paste0(j,"_v_",k),A=j,B=k)
          }else{data.frame(name=NA,A=NA,B=NA)}}))}))

    comparisons_sample <- do.call("rbind",list(
                                            comparisons_one_v_rest,
                                            comparisons_one_v_one))
    return(comparisons_sample)
}

comparisons_within <- do.call("rbind",lapply(celltypes,comparison_generator_within_celltypes))
comparisons_within <- comparisons_within[complete.cases(comparisons_within),]
row.names(comparisons_within)<-comparisons_within$name

comparisons<-comparisons_within
output_folder<-paste(project_data_directory,"04_dmr","dmr_within_celltype_across_group","dmr_out",sep="/")
dir.create(output_folder)
suffix="within_celltype_across_group"
dmr_out<-parallel::mclapply(row.names(comparisons),comparison_set_dmr,mc.cores=100)
dmr_out<-do.call("rbind",dmr_out)

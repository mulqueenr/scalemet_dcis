
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

## Run 500 or 250bp windows per celltype/group and per clone and lumhr

Then split in different ways to calculate DMR


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
library(RPhenograph)
library(matrixStats)
library(sparseMatrixStats)
library(irlba)
library(rtracklayer)

project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"

#read in object from directory
task_cpus=300
processing_folder="04_collapsed_bam_and_dmr"
wd=paste(sep="/",project_data_directory,processing_folder)
system(paste0("mkdir -p ",wd))
setwd(wd)

```
Functions
```R
#generate DMR input
write_binned_bigwigs<-function(celltype_tracks=celltype_tracks,
                                outdir=celltype_outdir, step=step,
                                i){
    #split bw into 4 files per cell type by methylation/average methylation
    out_dat<-celltype_tracks %>% distinct(chr,start,end, .keep_all=TRUE) %>% select(chr,start,end,i) 

    hg38_seq_info<-Seqinfo(genome="hg38")
    out_dat<-GRanges(out_dat[complete.cases(out_dat),]) #filter NA
    out_dat<-out_dat[out_dat@seqnames %in% hg38_seq_info@seqnames,] #filter chr
    out_dat<-resize(out_dat,width=step)
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
    seqlengths(out_dat_hypermet)<-as.data.frame(hg38_seq_info)[hg38_seq_info@seqnames %in% out_dat_hypermet@seqnames,]$seqlengths #filter by seqlengths

    out_dat_met_mid <- out_dat %>% 
                      as.data.frame() %>% 
                      filter(mcols(out_dat)$score <= mean_score & mcols(out_dat)$score > 50) %>% 
                      mutate(score=score-mean_score) %>% 
                      GRanges() 
    names(out_dat_met_mid@elementMetadata)<-"score"
    genome(out_dat_met_mid)<-"hg38"
    seqlengths(out_dat_met_mid)<-as.data.frame(hg38_seq_info)[hg38_seq_info@seqnames %in% out_dat_met_mid@seqnames,]$seqlengths #filter by seqlengths

    out_dat_met_low <- out_dat %>% 
                        as.data.frame() %>% 
                        filter(mcols(out_dat)$score <= 50 & mcols(out_dat)$score > 20) %>% 
                        mutate(score=score-mean_score) %>% 
                        GRanges() 
    names(out_dat_met_low@elementMetadata)<-"score"
    genome(out_dat_met_low)<-"hg38"
    seqlengths(out_dat_met_low)<-as.data.frame(hg38_seq_info)[hg38_seq_info@seqnames %in% out_dat_met_low@seqnames,]$seqlengths #filter by seqlengths

    out_dat_met_hypomet <- out_dat %>% 
                            as.data.frame() %>% filter(mcols(out_dat)$score <= 20) %>%
                            mutate(score=score-mean_score) %>% 
                            GRanges() 
    names(out_dat_met_hypomet@elementMetadata)<-"score"
    genome(out_dat_met_hypomet)<-"hg38"
    seqlengths(out_dat_met_hypomet)<-as.data.frame(hg38_seq_info)[hg38_seq_info@seqnames %in% out_dat_met_hypomet@seqnames,]$seqlengths #filter by seqlengths

    print(paste("Saving bedgraphs for...",i))
    export.bw(out_dat_hypermet,con=paste0(outdir,"/",i,".hypermet.bw"))
    export.bw(out_dat_met_mid,con=paste0(outdir,"/",i,".midmet.bw"))
    export.bw(out_dat_met_low,con=paste0(outdir,"/",i,".lowmet.bw"))
    export.bw(out_dat_met_hypomet,con=paste0(outdir,"/",i,".hypomet.bw"))
  #make a ucsc/igv trackhub
  col_out=c("hypermet"="#333333","midmet"="#444444","lowmet"="#bcbcbc","hypomet"="#FF13F0")
  
  #write out multiwig trackhub
  multiwig<-c(
      paste(c("track", paste0(i)),collapse=" "),
      paste(c("container", "multiWig"),collapse=" "),
      paste(c("shortLabel",  paste0(i)),collapse=" "),
      paste(c("longLabel",  paste0(i)),collapse=" "),
      paste(c("visibility", "full"),collapse=" "),
      paste(c("aggregate", "transparentOverlay"),collapse=" "),
      paste(c("showSubtrackColorOnUi", "on"),collapse=" "))
  track_lines<-lapply(names(col_out),function(met){
    return(c(
      paste(c("track", paste0(i,".",met,".bw")),collapse=" "),
      paste(c("shortLabel", paste0(i,".",met)),collapse=" "),
      paste(c("longLabel",  paste0(i,".",met)),collapse=" "),
      paste(c("parent",  paste0(i)),collapse=" "),
      paste(c("type", "bigWig"),collapse=" "),
      paste(c("group", i),collapse=" "),
      paste(c("autoscale", "group"),collapse=" "),
      paste(c("graphTypeDefault", "bar"),collapse=" "),
      paste(c("color", paste(col2rgb(col_out[met])[,1],collapse=","),collapse=" "))))
  })
  return(c(list(multiwig),track_lines))
  

}

generate_bigwig<-function(obj=immune,
                      suffix="immune",
                      groupBy="fine_cluster_phenograph",
                      outdir="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/03_fine_celltyping/",
                      step=500){
  bigwig_output_dir=paste0(outdir,"/bigwig_output_",suffix)
  system(paste("mkdir -p ", bigwig_output_dir))
  system(paste("mkdir -p ", paste0(bigwig_output_dir,"/","hg38")))

  obj@h5paths$path<-obj@h5paths$paths
  obj@h5paths$barcode<-row.names(obj@h5paths)

  #update to new smoothed windows
  celltype500bpwindows <- calcSmoothedWindows(obj, 
                                          type = "CG", 
                                          threads = 300,
                                          step = step, 
                                          smooth = 1,
                                          genome = "hg38",
                                          index = "chr_cg",
                                          groupBy = groupBy,
                                          returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                          returnPctMatrix = TRUE)
  saveRDS(celltype500bpwindows,file=paste0(bigwig_output_dir,"/","03.1.VMR_umap.",suffix,".fine_cluster.500bp_windows.rds"))
  obj@genomeMatrices[[paste0("cg_",suffix,"_cells_perc")]] <- celltype500bpwindows[["pct_matrix"]]
  #output tracks as bigwig
  groups<-colnames(celltype500bpwindows[["pct_matrix"]])[4:ncol(celltype500bpwindows[["pct_matrix"]])]
  groups<-groups[groups!="NA"]
  
  #make track hub.txt
  writeLines(c(
    paste(c("hub", groupBy),collapse=" "),
    paste(c("shortLabel", groupBy),collapse=" "),
    paste(c("longLabel", groupBy),collapse=" "),
    paste(c("genomesFile", "genomesFile.txt"),collapse=" "),
    paste(c("email ryan"),collapse=" ")),
  paste0(bigwig_output_dir,"/","hub.txt"))

  #make track genomesFile.txt
  writeLines(c(
    paste(c("genome", "hg38"),collapse=" "),
    paste(c("trackDb", "hg38/trackDb.txt"),collapse=" ")),
  paste0(bigwig_output_dir,"/","genomesFile.txt"))
  
  #write out bigwigs and trackhub data
  trackhub_dat<-lapply(groups,function(i) {
    write_binned_bigwigs(celltype_tracks=celltype500bpwindows[["pct_matrix"]],
                          outdir=paste0(bigwig_output_dir,"/","hg38"),
                          i=i,
                          step=step)})
  file_conn <- file(paste0(bigwig_output_dir,"/hg38/","trackDb.txt"), open = "a")
  for(i in trackhub_dat){writeLines(unlist(i), con = file_conn)}
  close(file_conn)

  return(obj)
}

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

## Run 500 or 250bp windows per celltype/group and per clone and lumhr

### Run per celltype x group
Running at 500bp and 250bp window sizes
```R
output_directory=paste0(wd,"/","dmr_celltype_group")
system(paste0("mkdir -p ",output_directory))

obj<-readRDS(file=paste(project_data_directory,"03_fine_celltyping","03_scaledcis.final_celltypes.amethyst.rds",sep="/"))
dat<-amethyst::subsetObject(obj,cells=row.names(obj@metadata[!(obj@metadata$celltype %in% c("stromal_unknown")) & !is.na(obj@metadata$celltype),]))

dat@metadata$celltype_group<-paste(dat@metadata$celltype,dat@metadata$Group,sep="_")

dat<-generate_bigwig(obj=dat,
                        suffix="dmr_celltype_group_500bp",
                        groupBy="celltype_group",
                        step=500,
                        outdir=getwd())

dat<-generate_bigwig(obj=dat,
                        suffix="dmr_celltype_group_250bp",
                        groupBy="celltype_group",
                        step=250,
                        outdir=getwd())                      
```

## Now running groupwise comparisons across cell types

First running on 500bp windows
```R
input_folder=paste0(wd,"/","bigwig_output_dmr_celltype_group_500bp")
suffix="dmr_celltype_group_500bp"

output_directory=paste0(input_folder,"/","dmr_celltype_group")
system(paste0("mkdir -p ",output_directory))

celltype500bpwindows<-readRDS(file=paste0(input_folder,"/","03.1.VMR_umap.",suffix,".fine_cluster.500bp_windows.rds"))

pct_mat<-celltype500bpwindows[["pct_matrix"]] 
sum_mat<-celltype500bpwindows[["sum_matrix"]] 


#do 1v1 and 1vrest for all celltype comparisons possible
#1vrest
comparisons_1vrest<-do.call("rbind",lapply(unique(dat@metadata$celltype),function(celltype){
  comparisons_set<-colnames(pct_mat)[grepl(colnames(pct_mat),pattern=celltype)]
  comparisons_set<-do.call("rbind",lapply(comparisons_set,function(i){
    data.frame(
    name=paste0(i,"_v_",celltype,"_rest"),
    A=i,
    B=paste(comparisons_set[grep(comparisons_set,invert=T,pattern=i)],collapse=","))
  }))
  }))

#1v1
comparisons_1v1<-do.call("rbind",lapply(unique(dat@metadata$celltype),function(celltype){
  comparisons_set<-colnames(pct_mat)[grepl(colnames(pct_mat),pattern=celltype)]
  comparisons_set<-do.call("rbind",lapply(comparisons_set,function(i){
    j_list=comparisons_set[grep(comparisons_set,invert=T,pattern=i)]
    do.call("rbind",lapply(j_list,function(j){
      data.frame(
      name=paste0(i,"_v_",j),
      A=i,
      B=j
      )}))
}))}))

comparisons<-rbind(comparisons_1vrest,comparisons_1v1)
row.names(comparisons)<-comparisons$name

#calculate DMRs for group-wise comparisons across cell types
#run DMR analysis per row in comparisons
dmr_out<-mclapply(row.names(comparisons), function(i) {
  celltype <- unlist(lapply(strsplit(i,"_"),"[",1))[1]
  print(paste("Running DMRs across cell types for:",i))
  comparisons_output_directory<-paste0(output_directory,"/",celltype,"/",i)
  system(paste0("mkdir -p ",paste0(output_directory,"/",celltype)))
  system(paste0("mkdir -p ",comparisons_output_directory))

  dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
          comparisons=comparisons[i,],
          nminTotal = 3, # Min number observations across all groups to include the region in calculations
          nminGroup = 3) # Min number observations across either members or nonmembers to include the region
  dmrs$comparison_name <- i
  dmrs$celltype <- celltype
  dmrs$groupA <-comparisons[i,]$A 
  dmrs$groupB <- comparisons[i,]$B
  saveRDS(dmrs,paste0(comparisons_output_directory,"/","04.1.",suffix,".",i,".dmr.rds"))

  print(paste("Filtering DMRs across clones for:",i))
  dmrs <- filterDMR(dmrs, 
              method = "BH", #bonferroni is too strict, lose a lot of significant results
              filter = FALSE, # If TRUE, removes insignificant results
              pThreshold = 0.05, # Maxmimum adjusted p value to allow if filter = TRUE
              logThreshold=0) # Minimum absolute value of the log2FC to allow if filter = TRUE
  dmrs<-dmrs %>% filter(padj<0.05)
  dmrs$comparison_name <- i
  dmrs <- collapseDMR(obj, 
                          dmrs, 
                          maxDist = 0, # Max allowable overlap between DMRs to be considered adjacent
                          minLength = 100, # Min length of collapsed DMR window to include in the output
                          reduce = T, # Reduce results to unique observations (recommended)
                          annotate = T) # Add column with overlapping gene names
  dmrs$comparison_name <- i
  dmrs$celltype <- unlist(lapply(strsplit(i,"_"),"[",1))[1]
  dmrs$groupA <-comparisons[i,]$A 
  dmrs$groupB <- comparisons[i,]$B

  saveRDS(dmrs,paste0(comparisons_output_directory,"/","04.1.",suffix,".",i,".dmr.filt.rds"))

  #output DMRs as bed files
  write.table(
    dmrs %>% filter(direction=="hypo") %>% select(chr,dmr_start,dmr_end) %>% as.data.frame(),
    col.names=F,row.names=F, sep="\t", quote=F,
    file=paste0(comparisons_output_directory,"/","04.1.",suffix,".",i,".dmr.filt.hypo.bed"))

  #output DMRs as bed files
  write.table(
    dmrs %>% filter(direction=="hyper") %>% select(chr,dmr_start,dmr_end) %>% as.data.frame(),
    col.names=F,row.names=F, sep="\t", quote=F,
    file=paste0(comparisons_output_directory,"/","04.1.",suffix,".",i,".dmr.filt.hyper.bed"))

  return(dmrs)
  },mc.cores=100)


dmr_out<-do.call("rbind",dmr_out)


```

Performing same comparisons using RNA data
```R

#compare dmr sites with rna marker genes overlap
library(Seurat)
library(GenomicRanges)


rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")
rna$celltype_group<-paste(rna$coarse_celltype,rna$Group,sep="_")
rna<-JoinLayers(rna)
input_folder=paste0(wd,"/","bigwig_output_dmr_celltype_group_500bp")
output_directory=paste0(input_folder,"/","dmr_celltype_group")

de_out<-lapply(row.names(comparisons), function(i) {
  celltype <- unlist(lapply(strsplit(i,"_"),"[",1))[1]
  print(paste("Running DE genes across cell types for:",i))
  comparisons_output_directory<-paste0(output_directory,"/",celltype,"/",i)

  #run paired comparisons that DMR underwent
  cell_A<-row.names(rna@meta.data)[rna@meta.data$celltype_group %in% unlist(strsplit(comparisons[i,]$A,","))]
  cell_B<-row.names(rna@meta.data)[rna@meta.data$celltype_group %in% unlist(strsplit(comparisons[i,]$B,","))]
  if(length(cell_A)>50 & length(cell_B)>50){
  print(table(rna@meta.data[cell_A,]$celltype_group))
  print(table(rna@meta.data[cell_B,]$celltype_group))
  rna_markers<-FindMarkers(rna@assays$RNA,cells.1=cell_A,cells.2=cell_B)
  saveRDS(rna_markers,file=paste0(comparisons_output_directory,"/","04.1.",suffix,".",i,".de.RNAmarkers.rds"))
  return(de_out)
  }
})

#also output bed file of significant RNA for plotting

```

Plotting of DMR count per comparison

```R

dmr_out<-do.call("rbind",lapply(list.files(path=output_directory,recursive=T,full.names=T,pattern="*.dmr.filt.rds"),readRDS))

#color assignment is fluor as cancer associated cell type, rest muted versions
celltype_col=c(
"pericyte"="#FF9900",
"fibroblast"="#FF0000",
"endothelial"="#FFFF66",
"myeloid"="#99FFFF",
"bcell"="#0099CC",
"tcell"="#99FF99",
"basal"="#990099",
"lumsec"="#CC0066",
"lumhr"="#FF00CC",
"cancer"="#00FF99")

mclapply(unique(dmr_out$celltype),function(i){
  print(paste("Plotting...",i))
  celltype_to_plot=i
  comparison_order_for_plotting=c(
    paste(celltype_to_plot,"HBCA","v",celltype_to_plot,"DCIS",sep="_"),
    paste(celltype_to_plot,"HBCA","v",celltype_to_plot,"Synchronous",sep="_"),
    paste(celltype_to_plot,"HBCA","v",celltype_to_plot,"IDC",sep="_"),
    paste(celltype_to_plot,"DCIS","v",celltype_to_plot,"Synchronous",sep="_"),
    paste(celltype_to_plot,"DCIS","v",celltype_to_plot,"IDC",sep="_"),
    paste(celltype_to_plot,"Synchronous","v",celltype_to_plot,"IDC",sep="_"),

    paste(celltype_to_plot,"HBCA","v",celltype_to_plot,"rest",sep="_"),
    paste(celltype_to_plot,"DCIS","v",celltype_to_plot,"rest",sep="_"),
    paste(celltype_to_plot,"Synchronous","v",celltype_to_plot,"rest",sep="_"),
    paste(celltype_to_plot,"IDC","v",celltype_to_plot,"rest",sep="_")
  )

  test_dat<-dmr_out %>% 
            filter(celltype==celltype_to_plot) %>% 
            filter(dmr_padj<0.05) %>% 
            filter(comparison_name %in% comparison_order_for_plotting) %>% 
            group_by(comparison_name,direction,.drop=FALSE) %>% 
            summarize(count=n(),dmr_kbp=sum(dmr_length)/1000) 

  test_dat$comparison_name <- factor(test_dat$comparison_name ,levels=rev(comparison_order_for_plotting))
  #plot count and size of DMRs
  plt1<-ggplot(test_dat %>% filter(direction=="hyper"),aes(x=comparison_name,y=count,fill=I(unname(celltype_col[celltype_to_plot]))))+
        geom_bar(stat="identity") + 
        geom_text(aes(label=paste0(dmr_kbp,"kbp")), vjust=0,color="yellow") +
        theme_minimal() + 
        coord_flip() 

  plt2<-ggplot(test_dat %>% filter(direction=="hypo"),aes(x=comparison_name,y=count,color=I(unname(celltype_col[celltype_to_plot]))))+
        geom_bar(stat="identity",fill="white") + 
        geom_text(aes(label=paste0(dmr_kbp,"kbp")), vjust=0.1,color="red") +
        theme_minimal() + 
        coord_flip()
  system(paste0("mkdir -p ",paste0(output_directory,"/","dmr_plots/")))
  ggsave(plt1|plt2,file=paste0(output_directory,"/","dmr_plots/","04.1.",suffix,".",celltype_to_plot,".dmr.counts.pdf"),height=20,width=20,units="in")

  #run GSEA on cell type
  system(paste0("mkdir -p ",output_directory,"/gsea_data"))
  for(j in comparison_order_for_plotting){
    gsea_across_sets(obj, 
                    collapsed_dmrs=dmr_out %>% filter(comparison_name == j),
                    sample_name=j, 
                    prefix=paste0(output_directory,"/gsea_data/","04.1.",suffix,".",celltype_to_plot))
  }

  tft_gsea<-do.call("rbind",lapply(
    list.files(path=paste0(output_directory,"/gsea_data/"),pattern=paste0("04.1.",suffix,".",celltype_to_plot,".*","TFT.rds"),full.names=T),
    readRDS))
  plot_gsea(gsea=tft_gsea,out_setname="TFT",prefix=paste0(output_directory,"/dmr_plots/","04.1.",suffix,".",celltype_to_plot))

  hallmark_gsea<-do.call("rbind",lapply(
    list.files(path=paste0(output_directory,"/gsea_data/"),pattern=paste0("04.1.",suffix,".",celltype_to_plot,".*","hallmark.rds"),full.names=T),
    readRDS))
  plot_gsea(gsea=hallmark_gsea,out_setname="hallmark",prefix=paste0(output_directory,"/dmr_plots/","04.1.",suffix,".",celltype_to_plot))

  position_gsea<-do.call("rbind",lapply(
    list.files(path=paste0(output_directory,"/gsea_data/"),pattern=paste0("04.1.",suffix,".",celltype_to_plot,".*","position.rds"),full.names=T),
    readRDS))
  plot_gsea(gsea=position_gsea,out_setname="position",prefix=paste0(output_directory,"/dmr_plots/","04.1.",suffix,".",celltype_to_plot))

  canceratlas_gsea<-do.call("rbind",lapply(
    list.files(path=paste0(output_directory,"/gsea_data/"),pattern=paste0("04.1.",suffix,".",celltype_to_plot,".*","3CA.rds"),full.names=T),
    readRDS))
  plot_gsea(gsea=canceratlas_gsea,out_setname="3CA",prefix=paste0(output_directory,"/dmr_plots/","04.1.",suffix,".",celltype_to_plot))
  },mc.cores=50)

#calculate DMR proximity to DE genes (or overlap with DE genes)
#DMR overlap with breast cancer gene list
#GSEA enrichment
#overlap with bulk methylation progression markers

#Run chromvar on hypomethylated 500bp windows
#with and without excluding CGI
```


```R
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
set.seed(2017)
library(JASPAR2024)
library(TFBSTools)
register(MulticoreParam(50, progressbar = TRUE))
library(BSgenome.Hsapiens.UCSC.hg38)
task_cpus=300

input_folder=paste0(wd,"/","bigwig_output_dmr_celltype_group_500bp")
suffix="dmr_celltype_group_500bp"
output_directory=paste0(input_folder,"/","dmr_celltype_group")
system(paste0("mkdir -p ",output_directory))

celltype500bpwindows<-readRDS(file=paste0(input_folder,"/","03.1.VMR_umap.",suffix,".fine_cluster.500bp_windows.rds"))

pct_mat<-celltype500bpwindows[["pct_matrix"]] 
pct_mat<-pct_mat[!duplicated(paste(pct_mat$chr,pct_mat$start,pct_mat$end,sep="_")),]
#set up binary counts matrix per window
binary_met<-pct_mat[,4:ncol(pct_mat)]
binary_met<-(binary_met>20)+0 #hypomethylated regions (<20% met) are 1, rest 0 or NA

windows_passing_filter<-which(rowSums(is.na(binary_met),na.rm=T)>1) #require at least 1 groups to have a hypomethylated region
binary_met<-binary_met[windows_passing_filter,]
ranges<-GRanges(pct_mat[windows_passing_filter,1:3])

binary_met[which(is.na(binary_met),arr.ind=T)]<-0 #set remaining NA to 0

#set up hypomet matrix
binary_met <- SummarizedExperiment(assays = list(counts = binary_met),rowRanges = sort(ranges))
ranges<-sort(ranges)

binary_met <- addGCBias(binary_met, genome = BSgenome.Hsapiens.UCSC.hg38)


#compute deviations for peaks
binary_met <- filterPeaks(binary_met,non_overlapping=FALSE)
bg <- getBackgroundPeaks(object = binary_met)


#set up motifs location calls across ranges
motifs <- getMatrixSet(x = JASPAR2024()@db, opts = list(collection = "CORE", species = "Homo sapiens"))

ranges_list<-split(binary_met@rowRanges, ntile(seq_along(binary_met@rowRanges),task_cpus))
motif_ix <- mclapply(1:length(ranges_list),function(i){
                        matchMotifs(motifs, ranges_list[[i]], genome = BSgenome.Hsapiens.UCSC.hg38)
                        },mc.cores=task_cpus)
motif_ix<-do.call("rbind",motif_ix)

#run chromvar
dev <- computeDeviations(object = binary_met, annotations = motif_ix, background_peaks = bg)
variability <- computeVariability(dev)

pdf("test.chromvar.variability.pdf")
plotVariability(variability, use_plotly = FALSE) 
dev.off()
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
  test_name<-setNames(nm=row.names(comp),1:nrow(comp))
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

 plt<-patchwork::wrap_plots(tft_plt,ncol=4,axes="collect")

  ggsave(tft_plt,
          file=paste0(prefix,".","GSEA_enrichment",sample_name,".dotplot.pdf"),
          width=40,height=length(unique(collapsed_dmrs$type))*5,limitsize = FALSE)

  plt<-patchwork::wrap_plots(list((tft_plt),(position_plt),(hallmark_plt),(cancercellatlas_plt)),ncol=4,axes="collect")

  ggsave(plt,
          file=paste0(prefix,".","GSEA_enrichment",sample_name,".dotplot.pdf"),
          width=40,height=length(unique(collapsed_dmrs$type))*5,limitsize = FALSE)

chromvar_on_dmr_sites<-function()

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
dat<-readRDS(file="08_scaledcis.final_celltype.amethyst.rds") #just saving because of the macrophage name correction

lumhr_hbca_cells <- row.names(dat@metadata[dat@metadata$Group=="HBCA" & dat@metadata$celltype=="lumhr",])

clones_passing_filter<-names(which(table(dat@metadata$cnv_clonename_500kb)>30))
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
dat<-readRDS(file="08_scaledcis.final_celltype.amethyst.rds") #just saving because of the macrophage name correction

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


```









4-6. Group comparisons

```R
#set up directories
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
dmr_output_dir<-paste0(project_data_directory,"/DMR_analysis/","/DMR_celltype_HBCAonly")
system(paste0("mkdir -p ",dmr_output_dir))
prefix=paste0(dmr_output_dir,"/","08_scaledcis.final_celltype.onevrest.HBCAonly_celltype")

#read in files
dat<-readRDS(file="08_scaledcis.final_celltype.amethyst.rds") #just saving because of the macrophage name correction
celltype500bpwindows<-readRDS(file=paste0(project_data_directory,"/merged_data/","08_scaledcis.final_celltype_by_group.500bp_windows.rds"))

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






<!--
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

-->

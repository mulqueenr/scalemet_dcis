Running subclustering on just immune cells for finer cell typing

1. using RNA based marker genes for hypomethylation over genes
2. overlapping dmrs and rna marker genes (fishers exact test)
3. ranking RNA marker genes and DMRs for better markers
4. integration with rliger

Output percent methylation as bedgraph to view in IGV over marker genes as well

Could rerun DMR per coarse cell type to see if that helps as well. 

# Prepare RNA marker genes per cell type

From processing/milestonev1_00.2_seurat_scrna_copykat.md seurat object.
<!--
RUN DURING IMMUNE SECTION
processing/milestonev1_01.2_amethyst_immune_fine_celltyping.md
```R

#compare dmr sites with rna marker genes overlap
library(Seurat)
library(GenomicRanges)

rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")

#just run on all cell types
table(rna$fine_celltype)
rna$celltype<-NA
rna$celltype<-rna$fine_celltype
Idents(rna)<-rna$celltype
table(Idents(rna))
rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)
rna <- ScaleData(rna)
rna <- JoinLayers(rna)
rna_markers<-FindAllMarkers(rna,assay="RNA",only.pos=TRUE)
saveRDS(rna_markers,file="/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.rna_markers.rds")
```
-->


# Read in methylation data and additional libraries

```R
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing

#source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
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
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)
system(paste0("mkdir -p ",project_data_directory,"/fine_celltyping"))
dat<-readRDS(file="05_scaledcis.coarse_clusters.amethyst.rds")
```

# stromal Cells

Wrapping each step into a function, so I can run it on stromal cell types as well.

## Subcluster on just stromal cells

Set up variables for output.
```R
prefix="stromal"
output_directory=paste0(project_data_directory,"/fine_celltyping/",prefix)
system(paste0("mkdir -p ",output_directory))
```

```R
celltype_umap<-function(obj=dat,prefix="allcells",dims=12,regressCov=TRUE,k_pheno=50,neigh=25,dist=1e-5,method="cosine",output_directory,window_name){
  print("Running IRLBA reduction...")
  obj@reductions[[paste(window_name,"irlba",sep="_")]] <- runIrlba(obj, genomeMatrices = c(window_name), dims = dims, replaceNA = c(0))

  if(regressCov){
    print("Running regression...")
  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- regressCovBias(obj, reduction = paste(window_name,"irlba",sep="_")) 
  } else {
      print("Skipping regression...")
  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- obj@reductions[[paste(window_name,"irlba",sep="_")]]
  }

  obj <- amethyst::runCluster(obj, k_phenograph = k_pheno, reduction = paste(window_name,"irlba_regressed",sep="_")) 

  print(paste("Running UMAP...",as.character(neigh),as.character(dist),as.character(method)))
  obj <- amethyst::runUmap(obj, neighbors = neigh, dist = dist, method = method, reduction = paste(window_name,"irlba_regressed",sep="_")) 

  outname=paste(prefix,"integrated_celltype",dims,as.character(regressCov),k_pheno,neigh,as.character(dist),method,sep="_")
  print(paste("Plotting...",outname))

  p1 <- dimFeature(obj, colorBy = sample, reduction = "umap") + ggtitle(paste(window_name,"Samples"))
  p2 <- dimFeature(obj, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red"),guide="colourbar") + ggtitle("Coverage distribution")
  p3 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
  p4 <- dimFeature(obj, colorBy = cluster_id, reduction = "umap") + ggtitle(paste(window_name,"Clusters"))
  p5 <- dimFeature(obj, colorBy = Group, reduction = "umap") + ggtitle(paste(window_name," Group"))
  p6<-ggplot()
  ggsave((p1|p2)/(p3|p4)/(p5|p6),file=paste0(output_directory,"/",outname,"_umap.pdf"),width=20,height=30)  
  return(obj)
}


cluster_subset<-function(
  dat=dat,
  broad_celltype=c("immune"), #note this is a list
  output_directory,
  window_name="coarse_cluster_dmr_sites",
  #umap function args
  prefix="immune",
  dims=12,
  regressCov=TRUE,
  k_pheno=50,
  neigh=25,
  dist=1e-5,
  method="cosine"){
  system(paste0("mkdir -p ",output_directory))
  print(paste("Subsetting object based on broad celltype:",broad_celltype))

  #subset amethyst object to broad celltype
  dat_sub<-subsetObject(dat,cells=row.names(dat@metadata[dat@metadata$broad_celltype %in% broad_celltype,]))

  #subset to windows with coverage
  dat_sub@genomeMatrices[[window_name]] <- dat_sub@genomeMatrices[[window_name]][rowSums(!is.na(dat_sub@genomeMatrices[[window_name]])) >= 45, ]
  est_dim<-dimEstimate(dat_sub, genomeMatrices = c(window_name), dims = c(20), threshold = 0.95)
  print(paste(est_dim,"estimated dimensions."))

  dat_sub<-celltype_umap(obj=dat_sub,
                        prefix=paste0("06_",prefix,"finecelltyping"),
                        dims=dims,
                        regressCov=regressCov,
                        k_pheno=k_pheno,
                        neigh=neigh,
                        dist=0.01,
                        method="cosine",
                        window_name=window_name,
                        output_directory=output_directory) 

  dat_sub@metadata$fine_cluster_id<-dat_sub@metadata$cluster_id
  return(dat_sub)
}


dat<-cluster_subset(
  dat=dat,
  broad_celltype=c("endothelial"), #note this is a list
  window_name="coarse_cluster_dmr_sites",
  prefix="endothelial",
  dims=14,
  regressCov=FALSE,
  k_pheno=200,
  neigh=10,
  dist=0.01,
  method="cosine",
  output_directory=output_directory)


#dat<-celltype_umap(obj=dat,prefix="06_stromal_finecelltyping",dims=16,regressCov=FALSE,k_pheno=15,neigh=6,dist=0.001,method="cosine") 

saveRDS(dat,file=paste0(output_directory,"/","06_scaledcis.",prefix,"_finecelltyping.amethyst.rds"))
dat<-readRDS(file=paste0(output_directory,"/","06_scaledcis.",prefix,"_finecelltyping.amethyst.rds"))

```


## Summarize new stromal clusters for marker plotting and calculate DMRs per final subclusters

```R


#bigwig output
bigwig_output<-function(
    obj,
    tracks="cg_cluster_tracks"){
    for(i in 4:ncol(obj@genomeMatrices[[tracks]])){
        cluster=names(obj@genomeMatrices[[tracks]])[i]
        out_bw<-as.data.frame(obj@genomeMatrices[[tracks]])
        out_bw<-out_bw[c("chr","start","end",cluster)]
        out_bw<-GRanges(out_bw[complete.cases(out_bw),])
        names(out_bw@elementMetadata)<-"score"
        out_bw<-out_bw[unique(findOverlaps(out_bw, type = "any", select = "first"))]
        out_bw <- resize(out_bw, width=1000, fix='start') #resize to avoid 1base overlap
        genome(out_bw)<-"hg38"
        hg38_seq_info<-Seqinfo(genome="hg38")
        seqlengths(out_bw)<-as.data.frame(hg38_seq_info)[hg38_seq_info@seqnames %in% out_bw@seqnames,]$seqlengths
        print(paste("Saving bigwig for...",cluster))
        export(out_bw,con=paste0(output_dir,"/",paste(tracks,cluster,"bw",sep=".")),format='bigWig')}
}
#updated test DMR function (correction by lauren in later releases of amethyst)
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

calculate_dmrs<-function(dat=dat,prefix=prefix,output_directory=output_directory){
  celltype500bpwindows <- calcSmoothedWindows(dat, 
                                          type = "CG", 
                                          threads = 100,
                                          step = 500, 
                                          smooth = 3,
                                          genome = "hg38",
                                          index = "chr_cg",
                                          groupBy = "fine_cluster_id",
                                          returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                          returnPctMatrix = TRUE)

  saveRDS(celltype500bpwindows,file=paste0(output_directory,"/","06_",prefix,"_finecelltyping.500bp_windows.rds"))
  
  dat@genomeMatrices[[paste0("cg_",prefix,"cells_perc")]] <- celltype500bpwindows[["pct_matrix"]]

  #output tracks as bigwig
  bigwig_output(obj=dat,tracks=paste0("cg_",prefix,"cells_perc"),output_directory=output_directory)

  saveRDS(dat,file=paste0(output_directory,"/","06_scaledcis.",prefix,"_finecelltyping.amethyst.rds"))


  pct_mat<-celltype500bpwindows[["pct_matrix"]] 
  sum_mat<-celltype500bpwindows[["sum_matrix"]] 

  dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
          eachVsAll = TRUE, # If TRUE, each group found in the sumMatrix will be tested against all others
          nminTotal = 3, # Min number observations across all groups to include the region in calculations
          nminGroup = 3) # Min number observations across either members or nonmembers to include the region
  saveRDS(dmrs,file=paste0(output_directory,"/","06_",prefix,"_finecelltyping",".500bp_dmrs.rds"))

  dmrs <- filterDMR(dmrs, 
              method = "bonferroni", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
              filter = FALSE, # If TRUE, removes insignificant results
              pThreshold = 0.05, # Maxmimum adjusted p value to allow if filter = TRUE
              logThreshold = 1.5) # Minimum absolute value of the log2FC to allow if filter = TRUE


  collapsed_dmrs <- collapseDMR(dat, 
                        dmrs, 
                          maxDist = 2000, # Max allowable overlap between DMRs to be considered adjacent
                          minLength = 500, # Min length of collapsed DMR window to include in the output
                          reduce = T, # Reduce results to unique observations (recommended)
                          annotate = T) # Add column with overlapping gene names


  saveRDS(collapsed_dmrs,file=paste0(output_directory,"/","06_",prefix,"_finecelltyping.500bp_dmrs.filt_collapsed.rds"))

  rename_dmr_output<-setNames(nm=1:length(unique(collapsed_dmrs$test)),gsub(colnames(sum_mat)[grepl(colnames(sum_mat),pattern="_t$")],pattern="_t",replacement=""))
  collapsed_dmrs$type <- rename_dmr_output[collapsed_dmrs$test]
  saveRDS(collapsed_dmrs,file=paste0(output_directory,"/","06_",prefix,"_finecelltyping.500bp_dmrs.filt_collapsed.rds"))

  #plot dmr counts per cluster
  pal <- colorRampPalette(colors = c("#F9AB60", "#E7576E", "#630661", "#B5DCA5"))
  COLS <- pal(length(unique(collapsed_dmrs$type)))

  plt<-ggplot(collapsed_dmrs |> dplyr::group_by(type, direction) |> dplyr::summarise(n = n()), 
      aes(y = type, x = n, fill = type)) + 
      geom_col() + 
      facet_grid(vars(direction), scales = "free_y") + 
      scale_fill_manual(values = COLS) + 
      theme_classic()
  ggsave(plt,file=paste0(output_directory,"/","06_",prefix,"_finecelltyping.500bp_dmrs.filt_collapsed.barplot.pdf"))
  return(collapsed_dmrs)

}

calculate_dmrs(dat=dat,prefix=prefix,output_directory=output_directory)
```

## DMR site overlap with RNA marker genes

Use RNA markers and DMR markers to define cell types.

```R
rna_markers<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.rna_markers.rds")
collapsed_dmrs<-readRDS(file=paste0(output_directory,"/","06_",prefix,"_finecelltyping.500bp_dmrs.filt_collapsed.rds"))

#perform fisher test on set
row_fisher <- function(i){
    df<-celltype_dmr_overlaps[i,]
    p11<-df["dmr_rnamarker_overlap"]
    p10<-df["dmr_cluster_total"]-df["dmr_rnamarker_overlap"]
    p01<-df["rna_celltype_total"]
    p00<-length(Features(rna))
    mat <- matrix(as.numeric(c(p11,p10,p01,p00)), ncol=2)
    f <- fisher.test(as.table(mat), alt="greater")
    return(f$p.value)
}

dmr_rna_marker_overlap<-function(rna_markers,collapsed_dmrs,prefix,output_directory){

  #overexpressed genes by cell type
  rna_markers_filt<-rna_markers %>% 
                      filter(avg_log2FC > 1) %>% 
                      filter(p_val_adj <0.05) %>% 
                      group_by(cluster)


  #filter collapsed dmrs by padj and extract gene names per type
  dmr_out<-collapsed_dmrs %>% 
                  filter(direction=="hypo") %>% 
                  filter(dmr_padj<0.05) %>% 
                  filter(abs(dmr_logFC) > 1.5) %>% 
                  filter(dmr_length<50000) %>% 
                  filter(!is.na(gene_names)) %>% 
                  as.data.frame()

  dmr_out$met_feature<-paste(dmr_out$chr,dmr_out$dmr_start,dmr_out$dmr_end,sep="_")

  #overlap with rna markers
  #expand dmr_out so each row is its own gene name
  #this will be used to summarize multiple DMRs to their genes
  dmrs_to_gene_names<-lapply(1:nrow(dmr_out), function(x) {
    temp<-dmr_out[x,]
    genes<-unlist(temp$gene_names %>% stringr::str_replace_all(" ","") %>% strsplit(split=","))
    return(cbind(genes=genes,type=temp$type,met_feature=temp$met_feature))
  })

  dmr_genes<-as.data.frame(do.call("rbind",dmrs_to_gene_names))
  colnames(dmr_genes)<-c("gene","cluster","met_feature")

  dmr_genes<-dmr_genes %>% filter(!duplicated(gene,cluster)) %>% group_by(cluster) 

  #count intersect between type in DMRs and cell type in RNA markers
  rna_markers_filt
  dmr_genes

  celltype_dmr_overlaps <- dmr_genes %>%
    inner_join(rna_markers_filt, by = "gene", relationship="many-to-many") %>%
    group_by(cluster.x,) %>%
    count(cluster.y, name = "overlap_count") %>%
    as.data.frame()

  colnames(celltype_dmr_overlaps)<-c("cluster","celltype","dmr_rnamarker_overlap")

  #add count of dmr features passing filter per cluster
  celltype_dmr_overlaps$dmr_cluster_total<-NA
  dmr_count<-dmr_genes %>% filter(!duplicated(met_feature)) %>% count() 
  dmr_count<-setNames(nm=dmr_count$cluster,dmr_count$n)
  celltype_dmr_overlaps$dmr_cluster_total<-NA
  celltype_dmr_overlaps$dmr_cluster_total<-dmr_count[celltype_dmr_overlaps$cluster]

  #add count of rna markers passing filter per cell type
  rna_marker_count<-rna_markers_filt %>% count() 
  rna_marker_count<-setNames(nm=rna_marker_count$cluster,rna_marker_count$n)
  celltype_dmr_overlaps$rna_celltype_total<-0
  celltype_dmr_overlaps$rna_celltype_total<-rna_marker_count[celltype_dmr_overlaps$celltype]

  fishers <- unlist(lapply(1:nrow(celltype_dmr_overlaps),row_fisher))
  celltype_dmr_overlaps$fishers.pval<-NA
  celltype_dmr_overlaps$fishers.pval<-fishers
  celltype_dmr_overlaps<-as.data.frame(celltype_dmr_overlaps)

  dat <- reshape2::dcast(celltype_dmr_overlaps, cluster~celltype,value.var="fishers.pval",fill=1) 
  dat_counts <- reshape2::dcast(celltype_dmr_overlaps, cluster~celltype,value.var="dmr_rnamarker_overlap",fill=0) 

  row.names(dat)<-dat$cluster
  row.names(dat_counts)<-dat_counts$cluster

  dmr_counts<-unique(celltype_dmr_overlaps[c("cluster","dmr_cluster_total")])
  dmr_counts<-setNames(nm=dmr_counts$cluster,dmr_counts$dmr_cluster_total)

  de_counts<-unique(celltype_dmr_overlaps[c("celltype","rna_celltype_total")])
  de_counts<-setNames(nm=de_counts$celltype,de_counts$rna_celltype_total)

  #negative log10 on dat
  dat<-dat[2:ncol(dat)]
  dat_counts<-dat_counts[2:ncol(dat_counts)]

  dat<- -(log10(dat))
  col_fun=colorRamp2(breaks=quantile(unlist(dat),probs=c(0,0.5,0.9)),colors=c("white","grey","#FF00FF"))

  #add annotation barplots on counts for rna and dmr
  row_ha = rowAnnotation(dmr_count = anno_barplot(dmr_counts[row.names(dat)]))
  column_ha = columnAnnotation(de_count = anno_barplot(de_counts[colnames(dat)]))

  #add dot size by overlap
  cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, 
              gp = gpar(col = "grey", fill = NA))
              grid.text(sprintf("%.1f", dat[i, j]), x = x, y = y)
      }
      
  plt<-Heatmap(dat,
              cell_fun=cell_fun,
              top_annotation = column_ha, 
              right_annotation = row_ha,
              col=col_fun,
              name="-log10p")

  #col=col_fun,)
  pdf(paste0(output_directory,"/","06_",prefix,"_finecelltyping.dmr_rnamarkergene.fisher.pdf"),width=20,height=20)
  print(plt)
  dev.off()


  #and compare group proportion per cluster
  dat_met<-readRDS(file=paste0(output_directory,"/","06_scaledcis.immune_finecelltyping.amethyst.rds"))
  prop_celltype<-table(dat_met@metadata$fine_cluster_id,dat_met@metadata$Group) %>% reshape2::melt()
  colnames(prop_celltype)<-c("cluster","group","count")
  plt<-ggplot(prop_celltype,aes(x=as.character(cluster),y=count,fill=group))+geom_bar(position="fill",stat="identity")+theme_minimal()
  ggsave(plt,file=paste0(output_directory,"/","06_",prefix,"finecelltyping.dmr_rnamarkergene.cellprop.pdf"))

}

rna_markers<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.rna_markers.rds")
collapsed_dmrs<-readRDS(file=paste0(output_directory,"/","06_",prefix,"_finecelltyping.500bp_dmrs.filt_collapsed.rds"))
dmr_rna_marker_overlap(rna_markers=rna_markers,collapsed_dmrs=collapsed_dmrs,prefix=prefix,output_directory=output_directory)
```



# Plotting marker genes over clusters

```R
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing

#source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
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
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)

dat<-readRDS(file=paste0(output_directory,"/","06_scaledcis.",prefix,"_finecelltyping.amethyst.rds"))

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

#https://www.nature.com/articles/s41588-024-01688-9/figures/2
cell_markers[["immune"]]<-c("PTPRC","CXCR4") #good!
cell_markers[["tcell"]]<-c("CD3D","CD3G","CD3E","IL7R") #good! 4,2,15,12,1,3,11,13,14

cell_markers[["tcell_CD4"]]<-c("SELL","TCF7","CD4","LEF1","CCR4","ITGB1","CD45RO") #not good
cell_markers[["tcell_CD8"]]<-c("CD8A","CD8B","GZMA","GZMB") #not good
cell_markers[["tcell_treg"]]<-c("FOXP3","CXCL13","IL2RA","IKZF2") #cd8b good


cell_markers[["myeloid"]]<-c("HLA-DRA","HLA-DPA1","CD74") #good! 6,9,5,17
cell_markers[["bcell"]]<-c("MS4A1","CD79A","CD19","CD79B","CD22") #good! 17,5

#want to find: mono, dc,  macro, neutrophil, mast, TAM (if possible)

cell_markers[["macro"]]<-c("CD68","FCER1G","TREM2","LYVE1","FOLR2","NRP1") #FCER1G good, TREM2 FOLR2 is TAM 14,8,6,9,16,7,17
cell_markers[["mono"]]<-c("CD14","CSF1R","FCGR3A","CX3CR1","ITGAM") #FCER1G good, TREM2 FOLR2 is TAM
cell_markers[["neutrophil"]]<-c("S100A8", "S100A9", "CXCL8","CSF3R","FCGR3B","CD177","MMP9")
cell_markers[["TAMs"]]<-c("CD163","CD206","MARCO","TREM2")
immune_markers<-cell_markers

cell_colors=c(
"bcell"="#65cbe4",
"plasma"="#7ecdc2",
"immune"="#8088c2",
"tcell"="#1d87c8",
"tcell_CD8"="#8088c2",
"tcell_CD4"="#8088c2",
"macro"="#8088c2"
)

#prepare cgi
cgisland="/data/rmulqueen/projects/scalebio_dcis/ref/cpgIslandExt.bed"
cgi<-rtracklayer::import(cgisland)
cgi<-as.data.frame(cgi)
colnames(cgi)<-c("chr","start","end","strand")


#rough ordering by what clustered together
order<-c(
  "4","2", "15",
  "12",
  "1","3","11",
  "13","14","8",
  "6",
  "9",
  "10",
  "16",
  "7",
  "17","5"
  )
   

mclapply(names(immune_markers),
function(celltype){
    genes<-unlist(immune_markers[celltype])
    genes<-genes[genes %in% dat@ref$gene_name]
    print(paste("Plotting:",celltype))
    print(paste("Genes to plot:",genes))
    plt<-histograModified(obj=dat, 
        baseline="mean",
        genes = unlist(immune_markers[celltype]),
        #colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
        matrix = paste0("cg_",prefix,"cells_perc"), arrowScale = .03, trackScale = .5,
        legend = F, cgisland=cgi,order=order) + ggtitle(celltype)
    ggsave(plt,
          file=paste0(output_directory,"/","06_",prefix,"finecelltyping.",celltype,".marker.pdf"),
          width=3*length(genes),
          height=ncol(dat@genomeMatrices[[paste0("cg_",prefix,"cells_perc")]])*1,
          limitsize=F)

},mc.cores=20)


```

Plot top marker sites based on DMR per group

```R

collapsed_dmrs<-readRDS(file=paste0(output_directory,"/","06_",prefix,"_finecelltyping.500bp_dmrs.filt_collapsed.rds"))

collapsed_dmrs %>% 
                filter(direction=="hypo") %>% 
                filter(dmr_padj<0.05) %>% 
                filter(abs(dmr_logFC) > 1.5) %>% 
                filter(dmr_length<50000) %>% 
                filter(gene_names!="NA") %>% 
                filter(!grepl(gene_names,pattern=",")) %>% 
                filter(type=="5") %>% select(gene_names) %>% slice_head(n=20)


collapsed_dmrs %>% filter(type=="9") %>% filter(direction=="hypo") %>% filter(dmr_padj<0.05)%>% slice_head(n=20) %>% select(gene_names)
#tcells https://cellxgene.cziscience.com/e/74520626-b0ba-4ee9-86b5-714649554def.cxg/

#myeloid https://cellxgene.cziscience.com/e/a6388a6f-6076-401b-9b30-7d4306a20035.cxg/

```


```R


order<-c(
"11",
"6",
"16",
"5",

"9"="tcell_nk", #maybe, def a tcell??
"1"="tcell_cd4","10"="tcell_cd4",
"2"="tcell_cd8","4"="tcell_cd8",
"6","8","3","7",
"5"="bcell_plasma"#cluster finer)
   

```

















# Use Liger to integrate DMRs and RNA markers for immune cells

```R
library(rliger)

dat<-readRDS(file="06_scaledcis.immune_finecelltyping.amethyst.rds")
collapsed_dmrs<-readRDS(file=paste0("06_immune_finecelltyping",".500bp_dmrs.filt_collapsed.rds"))

#select dmrs for integration (same way we selected for dmr and gene overlap)
dmrs_for_integration <- collapsed_dmrs %>% 
                          filter(direction=="hypo") %>%
                          filter(dmr_padj < 0.05) %>%
                          filter(abs(dmr_logFC) > 1.5) %>% 
                          filter(dmr_length<50000) %>% 
                          filter(gene_names != "NA") %>%
                          group_by(type) 
table(dmrs_for_integration$type)

genes_for_integration <- dmrs_for_integration %>% 
                          as.data.frame() %>%
                          select(gene_names) %>% 
                          unlist %>%
                          strsplit(split=",") %>% 
                          unlist() %>% 
                          stringr::str_replace_all(" ","") %>% 
                          unique() %>%
                          sort()

#reduce overlapping dmrs regardless of cluster
dmr_out<-dmrs_for_integration %>% 
        as.data.frame() %>% 
        select(chr,dmr_start,dmr_end) %>% 
        distinct(chr,dmr_start,dmr_end) %>% 
        makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>%
        reduce()

dmr_out<-dmr_out[width(dmr_out)<50000,]
summary(width(dmr_out))

#assign genes to merged dmr_out features
reduced_overlap<-findOverlaps(dmr_out,
            makeGRangesFromDataFrame(dmrs_for_integration,keep.extra.columns=TRUE))

#add list of gene names and store original dmr sites in bed
dmr_out$gene_names<-NA
dmr_out$gene_names<-mclapply(unique(queryHits(reduced_overlap)), function(x){
    genes=dmrs_for_integration[subjectHits(reduced_overlap[queryHits(reduced_overlap)==x,]),]$gene_names
    genes <- genes %>% strsplit(split=",") %>% unlist() %>% stringr::str_replace_all(" ","") %>% unique() %>% paste(collapse=",")
    genes<-paste(unlist(genes),sep=",")
    dmr_out[x,]$gene_names<-genes
},mc.cores=100)
dmr_out$met_feature<-paste(seqnames(dmr_out),start(dmr_out),end(dmr_out),sep="_")

#output bed format for making windows
dmr_bed<-data.frame(chr=seqnames(dmr_out),start=start(dmr_out),end=end(dmr_out))

#create new matrix from DMR sites for refined clustering
dat@genomeMatrices[["immune_cluster_dmr_sites"]] <- makeWindows(dat, 
                                                     bed = dmr_bed,
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 50, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 
saveRDS(dat,file="06_scaledcis.immune_finecelltyping.amethyst.rds")

#add filter to DMRs with atleast 10% coverage per cell
met_features<-which(rowSums(!is.na(dat@genomeMatrices[["immune_cluster_dmr_sites"]]))>(ncol(dat@genomeMatrices[["immune_cluster_dmr_sites"]])*0.1))

#expand dmr_out so each row is its own gene name
#this will be used to summarize multiple DMRs to their genes
dmrs_to_gene_names<-lapply(1:length(dmr_out), function(x) {
  temp<-dmr_out[x,]
  genes<-unlist(temp$gene_names %>% unlist() %>% strsplit(split=",")) 
  return(cbind(genes=genes,met_feature=temp$met_feature))
})

dmrs_to_gene_key<-as.data.frame(do.call("rbind",dmrs_to_gene_names))

#set new matrix with average DMR score per gene using the key above
gene_methylation<-mclapply(unique(dmrs_to_gene_key$genes), function(i){
  gene_name=i
  met_features<-dmrs_to_gene_key[dmrs_to_gene_key$genes==gene_name,2]
  return(colMeans(dat@genomeMatrices[["immune_cluster_dmr_sites"]][row.names(dat@genomeMatrices[["immune_cluster_dmr_sites"]]) %in% met_features,],na.rm=TRUE))
},mc.cores=100)

gene_met<-as.data.frame(do.call("rbind",gene_methylation))
row.names(gene_met)<-unique(dmrs_to_gene_key$genes)

# Create from raw score matrices MET
met_dat <- Matrix(as.matrix(gene_met), sparse = TRUE) # Convert the dense matrix to a sparse dgCMatrix

#set NA values to 0 (following amethyst example)
met_dat[which(is.na(met_dat),arr.ind=TRUE)]<-0

system(paste0("mkdir -p ",project_data_directory,"/integration_immune/"))
saveRDS(met_dat,file=paste0(project_data_directory,"/integration_immune/","scaledcis.gene_dmr_methylation.immune.mat.rds"))

table(rna_markers_filt$cluster)

# RNA markers continued from previous section
#get raw RNA values on filtered methylation genes
rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")
rna<-subset(rna,coarse_celltype %in% c("bcell","myeloid","plasma","tcell"))
Idents(rna)<-rna$fine_celltype
rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)
rna <- ScaleData(rna)
rna <- JoinLayers(rna)

#overexpressed genes by fine cell type
rna_markers<-FindAllMarkers(rna,assay="RNA",only.pos=TRUE)
rna_markers_filt<-rna_markers %>% 
                    filter(avg_log2FC > 3) %>% 
                    filter(p_val_adj <0.01) %>% 
                    group_by(cluster)


#filter met_dat to only genes in rna markers
met_dat<-met_dat[row.names(met_dat) %in% unique(rna_markers_filt$gene),]


rna_dat<-LayerData(rna,features=row.names(met_dat),assay="RNA",layer="counts")
saveRDS(rna_dat,file=paste0(project_data_directory,"/integration_immune/","scaledcis.gene_rna.mat.immune.rds"))

rna$broad_celltype<-rna$fine_celltype
rna$fine_cluster_id <- "NA"
rna$mcg_pct <- "NA"

#create shared metadata, add coarse_cluster_id to meta
met_dat<-met_dat[,order(colnames(met_dat))]
meta<-rbind(dat@metadata[colnames(met_dat),c("sample","broad_celltype","fine_cluster_id","mcg_pct")],rna@meta.data[c("sample","broad_celltype","fine_cluster_id","mcg_pct")])

summary(unlist(as.data.frame(met_dat)))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-1.00000  0.00000  0.00000  0.06306  0.00000  1.00000 

summary(unlist(as.data.frame(rna_dat)))
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#    0.000     0.000     0.000     0.374     0.000 15099.000 

met_dat<-met_dat[row.names(met_dat)%in% row.names(rna_dat),]

#make liger object
met.liger <- createLiger(rawData=list(met=met_dat,rna=rna_dat), 
                        modal=c("meth","rna"),cellMeta=meta,
                        removeMissing=TRUE)

#now follow liger for integration
rna.met <- met.liger %>% 
            rliger::normalize() %>%
            rliger::selectGenes(useDatasets = "rna",nGenes=nrow(met_dat)) %>%
            rliger::scaleNotCenter()



#try https://welch-lab.github.io/liger/reference/runUINMF.html on unshared features as well

summary(unlist(as.data.frame(rna.met@datasets$met@scaleData)))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  1.0000  1.0000  0.9002  1.0000  2.0000 

summary(unlist(as.data.frame(rna.met@datasets$rna@scaleData)))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.000   0.000   0.000   0.116   0.000 258.780 

rna.met <- rliger::runIntegration(rna.met, k = 20)
rna.met <- rliger::quantileNorm(rna.met)
rna.met <- runCluster(rna.met)
rna.met <- runUMAP(rna.met)

plt1<-plotDatasetDimRed(rna.met,splitBy="dataset")
ggsave(wrap_plots(plt1),file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","dataset",".integrated.umap.pdf"),width=20,height=10)

plt2<-plotClusterDimRed(rna.met,splitBy="dataset","broad_celltype")
ggsave(wrap_plots(plt2),file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","celltype",".integrated.umap.pdf"),width=20,height=10)

plt3<-plotClusterDimRed(rna.met,splitBy="dataset","sample")
ggsave(wrap_plots(plt3),file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","sample",".integrated.umap.pdf"),width=20,height=10)

plt4<-plotClusterDimRed(rna.met,splitBy="dataset","fine_cluster_id")
ggsave(wrap_plots(plt4),file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","fine_cluster_id",".integrated.umap.pdf"),width=20,height=10)

#plt5<-plotClusterDimRed(rna.met,slot="cellMeta","mcg_pct")
#ggsave(wrap_plots(plt5),file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","mcg",".integrated.umap.pdf"),width=20,height=10)

umap_dims<-rna.met@dimReds$UMAP
colnames(umap_dims)<-c("X","Y")
umap_dims<-umap_dims[paste0("met_",colnames(met_dat)),]
plt<-ggplot(umap_dims,aes(x=X,y=Y,color=rna.met@cellMeta[row.names(umap_dims),]$fine_cluster_id))+geom_point()+theme_void()
ggsave(plt,file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","fine_cluster_id",".integrated.umap.pdf"),width=20,height=10)

plots <- plotGeneDimRed(rna.met, c("PDGFB", "IL3RA"), splitBy = "dataset",
                        titles = c(names(rna.met), names(rna.met)))
plt<-cowplot::plot_grid(plotlist = plots, nrow = 2)
ggsave(plt,file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","fine_cluster_id",".integrated.umap.pdf"),width=20,height=10)

saveRDS(rna.met,file=paste0(project_data_directory,"/integration_immune/","immune.met_rna.integrated.liger.rds"))

rna.met@cellMeta[rna.met@cellMeta$dataset=="met",]$fine_cluster_id == dat@metadata[colnames(met_dat),c("fine_cluster_id")]
gsub("met_","",row.names(rna.met@cellMeta[rna.met@cellMeta$dataset=="met",])) == row.names(dat@metadata[colnames(met_dat),])
colnames(met_dat)==gsub("met_","",row.names(rna.met@cellMeta[rna.met@cellMeta$dataset=="met",])) 
```





Promoter regions on marker genes
Also doesnt work particularly well.

```R


rna_dat<-LayerData(rna,features=unique(rna_markers_filt$gene),assay="RNA",layer="counts")
saveRDS(rna_dat,file=paste0(project_data_directory,"/integration_immune/","scaledcis.gene_rna.mat.immune.rds"))

rna$broad_celltype<-rna$fine_celltype
rna$fine_cluster_id <- "NA"
rna$mcg_pct <- "NA"

#create new matrix from DMR sites for refined clustering
promoter_met <- makeWindows(dat, 
                    genes = unique(rna_markers_filt$gene),
                    promoter=TRUE,
                    type = "CG", 
                    metric = "score", 
                    threads = 50, 
                    index = "chr_cg", 
                    nmin = 2) 

#get the perc matrix before progressing

promoter_met2<-promoter_met[which(rowSums(!is.na(promoter_met))>(ncol(promoter_met)*0.1)),]

rna_dat<-rna_dat[row.names(promoter_met2),]

promoter_met <- Matrix(as.matrix(promoter_met2), sparse = TRUE) # Convert the dense matrix to a sparse dgCMatrix
promoter_met[which(is.na(promoter_met),arr.ind=TRUE)]<-0

#make liger object
met.liger <- createLiger(rawData=list(met=promoter_met,rna=rna_dat), 
                        modal=c("meth","rna"),cellMeta=meta,
                        removeMissing=TRUE)

#now follow liger for integration
rna.met <- met.liger %>% 
            rliger::normalize() %>%
            rliger::selectGenes(useDatasets = "rna") %>%
            rliger::scaleNotCenter()

rna.met <- rliger::runIntegration(rna.met, k = 20)
rna.met <- rliger::quantileNorm(rna.met)
rna.met <- runCluster(rna.met)
rna.met <- runUMAP(rna.met)

plt1<-plotDatasetDimRed(rna.met,splitBy="dataset")
ggsave(wrap_plots(plt1),file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","dataset",".integrated.umap.pdf"),width=20,height=10)

plt2<-plotClusterDimRed(rna.met,splitBy="dataset","broad_celltype")
ggsave(wrap_plots(plt2),file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","celltype",".integrated.umap.pdf"),width=20,height=10)

plt3<-plotClusterDimRed(rna.met,splitBy="dataset","sample")
ggsave(wrap_plots(plt3),file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","sample",".integrated.umap.pdf"),width=20,height=10)

plt4<-plotClusterDimRed(rna.met,splitBy="dataset","fine_cluster_id")
ggsave(wrap_plots(plt4),file=paste0(project_data_directory,"/integration_immune/","scaledcis.met_rna.","fine_cluster_id",".integrated.umap.pdf"),width=20,height=10)


```
Running subclustering on just immune cells for finer cell typing

1. subset out major cell types
2. perform reclustering on dmrs for just those cells
3. generate new subclusters and summarize DMRs per subcluster
4. overlap hypomethylated DMRs with RNA marker genes
5. generate methylation plots over canonical marker genes

# Prepare RNA marker genes per cell type

From processing/milestonev1_00.2_seurat_scrna_copykat.md seurat object.
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

#treg vs cd4
Idents(rna)<-rna$fine_celltype
rna_markers<-FindMarkers(rna,assay="RNA",ident.1="tcell_nk",ident.2=NULL,only.pos=TRUE)
rna_markers %>% mutate(gene=row.names(rna_markers)) %>% filter(p_val_adj<0.05) %>% filter(avg_log2FC>2) %>% head(n=20) %>% select(gene)

Idents(rna)<-rna$coarse_celltype
rna_markers<-FindMarkers(rna,assay="RNA",ident.1="bcell",ident.2="plasma",only.pos=TRUE)
rna_markers %>% mutate(gene=row.names(rna_markers)) %>% filter(p_val_adj<0.05) %>% filter(avg_log2FC>2) %>% head(n=20) %>% select(gene)

rna_markers<-FindAllMarkers(rna,assay="RNA",only.pos=TRUE)
saveRDS(rna_markers,file="/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.rna_markers.rds")

#b plasma vs mem
#then tackle myeloid
rna_markers<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.rna_markers.rds")
rna_markers %>% filter(p_val_adj<0.05) %>% filter(avg_log2FC>1) %>% head(n=20) %>% select(gene)

rna_markers %>% filter(p_val_adj<0.05) %>% filter(avg_log2FC>1) %>% filter(cluster=="tcell_cd4") %>% head(n=20) %>% select(gene)

```


# Read in methylation data and additional libraries
From processing/milestonev1_01_amethyst_coarse_celltyping.md
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

task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)
system(paste0("mkdir -p ",project_data_directory,"/fine_celltyping"))


cell_colors=c(
"basal"="#844c9d",
"lumhr"="#e23e96",
"cancer"="#39FF14",
"lumsec"="#ff6498",
"fibro"="#f58e90",
"endo"="#5bbb5a",
"myeloid"="#8088c2",
"tcell"="#1d87c8",
"bcell"="#65cbe4")

```

## Subcluster on just immune cells
Wrapping each step into a function, so I can run it on stromal cell types as well.

```R

celltype_umap<-function(obj=dat,prefix="allcells",dims=12,pc_use,regressCov=TRUE,regressCG=FALSE,k_pheno=50,k_umap=50,neigh=25,dist=1e-5,method="cosine",output_directory,window_name,cluster_on_umap=FALSE){
  print("Running IRLBA reduction...")
  obj@reductions[[paste(window_name,"irlba",sep="_")]] <- runIrlba(obj, genomeMatrices = c(window_name), dims = dims, replaceNA = c(0))

  if(regressCov){
    print("Running regression...")
  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- regressCovBias(obj, reduction = paste(window_name,"irlba",sep="_")) 
  } else {
      print("Skipping regression...")
  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- obj@reductions[[paste(window_name,"irlba",sep="_")]]
  }

  if(regressCov & regressCG){
      print("Running regression on CG percent met...")
      obj@metadata$cov<-obj@metadata$mcg_pct  #add cg percent to meta for regression, then revert back
      obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- regressCovBias(obj, reduction = paste(window_name,"irlba_regressed",sep="_"),method="glm") 
      obj@metadata <- meta
  } else if(regressCG){
      print("Running regression on CG percent met...")
      obj@metadata$cov<-obj@metadata$mcg_pct  #add cg percent to meta for regression, then revert back
      obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- regressCovBias(obj, reduction = paste(window_name,"irlba",sep="_"),method="glm") 
      obj@metadata <- meta
  } else {
      print("Skipping regression on CG...")
  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- obj@reductions[[paste(window_name,"irlba",sep="_")]]
  }

  #filter irlba output to just pcs to use
  print(paste("Using IRLBA PC:", min(pc_use), "to", max(pc_use)))
  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]]<-obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]][,pc_use]

  obj <- amethyst::runCluster(obj, k_phenograph = k_pheno, reduction = paste(window_name,"irlba_regressed",sep="_")) 

  print(paste("Running UMAP...",as.character(neigh),as.character(dist),as.character(method)))
  obj <- amethyst::runUmap(obj, neighbors = neigh, dist = dist, method = method, reduction = paste(window_name,"irlba_regressed",sep="_")) 

  print("Clustering on umap.")
  if(cluster_on_umap){
  umap_clus<-obj@metadata%>% select(c("umap_x","umap_y"))
  umap_clus<-Rphenograph::Rphenograph(umap_clus,k=k_pheno)
  obj@metadata$cluster_id<-paste0(prefix,"_",igraph::membership(umap_clus[[2]]))
  }

  outname=paste(prefix,"integrated_celltype",dims,as.character(regressCov),k_pheno,neigh,as.character(dist),method,sep="_")
  print(paste("Plotting...",outname))

  p1 <- dimFeature(obj, colorBy = sample, reduction = "umap",pointSize=3) + ggtitle(paste(window_name,"Samples"))
  p2 <- dimFeature(obj, colorBy = log10(cov), pointSize = 3) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red"),guide="colourbar") + ggtitle("Coverage distribution")
  p3 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 3) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
  p4 <- dimFeature(obj, colorBy = cluster_id, reduction = "umap",pointSize=3) + ggtitle(paste(window_name,"Clusters"))
  p5 <- dimFeature(obj, colorBy = Group, reduction = "umap",pointSize=3) + ggtitle(paste(window_name," Group"))
  p6 <- dimFeature(obj, colorBy = broad_celltype, reduction = "umap",pointSize=3) + ggtitle(paste(window_name,"broad_celltype"))
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
  regressCG=FALSE,
  k_pheno=50,
  pc_use=NA,
  k_umap=50,
  neigh=25,
  var_features=NA,
  perc_cell_cov_per_window=0.01, #required window coverage for clustering (1%)
  dist=1e-5,
  method="cosine",cluster_on_umap=FALSE){

  print(paste("Making fine celltyping directory:",output_directory))
  system(paste0("mkdir -p ",output_directory))
  print(paste("Subsetting object based on broad celltype:",broad_celltype))

  #subset amethyst object to broad celltype
  dat_sub<-subsetObject(dat,cells=row.names(dat@metadata[dat@metadata$broad_celltype %in% broad_celltype,]))

  if(any(is.na(pc_use))){
    pc_use=1:max(dims)
  }
  #subset to windows with 5% coverage
  req_cov<-as.integer(nrow(dat_sub@metadata)*perc_cell_cov_per_window)
  print(paste("Requiring",perc_cell_cov_per_window,"of cells for window coverage filter:",as.character(req_cov),"cells."))
  dat_sub@genomeMatrices[[window_name]] <- dat_sub@genomeMatrices[[window_name]][rowSums(!is.na(dat_sub@genomeMatrices[[window_name]])) >= req_cov, ]


  if(!is.na(var_features)){
  print(paste("Slicing to top variance of features:",var_features))
  variable_windows <- rowVars(as.matrix(dat_sub@genomeMatrices[[window_name]]),na.rm=T)
  top_var<-variable_windows %>% order(decreasing=TRUE) %>% head(n=var_features)
  dat_sub@genomeMatrices[[window_name]] <- dat_sub@genomeMatrices[[window_name]][top_var,]
  }

  print(paste("Matrix windows passing filter:",nrow(dat_sub@genomeMatrices[[window_name]])))

  #estimate dims
  est_dim<-dimEstimate(dat_sub, genomeMatrices = c(window_name), dims = c(20), threshold = 0.95)
  print(paste(est_dim,"estimated dimensions."))

  dat_sub<-celltype_umap(obj=dat_sub,
                        prefix=paste0("06_",prefix,".finecelltyping"),
                        dims=dims,
                        regressCov=regressCov,
                        k_pheno=k_pheno,
                        neigh=neigh,
                        dist=dist,
                        method=method,
                        pc_use=pc_use,
                        window_name=window_name,
                        output_directory=output_directory,
                        cluster_on_umap=cluster_on_umap) 

  dat_sub@metadata$fine_cluster_id<-paste(prefix,dat_sub@metadata$cluster_id,sep="_")
  return(dat_sub)
}

bigwig_output<-function(
                        obj,
                        tracks="cg_cluster_tracks",
                        prefix,
                        output_directory){
            hg38_seq_info<-Seqinfo(genome="hg38")

    #https://genome.ucsc.edu/FAQ/FAQformat.html#format1
    for(i in 4:ncol(obj@genomeMatrices[[tracks]])){
        #set genomic ranges
        cluster=names(obj@genomeMatrices[[tracks]])[i]
        out_bw<-as.data.frame(obj@genomeMatrices[[tracks]])
        out_bw<-out_bw[c("chr","start","end",cluster)]
        colnames(out_bw)<-c("chrom","chromStart","chromEnd","score")
        out_bw<-GRanges(out_bw[complete.cases(out_bw),])
        genome(out_bw)<-"hg38"
        seqlengths(out_bw)<-as.data.frame(hg38_seq_info)[hg38_seq_info@seqnames %in% out_bw@seqnames,]$seqlengths
        out_bw<- resize(out_bw, width=499, fix='start') #resize to avoid 1base overlap
        out_bw<-trim(out_bw)

        #track_line <- new("GraphTrackLine",
        #                  autoScale=FALSE,
        #                  color=col2rgb(color)[,1],
        #                  graphType="bar",viewLimits=c(0,100),
        #                  name = cluster)

        print(paste("Saving bigwig for...",cluster))
        rtracklayer::export(object=out_bw,
                          #trackLine=track_line,
                          con=paste0(output_directory,"/","06_",prefix,".",paste(tracks,cluster,"bw",sep=".")))
      }
}

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

calculate_dmrs<-function(dat=dat,prefix=prefix,output_directory=output_directory,split_by_group=TRUE,min_cells=10){

  groupBy="fine_cluster_id"
  if(split_by_group){
  dat@metadata$fine_cluster_id_group<-paste(dat@metadata$fine_cluster_id,dat@metadata$Group,sep="_")
  groupBy="fine_cluster_id_group"
  }

  #filter out groups that are too small
  groups_passing_filter <- names(table(dat@metadata[groupBy]))[table(dat@metadata[groupBy])>min_cells]

  dat<-subsetObject(dat,
                  cells=row.names(dat@metadata)[dat@metadata[groupBy] %in% groups_passing_filter])

  celltype500bpwindows <- calcSmoothedWindows(dat, 
                                          type = "CG", 
                                          threads = 200,
                                          step = 500, 
                                          smooth = 3,
                                          genome = "hg38",
                                          index = "chr_cg",
                                          groupBy = groupBy,
                                          returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                          returnPctMatrix = TRUE)
  
  saveRDS(celltype500bpwindows,file=paste0(output_directory,"/","06_",prefix,"_finecelltyping.500bp_windows.rds"))
  dat@genomeMatrices[[paste0("cg_",prefix,"_cells_perc")]] <- celltype500bpwindows[["pct_matrix"]]

  #output tracks as bigBed
  bigwig_output(obj=dat,
                tracks=paste0("cg_",prefix,"_cells_perc"),
                output_directory=output_directory,
                prefix=prefix)
  #save subset amethyst file
  print(paste("Saving new amethyst file:",paste("06_scaledcis.",prefix,"_finecelltyping.amethyst.rds")))
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

dmr_rna_marker_overlap<-function(dat,rna_markers,collapsed_dmrs,prefix,output_directory){

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

  #make nested lists for overlap
  dmr_list<-list()
  for(x in unique(dmr_genes$cluster)){
      genes<-dmr_genes[dmr_genes$cluster==x,]$gene
      dmr_list[[x]]<-genes[!duplicated(genes)]
    }
  
  rna_markers_list<-list()
  for(x in unique(rna_markers_filt$cluster)){
      genes<-rna_markers_filt[rna_markers_filt$cluster==x,]$gene
      rna_markers_list[[x]]<-genes[!duplicated(genes)]
    }
  #count intersect between type in DMRs and cell type in RNA markers

  gene_overlap<-newGOM(dmr_list,
                  rna_markers_list,
                  genome.size=length(unique(dat@ref$gene_name)))
  
  saveRDS(gene_overlap,paste0(output_directory,"/","06_",prefix,"_finecelltyping.dmr_rnamarkergene.overlap.rds"))

  gene_overlap_pval <- getMatrix(gene_overlap, name="pval")
  gene_overlap_count <- getMatrix(gene_overlap, name="intersection")

  dmr_counts<-unlist(lapply(dmr_list,length))
  de_counts<-unlist(lapply(rna_markers_list,length))

  gene_overlap_pval<- -(log10(gene_overlap_pval))
  col_fun=colorRamp2(breaks=quantile(unlist(gene_overlap_pval),probs=c(0,0.5,0.9)),colors=c("white","grey","#FF00FF"))

  #add annotation barplots on counts for rna and dmr
  row_ha = rowAnnotation(dmr_count = anno_barplot(dmr_counts[row.names(gene_overlap_pval)]))
  column_ha = columnAnnotation(de_count = anno_barplot(de_counts[colnames(gene_overlap_pval)]))
      
  #add dot size by overlap
  cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, 
              gp = gpar(col = "grey", fill = NA))
              grid.text(gene_overlap_count[i,j], x = x, y = y)
      }

  plt<-Heatmap(gene_overlap_pval,
              cell_fun=cell_fun,
              top_annotation = column_ha, 
              right_annotation = row_ha,
              col=col_fun,
              name="-log10p")

  #col=col_fun,)
  pdf(paste0(output_directory,"/","06_",prefix,"_finecelltyping.dmr_rnamarkergene.fisher.pdf"),width=30,height=20)
  print(plt)
  dev.off()

  #and compare group proportion per cluster
  prop_celltype<-table(dat@metadata$fine_cluster_id,dat@metadata$Group) %>% reshape2::melt()
  colnames(prop_celltype)<-c("cluster","group","count")
  plt<-ggplot(prop_celltype,aes(x=as.character(cluster),y=count,fill=group))+geom_bar(position="stack",stat="identity")+theme_minimal()+theme(axis.text.y = element_text(angle = 90, hjust = 1))

  ggsave(plt,file=paste0(output_directory,"/","06_",prefix,".finecelltyping.dmr_rnamarkergene.cellprop.pdf"))

}

find_cluster_markers<-function(dat,celltype500bp_windows,comp){
  #comparisons: If eachVsAll is not desired, provide a data frame
   #       describing which tests to run. The data.frame should have
   #       three columns with rows describing conditions of each test.
   #       "name" determines the name of the test in the output; "A"
   #       lists group members, and "B" lists group nonmembers.

  pct_mat<-celltype500bpwindows[["pct_matrix"]] 
  sum_mat<-celltype500bpwindows[["sum_matrix"]] 

  dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
          comparisons = comp, # If TRUE, each group found in the sumMatrix will be tested against all others
          nminTotal = 3, # Min number observations across all groups to include the region in calculations
          nminGroup = 3) # Min number observations across either members or nonmembers to include the region

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

  rename_dmr_output<-setNames(nm=1:length(unique(collapsed_dmrs$test)),gsub(colnames(sum_mat)[grepl(colnames(sum_mat),pattern="_t$")],pattern="_t",replacement=""))
  collapsed_dmrs$type <- rename_dmr_output[collapsed_dmrs$test]
  return(collapsed_dmrs)
}

```

#Recluster at windows with differing resolutions

```R

set.seed(111)
options(future.globals.maxSize= 500000*1024^2) #80gb limit for parallelizing
library(amethyst)
library(data.table)

project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)

#Set up variables for output.
dat<-readRDS(file="07_scaledcis.cnv_clones.amethyst.rds")

#filter out contigs
x<-50 #5, 10running, 25running, 50running, 100done

chr_length<-read.table("/data/rmulqueen/projects/scalebio_dcis/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/star/chrNameLength.txt")
colnames(chr_length)<-c("chr","end")
chr_length$start<-1
chr_length<-chr_length[chr_length$chr %in%  paste0("chr",c(1:22,"X")),]
chr_length$chr<-factor(chr_length$chr,levels=paste0("chr",c(1:22,"X")))
chr_length<-chr_length[order(chr_length$chr),]
chr_length<-GRanges(
                seqnames=paste0("chr",c(1:22,"X")),
                ranges=IRanges(chr_length$start,chr_length$end),
                seqlengths=setNames(nm=chr_length$chr,chr_length$end))

win<-unlist(tile(GRanges(chr_length),width=x*1000))
bed<-data.frame(chr=seqnames(win),start=start(win),end=end(win))

print(paste("Generating win size",x,"kb"))
win_out <- makeWindows(dat, 
                    bed = bed,
                    type = "CG", 
                    metric = "score", 
                    threads = 100, 
                    index = "chr_cg", 
                    nmin = 2) 
saveRDS(win_out,file=paste0("07_scaledcis.",x,"kb.windows.rds"))

#at 100kb, windows have ~50% cells with coverage
#at 50kb, windows have ~20% cells with coverage
#at 25kb, windows have ~6% cells with coverage

summary(rowSums(!is.na(win_out))/dim(win_out)[1])
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0000171 0.1723718 0.1914405 0.1882263 0.2131647 0.4055990 

x=50 #was 25, for 25kb, 0.3% cell cov is a good cutoff, for 10, probably 15-20% is good?
win_in<-readRDS(file=paste0("07_scaledcis.",x,"kb.windows.rds"))
prefix=paste0("all_cells","_",as.character(x),"kb")
output_directory=paste0(project_data_directory,"/fine_celltyping/",prefix)
system(paste("mkdir -p",output_directory))
dat@genomeMatrices[["cg_win_score"]]<-win_in
saveRDS(dat,file="07_scaledcis.cnv_clones.amethyst.rds")


#### All cells, final setup for all cell umap
#clustering on umap, to assign broad cell types, regressing cov, using 12 dims, of 50kb windows, filtering to 20% coverage  taking top 10000 featuers (by var)
dat_sub<-subsetObject(dat,cells=row.names(dat@metadata)[dat@metadata$mcg_pct>65])
dat_sub<-cluster_subset(
  dat=dat_sub,
  broad_celltype=unique(dat_sub@metadata$broad_celltype), #note this is a list, doing all cells together here
  window_name="cg_win_score", 
  prefix=prefix,
  dims=12, 
  regressCov=TRUE, 
  regressCG=FALSE,
  var_features=10000,
  k_pheno=200,
  k_umap=20,
  neigh=20, 
  dist=0.001,
  perc_cell_cov_per_window=0.2,
  method="cosine",
  output_directory=output_directory,
  cluster_on_umap=TRUE)

scale(table(dat_sub@metadata$broad_celltype,dat_sub@metadata$fine_cluster_id))

#reassign broad celltypes based on broad celltyping clustering
dat_sub@metadata$broad_celltype<-"lumhr"
dat_sub@metadata[dat_sub@metadata$fine_cluster_id %in% c("all_cells_50kb_06_all_cells_50kb.finecelltyping_7","all_cells_50kb_06_all_cells_50kb.finecelltyping_9"),]$broad_celltype<-"basal"
dat_sub@metadata[dat_sub@metadata$fine_cluster_id %in% c("all_cells_50kb_06_all_cells_50kb.finecelltyping_2","all_cells_50kb_06_all_cells_50kb.finecelltyping_5"),]$broad_celltype<-"lumsec"

dat_sub@metadata[dat_sub@metadata$fine_cluster_id %in% 
c("all_cells_50kb_06_all_cells_50kb.finecelltyping_1",
"all_cells_50kb_06_all_cells_50kb.finecelltyping_12",
"all_cells_50kb_06_all_cells_50kb.finecelltyping_10",
"all_cells_50kb_06_all_cells_50kb.finecelltyping_11"),]$broad_celltype<-"stromal"

dat_sub@metadata[dat_sub@metadata$fine_cluster_id %in% 
c("all_cells_50kb_06_all_cells_50kb.finecelltyping_15",
"all_cells_50kb_06_all_cells_50kb.finecelltyping_19",
"all_cells_50kb_06_all_cells_50kb.finecelltyping_21",
"all_cells_50kb_06_all_cells_50kb.finecelltyping_3",
"all_cells_50kb_06_all_cells_50kb.finecelltyping_20"),]$broad_celltype<-"immune"
dat_sub@metadata[which(!endsWith(dat@metadata$cnv_clonename!="diploid") & !is.na(dat@metadata$cnv_clonename)),]$broad_celltype<-"cancer"
#persisting NA values from CNV calls just have too low read count


saveRDS(dat@genomeMatrices[["cg_win_score"]],file="07_scaledcis.filtered_50kbwin.rds")


#save filtered windows for integration with RNA

#save the updated broad celltypes
dat_sub@genomeMatrices[["cg_win_score"]]<-win_in

saveRDS(dat_sub,file="07_scaledcis.cnv_clones.amethyst.rds")
dat<-dat_sub
```

Immune cells

```R
### Immune
dat_sub<-subsetObject(dat,cells=row.names(dat@metadata)[dat@metadata$mcg_pct>65])
prefix=paste0("immune","_",as.character(x),"kb")
output_directory=paste0(project_data_directory,"/fine_celltyping/",prefix)
dat_sub<-cluster_subset(
  dat=dat_sub,
  broad_celltype=c("immune"), #note this is a list, doing all cells together here
  window_name="cg_win_score", #coarse_cluster_dmr_sites orinitial_cluster_5kb_win for 5kb require more coverage
  prefix=prefix,
  dims=12, #i think use at least 14 dims
  regressCov=TRUE, 
  regressCG=FALSE,
  k_pheno=45,
  k_umap=50,
  neigh=80, 
  dist=0.001,
  perc_cell_cov_per_window=0.3,
  method="cosine",
  output_directory=output_directory,
  cluster_on_umap=FALSE)
saveRDS(dat_sub,file=paste0(output_directory,"/","06_scaledcis.",prefix,"_finecelltyping.amethyst.rds"))
dat_sub<-readRDS(file=paste0(output_directory,"/","06_scaledcis.",prefix,"_finecelltyping.amethyst.rds"))


#calculate dmrs per fine celltype grouping
collapsed_dmrs<-calculate_dmrs(dat=dat_sub,split_by_group=TRUE,
                prefix=prefix,
                output_directory=output_directory, min_cells=20)

```

Stromal

```R
dat_sub<-subsetObject(dat,cells=row.names(dat@metadata)[dat@metadata$mcg_pct>65])
prefix=paste0("stromal","_",as.character(x),"kb")
output_directory=paste0(project_data_directory,"/fine_celltyping/",prefix)
dat_sub<-cluster_subset(
  dat=dat_sub,
  broad_celltype=c("stromal"), #note this is a list, doing all cells together here
  window_name="cg_win_score", #coarse_cluster_dmr_sites orinitial_cluster_5kb_win for 5kb require more coverage
  prefix=prefix,
  dims=12, #i think use at least 14 dims
  regressCov=TRUE, 
  regressCG=FALSE,
  k_pheno=45,
  k_umap=50,
  neigh=80, 
  dist=0.001,
  perc_cell_cov_per_window=0.3,
  method="cosine",
  output_directory=output_directory,
  cluster_on_umap=FALSE)
saveRDS(dat_sub,file=paste0(output_directory,"/","06_scaledcis.",prefix,"_finecelltyping.amethyst.rds"))

#calculate dmrs per fine celltype grouping
collapsed_dmrs<-calculate_dmrs(dat=dat_sub,split_by_group=TRUE,
                prefix=prefix,
                output_directory=output_directory)

```

Epithelial

```R
dat_sub<-subsetObject(dat,cells=row.names(dat@metadata)[dat@metadata$mcg_pct>65])
prefix=paste0("epithelial","_",as.character(x),"kb")
output_directory=paste0(project_data_directory,"/fine_celltyping/",prefix)
dat_sub<-cluster_subset(
  dat=dat_sub,
  broad_celltype=c("lumhr","lumsec","basal"), #note this is a list, doing all cells together here
  window_name="cg_win_score", #coarse_cluster_dmr_sites orinitial_cluster_5kb_win for 5kb require more coverage
  prefix=prefix,
  dims=10, #i think use at least 14 dims
  regressCov=FALSE, 
  regressCG=FALSE,
  k_pheno=200,
  k_umap=50,
  neigh=50, 
  dist=0.001,
  perc_cell_cov_per_window=0.3,
  method="cosine",
  output_directory=output_directory,
  cluster_on_umap=FALSE)

saveRDS(dat_sub,file=paste0(output_directory,"/","06_scaledcis.",prefix,"_finecelltyping.amethyst.rds"))

#calculate dmrs per fine celltype grouping
collapsed_dmrs<-calculate_dmrs(dat=dat_sub,split_by_group=TRUE,
                prefix=prefix,
                output_directory=output_directory)

```







#trying reclustering on initial dmrs for cell type
#filter on sig, and width
collapsed_dmrs<-readRDS(collapsed_dmrs,file=paste0(output_directory,"/","06_",prefix,"_finecelltyping.500bp_dmrs.filt_collapsed.rds"))
filt_init_dmrs <- collapsed_dmrs %>% filter(dmr_padj<0.05) %>% GRanges() 
filt_init_dmrs$width <- width(filt_init_dmrs)
summary(filt_init_dmrs[filt_init_dmrs$width<50000,]$width)
filt_init_dmrs <- filt_init_dmrs %>% as.data.frame() %>% filter(width<50000) %>% select(seqnames,start,end) %>% GRanges() %>% reduce() %>% as.data.frame()
dmr_bed<-as.data.frame(filt_init_dmrs)[,1:3]

window_name<-paste0("cg_",prefix,"_dmrs")
dat_sub@genomeMatrices[[window_name]] <- makeWindows(dat_sub, 
                                                     bed=dmr_bed,
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 200, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

#define subclusters
dat_sub<-cluster_subset(
        dat=dat_sub,
        broad_celltype=c("tcell","bcell","myeloid"), #note this is a list and can include multiple cell types
        window_name=window_name,
        prefix=paste0(prefix,".dmr"),
        dims=10,
        regressCov=FALSE,
        k_pheno=120,
        k_umap=50,
        neigh=8,
        dist=0.01,
        perc_cell_cov_per_window=0.2,
        method="cosine",
        output_directory=output_directory)

saveRDS(dat_sub,file=paste0(output_directory,"/","06_scaledcis.",prefix,"_finecelltyping.amethyst.rds"))
dat_sub<-readRDS(file=paste0(output_directory,"/","06_scaledcis.",prefix,"_finecelltyping.amethyst.rds"))

#calculate dmrs per fine celltype grouping
collapsed_dmrs<-calculate_dmrs(dat=dat_sub,
                prefix=prefix,
                output_directory=output_directory)

dat_immune<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/fine_celltyping/immune/06_scaledcis.immune_finecelltyping.amethyst.rds")
immune_celltype<-c(
            "06_immune.dmr.finecelltyping_6"="macrophage",
            "06_immune.dmr.finecelltyping_4"="monocyte",
            "06_immune.dmr.finecelltyping_13"="TAM",
            "06_immune.dmr.finecelltyping_15"="monocyte",
            "06_immune.dmr.finecelltyping_7"="DC",
            "06_immune.dmr.finecelltyping_2"="DC",
            "06_immune.dmr.finecelltyping_11"="TAM_2",
            "06_immune.dmr.finecelltyping_8"="bcell",
            "06_immune.dmr.finecelltyping_14"="bcell",
            "06_immune.dmr.finecelltyping_12"="tcell_cd4",
            "06_immune.dmr.finecelltyping_16"="tcell_cd4",
            "06_immune.dmr.finecelltyping_3"="tcell_cd4",
            "06_immune.dmr.finecelltyping_1"="tcell_cd8",
            "06_immune.dmr.finecelltyping_9"="tcell_cd8",
            "06_immune.dmr.finecelltyping_5"="nk_tnk",
            "06_immune.dmr.finecelltyping_10"="tcell_cd8_2")
dat_immune@metadata$celltype<-immune_celltype[dat_immune@metadata$cluster_id]
saveRDS(dat_immune,file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/fine_celltyping/immune/06_scaledcis.immune_finecelltyping.amethyst.rds")

```

### Some additional checking scripts

```R

dmr_rna_marker_overlap(dat=dat_sub,rna_markers,collapsed_dmrs,prefix,output_directory)

## DMR site overlap with RNA marker genes 
rna_markers<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.rna_markers.rds")
collapsed_dmrs<-readRDS(collapsed_dmrs,file=paste0(output_directory,"/","06_",prefix,"_finecelltyping.500bp_dmrs.filt_collapsed.rds"))

collapsed_dmrs %>% filter(direction=="hypo") %>% filter(dmr_padj<0.01) %>% filter(!is.na(gene_names)) %>% filter(type=="immune.dmr_06_immune.dmr.finecelltyping_10") %>% head(n=20) %>% select(gene_names)

bed_out<-as.data.frame(collapsed_dmrs)[collapsed_dmrs$dmr_padj<0.05,]
bed_out$strand<-ifelse(bed_out$direction=="hyper","+","-")
bed_out<-bed_out[,c("chr","dmr_start","dmr_end","strand")]
colnames(bed_out)<-c("chr","start","end","strand")
bed_out<-bed_out[!duplicated(bed_out),]
write.table(bed_out,row.names=F,col.names=F,quote=F,file=paste0(output_directory,"/","06_",prefix,"_finecelltyping.500bp_dmrs.filt_collapsed.bed"))

#one v group comparisons
celltype500bpwindows<-readRDS(file=paste0(output_directory,"/","06_",prefix,"_finecelltyping.500bp_windows.rds"))

#pairwise comparisons to check clusters
test_cluster<-"immune.dmr_06_immune.dmr.finecelltyping_4"
comparison_to<-c("immune.dmr_06_immune.dmr.finecelltyping_15")

comparison_to<-colnames(celltype500bpwindows[["pct_matrix"]])[4:ncol(celltype500bpwindows[["pct_matrix"]])]
comparison_to<-comparison_to[!(comparison_to)==test_cluster]
comparison_to<-c("immune.dmr_06_immune.dmr.finecelltyping_13","immune.dmr_06_immune.dmr.finecelltyping_15","immune.dmr_06_immune.dmr.finecelltyping_4")


comp<-cbind("name"="test",
            "A"=test_cluster,
            "B"=paste(comparison_to,collapse=","))

collapsed_dmrs<-find_cluster_markers(dat=dat_sub,celltype500bp_windows,comp)


collapsed_dmrs %>% 
    filter(direction=="hypo") %>% 
    filter(dmr_length<20000) %>% 
    filter(dmr_padj<0.05) %>% 
    filter(!(gene_names=="NA")) %>% 
    arrange(dmr_logFC) %>% 
    as.data.frame() %>% 
    slice_min(dmr_logFC,n=300) %>% 
    select(gene_names) %>% 
    paste(collapse=",")


#comparisons: If eachVsAll is not desired, provide a data frame
  #       describing which tests to run. The data.frame should have
  #       three columns with rows describing conditions of each test.
  #       "name" determines the name of the test in the output; "A"
  #       lists group members, and "B" lists group nonmembers.

```

## stromal

```R
#Set up variables for output.
prefix="fibroblast"
output_directory=paste0(project_data_directory,"/fine_celltyping/",prefix)
dat<-readRDS(file="05_scaledcis.coarse_clusters.amethyst.rds")

dat_sub<-subsetObject(dat,cells=row.names(dat@metadata)[dat@metadata$mcg_pct>65])

#define subclusters
dat_sub<-cluster_subset(
  dat=dat_sub,
  broad_celltype=c("fibroblast"), #note this is a list
  window_name="coarse_cluster_dmr_sites",
  prefix=prefix,
  dims=18,
  regressCov=TRUE,
  regressCG=TRUE,
  k_pheno=75,
  k_umap=8,
  neigh=10,
  dist=0.01,
  perc_cell_cov_per_window=0.05,
  method="cosine",
  output_directory=output_directory)


saveRDS(dat_sub,file=paste0(output_directory,"/","06_scaledcis.",prefix,"_finecelltyping.amethyst.rds"))

dat_sub<-readRDS(file=paste0(output_directory,"/","06_scaledcis.",prefix,"_finecelltyping.amethyst.rds"))

#calculate dmrs per fine celltype grouping
collapsed_dmrs<-calculate_dmrs(dat=dat_sub,
                prefix=prefix,
                output_directory=output_directory)

collapsed_dmrs<-readRDS(file=paste0(output_directory,"/","06_",prefix,"_finecelltyping.500bp_dmrs.filt_collapsed.rds"))
#trying reclustering on initial dmrs for cell type
#filter on sig, and width
filt_init_dmrs <- collapsed_dmrs %>% filter(dmr_padj<0.05) %>% GRanges() 
filt_init_dmrs$width <- width(filt_init_dmrs)
summary(filt_init_dmrs[filt_init_dmrs$width<50000,]$width)
filt_init_dmrs <- filt_init_dmrs %>% as.data.frame() %>% filter(width<50000) %>% select(seqnames,start,end) %>% GRanges() %>% reduce() %>% as.data.frame()
dmr_bed<-as.data.frame(filt_init_dmrs)[,1:3]

window_name<-paste0("cg_",prefix,"_dmrs")
dat_sub@genomeMatrices[[window_name]] <- makeWindows(dat_sub, 
                                                     bed=dmr_bed,
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 200, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

#define subclusters
dat_sub2<-cluster_subset(
        dat=dat_sub,
        broad_celltype=c("endothelial"), #note this is a list and can include multiple cell types
        window_name=window_name,
        prefix=paste0(prefix,".dmr"),
        dims=18,
        regressCov=FALSE,
        k_pheno=120,
        k_umap=20,
        neigh=6,
        dist=0.01,
        perc_cell_cov_per_window=0.2,
        method="cosine",
        output_directory=output_directory)


saveRDS(dat_sub2,file=paste0(output_directory,"/","06_scaledcis.",prefix,"_finecelltyping.amethyst.rds"))
dat_sub<-readRDS(file=paste0(output_directory,"/","06_scaledcis.",prefix,"_finecelltyping.amethyst.rds"))

#calculate dmrs per fine celltype grouping
collapsed_dmrs<-calculate_dmrs(dat=dat_sub,
                prefix=prefix,
                output_directory=output_directory)

## DMR site overlap with RNA marker genes 
rna_markers<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.rna_markers.rds")
collapsed_dmrs<-readRDS(collapsed_dmrs,file=paste0(output_directory,"/","06_",prefix,"_finecelltyping.500bp_dmrs.filt_collapsed.rds"))
dmr_rna_marker_overlap(dat=dat_sub,rna_markers,collapsed_dmrs,prefix,output_directory)


dat_stromal<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/fine_celltyping/stromal/06_scaledcis.stromal_finecelltyping.amethyst.rds")

stromal_celltype<-c(
            "06_stromal.dmr.finecelltyping_9"="endothelial",
            "06_stromal.dmr.finecelltyping_3"="endothelial",
            "06_stromal.dmr.finecelltyping_1"="TEC",
            "06_stromal.dmr.finecelltyping_13"="fibroblast",
            "06_stromal.dmr.finecelltyping_11"="fibroblast",
            "06_stromal.dmr.finecelltyping_2"="fibroblast",
            "06_stromal.dmr.finecelltyping_4"="fibroblast",
            "06_stromal.dmr.finecelltyping_10"="fibroblast",
            "06_stromal.dmr.finecelltyping_6"="fibroblast",
            "06_stromal.dmr.finecelltyping_12"="fibroblast",
            "06_stromal.dmr.finecelltyping_7"="pericyte_VSMC",
            "06_stromal.dmr.finecelltyping_8"="CAF",
            "06_stromal.dmr.finecelltyping_5"="CAF")

dat_stromal@metadata$celltype<-stromal_celltype[dat_stromal@metadata$cluster_id]
saveRDS(dat_stromal,file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/fine_celltyping/stromal/06_scaledcis.stromal_finecelltyping.amethyst.rds")

```

From here I looked at bigwig files in IGV to assess genes quickly. For most genes that worked, it was clear from promoter methylation. I tried to rely mostly on lineage markers and cannonical genes. But I also used sn/sc RNA data sets and some markers of methylation changes from monocyte to macrophage differentiation.

For detailed loci and cell type calling:
See /Volumes/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/fine_celltyping/fine_celltyping_markers.xlsx 

Assign cell types:

```R
dat<-readRDS(file="05_scaledcis.coarse_clusters.amethyst.rds")
dat_immune<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/fine_celltyping/immune/06_scaledcis.immune_finecelltyping.amethyst.rds")
dat_stromal<-readRDS(file="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/fine_celltyping/stromal/06_scaledcis.stromal_finecelltyping.amethyst.rds")

dat@metadata$celltype<-dat@metadata$broad_celltype
dat@metadata$immune_subcluster<-NA
dat@metadata$stromal_subcluster<-NA
dat@metadata[row.names(dat_immune@metadata),]$celltype<-dat_immune@metadata$celltype
dat@metadata[row.names(dat_stromal@metadata),]$celltype<-dat_stromal@metadata$celltype

dat@metadata[row.names(dat_immune@metadata),]$immune_subcluster<-dat_immune@metadata$cluster_id
dat@metadata[row.names(dat_stromal@metadata),]$stromal_subcluster<-dat_stromal@metadata$cluster_id

saveRDS(dat,file="06_scaledcis.celltype.amethyst.rds")



```

Next we will add in cancer cell identities using CNV data.


Use RNA markers and DMR markers to define cell types.

Using library GeneOverlap because I always get turned around when making contingency tables.


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




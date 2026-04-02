# Final cell typing

Cell type assignment by order of importance:
1. CNVs identified by copy number clones assigned cancer in scalemet_dcis/processing/milestonev1_03.2_copykit_cnv_calling.md
2. Cluster based on DMRs in scalemet_dcis/processing/milestonev1_04.1_liger_scrna_met_integration.md assigned cell types after integration
3. Assign cell types for clusters with clear (>75% of cell type consistently assigned)
4. Identify clusters which group together on umap at lower cluster resolution (assign celltype by shared cluster)
5. Confirm cell types by RNA based marker genes and hypomethylation over gene body

Use integration to finalize cell typing. Note some cell types were mislabelled based on small set of marker genes originally used. 
Using integration data for final classifications.

```R
BiocManager::install("HDF5Array")
install.packages('leidenAlg')
install.packages('rliger')
```


```R
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

set.seed(111)
options(future.globals.maxSize= 200000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)

obj<-readRDS(file="07_scaledcis.cnv_clones.amethyst.rds")

#color assignment is fluor as cancer associated cell type, rest muted versions
celltype_col=c(
"cancer"="#ff00ff",
"basal"="#92278F",
"lumsec"="#AD8CFF",
"lumhr"="#FF6AD5",

"TEC"="#ff8000",
"endothelial"="#ffab5f",

"CAF"="#ff2222",
"pericyte_VSMC"="#FA7876",
"fibroblast"="#9b1c31",

"TAM"="#ccff00",
"TAM2"="#00ff80",
"monocyte"="#98d3b9",
"macrophage"="#00af5f",
"DC"="#008080",

"nk_tnk"="#00ffff",
"tcell_cd4"="#00bae5",
"tcell_cd8"="#1800ff",
"tcell_cd8_2"="#0016b7",
"bcell"="#87ceeb")

```

## Subcluster on just immune cells
Using the same wrapping functions for final clustering as I used for fine celltyping of immune and stromal
processing/milestonev1_01.2_amethyst_fine_celltyping.md

```R

celltype_umap<-function(obj=dat,prefix="allcells",dims=12,regressCov=TRUE,cluster_on_umap=TRUE,k_pheno=50,k_umap=50,neigh=25,dist=1e-5,method="cosine",output_directory,window_name){
  print("Running IRLBA reduction...")
  obj@reductions[[paste(window_name,"irlba",sep="_")]] <- runIrlba(obj, genomeMatrices = c(window_name), dims = dims, replaceNA = c(0))

  if(regressCov){
    print("Running regression...")
  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- regressCovBias(obj, reduction = paste(window_name,"irlba",sep="_")) 
  } else {
      print("Skipping regression...")
  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- obj@reductions[[paste(window_name,"irlba",sep="_")]]
  }

  obj <- amethyst::runCluster(obj, k_phenograph = k_umap, reduction = paste(window_name,"irlba_regressed",sep="_")) 

  print(paste("Running UMAP...",as.character(neigh),as.character(dist),as.character(method)))
  obj <- amethyst::runUmap(obj, neighbors = neigh, dist = dist, method = method, reduction = paste(window_name,"irlba_regressed",sep="_")) 

  if(cluster_on_umap){
      print("Clustering on umap. Overwritting clustering on IRLBA")
      umap_clus<-obj@metadata%>% select(c("umap_x","umap_y"))
      umap_clus<-Rphenograph::Rphenograph(umap_clus,k=k_pheno)
      obj@metadata$cluster_id<-paste0(prefix,"_",igraph::membership(umap_clus[[2]]))
  } else {
      print("Clustering on IRLBA. Overwritting current clustering on IRLBA")
      irlba_clus<-obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]][row.names(obj@metadata),]
      irlba_clus<-Rphenograph::Rphenograph(irlba_clus,k=k_pheno)
      obj@metadata$cluster_id<-paste0(prefix,"_",igraph::membership(irlba_clus[[2]]))
  }

  outname=paste(prefix,"integrated_celltype",dims,as.character(regressCov),k_pheno,neigh,as.character(dist),method,sep="_")
  print(paste("Plotting...",outname))

  p1 <- dimFeature(obj, colorBy = cluster_id, reduction = "umap",pointSize=1) + ggtitle(paste(window_name,"Cluster"))
  p2 <- dimFeature(obj, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red"),guide="colourbar") + ggtitle("Coverage distribution")
  p3 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
  p4 <- dimFeature(obj, colorBy = celltype, reduction = "umap",pointSize=1) + ggtitle(paste(window_name,"Clusters")) + scale_color_manual(values=celltype_col)
  p5 <- dimFeature(obj, colorBy = Group, reduction = "umap",pointSize=1) + ggtitle(paste(window_name," Group"))
  p6 <- dimFeature(obj, colorBy = Sample, reduction = "umap",pointSize=1) + ggtitle(paste(window_name," Samples"))

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
  k_umap=50,
  neigh=25,
  perc_cell_cov_per_window=0.01, #required window coverage for clustering (1%)
  dist=1e-5,
  method="cosine"){

  print(paste("Making fine celltyping directory:",output_directory))
  system(paste0("mkdir -p ",output_directory))
  print(paste("Subsetting object based on broad celltype:",broad_celltype))

  #subset amethyst object to broad celltype
  dat_sub<-subsetObject(dat,cells=row.names(dat@metadata[dat@metadata$broad_celltype %in% broad_celltype,]))

  #subset to windows with 5% coverage
  req_cov<-as.integer(nrow(dat_sub@metadata)*perc_cell_cov_per_window)
  print(paste("Requiring",perc_cell_cov_per_window,"of cells for window coverage filter:",as.character(req_cov),"cells."))
  dat_sub@genomeMatrices[[window_name]] <- dat_sub@genomeMatrices[[window_name]][rowSums(!is.na(dat_sub@genomeMatrices[[window_name]])) >= req_cov, ]
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
                        window_name=window_name,
                        output_directory=output_directory) 

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

calculate_dmrs<-function(dat=dat,prefix=prefix,groupBy="fine_cluster_id",output_directory=output_directory){
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

  saveRDS(celltype500bpwindows,file=paste0(output_directory,"/","06_",prefix,"_celltype.500bp_windows.rds"))
  dat@genomeMatrices[[paste0("cg_",prefix,"_cells_perc")]] <- celltype500bpwindows[["pct_matrix"]]

  #output tracks as bigBed
  bigwig_output(obj=dat,
                tracks=paste0("cg_",prefix,"_cells_perc"),
                output_directory=output_directory,
                prefix=prefix)
  #save subset amethyst file
  print(paste("Saving new amethyst file:",paste("06_scaledcis.",prefix,"_celltype.amethyst.rds")))
  saveRDS(dat,file=paste0(output_directory,"/","06_scaledcis.",prefix,"_celltype.amethyst.rds"))

  pct_mat<-celltype500bpwindows[["pct_matrix"]] 
  sum_mat<-celltype500bpwindows[["sum_matrix"]] 

  dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
          eachVsAll = TRUE, # If TRUE, each group found in the sumMatrix will be tested against all others
          nminTotal = 3, # Min number observations across all groups to include the region in calculations
          nminGroup = 3) # Min number observations across either members or nonmembers to include the region
  saveRDS(dmrs,file=paste0(output_directory,"/","06_",prefix,"_celltype",".500bp_dmrs.rds"))

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
  saveRDS(collapsed_dmrs,file=paste0(output_directory,"/","06_",prefix,"_celltype.500bp_dmrs.filt_collapsed.rds"))

  rename_dmr_output<-setNames(nm=1:length(unique(collapsed_dmrs$test)),gsub(colnames(sum_mat)[grepl(colnames(sum_mat),pattern="_t$")],pattern="_t",replacement=""))
  collapsed_dmrs$type <- rename_dmr_output[collapsed_dmrs$test]
  saveRDS(collapsed_dmrs,file=paste0(output_directory,"/","06_",prefix,"_celltype.500bp_dmrs.filt_collapsed.rds"))

  #plot dmr counts per cluster
  pal <- colorRampPalette(colors = c("#F9AB60", "#E7576E", "#630661", "#B5DCA5"))
  COLS <- pal(length(unique(collapsed_dmrs$type)))

  plt<-ggplot(collapsed_dmrs |> dplyr::group_by(type, direction) |> dplyr::summarise(n = n()), 
      aes(y = type, x = n, fill = type)) + 
      geom_col() + 
      facet_grid(vars(direction), scales = "free_y") + 
      scale_fill_manual(values = COLS) + 
      theme_classic()
  ggsave(plt,file=paste0(output_directory,"/","06_",prefix,"_celltype.500bp_dmrs.filt_collapsed.barplot.pdf"))
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

#################################
#Get DMR per Integrated Celltypes
#################################


```R
#celltype
dmr_outdir<-"/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/DMR_analysis"
dmr_celltype_outdir=paste(sep="/",dmr_outdir,"celltype")
system(paste("mkdir -p", dmr_celltype_outdir))


celltype500bpwindows <- calcSmoothedWindows(obj, 
                                        type = "CG", 
                                        threads = 100,
                                        step = 500, # change to 500 for real data unless you have really low coverage
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "celltype",
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)
saveRDS(celltype500bpwindows,file=paste0(dmr_celltype_outdir,"/","dmr_analysis.celltype.500bp_windows.rds"))
celltype500bpwindows<-readRDS(file=paste0(dmr_celltype_outdir,"/","dmr_analysis.celltype.500bp_windows.rds"))

pct_mat<-celltype500bpwindows[["pct_matrix"]] 
sum_mat<-celltype500bpwindows[["sum_matrix"]] 

#save clone object for future genome track plotting
obj@genomeMatrices[["cg_celltype_tracks"]] <- pct_mat #load it into amethyst object for plotting

dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
        eachVsAll = TRUE, # If TRUE, each group found in the sumMatrix will be tested against all others
        nminTotal = 3, # Min number observations across all groups to include the region in calculations
        nminGroup = 3) # Min number observations across either members or nonmembers to include the region

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



rename_dmr_output<-setNames(nm=1:length(unique(collapsed_dmrs$test)),gsub(colnames(sum_mat)[grepl(colnames(sum_mat),pattern="_t$")],pattern="_t",replacement=""))
collapsed_dmrs$type <- rename_dmr_output[collapsed_dmrs$test]
saveRDS(collapsed_dmrs,file=paste0(dmr_celltype_outdir,"/","celltype_allcells",".dmr_filt_collapse.rds"))


#clear that plasma and epithelial is too few cell count
# rename them when happy with umap projections

#filtering dmrs
#use only hypo dmrs, take top 5000 per celltype, filter by length and merge any that overlap (across cell types)
dmrs <- collapsed_dmrs %>% 
  filter(direction=="hypo") %>% 
  filter(dmr_padj<0.05) %>% 
  filter(dmr_length < 20000) %>% 
  group_by(type) %>% 
  #slice_min(n=20000, dmr_logFC) %>%
  filter(abs(dmr_logFC)>1.5) %>%
  as.data.frame() %>%
  GenomicRanges::makeGRangesFromDataFrame() %>%   GenomicRanges::reduce()

summary(width(dmrs))
length(dmrs)
dmr_bed<-data.frame(chr=seqnames(dmrs),start=start(dmrs),end=end(dmrs))


#create new matrix from DMR sites for refined clustering
window_name="celltype_dmr_sites"


obj@genomeMatrices[[window_name]] <- makeWindows(obj, 
                                                     bed = dmr_bed,
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 50, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

obj@genomeMatrices[[window_name]] <- obj@genomeMatrices[[window_name]][rowSums(!is.na(obj@genomeMatrices[[window_name]])) >= 45, ]
est_dim<-dimEstimate(obj, genomeMatrices = c(window_name), dims = c(20), threshold = 0.95)
print(est_dim)

saveRDS(obj,file="06_scaledcis.celltype.amethyst.rds")

#plotting all cells together, also overclustering a bit just to clean up cell types potentially
#if a cluster is >90% one celltype others are assigned the same

obj<-celltype_umap(obj=obj,
              prefix="06_scaledcis.celltype",
              dims=13, #13 for b cell separation
              regressCov=FALSE,
              k_pheno=10, #100 (for reassignment of celltype per cluster)
              k_umap=20, #20
              cluster_on_umap=FALSE, #on irlba looks better and tracks with cell types better
              neigh=15, #50
              dist=0.01,
              method="cosine",
              output_directory=wd,
              window_name)

saveRDS(obj,file="06_scaledcis.celltype.amethyst.rds")


celltype_per_cluster<-table(obj@metadata$celltype,obj@metadata$cluster_id)
celltype_per_cluster<-t(t(celltype_per_cluster)/colSums(celltype_per_cluster))
col_fun = colorRamp2(c(0, 2), c("white", "red"))
plt<-Heatmap(celltype_per_cluster,col=col_fun)
pdf(paste0("06_scaledcis.celltype",".cluster.heatmap.pdf"))
plt
dev.off()

celltype_per_cluster<-table(obj@metadata$celltype,obj@metadata$Group)
celltype_per_cluster<-t(t(celltype_per_cluster)/colSums(celltype_per_cluster))
col_fun = colorRamp2(c(0, 2), c("white", "red"))
plt<-Heatmap(celltype_per_cluster,col=col_fun)
pdf(paste0("06_scaledcis.celltype",".group.heatmap.pdf"))
plt
dev.off()

#correct cluster assignment and run on celltype again?
#the proportion of cells in other clusters is low, or it makes sense, so i'm trusting the fine celltyping output.
collapsed_dmrs<-calculate_dmrs(dat=obj,groupBy="celltype",prefix="06_scaledcis.celltype",output_directory=wd)

#a little bonus object cleanup
#clean up unneeded tracks (they are in previous amethyst objects)
obj@genomeMatrices[["cg_100k_score"]]<-NULL
obj@genomeMatrices[["initial_cluster_5kb_win"]]<-NULL
obj@genomeMatrices[["cg_5kbp_initial_clusters_tracks"]]<-NULL
obj@genomeMatrices[["5kbp_initial_clusters_dmr_sites"]]<-NULL
obj@genomeMatrices[["window_name"]]<-NULL
saveRDS(obj,file="06_scaledcis.celltype.amethyst.rds")

```

Output bedgraph tracks with colors to match cell type marker plots. 
Save as IGV session for sharing.

```R

#output DMRs as bed file as well
obj@genomeMatrices$cg_celltype_tracks

#output as bed with RGB and value, convert bed to bigbed?
#bigwig output of 500bp resolution tracks along the genome per celltype


#split bw into 4 files per cell type by methylation/average methylation
for(i in colnames(obj@genomeMatrices$cg_celltype_tracks)[4:ncol(obj@genomeMatrices$cg_celltype_tracks)]){
    out_dat<-obj@genomeMatrices$cg_celltype_tracks %>% select(chr,start,end,i) 

    hg38_seq_info<-Seqinfo(genome="hg38")
    out_dat<-GRanges(out_dat[complete.cases(out_dat),]) #filter NA
    out_dat<-out_dat[out_dat@seqnames %in% hg38_seq_info@seqnames,] #filter chr
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
    rtracklayer::export(out_dat_hypermet,con=paste(i,"hypermet","bw",sep="."))
    rtracklayer::export(out_dat_met_mid,con=paste(i,"midmet","bw",sep="."))
    rtracklayer::export(out_dat_met_low,con=paste(i,"lowmet","bw",sep="."))
    rtracklayer::export(out_dat_met_hypomet,con=paste(i,"hypomet","bw",sep="."))
    print(paste(i,celltype_col[i]))
    }
```


```R
obj<-readRDS(file="06_scaledcis.celltype.amethyst.rds")
#run 3d plotting just for fun
dim3d<-uwot::umap(
  X=obj@reductions[["celltype_dmr_sites_irlba_regressed"]],
  n_neighbors = 15,
  n_components = 3,
  metric = "cosine",
  seed = 123,
  n_threads=50,
)


dim3d<-as.data.frame(dim3d)
colnames(dim3d)<-c("X","Y","Z")
dim3d$cellname<-row.names(dim3d)
dim3d$celltype<-as.data.frame(obj@metadata[row.names(dim3d),])$celltype
dim3d$hex_color<-paste0(celltype_col[dim3d$celltype]
dim3d$r_col<-unlist(col2rgb(dim3d$hex_color)["red",])
dim3d$g_col<-unlist(col2rgb(dim3d$hex_color)["green",])
dim3d$b_col<-unlist(col2rgb(dim3d$hex_color)["blue",])

write.table(dim3d,col.names=T,file="06_scaledcis.celltype.3dumap.csv",sep=",",row.names=F)

```

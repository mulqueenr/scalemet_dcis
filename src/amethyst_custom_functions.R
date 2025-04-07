#run in amethyst.sif
#singularity shell \
#--bind ~/projects/ \
#--bind /volumes/seq/projects/metACT \
#~/singularity/amethyst.sif

library(amethyst)
library(rhdf5)
library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(plyr)
library(future)
library(furrr)
library(purrr)
library(cowplot)
library(pheatmap)
library(plyr)
library(parallel)
library(rtracklayer)
library(tidyverse)
library(colorspace)
library(GeneNMF) #new from here down
library(magrittr)
library(universalmotif)
library(ape)
library(ggtree)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(parallel)
library(chromVAR)
library(SummarizedExperiment)
library(ggseqlogo)
library(patchwork)

#source("~/projects/metact/src/amethyst_custom_functions.R") to load in

#promoter extension
extend <- function(x, upstream=0, downstream=0)     
{
    if (any(strand(x) == "*"))
        warning("'*' ranges were treated as '+'")
    on_plus <- strand(x) == "+" | strand(x) == "*"
    new_start <- start(x) - ifelse(on_plus, upstream, downstream)
    new_end <- end(x) + ifelse(on_plus, downstream, upstream)
    ranges(x) <- IRanges(new_start, new_end)
    trim(x)
}

#read in metadata, if else for alpha vs new file
setup_amethyst_metadata<-function(
    in_dir="/volumes/USR2/Ryan/projects/metact/240526_RMMM_scalebio_dcis/transfer_dat/", #sample directory
    sample_list=c("hbca-83l","hbca-16r","DCIS-41T","DCIS-66T"), #list of sample names
    runid=1,
    cg_cov_filter=10000
    ){

    if(endsWith(in_dir,suffix="/transfer_dat/")){ #processing of earlier scalebio pipeline output (alpha kit), to be deprocated in future
        samp_dir=list.files(path=paste0(in_dir,"/report/"))

        met<-lapply(samp_dir,function(i){
            samp_met<-read.csv(paste0(in_dir,"/report/",i,"/csv/",i,".passingCellsMapMethylStats.csv"))
            row.names(samp_met)<-samp_met$BC
            return(samp_met)
        })
        metadat<-do.call("rbind",met)
        metadat$run<-runid
        metadat<-metadat[metadat$CG_Cov>cg_cov_filter,]
        metadat<-metadat[toupper(metadat$sampleName) %in% toupper(sample_list),]
        colnames(metadat)<-c("cell_id","cov","cg_cov","mcg_pct","ch_cov","mch_pct","sample","total_reads",
        "passing_reads","unique_reads","mito_reads","percent","pct_uniq_pass","pct_pass_total","pct_pass_cell",
        "pass","threshold","tgmt","tgmt_well","i5","i5_well","i7","i7_well","cg_total_ratio","run") #matching to newer metadata outputs

    } else {
        samp_dir=list.files(path=in_dir,pattern="*csv",full.names=T)
        met<-lapply(samp_dir,function(i){
            samp_met<-read.csv(i)
            row.names(samp_met)<-samp_met$cell_id
            return(samp_met)
        })
        metadat<-do.call("rbind",met)
        metadat$run<-runid
        metadat<-metadat[metadat$cg_cov>cg_cov_filter,]
        metadat<-metadat[complete.cases(metadat),]
        metadat<-metadat[toupper(metadat$sample) %in% toupper(sample_list),]
    }
    return(metadat)
    }

correct_h5_cellnames<-function(h5,runid){
    print(paste("Correcting names for...", basename(h5)))
    run_number=as.character(runid)
    h5list = h5ls(h5)
    for (i in 1:nrow(h5list)){
        tryCatch({if(endsWith(h5list[i,"group"],"CG")){
            if( !(paste0(h5list[i,"name"],"_",run_number) %in% h5list$name) &&
              !endsWith(h5list[i,"name"],run_number)){
                celldat<-h5read(h5,paste0("CG/",h5list[i,"name"]))
                h5write(celldat, file=h5, name=paste0("CG/",h5list[i,"name"],"_",run_number))
            }
        }
        },error =function(e) { cat("Proceeding past line",i,"for",basename(h5),"\n")} )
    }
 }


#clustering function, supply a bed file or arguments for stepsize windows
cluster_by_windows<-function(
    obj,
    window_name,
    prefix,
    stepsize=NULL,
    bed=NULL,
    metric="score",
    threads=200,
    neighbors=50,
    distance=0.2,
    genes=NULL,
    promoter=FALSE,
    overwrite_windows=TRUE,
    k_phenograph=175){
    if(!(window_name %in% names(obj@genomeMatrices)) || overwrite_windows){
    print(paste("Making window summaries for ",window_name))
    obj@genomeMatrices[[window_name]] <- makeWindows(obj,
                                                      stepsize = stepsize, 
                                                      type = "CG", 
                                                      metric = metric, 
                                                      bed = bed, genes=genes, promoter=promoter,
                                                      threads = threads, 
                                                      index = "chr_cg", 
                                                      nmin = 2) 
    }
  print(paste("Estimating dimensions..."))                                           
  #filter windows by cell coverage
  obj@genomeMatrices[[window_name]] <- obj@genomeMatrices[[window_name]][rowSums(!is.na(obj@genomeMatrices[[window_name]])) >= 45, ]
  est_dim<-dimEstimate(obj, genomeMatrices = c(window_name), dims = c(10), threshold = 0.95)
  print(est_dim)
  set.seed(123)
  print("Running IRLBA reduction...")
  obj@reductions[[paste(window_name,"irlba",sep="_")]] <- runIrlba(obj, genomeMatrices = c(window_name), dims = est_dim, replaceNA = c(0))
  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- regressCovBias(obj, reduction = paste(window_name,"irlba",sep="_")) # Optional; helps reduce coverage bias
  print("Clustering on coverage regressed reduction...")

  obj <- runCluster(obj, k_phenograph = k_phenograph, reduction = paste(window_name,"irlba_regressed",sep="_")) # consider increasing k_phenograph to 50 for larger datasets
  obj <- runUmap(obj, neighbors = neighbors, dist = distance, method = "euclidean", reduction = paste(window_name,"irlba_regressed",sep="_")) 
    
  print("Plotting...")

  p1 <- dimFeature(obj, colorBy = cluster_id, reduction = "umap") + ggtitle(paste(window_name,"Clusters"))
  p2 <- dimFeature(obj, colorBy = sample, reduction = "umap") + ggtitle(paste(window_name,"Samples"))
  p3 <- dimFeature(obj, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Coverage distribution")
  p4 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
  
  plt<-plot_grid(p1, p2,p3, p4,ncol=2)
  ggsave(plt,file=paste(prefix,window_name,"umap.pdf",sep="."))     
  return(obj)                                             
}



###### 1kb window generation, and DMR Calculation ####
dmr_and_1kb_window_gen<-function(
    obj=hbca,
    prefix="hbca",
    groupBy="cluster_id"){
    cluster1kbwindows <- calcSmoothedWindows(obj, 
                                            type = "CG", 
                                            threads = 300,
                                            step = 1000,
                                            smooth = 3,
                                            index = "chr_cg",
                                            groupBy = groupBy, 
                                            returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                            returnPctMatrix = TRUE)
    obj@genomeMatrices[[paste0("cg_",groupBy,"_tracks")]] <- cluster1kbwindows[["pct_matrix"]]

    pal=c("#E5E6E4","#CFD2CD","#A6A2A2","#847577","#6E44FF")
    dmrs<-testDMR(cluster1kbwindows[["sum_matrix"]], eachVsAll = TRUE, nminTotal = 0, nminGroup = 0) 
    dmrs2<-filterDMR(dmrs, method = "bonferroni", filter = TRUE,pThreshold=0.05,logThreshold=0.5) #add additional columns direction column
    celltype_test<-setNames(nm=unique(dmrs2$test), gsub("_c$","",colnames(cluster1kbwindows[["sum_matrix"]])[grepl(pattern="_c$",colnames(cluster1kbwindows[["sum_matrix"]]))]))
    dmrs2$celltype<-celltype_test[dmrs2$test]
    saveRDS(dmrs2,file=paste0(prefix,".dmr.",groupBy,".rds"))

    collapsed_dmrs <- collapseDMR(obj, dmrs2, maxDist = 4000, minLength = 1000, reduce = T, annotate = T) 
    collapsed_dmrs$celltype<-celltype_test[collapsed_dmrs$test]
    saveRDS(collapsed_dmrs,file=paste0(prefix,".dmr.",groupBy,".collapsed.rds"))

    plt<-ggplot(collapsed_dmrs |> dplyr::group_by(celltype, direction) |> dplyr::summarise(n = dplyr::n()), 
       aes(y = celltype, x = n, fill = celltype)) + geom_col() + 
        facet_grid(vars(direction), scales = "free_y") + 
        scale_fill_manual(values = makePalette(option = 7, n = length(unique(collapsed_dmrs$celltype)) ) ) + 
        theme_classic()
    
    ggsave(plt,file=paste0(prefix,".met_per_dmr.",groupBy,".pdf"))
    median(collapsed_dmrs$dmr_end-collapsed_dmrs$dmr_start)

    top_dmrs <- collapsed_dmrs |> 
    dplyr::group_by(celltype, direction) |> 
    dplyr::arrange(dmr_padj, .by_group = TRUE) |> dplyr::mutate(rank_padj = 1:dplyr::n()) |>
    dplyr::arrange(desc(abs(dmr_logFC)), .by_group = TRUE) |> dplyr::mutate(rank_logFC = 1:dplyr::n()) |>
    rowwise() |> dplyr::mutate(total_rank = sum(rank_padj, rank_logFC)) |> 
    group_by(celltype, direction) |> slice_min(n = 20, order_by = total_rank) |>
    dplyr::mutate(location = paste0(chr, "_", (dmr_start - 2000), "_", (dmr_end + 2000))) |> dplyr::arrange(direction)
    write.table(top_dmrs,file=paste0(prefix,".cluster_dmr.",groupBy,".tsv"))
    return(obj)
}


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
        export(out_bw,con=paste(tracks,cluster,"bw",sep="."),format='bigWig')}
}


#GSEA of clone DMRs
gsea_enrichment<-function(prefix,dmrs,gene_universe,category="C3",subcategory="TFT:GTRD",out_setname="TFT"){
  top_p_gsea <- do.call("rbind",
    lapply(unique(dmrs$celltype), 
    function(i) {
    #gene set
    gene_list<-dmrs |> dplyr::filter(celltype==i) |> pull(gene_names)
    gene_list<-str_replace_all(unlist(lapply(strsplit(gene_list, ","),unlist)), " ", "") 
    out<-runGSEA(gene_list, universe=gene_universe, category = category,subcategory=subcategory)
    out$celltype<-i
    return(out)
    }
    ))
pltdat<-top_p_gsea %>% group_by(celltype) %>% slice_max(order_by = -padj, n = 5)
plt<-ggplot(pltdat,aes(x=celltype,y=pathway))+
geom_point(aes(size = -log10(padj), fill = overlap/size), shape=21)+
theme_minimal() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(plt,file=paste0(prefix,"_GSEA","_",out_setname,".pdf"),width=20)
saveRDS(top_p_gsea,file=paste0(prefix,"_GSEA","_",out_setname,".rds"))
}

clone_dmr<-function(obj=obj,sample=c("DCIS-41T"),prefix="DCIS-41T",k_phenograph=50){
  dcis<-subsetObject(obj, cells = row.names(obj@metadata[obj@metadata$sample %in% sample,]))
  #recluster just 41 on 100kb windows
  dcis<-cluster_by_windows(obj=dcis,
                          prefix=prefix,
                          window_name="cg_100k_score",
                          stepsize=100000,
                          threads=200,
                          neighbors=50,
                          distance=0.05,
                          overwrite_windows=FALSE,
                          k_phenograph=200)
  p1 <- dimFeature(dcis, colorBy = subclones, reduction = "umap") + ggtitle(paste(prefix, "Subclones"))
  ggsave(p1,file=paste(prefix,"subclones_umap.pdf",sep="."))     
  dcis<-dmr_and_1kb_window_gen(obj=dcis,prefix=prefix,groupBy="subclones")
  dcis<-dmr_and_1kb_window_gen(obj=dcis,prefix=prefix,groupBy="cluster_id")
  saveRDS(dcis,paste0(prefix,".amethyst.rds"))
  #GSEA processing to make some sense of DMRS
  dmrs<-readRDS(paste0(prefix,".dmr.subclones.collapsed.rds"))
  dmrs <- dmrs |> dplyr::filter(direction=="hypo") |> dplyr::filter(gene_names != "NA") 

  #set up gene universe
  gene_universe<-dmrs |> pull(gene_names)
  gene_universe<-str_replace_all(unlist(lapply(strsplit(gene_universe, ","),unlist)), " ", "") #flatten and remove 
  gene_universe<-gene_universe[!duplicated(gene_universe)]

  #run gsea enrichment at gene level on different sets
  gsea_enrichment(prefix=prefix,dmrs=dmrs,gene_universe=gene_universe,
  category="C3",subcategory="TFT:GTRD",out_setname="TFT") #find enrichment in tft (transcription factor targets)

  gsea_enrichment(prefix=prefix,dmrs=dmrs,gene_universe=gene_universe,
  category="C1",subcategory=NULL,out_setname="position") #find enrichment in c1 signatures (positional)

  gsea_enrichment(prefix=prefix,dmrs=dmrs,gene_universe=gene_universe,
  category="H",subcategory=NULL,out_setname="hallmark") #find cancer hallmark signatures

}


methyltree_output<-function(obj=obj,prefix="DCIS-41T",sample="DCIS-41T",filt_min_pct=10,filt_max_pct=80,threads=1){
        dcis<-subsetObject(obj, cells = row.names(obj@metadata[obj@metadata$sample %in% sample,]))
        dcis@metadata$methyltree_group<-"all"
        #make 500bp windows with methylation percentages
        methyltreewindows <- calcSmoothedWindows(dcis, 
                                            type = "CG", 
                                            threads = threads,
                                            step = 500,
                                            smooth = 1,
                                            index = "chr_cg",
                                            groupBy = "methyltree_group", 
                                            returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                            returnPctMatrix = TRUE)
        print(paste("Starting number of windows:",as.character(nrow(methyltreewindows[["pct_matrix"]]))))
        methyltreewindows[["pct_matrix"]]<-methyltreewindows[["pct_matrix"]][methyltreewindows[["pct_matrix"]]$all>=filt_min_pct & methyltreewindows[["pct_matrix"]]$all<=filt_max_pct,]
        #filter to windows to middling methylation values
        print(paste("Filtering by m0 >=", as.character(filt_min_pct), "m1 <=", as.character(filt_max_pct),as.character(nrow(methyltreewindows[["pct_matrix"]]))))
        #merge windows that are touching
        methyltreewindows<-reduce(GenomicRanges::makeGRangesFromDataFrame(methyltreewindows[["pct_matrix"]]))
        print(paste("Filtered window count:",as.character(nrow(as.data.frame((methyltreewindows))))))
        print(paste("Filtered window average width:",as.character(mean(width(methyltreewindows)))))
        print(paste("Total genome covered:",as.character(sum(width(methyltreewindows))/1000000),"Mbp"))
        #make a merged windows percentile matrix per cell for methyltree

        methyltreeoutput<-makeWindows(dcis,
                    type = "CG", 
                    metric = "percent", 
                    bed = as.data.frame(methyltreewindows,col.names=NULL)[,1:3],
                    threads = threads, 
                    index = "chr_cg", 
                    nmin = 1) 
      print(paste("Mean percent cells covered per window:",
            mean((rowSums(!is.na(methyltreeoutput))/ncol(methyltreeoutput))*100)))
      print("Filtering to windows with >10% of cells with coverage")
      methyltreeoutput<-methyltreeoutput[(rowSums(!is.na(methyltreeoutput))/ncol(methyltreeoutput)*100)>=10,]
      methyltreewindows<-data.frame(do.call("rbind",strsplit(row.names(methyltreeoutput),"_")))
      colnames(methyltreewindows)<-c("chr","start","end")
      methyltreewindows<-GenomicRanges::makeGRangesFromDataFrame(methyltreewindows)
      print(paste("Final Filtered window count:",as.character(nrow(as.data.frame((methyltreewindows))))))
      print(paste("Final Filtered window average width:",as.character(mean(width(methyltreewindows)))))
      print(paste("Final Total genome covered:",as.character(sum(width(methyltreewindows))/1000000),"Mbp"))
      methyltreeoutput<-makeWindows(dcis,
                  type = "CG", 
                  metric = "percent", 
                 bed = as.data.frame(methyltreewindows,col.names=NULL)[,1:3],
                  threads = threads, 
                  index = "chr_cg", 
                  nmin = 1) 

      methyltreeoutput$genomic_region_id<-row.names(methyltreeoutput)
      methyltreeoutput <- methyltreeoutput |> 
          pivot_longer(
          cols = !genomic_region_id, 
          names_to = "cell_id",
          values_to = "value",
          values_drop_na = TRUE)
    #make a metadata sheet with cluster info
    out_metadata<-dcis@metadata[,c("pass","cluster_id","cg_cov","mcg_pct","subclones")]
    colnames(out_metadata)<-c("HQ","celltype","nCG","met_rate","large_clone_id") #match names
    out_metadata$sample<-row.names(out_metadata) #sample (cell) names
    out_metadata$met_rate<-out_metadata$met_rate/100 #percentage to rate

    if(file.exists(paste0(prefix,"_methyltree_input.h5"))){
        system(paste0("rm -rf ",prefix,"_methyltree_input.h5"))
    }
      h5createFile(file=paste0(prefix,"_methyltree_input.h5"))
      h5write(methyltreeoutput,file=paste0(prefix,"_methyltree_input.h5"),name="data")
      h5write(out_metadata,file=paste0(prefix,"_methyltree_input.h5"),name="metadata")
    }



##############CHROMVAR BLOCK################

#chromvar preparation for all cells
chromvar_met_per_cell<-function(obj,stepsize=500,threads=50,percent_cell_coverage=2.5){
  #run window scoring for all cells
  chromvar_windows <- makeWindows(obj,
                                  stepsize = stepsize, 
                                  type = "CG", 
                                  metric = "score", 
                                  threads = threads, 
                                  index = "chr_cg", 
                                  nmin = 2) 
  #binarized based on hypomethylation score (1 for hypo, 0 for hyper)
  counts<-ifelse(chromvar_windows<0.25,1,0) #score ranges from -1 to 1, using 0.25 for cutoff
  #require 2.5% cell coverage for windows, this is kinda on par with ATAC data, kinda an arbitrary cutoff. but mostly to decrease computational time on which windows we scan for motifs (5% is 10k window, 2% is 909k windows)
  counts<-counts[rowSums(!is.na(counts))>=(ncol(counts)/100)*percent_cell_coverage,] 
  counts[is.na(counts)]<-0
  #remove purely hypermethylated windows
  counts<-counts[rowSums(counts)>1,] 
  #dim(counts)
  return(counts)}

chromvar_met_per_cluster<-function(obj,stepsize=500,threads=50,mincov=30,percent_cell_coverage=50,groupBy="cluster_id"){
  #run window scoring for all cells
  clusterwindows <- calcSmoothedWindows(obj, 
              type = "CG", 
              threads = threads,
              step = 500,
              smooth = 1,
              index = "chr_cg",
              groupBy = groupBy, 
              returnSumMatrix = TRUE, 
              returnPctMatrix = TRUE)
  cov_mat<-clusterwindows[["sum_matrix"]]
  cov_mat_col<-colnames(cov_mat)[which(grepl(colnames(cov_mat),pattern="_t$"))]
  cov_mat<-as.data.frame(cov_mat)[,colnames(cov_mat) %in% cov_mat_col]
  pct_mat<-as.data.frame(clusterwindows[["pct_matrix"]][,4:ncol(clusterwindows[["pct_matrix"]])])
  row.names(pct_mat)<-paste(clusterwindows[["pct_matrix"]]$chr,clusterwindows[["pct_matrix"]]$start,clusterwindows[["pct_matrix"]]$end,sep="_")
  pct_mat[which(cov_mat<min_cov,arr.in=T)]<-NA
  #binarized based on hypomethylation score (1 for hypo, 0 for hyper)
  counts<-ifelse(pct_mat<0.25,1,0) #score ranges from -1 to 1, using 0.25 for cutoff
  #require 50% cell coverage for windows
  idx<-which(rowSums(!is.na(counts))>=(ncol(counts)/100)*percent_cell_coverage)
  counts<-counts[idx,] 
  counts[is.na(counts)]<-0
  #remove purely hypermethylated windows
  counts<-counts[rowSums(counts)>1,] 
  #dim(counts)
  return(counts)}

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


###############################################################################################################
########################################HYPOMETHYLATION PLOT###################################################
###############################################################################################################


col_celltype=c(
"Tcell_1"="#00ffff",
"Tcell_2"="#3399cc",
"Tcell_3"="#0099ff",
"Myeloid"="#66cc99",
"Endo"="#ffcc00",
"Fibro"="#cc3366",
"LumHRCancer_1"="#ff3399",
"LumHRCancer_2"="#993399",
"LumHRCancer_3"="#cc3399",
"LumHRCancer_4"="#ff6666",
"LumHRNormal"="#cc3399",
"LumSec_1"="#99cc99",
"LumSec_2"="#66cc66",
"Basal_1"="#993399",
"Basal_2"="#9966ff")



# #set colors for each facet and use greyscaling for higher met levels
# #https://stackoverflow.com/questions/33221794/separate-palettes-for-facets-in-ggplot-facet-grid
# #histograM <- function(obj,

hypoM<-function(obj=obj,genes=c("KIT"),prefix="hbca_dcis",matrix = "cg_celltype_tracks",colors=col_celltype,
                trackOverhang = 3000,arrowOverhang = 0,ncol = length(genes),
                invert_met = TRUE,legend = TRUE,removeNA = TRUE,
                width = 500,trackScale = 1.5,colorMax = 100,
                cgisland="/volumes/USR2/Ryan/projects/metact/ref/cpgIslandExt.bed"){
    #prepare cgi
    cgi<-rtracklayer::import(cgisland)
    cgi<-as.data.frame(cgi)
    colnames(cgi)<-c("chr","start","end","strand")
    #use named list for colors
    if (!is.null(colors)) {
    pal <- colors
    } else {
    pal <- rev(c("#FF0082", "#dbdbdb", "#cccccc", "#999999"))
    }
    #req ref
    if (is.null(obj@ref)) {
    stop("Please make sure a genome annotation file has been added to the obj@ref slot with makeRef().")
    }
    p <- vector("list", length(genes)) # empty plot list
    for (i in 1:length(genes)) { #make a plot per gene
        print(paste("Plotting gene:",genes[i]))
        ref <- obj@ref |> dplyr::filter(gene_name == genes[i])
        aggregated <- obj@genomeMatrices[[matrix]]

        toplot <- aggregated[c((aggregated$chr == ref$seqid[ref$type == "gene"] &
                                    aggregated$start > (ref$start[ref$type == "gene"] - trackOverhang) &
                                    aggregated$end < (ref$end[ref$type == "gene"] + trackOverhang))), ]
        cgi <- cgi[cgi$chr == as.character(toplot$chr[1]) &
                            cgi$start > min(toplot$start) &
                            cgi$end < max(toplot$end),]
                                
        ngroups <- ncol(toplot) - 3 #subtract chr start end
        trackHeight <- ngroups * trackScale
        toplot <- tidyr::pivot_longer(toplot, cols = c(4:ncol(toplot)), names_to = "group", values_to = "pct_m") |> 
        dplyr::rowwise() |> 
        dplyr::mutate(middle = mean(c(start, end), na.rm = TRUE))
        toplot <- toplot |> dplyr::filter(!is.na(pct_m))

        if (invert_met){
            toplot$pct_m<-toplot$pct_m-100
            colorMax=-100
        }
        #initialize first row
        p[[i]] <- vector("list", length(colors)+1) # empty plot list
        #plot group tracks
        xmin_pos<-NA
        xmax_pos<-NA
        for (j in names(colors)){
            toplot_sub<-toplot[toplot$group==j,]
            j_panel=which(j==names(colors))+1
            p[[i]][[j_panel]] <- ggplot2::ggplot() + 
                ggplot2::geom_col(data = toplot_sub, 
                ggplot2::aes(x = middle, y = pct_m, fill = pct_m), width = width) + 
                scale_fill_gradient2(high="black",mid="grey",low=colors[j],midpoint=-60) + 
                ylab(j) +
                theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_blank(),
                axis.text.y=element_text(size=3),
                axis.title.y=element_text(angle=90, vjust=1,size=4,color=colors[j]),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank())+
                scale_y_continuous(limits = c(-100,0),labels=rev(c("100", "75","50","25","0")))+
                theme(legend.position="none")
            xmin_pos=min(toplot_sub$middle)-width/2
            xmax_pos=max(toplot_sub$middle)-width/2
            }
        #first track is annot
        p[[i]][[1]] <-ggplot()+
        ggplot2::geom_rect(fill = "pink",alpha=0.8, data = ref |> dplyr::filter(type == "gene") |>
                            dplyr::mutate(promoter_start = ifelse(strand == "+", (start - 1500), (end + 1500)),
                            promoter_end = ifelse(strand == "+", (promoter_start+3000), (promoter_start-3000))),
                            ggplot2::aes(xmin = promoter_start, xmax = promoter_end, ymin = 0, ymax = trackHeight/4)) +
        ggplot2::geom_rect(alpha=0.8,fill = "black", data = ref |> dplyr::filter(type == "exon"),ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = trackHeight/4)) +
        ggplot2::geom_rect(alpha=0.8,fill = "green", data = cgi, ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = trackHeight/4)) +
        ggplot2::geom_segment(data = ref, color="gray",
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
                            title=element_text(size=4))+
                            theme(legend.position="none") +xlim(c(xmin_pos,xmax_pos))
    ggsave(wrap_plots(p[[i]])+plot_layout(ncol=1),file=paste(prefix,genes[i],"hypoM.pdf",sep="."),limitsize=F,height=10)
    }
}

#hypoM(obj=obj,genes=c("KIT","COL1A1","ESR1"),prefix="hbca_dcis.celltype")


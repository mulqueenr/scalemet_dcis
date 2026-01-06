#run in amethyst.sif
#singularity shell \
#--bind ~/projects/ \
#--bind /volumes/seq/projects/metACT \
#~/singularity/amethyst.sif

library(amethyst)
library(rhdf5)
library(data.table)
library(ggplot2)
library(patchwork)
library(tibble)
library(tidyr)
library(plyr); library(dplyr)
library(future)
library(furrr)
library(purrr)
library(cowplot)
library(pheatmap)
library(parallel)
library(rtracklayer)
library(gridExtra)
library(tidyverse)
library(colorspace)
library(GeneNMF) #new from here down
library(magrittr)
library(universalmotif)
library(ape)
library(ggtree)
library(BSgenome.Hsapiens.UCSC.hg38)
#library(motifmatchr)
#library(JASPAR2020)
#library(TFBSTools)
#library(parallel)
#library(chromVAR)
#library(SummarizedExperiment)
#library(ggseqlogo)
#library(patchwork)


#source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in

cluster_by_windows<-function(obj,
    window_name,
    metric="score",
    threads=100,
    neighbors=50,
    est_dim=10,
    k_pheno=15,
    dist=0.1,
    outname="clustering"){
  print(paste("Estimating dimensions..."))                                           
  #filter windows by cell coverage
  obj@genomeMatrices[[window_name]] <- obj@genomeMatrices[[window_name]][rowSums(!is.na(obj@genomeMatrices[[window_name]])) >= 45, ]
  est_dim<-dimEstimate(obj, genomeMatrices = c(window_name), dims = c(20), threshold = 0.95)
  print(est_dim)

  print("Running IRLBA reduction...")
  obj@reductions[[paste(window_name,"irlba",sep="_")]] <- runIrlba(obj, genomeMatrices = c(window_name), dims = est_dim, replaceNA = c(0))

  obj@reductions[[paste(window_name,"irlba_regressed",sep="_")]] <- regressCovBias(obj, reduction = paste(window_name,"irlba",sep="_")) 
  # Optional; helps reduce coverage bias
  print("Clustering on coverage regressed reduction...")

  #obj <- runCluster(obj, k_phenograph = 25, reduction = paste(window_name,"irlba_regressed",sep="_")) 
    obj <- runCluster(obj, k_phenograph = k_pheno, reduction = paste(window_name,"irlba_regressed",sep="_")) 
  # consider increasing k_phenograph to 50 for larger datasets

  obj <- runUmap(obj, neighbors = neighbors, dist = dist, method = "euclidean", reduction = paste(window_name,"irlba_regressed",sep="_")) 
  
  print("Plotting...")
  p1 <- dimFeature(obj, colorBy = cluster_id, reduction = "umap") + ggtitle(paste(window_name,"Clusters"))
  p2 <- dimFeature(obj, colorBy = sample, reduction = "umap") + ggtitle(paste(window_name,"Samples"))
  p3 <- dimFeature(obj, colorBy = log10(cov), pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red"),guide="colourbar") + ggtitle("Coverage distribution")
  p4 <- dimFeature(obj, colorBy = mcg_pct, pointSize = 1) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + ggtitle("Global %mCG distribution")
  plt<-plot_grid(p1, p2,p3, p4,ncol=2)
  ggsave(plt,file=paste0(outname,"_umap.pdf"),width=20,height=20)  

  return(obj)                                             
}

dmr_and_1kb_window_gen<-function(
    obj=dat,
    prefix="nakshatri_peaks",
    groupBy="cluster_id",
    threads=20,
    step=1000){
    cluster1kbwindows <- calcSmoothedWindows(obj, 
                                            type = "CG", 
                                            threads = threads,
                                            step = step,
                                            smooth = 3,
                                            index = "chr_cg",
                                            groupBy = groupBy, 
                                            returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                            returnPctMatrix = TRUE)
    obj@genomeMatrices[[paste0("cg_",groupBy,"_tracks")]] <- cluster1kbwindows[["pct_matrix"]]

    pal=c("#E5E6E4","#CFD2CD","#A6A2A2","#847577","#6E44FF")
    dmrs<-testDMR(cluster1kbwindows[["sum_matrix"]], eachVsAll = TRUE, nminTotal = 0, nminGroup = 0) 
    dmrs2<-filterDMR(dmrs, method = "bonferroni", filter = TRUE, pThreshold=0.05,logThreshold=0.5) #add additional columns direction column
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

################################################
################FINE CELLTYPING FUNCTIONS######
################################################
################################################

plot_histogram_page<-function(celltype){
    plt<-histograModified(obj, 
        baseline="mean",
        genes = unlist(cell_markers[celltype]),
        colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
        matrix = "cg_cluster_id_tracks", arrowScale = .03, trackScale = .5,
        legend = F, cgisland=cgi)
    return(plt)
}

plot_dmr_RNAmarker_overlap<-function(dmr,rna_markers,celltype,plot_value="Jaccard",clus_order){
    #plot_value can be one of pval, odds.ratio, intersection, union, Jaccard
    #flow through:
    # - filter DMRs to HYPO
    # - filter RNA markers to overexpressed
    # - table dmr to marker count, and plot odds ratio
    # - plot gene body methylation hypoM plots
    output_directory=paste0(project_data_directory,"/fine_celltyping/",celltype)

    #PLOT HYPO
    #filter to hypo DMRs (more likely to increase expression)
    dmr_hypo<- dmr %>% group_by(test) %>% filter(direction=="hypo") 
    rna_markers_up <- rna_markers %>% filter(avg_log2FC>0) #look only for increased expression markers
    rna_markers_up$gene_name<-row.names(rna_markers_up)
    dmr_gene_lists<-split(row.names(dmr_hypo),dmr_hypo$test)
    dmr_gene_lists<-lapply(dmr_gene_lists,function(i) {
            unlist(
                lapply(unlist(strsplit(dmr_hypo[unlist(dmr_gene_lists[i]),]$gene_names,",")),
                function(x) {gsub("[[:space:]]", "", x)})
            )
        })
    rna_markers_lists<-split(rna_markers_up$gene,rna_markers_up$cluster)
    gom.obj<-newGOM(dmr_gene_lists,
                rna_markers_lists,
                genome.size=length(unique(c(unlist(dmr_gene_lists),unlist(rna_markers_lists))))
    )
    or_mat<-getMatrix(gom.obj, plot_value)

    dmr_count<-dmr_hypo %>% 
                    group_by(test,.drop=FALSE) %>% 
                    tally()
    dmr_count<-setNames(dmr_count$n,nm=dmr_count$test)
    dmr_count=dmr_count[clus_order]

    marker_count<-rna_markers_up %>% 
                    group_by(cluster,.drop=FALSE) %>% 
                    tally() 
    marker_count<-setNames(marker_count$n,nm=marker_count$cluster)


    print("Plotting hypo overlap of DMRs...")
    pdf(paste0(output_directory,"/","04_scaledcis.",broad_celltype,".dmr_marker_overlap.heatmap.pdf"))
    print(
        ComplexHeatmap::Heatmap(or_mat,
            row_order=clus_order,cluster_rows=FALSE,cluster_columns=FALSE,
            bottom_annotation=HeatmapAnnotation(bar=anno_barplot(marker_count)),
            show_row_names = TRUE, show_column_names = TRUE, row_names_side="right",column_names_side="bottom",
            column_title = "Hypomethylated and Overexpressed",
            left_annotation=rowAnnotation(bar=anno_barplot(dmr_count)),
            name=plot_value,
            cell_fun = function(j, i, x, y, width, height, fill) {
                                grid.text(format(round(or_mat,digits=2),nsmall=2)[i, j], x, y, gp = gpar(fontsize = 10))
                                }
            ))
    #PLOT HYPER
    #filter to hypo DMRs (more likely to increase expression)
    dmr_hyper<- dmr %>% group_by(test) %>% filter(direction=="hyper") 
    rna_markers_down <- rna_markers %>% filter(avg_log2FC<0) #look only for increased expression markers
    rna_markers_down$gene_name<-row.names(rna_markers_down)
    dmr_gene_lists<-split(row.names(dmr_hyper),dmr_hyper$test)
    dmr_gene_lists<-lapply(dmr_gene_lists,function(i) {
            unlist(
                lapply(unlist(strsplit(dmr_hyper[unlist(dmr_gene_lists[i]),]$gene_names,",")),
                function(x) {gsub("[[:space:]]", "", x)})
            )
        })
    rna_markers_lists<-split(rna_markers_down$gene,rna_markers_down$cluster)
    gom.obj<-newGOM(dmr_gene_lists,
                rna_markers_lists,
                genome.size=length(unique(c(unlist(dmr_gene_lists),unlist(rna_markers_lists))))
    )
    or_mat<-getMatrix(gom.obj, plot_value)
    
    dmr_count<-dmr_hyper %>% 
                    group_by(test,.drop=FALSE) %>% 
                    tally()
    dmr_count<-setNames(dmr_count$n,nm=dmr_count$test)
    dmr_count=dmr_count[clus_order]

    marker_count<-rna_markers_down %>% 
                    group_by(cluster,.drop=FALSE) %>% 
                    tally() 
    marker_count<-setNames(marker_count$n,nm=marker_count$cluster)

    
    print(
        ComplexHeatmap::Heatmap(or_mat,
            bottom_annotation=HeatmapAnnotation(bar=anno_barplot(marker_count)),
            row_order=clus_order,cluster_rows=FALSE,cluster_columns=FALSE,
            show_row_names = TRUE, show_column_names = TRUE, row_names_side="right",column_names_side="bottom",
            column_title = "Hypermethylated and Underexpressed",
            left_annotation=rowAnnotation(bar=anno_barplot(dmr_count)),
            name=plot_value,
            cell_fun = function(j, i, x, y, width, height, fill) {
                                grid.text(format(round(or_mat,digits=2),nsmall=2)[i, j], x, y, gp = gpar(fontsize = 10))
                                }
            ))
    
    dev.off()
}


plot_dmr_RNAMarker_hypoM<-function(dat=dat_fine,dmr,order=NULL,rna_markers,broad_celltype,gtf_file="/container_ref/gencode.v43.annotation.gtf.gz"){

    rna_markers$gene_name<-row.names(rna_markers)
    #subset markers by size for plotting
    #use gtf file to set gene sizes

    met_clusters<-length(unique(dmr$group))

    print("Reading in GTF file to append to markers...")
    gtf <- rtracklayer::readGFF(gtf_file)
        gtf<- gtf %>% 
        filter(type=="gene" & gene_type %in% c("lncRNA","protein_coding")) %>% 
        filter(gene_name %in% rna_markers$gene_name) %>% 
        mutate(gene_length=abs(end-start))

    print("Merging GTF data and marker gene information...")
    rna_markers_merge<-merge(rna_markers,gtf,by="gene_name")
    #overlap hypo DMR and RNA markers per fine celltype
    #limit markers to those that will plot well for visual inspection (i.e. correct size)
    #plot celltype markers for fine celltyping

    print("Filtering marker genes to those that are between 5 and 25kbp")
    rna_markers_merge<-rna_markers_merge %>% 
    filter(gene_length>5000) %>% 
    filter(gene_length<50000)

    print(paste(nrow(rna_markers_merge),"RNA markers remaining..."))
    #so now we have a list of RNA markers that are set to size thresholds
    #pick based on significant DMRs (across all clusters) per fine_celltype rna marker

    dmr<-dmr %>% filter(direction=="hypo" & dmr_logFC< -1 & dmr_length > 2000)
    print(paste(nrow(dmr),"DMR sites that are hypomethylated, log2FC >2, and >2kbp"))
    #dmr gene names is a nested list (dmr associated with multiple genes, so unlist them) and remove white space
    dmr_gene_names<-unlist(
        lapply(
            unlist(
                strsplit(
                    unlist(dmr$gene_names),",")),
                function(x) gsub("[[:space:]]", "", x)
            )
        )
    genes_to_plot<-rna_markers_merge %>% 
                    group_by(cluster) %>% 
                    filter(gene %in% dmr_gene_names) %>% 
                    slice_min(n=10,order=p_val_adj) %>%
                    select(gene_name)
    print(paste("Genes passing all filters for plotting:",nrow(genes_to_plot)))

    #make a named list of fine cell types and genes for plotting
    plotting_list<-setNames(nm=unique(genes_to_plot$cluster),
        lapply(unique(genes_to_plot$cluster), function(x){
            return(genes_to_plot[genes_to_plot$cluster==x,]$gene_name)
        }) )
    
    print(paste("Making genome plots..."))
    plt_list<-lapply(names(plotting_list),
        function(celltype){
                genes<-unlist(plotting_list[celltype])
                genes<-genes[genes %in% dat_fine@ref$gene_name]
                print(paste("Plotting:",celltype))
                print(paste("Genes to plot:",genes))
                plt<-histograModified(obj=dat_fine, 
                    baseline="mean",
                    genes = unlist(plotting_list[celltype]),
                    order=order,
                    matrix = paste0("cg_",broad_celltype,".clusters_tracks"), arrowScale = .03, trackScale = .5,
                    legend = F, cgisland=cgi) + ggtitle(celltype)
                ggsave(plt,
                file=paste(broad_celltype,"clusters",celltype,"fine_celltype.markers.pdf",sep="."),
                width=length(genes)*5,
                height=length(met_clusters)*20,
                limitsize=F)
            })

}


calc_overlap_stats<-function(i){
    i_name<-names(dmr_granges)[i]
    out<-do.call("rbind",
        lapply(1:length(atac_markers_granges),
            function(j){
                j_name <- names(atac_markers_granges)[j]
                sizeA <- length(dmr_granges[[i]])
                sizeB <- length(atac_markers_granges[[j]])
                intersection <- GenomicRanges::countOverlaps(dmr_granges[i],atac_markers_granges[[j]])
                union <- sizeA+sizeB - intersection
                cont.tbl <- matrix(c(union, 
                                            sizeA - intersection, 
                                            sizeB - intersection, 
                                            intersection), 
                                            ncol=2)
                rownames(cont.tbl) <- c('notB', 'inB')
                colnames(cont.tbl) <- c('notA', 'inA')

                # Perform Fisher's exact test.
                res.fisher <- try(fisher.test(cont.tbl, alternative='greater'), silent=TRUE)

                if(is.list(res.fisher)) {
                    odds.ratio <- setNames(res.fisher$estimate, NULL)
                    pval <- res.fisher$p.value
                    } else {
                    odds.ratio <- .0
                    pval <- 1.
                    }
                # Calculate Jaccard index. how many DMRs overlap out of all DMRs
                Jaccard_val <- ifelse(union == 0, 0, 
                                            intersection / 
                                                sizeA
                )

                return(c(i_name,j_name,as.numeric(odds.ratio),as.numeric(pval),as.numeric(Jaccard_val)))
            }))
    return(out)
}

################################################
################################################
################################################
################################################

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

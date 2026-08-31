
Plotting differential methylation results
Running DMR comparisons across different sets.

1. Cell type level comparison

2. Cell type level comparison (HBCA only)

3. Cell type x group comparison per cell type

Load libraries

```R
library(Seurat)
library(SeuratExtend)

library(amethyst)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(Matrix)

library(GeneOverlap)
library(fgsea)
library(msigdbr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(annotatr)
library(GenomeInfoDb)

library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(cowplot)
library(patchwork)

library(parallel)
options(future.globals.maxSize= 500000*1024^2) #80gb limit for parallelizing
options(scipen = 0)
set.seed(111)
```

Set up environment
```R
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"

#read in object from directory
task_cpus=300
processing_folder="04_dmr"
wd=paste(sep="/",project_data_directory,processing_folder)
system(paste0("mkdir -p ",wd))
setwd(wd)

#set up output directory
output_directory=paste0(wd,"/","dmr_celltype_group")
system(paste0("mkdir -p ",output_directory))

#read RNA
rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")
rna<-subset(rna,cells=row.names(rna@meta.data)[!(rna$coarse_celltype %in% c("suspected_doublet"))])
Idents(rna)<-rna$coarse_celltype

#read methylation
obj<-readRDS(file=paste(project_data_directory,"03_fine_celltyping","03_scaledcis.final_celltypes.amethyst.rds",sep="/"))

#read GTF for RNA annotation
gtf<-import("/data/rmulqueen/projects/scalebio_dcis/ref/gencode.v43.annotation.gtf.gz")
gtf$gene<-gtf$gene_name
```

Functions
```R

#gsea of top DMR genes (hypomet)
gsea_enrichment<-function(dmrs,species="human",
                          category="C3",
                          subcategory="TFT:GTRD",
                          out_setname="TFT",
                          prefix=prefix,
                          obj=obj,
                          sample_name=sample_name,
                          output_gsea_directory){

  print(paste("Calculating DMR overlap with:",out_setname))

  pathwaysDF <- msigdbr(species=species, 
                        collection=category, 
                        subcollection = subcategory)

  #limit pathways to genes in our data
  pathwaysDF <- pathwaysDF[pathwaysDF$ensembl_gene %in% unlist(lapply(strsplit(unique(obj@ref$gene_id),"[.]"),"[",1)),]
  
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
  saveRDS(fgseaRes,file=paste0(output_gsea_directory,"/",prefix,".GSEA_enrichment.",sample_name,".",out_setname,".rds"))
}

gsea_across_sets<-function(obj, 
                    collapsed_dmrs,
                    sample_name, 
                    prefix,
                    output_gsea_directory){
  print(paste("Loading DMRS for sample:",sample_name))

  #run gsea enrichment on different sets
  print(paste("Calculating TF Binding Enrichment"))
  tft_gsea<-gsea_enrichment(species="human",
                          category="C3",
                          subcategory="TFT:GTRD",
                          out_setname="TFT",
                          prefix=prefix,
                          sample_name=sample_name,
                          dmrs=collapsed_dmrs,obj=obj,output_gsea_directory=output_gsea_directory)


  print(paste("Calculating Position Enrichment"))
  position_gsea<-gsea_enrichment(species="human",
              category="C1",
              subcategory=NULL,
              out_setname="position",
              prefix=prefix,
              sample_name=sample_name,
              dmrs=collapsed_dmrs,obj=obj,output_gsea_directory=output_gsea_directory)

  print(paste("Calculating Hallmark Enrichment"))
  hallmark_gsea<-gsea_enrichment(species="human",
              category="H",
              subcategory=NULL,
              out_setname="hallmark",              
              prefix=prefix,
              sample_name=sample_name,
              dmrs=collapsed_dmrs,obj=obj,output_gsea_directory=output_gsea_directory)

  print(paste("Calculating Cancer Cell Atlas Enrichment"))
  cancercellatlas_gsea<-gsea_enrichment(species="human",
              category="C4",
              subcategory="3CA",
              out_setname="3CA",
              prefix=prefix,
              sample_name=sample_name,
              dmrs=collapsed_dmrs,obj=obj,output_gsea_directory=output_gsea_directory)

  print(paste("Finished sample:",sample_name))
}

plot_gsea<-function(gsea=gsea,
                    out_setname="hallmark",
                    output_directory,
                    prefix=prefix,
                    order_by=names(celltype_col)){

  #get top hits per group
  features_to_keep_hypo <- gsea %>% filter(padj<0.05) %>% group_by(test) %>% slice_min(NES,n=10) %>% dplyr::select(names)
  features_to_keep_hyper <- gsea %>% filter(padj<0.05) %>% group_by(test) %>% slice_max(NES,n=10) %>% dplyr::select(names)
  features_to_keep<-unique(unlist(features_to_keep_hypo$names,features_to_keep_hyper$names))

  if(length(features_to_keep)>1){
    if(out_setname!="position"){
      #cluster by features
      gsea_nes <- gsea %>% 
                  filter(names %in% features_to_keep) %>% 
                  tidyr::pivot_wider(names_from=test, id_cols=names, values_from=NES) %>% 
                  as.data.frame()
      row.names(gsea_nes)<-gsea_nes$names
      gsea_nes<-gsea_nes[,2:ncol(gsea_nes)]
      clus<-dist(gsea_nes, method = "euclidean") %>%  
            hclust( method = "complete")
      row_order <- row.names(gsea_nes)[clus$order ]
      gsea_plot<-gsea[gsea$names %in% row_order,]
      gsea_plot$names<-factor(gsea_plot$names,levels=row_order)
      gsea_plot$celltype<-unlist(lapply(strsplit(gsea_plot$test,"_"),"[",1))

      gsea_plot$order_by<-factor(gsea_plot$test,levels=order_by)

      # Simple dot plot of GSEA results
      plt<-ggplot(gsea_plot, aes(x = order_by, y = names)) +
        geom_point(aes(size = -log10(padj), color = NES)) +
        scale_radius(range = c(0, 10), name = "-log10(p-value)") +
        scale_color_gradient2(low = "#b84d9c", mid = "white", high = "#666666", midpoint = 0, name = "NES") +
        theme_bw() + facet_wrap (~celltype,scales = "free_x",nrow=1) + 
        theme_minimal() + theme(
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.minor = element_blank()
        ) 

        ggsave(plt,file=paste0(output_directory,"/",prefix,".",out_setname,".NES.heatmap.pdf"),width=50,height=10,limitsize=FALSE)
        

  } else { #sort by position in genome for visualization sake
    gsea_nes <- gsea %>% filter(names %in% features_to_keep) %>% tidyr::pivot_wider(names_from=test, id_cols=names, values_from=NES)  %>% as.data.frame()
    row.names(gsea_nes)<-gsea_nes$names

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

    gsea<-as.data.frame(gsea)

    gsea_plot<-gsea[gsea$names %in% paste0("chr",row.names(row_order)),]


    #gsea_plot<-tidyr::pivot_longer(gsea_plot,cols=colnames(gsea_plot)[endsWith(colnames(gsea_plot),suffix="_v_rest")],names_to="celltype")
    gsea_plot$names<-factor(gsea_plot$names,levels=paste0("chr",row.names(row_order)))
    gsea_plot$celltype<-unlist(lapply(strsplit(gsea_plot$test,"_"),"[",1))

    gsea_plot$order_by<-factor(gsea_plot$test,levels=order_by)
    
    # Simple dot plot of GSEA results
    plt<-ggplot(gsea_plot, aes(x = order_by, y = names)) +
      geom_point(aes(size = -log10(padj), color = NES)) +
      scale_radius(range = c(0, 10), name = "-log10(p-value)") +
      scale_color_gradient2(low = "#b84d9c", mid = "white", high = "#666666", midpoint = 0, name = "NES") +
      theme_bw() + facet_wrap (~celltype,scales = "free_y",ncol=1) + 
      theme_minimal() + theme(
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.minor = element_blank()
        ) +coord_flip()

    ggsave(plt,file=paste0(output_directory,"/",prefix,".",out_setname,".NES.heatmap.pdf"),width=50,height=50,limitsize=FALSE)

    
    } 
  } else { print(paste("No significant findings in",out_setname))}
}


calculate_rna_pathway<-function(rna,
  out_name="hallmark",
  cores=50,
  genesets=SeuratExtendData::Genesets_data$human$GSEA[["hallmark gene sets"]],
  output_directory=output_plot_directory,
  prefix=prefix){

  #names(SeuratExtendData::Genesets_data$human$GSEA)
  print(paste("Calculating",out_name,"Enrichment"))
  rna <- GeneSetAnalysis(subset(rna,downsample=100), 
                          title=out_name,
                          nCores=cores, 
                          genesets = genesets)

  matr <- rna@misc[["AUCell"]][[out_name]]
  saveRDS(matr,file=paste0(output_directory,"/",prefix,".",out_name,".AUC.gsea.rds"))
  celltype_matr<-matr %>% 
                t() %>% 
                as.data.frame() %>% 
                dplyr::group_by(rna$coarse_celltype[colnames(matr)]) %>% 
                summarise(across(everything(), mean)) %>% 
                as.data.frame()

  celltype_matr_sd<-matr %>% 
                t() %>% 
                as.data.frame() %>% 
                dplyr::group_by(rna$coarse_celltype[colnames(matr)]) %>% 
                summarise(across(everything(), sd)) %>% 
                as.data.frame()
  row.names(celltype_matr_sd)<-celltype_matr_sd[,1]
  row.names(celltype_matr)<-celltype_matr[,1]
  celltype_matr<-as.data.frame(t(celltype_matr[,2:ncol(celltype_matr)]))
  celltype_matr<-data.frame(celltype_matr)

  if(out_name!="position"){
  #clustering for order
    clus<-dist(celltype_matr, method = "euclidean") %>%  
          hclust( method = "complete")
    row_order <- row.names(celltype_matr)[clus$order ]
    #gsea_plot<-tidyr::pivot_longer(gsea_plot,cols=colnames(gsea_plot)[endsWith(colnames(gsea_plot),suffix="_v_rest")],names_to="celltype")
    celltype_matr_plot<-celltype_matr %>% 
                        tibble::rownames_to_column(var = "pathway") %>% 
                        tidyr::pivot_longer(cols = -pathway,
                                          names_to = "celltype",
                                          values_to = "auc_score")

    #add sd value in
  celltype_matr_plot$sd<-unlist(lapply(1:nrow(celltype_matr_plot),function(x){
    celltype_matr_sd[celltype_matr_plot[x,]$celltype,celltype_matr_plot[x,]$pathway]
  }))
  celltype_matr_plot$pathway<-factor(celltype_matr_plot$pathway,levels=row_order)
  celltype_matr_plot$order_by<-factor(celltype_matr_plot$celltype,levels=order_by)

  # Simple dot plot of GSEA results
  plt<-ggplot(celltype_matr_plot, aes(x = order_by, y = pathway)) +
  geom_point(aes(size = sd, color = auc_score)) +
  scale_radius(range = c(0, 10), name = "AUC Score") +
  scale_color_gradientn(colors=c("white","grey","orange"), name = "NES",breaks=c(0,0.5,1)) +
  theme_bw() + facet_wrap (~celltype,scales = "free_x",nrow=1) + 
  theme_minimal() + theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.minor = element_blank()
    )

  ggsave(plt,file=paste0(output_directory,"/",prefix,".",out_name,".AUC.heatmap.pdf"),width=50,height=50,limitsize=FALSE)
  } else {

  #for plotting position in order
  row.names(celltype_matr)<-gsub(row.names(celltype_matr),pattern="chr",replacement="")
  row_order_chr<-unlist(lapply(strsplit(row.names(celltype_matr),split="p|q"),"[",1))
  row_order_arm<-gsub("\\d", "", row.names(celltype_matr))
  row_order_band<-unlist(lapply(strsplit(row.names(celltype_matr),split="p|q"),"[",2))
  row_order<-data.frame(chr=row_order_chr,arm=row_order_arm,band=row_order_band)
  row_order$chr<-factor(row_order$chr,c(1:22,"X"))
  row_order$arm<-factor(row_order$arm,c("p","q","Xp","Xq"))
  row.names(row_order)<-row.names(celltype_matr)
  row_order<-row_order %>% arrange(chr,arm,band) 
  celltype_matr<-as.data.frame(celltype_matr)
  celltype_matr_plot<-celltype_matr[row.names(celltype_matr) %in% row.names(row_order),]

  #gsea_plot<-tidyr::pivot_longer(gsea_plot,cols=colnames(gsea_plot)[endsWith(colnames(gsea_plot),suffix="_v_rest")],names_to="celltype")
  celltype_matr_plot<-celltype_matr_plot %>% 
                tibble::rownames_to_column(var = "pathway") %>% 
                tidyr::pivot_longer(cols = -pathway,
                                  names_to = "celltype",
                                  values_to = "auc_score")
  celltype_matr_plot$sd<-as.numeric(unlist(lapply(1:nrow(celltype_matr_plot),function(x){
    celltype_matr_sd[celltype_matr_plot[x,]$celltype,celltype_matr_plot[x,]$pathway]
  })))
  celltype_matr_plot$pathway<-factor(celltype_matr_plot$pathway,levels=row.names(row_order))
  celltype_matr_plot$order_by<-factor(celltype_matr_plot$celltype,levels=order_by)

  # Simple dot plot of GSEA results
  plt<-ggplot(celltype_matr_plot, aes(x = order_by, y = pathway)) +
    geom_point(aes(size = sd, color = auc_score)) +
    scale_radius(range = c(0, 10), name = "AUC Score") +
    scale_color_gradientn(colors=c("white","grey","orange"), name = "NES",breaks=c(0,0.5,1)) +
    theme_bw() + facet_wrap (~celltype,scales = "free_y",ncol=1) + 
    theme_minimal() + theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.minor = element_blank()
      ) +coord_flip()

  ggsave(plt,file=paste0(output_directory,"/",prefix,".",out_name,".AUC.heatmap.pdf"),width=50,height=50,limitsize=FALSE)
  }
}


plot_marker_dmr_overlap<-function(markers=markers,
                                  dmrs=dmrs,
                                  prefix="celltype_v_rest",
                                  output_plots_folder=output_plots_folder,
                                  order_by=names(celltype_col)){

  #get overlap of genes and dmrs per celltype (need to add genomic ranges to RNA markers)
  print("Using GTF file to annotation gene locations...")
  gtf_markers<- gtf %>% 
    as.data.frame() %>% 
    filter(type=="gene") %>% 
    filter(gene_name %in% markers$gene) 

  ensembl_vector <- mapIds(
    x = org.Hs.eg.db, 
    keys = markers$gene, 
    column = "ENSEMBL", 
    keytype = "SYMBOL", 
    multiVals = "first" )
  
  markers$gene_id<-ensembl_vector[markers$gene]
  markers$direction_de<-ifelse(markers$avg_log2FC>0,"upregulated","downregulated")
  markers[is.na(markers$gene_id) & startsWith(markers$gene,prefix="ENSG00"),]$gene_id<-markers[is.na(markers$gene_id) & startsWith(markers$gene,prefix="ENSG00"),]$gene
  markers[is.na(markers$gene_id),] #some antisense genes are still lost
  gtf_markers$gene_id<-unlist(lapply(strsplit(gtf_markers$gene_id,"[.]"),'[',1))
  markers2<-left_join(markers,gtf_markers,"gene_id",relationship="many-to-many")
  markers2<-markers2[!is.na(markers2$seqnames),]

  print("Converting DMRs and DE output to GenomicRanges...")

  markers_grange<-GRanges(markers2)
  dmr_grange<-GRanges(dmr)

  #get count of dmr overlaps per marker genes
  #add 5kb upstream of promoter
  print("Running DMR and DE overlap...")

  markers_grange_overlap<-list()
  dmr_grange_overlap<-list()

  for(comp in as.character(unique(markers$comparison_name))){
      print(paste("Preparing overlap of...",comp))
      #prepare RNA differences
      markers_grange_comp<-markers_grange %>% 
                          as.data.frame() %>% 
                          filter(comparison_name==comp) %>% 
                          GRanges()
      promoter_ranges <- promoters(markers_grange_comp, upstream = 5000, downstream = 0)
      
      #prepare dmr
      dmr_grange_comp<-dmr_grange %>% 
                        as.data.frame() %>% 
                        filter(dmr_grange$comparison_name==comp) %>% 
                        GRanges()

      dmr_grange_comp$dmr_name<-paste(seqnames(dmr_grange_comp),
                                      start(dmr_grange_comp),
                                      end(dmr_grange_comp),sep="_")

      #dmr_grange_comp <- split(dmr_grange_comp,~direction)
      
      overlaps_genebody <- findOverlaps(query = dmr_grange_comp, 
        subject = markers_grange_comp)
  
      overlaps_promoter <- findOverlaps(query = dmr_grange_comp, 
        subject = promoter_ranges)

      # and add collapsed list of genes markers to DMRs (promoter)
      dmr_grange_comp$comparison_de_promoter_overlap<-NA
      dmr_grange_comp[overlaps_promoter@from,]$comparison_de_promoter_overlap<-markers_grange_comp[overlaps_promoter@to,]$gene.x
      
      #add collapsed list of DMRs as metadata column into markers_grange_celltype
      markers_grange_comp$comparison_dmr_hypo_promoter_overlap<-NA
      markers_grange_comp$comparison_dmr_hyper_promoter_overlap<-NA

      if(length(overlaps_promoter)>1){
        hypo_overlap<-dmr_grange_comp[dmr_grange_comp$direction=="hypo" & dmr_grange_comp$comparison_de_promoter_overlap %in% markers_grange_comp$gene.x,]$comparison_de_promoter_overlap
        markers_grange_comp[markers_grange_comp$gene.x %in% hypo_overlap,]$comparison_dmr_hypo_promoter_overlap<-1

        hyper_overlap<-dmr_grange_comp[dmr_grange_comp$direction=="hyper" & dmr_grange_comp$comparison_de_promoter_overlap %in% markers_grange_comp$gene.x,]$comparison_de_promoter_overlap
        markers_grange_comp[markers_grange_comp$gene.x %in% hyper_overlap,]$comparison_dmr_hyper_promoter_overlap<-1
      }

      # and add collapsed list of genes markers to DMRs  (genebody)
      dmr_grange_comp$comparison_de_genebody_overlap<-NA
      dmr_grange_comp[overlaps_genebody@from,]$comparison_de_genebody_overlap<-markers_grange_comp[overlaps_genebody@to,]$gene.x
      
      #add collapsed list of DMRs as metadata column into markers_grange_celltype
      markers_grange_comp$comparison_dmr_hypo_genebody_overlap<-NA
      markers_grange_comp$comparison_dmr_hyper_genebody_overlap<-NA

      if(length(overlaps_genebody)>1){
        hypo_overlap<-dmr_grange_comp[dmr_grange_comp$direction=="hypo" & dmr_grange_comp$comparison_de_genebody_overlap %in% markers_grange_comp$gene.x,]$comparison_de_genebody_overlap
        markers_grange_comp[markers_grange_comp$gene.x %in% hypo_overlap,]$comparison_dmr_hypo_genebody_overlap<-1

        hyper_overlap<-dmr_grange_comp[dmr_grange_comp$direction=="hyper" & dmr_grange_comp$comparison_de_genebody_overlap %in% markers_grange_comp$gene.x,]$comparison_de_genebody_overlap
        markers_grange_comp[markers_grange_comp$gene.x %in% hyper_overlap,]$comparison_dmr_hyper_genebody_overlap<-1
      }

      markers_grange_overlap<-list(markers_grange_overlap,markers_grange_comp)
      dmr_grange_overlap<-list(dmr_grange_overlap,dmr_grange_comp)

    }

    dmr_grange_overlap<-do.call("c",unlist(dmr_grange_overlap))
    markers_grange_overlap<-do.call("c",unlist(markers_grange_overlap))
    
    #sanity check
    dmr_grange_overlap %>% as.data.frame() %>% 
    group_by(direction,comparison_name) %>% 
    summarize(total=n(),overlap=sum(!is.na(comparison_de_promoter_overlap)))

    markers_grange_overlap %>% as.data.frame() %>% 
    group_by(comparison_name) %>% 
    summarize(total=n(),overlap=sum(!is.na(comparison_dmr_hypo_genebody_overlap)))

    print("Saving DMR and DE overlaps...")

    saveRDS(dmr_grange_overlap,paste0(output_plots_folder,"/","04_DMR_with_DE_overlaps.celltype_v_rest.rds"))
    saveRDS(markers_grange_overlap,paste0(output_plots_folder,"/","04_DE_with_DMR_overlaps.celltype_v_rest.rds"))

    marker_overlap_counts<-markers_grange_overlap %>% 
      as.data.frame() %>% 
      group_by(comparison_name,direction_de) %>% 
      summarize(
        total_gene_count=n(),
        hypo_genebody_overlap=sum(!is.na(comparison_dmr_hypo_genebody_overlap)),
        hyper_genebody_overlap=sum(!is.na(comparison_dmr_hyper_genebody_overlap)),
        hypo_promoter_overlap=sum(!is.na(comparison_dmr_hypo_promoter_overlap)),
        hyper_promoter_overlap=sum(!is.na(comparison_dmr_hyper_promoter_overlap)),
        hypo_promoter_or_genebody_overlap=sum(!is.na(comparison_dmr_hypo_promoter_overlap) | !is.na(comparison_dmr_hypo_genebody_overlap))) %>%
        as.data.frame()

    dmr_overlap_counts<-dmr_grange_overlap %>% 
      as.data.frame() %>% 
      group_by(comparison_name,direction) %>% 
      summarize(
        total_dmr_count=n(),
        genebody_overlap=sum(!is.na(comparison_de_genebody_overlap)),
        promoter_overlap=sum(!is.na(comparison_de_promoter_overlap)),
        genebody_or_promoter_overlap=sum(!is.na(comparison_de_promoter_overlap) | !is.na(comparison_de_genebody_overlap))) %>%
        as.data.frame()

    saveRDS(marker_overlap_counts,paste0(output_plots_folder,"/","04_DE_genes_with_DMR_overlaps.",comparison_set,".summary.counts.rds"))
    saveRDS(dmr_overlap_counts,paste0(output_plots_folder,"/","04_DMRs_with_DE_overlaps.",comparison_set,".summary.counts.rds"))

    marker_overlap_counts$celltype<-unlist(lapply(strsplit(as.character(marker_overlap_counts$comparison_name),"_"),"[",1))
    marker_overlap_counts$order_by<-factor(marker_overlap_counts$comparison_name,levels=order_by)

    dmr_overlap_counts$celltype<-unlist(lapply(strsplit(as.character(dmr_overlap_counts$comparison_name),"_"),"[",1))
    dmr_overlap_counts$order_by<-factor(dmr_overlap_counts$comparison_name,levels=order_by)

    print("Making Summary plots of overlaps...")
    #promoter hypo overlap
    overlap_counts_promoter <-marker_overlap_counts %>% 
      filter(direction_de=="upregulated") %>% 
      dplyr::select(comparison_name,celltype,total_gene_count,hypo_promoter_overlap)
    overlap_counts_promoter$hypo_promoter_nonoverlap <-overlap_counts_promoter$total_gene_count - overlap_counts_promoter$hypo_promoter_overlap
    overlap_counts_promoter <-tidyr::pivot_longer(overlap_counts_promoter,cols=c("hypo_promoter_nonoverlap","hypo_promoter_overlap"))
    overlap_counts_promoter$percent<-unlist(lapply(1:nrow(overlap_counts_promoter),function(x) 
      ifelse(overlap_counts_promoter[x,]$name=="hypo_promoter_overlap",
            round((overlap_counts_promoter[x,]$value/overlap_counts_promoter[x,]$total_gene_count)*100,2),
            NA)))
    overlap_counts_promoter$comparison_name<-factor(overlap_counts_promoter$comparison_name,levels=order_by)

    plt1<-ggplot(overlap_counts_promoter,aes(x=comparison_name,y=value,fill=name))+
          geom_col()+ geom_text(aes(label = percent), vjust = -0.5, size = 5) +
          facet_wrap(~comparison_name,nrow=1,scales="free_x")+
          ylab("Promoter hypo DMR overlap for DE upregulated genes")+  
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    #genebody hypo overlap

    overlap_counts_genebody<-overlap_counts %>% 
      filter(direction_de=="upregulated") %>% 
      dplyr::select(comparison_name,celltype,total_gene_count,hypo_genebody_overlap)
      
    overlap_counts_genebody$hypo_genebody_nonoverlap<-overlap_counts_genebody$total_gene_count - overlap_counts_genebody$hypo_genebody_overlap
    overlap_counts_genebody<-tidyr::pivot_longer(overlap_counts_genebody,cols=c("hypo_genebody_nonoverlap","hypo_genebody_overlap"))
    overlap_counts_genebody$percent<-unlist(lapply(1:nrow(overlap_counts_genebody),function(x) 
      ifelse(overlap_counts_genebody[x,]$name=="hypo_genebody_overlap",
            round((overlap_counts_genebody[x,]$value/overlap_counts_genebody[x,]$total_gene_count)*100,2),
            NA)))

    overlap_counts_genebody$comparison_name<-factor(overlap_counts_genebody$comparison_name,levels=order_by)

    plt2<-ggplot(overlap_counts_genebody,aes(x=comparison_name,y=value,fill=name))+
          geom_col()+ geom_text(aes(label = percent), vjust = -0.5, size = 5) +
          facet_wrap(~comparison_name,nrow=1,scales="free_x")+
          ylab("Gene body hypo DMR overlap for DE upregulated genes")+  
          theme(axis.text.x = element_text(angle = 90, hjust = 1))

    #hypo overlap on gene body OR promoter
    overlap_counts_genebody_or_promoter<-overlap_counts %>% 
      filter(direction_de=="upregulated") %>% 
      dplyr::select(comparison_name,celltype,total_gene_count,hypo_promoter_or_genebody_overlap)
    
    overlap_counts_genebody_or_promoter$hypo_promoter_or_genebody_nonoverlap<-overlap_counts_genebody_or_promoter$total_gene_count - overlap_counts_genebody_or_promoter$hypo_promoter_or_genebody_overlap
    overlap_counts_genebody_or_promoter<-tidyr::pivot_longer(overlap_counts_genebody_or_promoter,cols=c("hypo_promoter_or_genebody_nonoverlap","hypo_promoter_or_genebody_overlap"))

    overlap_counts_genebody_or_promoter$percent<-unlist(lapply(1:nrow(overlap_counts_genebody_or_promoter),function(x) 
      ifelse(overlap_counts_genebody_or_promoter[x,]$name=="hypo_promoter_or_genebody_overlap",
            round((overlap_counts_genebody_or_promoter[x,]$value/overlap_counts_genebody_or_promoter[x,]$total_gene_count)*100,2),
            NA)))
  overlap_counts_genebody_or_promoter$comparison_name<-factor(overlap_counts_genebody_or_promoter$comparison_name,levels=order_by)

    plt3<-ggplot(overlap_counts_genebody_or_promoter,aes(x=comparison_name,y=value,fill=name))+
          geom_col()+facet_wrap(~comparison_name,nrow=1,scales="free_x")+
          geom_text(aes(label = percent), vjust = -0.5, size = 5) +
          ylab("Gene body or promoter hypo DMR overlap for DE upregulated genes")+  
          theme(axis.text.x = element_text(angle = 90, hjust = 1))


    print("Hypomethylated DMRs that overlap promoters/genebodies")
    print(hypodmrs_that_overlap_promoters)
    setNames(nm=paste(dmr_overlap_counts$comparison_name,dmr_overlap_counts$direction),round((dmr_overlap_counts$genebody_or_promoter_overlap/dmr_overlap_counts$total_dmr_count)*100,2))
  #basal_v_rest hyper        basal_v_rest hypo       bcell_v_rest hyper 
  #                   6.24                     9.28                     2.52 
  #      bcell_v_rest hypo      cancer_v_rest hyper       cancer_v_rest hypo 
  #                   3.67                     3.85                     3.73 
  #endothelial_v_rest hyper  endothelial_v_rest hypo  fibroblast_v_rest hyper 
  #                    8.02                    11.00                    11.18 
  #  fibroblast_v_rest hypo       lumhr_v_rest hyper        lumhr_v_rest hypo 
  #                   10.81                     4.33                     5.18 
  #     lumsec_v_rest hyper       lumsec_v_rest hypo     myeloid_v_rest hyper 
  #                    3.86                     3.08                     4.11 
  #     myeloid_v_rest hypo    pericyte_v_rest hyper     pericyte_v_rest hypo 
  #                    9.75                     6.48                    11.31 
  #      tcell_v_rest hyper        tcell_v_rest hypo 
  #                    4.64                     5.40 

    dmr_overlap_counts$genebody_or_promoter_nonoverlap<-dmr_overlap_counts$total_dmr_count - dmr_overlap_counts$genebody_or_promoter_overlap

    dmr_overlap_counts_genebody_or_promoter<-tidyr::pivot_longer(dmr_overlap_counts,cols=c("genebody_or_promoter_nonoverlap","genebody_or_promoter_overlap"))

  dmr_overlap_counts_genebody_or_promoter$percent<-unlist(lapply(1:nrow(dmr_overlap_counts_genebody_or_promoter),function(x) 
      ifelse(dmr_overlap_counts_genebody_or_promoter[x,]$name=="genebody_or_promoter_overlap",
            round((dmr_overlap_counts_genebody_or_promoter[x,]$value/dmr_overlap_counts_genebody_or_promoter[x,]$total_dmr_count)*100,2),
            NA)))
  dmr_overlap_counts_genebody_or_promoter$comparison_name<-factor(dmr_overlap_counts_genebody_or_promoter$comparison_name,levels=order_by)
    plt4<-ggplot(dmr_overlap_counts_genebody_or_promoter,aes(x=comparison_name,y=value,fill=name))+
            geom_col()+
            geom_text(aes(label = percent), vjust = -0.5, size = 5) +
            facet_wrap(direction~comparison_name,nrow=1,scales="free_x")+
            ylab("hypo DMRs that overlap a promoter for DE upregulated genes")+  
            theme(axis.text.x = element_text(angle = 90, hjust = 1))

  ggsave(plt1/plt2/plt3/plt4,file=paste0(output_plots_folder,"/",comparison_set,".de_gene_dmr_hypo_overlap.pdf"),height=50,width=50,limitsize=FALSE)

    return(dmr_gene_overlap)

}

write_dmr_bed<-function(dmr=dmr,
                        prefix="celltype_v_rest",
                        output_folder=output_plots_folder){
  dmr_split<-dmr %>% split(~comparison_name)
  lapply(dmr_split,function(dmr_out){
    outname=dmr_out$comparison_name[1]
    #rename features
    mcols(dmr_out)$name <- paste(seqnames(dmr_out),start(dmr_out),end(dmr_out),sep="_")
    #convert score to integer
    mcols(dmr_out)$score <- ifelse(dmr_out$dmr_logFC<0,1,2)

    #set color for hypo and hyper
    mcols(dmr_out)$itemRgb <- ifelse(dmr_out$dmr_logFC<0,"red","gray")
    mcols(dmr_out)$thickStart <- start(dmr_out)
    mcols(dmr_out)$thickEnd <- end(dmr_out)

    mcols(dmr_out)<-mcols(dmr_out)[c("name","score","thickStart","thickEnd","itemRgb")]
    export.bed(dmr_out, con = paste0(output_folder,"/",outname,".",prefix,".dmr.bed"))
  })
}

write_de_bed<-function(markers=markers,
                        prefix="celltype_v_rest",
                        output_folder=output_plots_folder){
  

  #get overlap of genes and dmrs per celltype (need to add genomic ranges to RNA markers)
  print("Using GTF file to annotation gene locations...")
  gtf_markers<- gtf %>% 
    as.data.frame() %>% 
    filter(type=="gene") %>% 
    filter(gene_name %in% markers$gene) 

  ensembl_vector <- mapIds(
    x = org.Hs.eg.db, 
    keys = markers$gene, 
    column = "ENSEMBL", 
    keytype = "SYMBOL", 
    multiVals = "first" 
  )
  
  markers$gene_id<-ensembl_vector[markers$gene]

  markers$direction_de<-ifelse(markers$avg_log2FC>0,"upregulated","downregulated")

  markers[is.na(markers$gene_id) & startsWith(markers$gene,prefix="ENSG00"),]$gene_id<-markers[is.na(markers$gene_id) & startsWith(markers$gene,prefix="ENSG00"),]$gene

  markers[is.na(markers$gene_id),] #some antisense genes are still lost
  gtf_markers$gene_id<-unlist(lapply(strsplit(gtf_markers$gene_id,"[.]"),'[',1))

  markers2<-left_join(markers,gtf_markers,"gene_id",multiple="first")
  markers2<-markers2[!is.na(markers2$seqnames),]

  print("Converting DMRs and DE output to GenomicRanges...")

  markers_grange<-GRanges(markers2)
  markers_split<-markers_grange %>% split(~comparison_name)

  lapply(markers_split,function(marker_out){
    outname=marker_out$comparison_name[1]
    #rename features
    mcols(marker_out)$name <- mcols(marker_out)$gene_name
    #convert score to integer
    mcols(marker_out)$score <- ifelse(marker_out$avg_log2FC<0,1,2)

    #set color for hypo and hyper
    mcols(marker_out)$itemRgb <- ifelse(marker_out$avg_log2FC<0,"gray","orange")
    mcols(marker_out)$thickStart <- start(marker_out)
    mcols(marker_out)$thickEnd <- end(marker_out)

    mcols(marker_out)<-mcols(marker_out)[c("name","score","thickStart","thickEnd","itemRgb")]
    export.bed(marker_out, con = paste0(output_folder,"/",outname,".",prefix,".de.bed"))
  })
}

annotate_dmrs<-function(dmr=dmr,
                        prefix="celltype_v_rest",
                        output_folder=output_plots_folder){

  hg38_lengths <- seqlengths(Seqinfo(genome = "hg38"))
  dmr<-GRanges(dmr)
  seqinfo(dmr) <- Seqinfo(
    seqnames = levels(seqnames(dmr)),
    seqlengths = hg38_lengths[levels(seqnames(dmr))],
    isCircular = rep(FALSE,length(levels(seqnames(dmr)))),
    genome = "hg38")

  dm_random_regions = randomize_regions(
      regions = dmr,
      allow.overlaps = TRUE,
      per.chromosome = TRUE)

  #annotate regions of dmrs
  print("Annotating regions respective of genes....")

    feature_priority <- c(
      "hg38_genes_1to5kb",
      "hg38_genes_promoters",
      "hg38_genes_5UTRs",
      "hg38_genes_3UTRs",
      "hg38_genes_exons",
      "hg38_genes_introns",
      "hg38_genes_intergenic"
    )

    annotations = build_annotations(genome = 'hg38', annotations = feature_priority)

    annotated_gr <- annotate_regions(
      regions = dmr,
      annotations = annotations)

    #have to assign features since it allows for overlap
    print("Assigning features by hierarchy to single type....")
    dm_hierarchical_df <- annotated_gr %>% as.data.frame() %>%
      group_by(seqnames, start, end, comparison_name, direction) %>%
      mutate(annot.type = factor(annot.type, levels = feature_priority)) %>%
      arrange(annot.type, .by_group = TRUE) %>%
      slice(1) %>%
      ungroup()

    dm_catsum = summarize_categorical(
      annotated_regions = GRanges(dm_hierarchical_df),
      by = c("comparison_name",'annot.type',"direction"),
      quiet = FALSE) %>% as.data.frame()

    dmr_totals <- dmr %>% as.data.frame() %>% group_by(comparison_name,direction) %>% summarize(total=n()) %>% as.data.frame()
    dmr_outside_annotation<-lapply(1:nrow(dmr_totals), function(x) {
      sum_in_annot<-sum(dm_catsum[dm_catsum$comparison_name == dmr_totals$comparison_name[x] & dm_catsum$direction == dmr_totals$direction[x],]$n)
      return(data.frame(comparison_name=dmr_totals[x,]$comparison_name,
      direction=dmr_totals[x,]$direction,
      dmr_outside_annotation=dmr_totals[x,]$total-sum_in_annot))
      })
    
    dmr_outside_annotation<-do.call("rbind",dmr_outside_annotation)
    
    print("DMR sites that fall outside of annotation (sanity check)...")
    print(dmr_outside_annotation$dmr_outside_annotation)

    dm_catsum<-dm_catsum %>% 
              dplyr::group_by(comparison_name,direction) %>% 
              mutate(percent=round((n/sum(n))*100,2)) %>% 
              as.data.frame() 
    dm_catsum$celltype<-unlist(lapply(strsplit(dm_catsum$comparison_name,"_"),"[",1))
    dm_catsum$celltype<-factor(dm_catsum$celltype,levels=names(celltype_col))
    dm_catsum$annot.type<-factor(dm_catsum$annot.type,levels=feature_priority)

    print("Plotting randomized set for comparison (and statistical test)...")
    # Annotate the random regions using the same annotations as above
    dm_random_annotated = annotate_regions(
        regions = dm_random_regions,
        annotations = annotations,
        ignore.strand = TRUE,
        quiet = FALSE)

    print("Perform same hierarchical assignment for random features...")
    dm_hierarchical_df_random <- dm_random_annotated %>% as.data.frame() %>%
        group_by(seqnames, start, end) %>%
        mutate(annot.type = factor(annot.type, levels = feature_priority)) %>%
        arrange(annot.type, .by_group = TRUE) %>%
        slice(1) %>%
        ungroup()

    dm_catsum_rando = summarize_categorical(
      annotated_regions = GRanges(dm_hierarchical_df_random),
      by = c('annot.type'),
      quiet = FALSE) %>% as.data.frame()
    dm_catsum_rando$percent<-round(dm_catsum_rando$n/sum(dm_catsum_rando$n),2)
    dm_catsum_rando$annot.type<-factor(dm_catsum_rando$annot.type,levels=feature_priority)

    plt<-ggplot(dm_catsum,aes(x=celltype,y=n,fill=annot.type))+
    geom_bar(stat="identity")+#geom_text(aes(label = percent),  size = 5)+
    facet_wrap(celltype~direction,scale="free",ncol=2)
    
    plt2<-ggplot(dm_catsum_rando,aes(x="random",y=n,fill=annot.type))+
    geom_bar(stat="identity")#+geom_text(aes(label = percent),size = 5)

    ggsave(plt|plt2,file=paste0(output_folder,"/",prefix,".","dmr.gene_feature.annotation.pdf"),height=20,width=20)
    saveRDS(dm_catsum_rando,file=paste0(output_folder,"/",prefix,".","dmr.gene_feature.annotation.randomset.rds"))
    saveRDS(dm_catsum,file=paste0(output_folder,"/",prefix,".","dmr.gene_feature.annotation.counts.rds"))


  #annotate regions of dmrs
  print("Annotating regions respective of CG Islands....")

    feature_priority <- c(
        'hg38_cpg_islands',
        'hg38_cpg_shores',
        'hg38_cpg_shelves',
        'hg38_cpg_inter')

    annotations = build_annotations(genome = 'hg38', annotations = feature_priority)

    annotated_gr <- annotate_regions(
      regions = dmr,
      annotations = annotations)

    #have to assign features since it allows for overlap
    print("Assigning features by hierarchy to single type....")
    dm_hierarchical_df <- annotated_gr %>% as.data.frame() %>%
      group_by(seqnames, start, end, comparison_name, direction) %>%
      mutate(annot.type = factor(annot.type, levels = feature_priority)) %>%
      arrange(annot.type, .by_group = TRUE) %>%
      slice(1) %>%
      ungroup()

    dm_catsum = summarize_categorical(
      annotated_regions = GRanges(dm_hierarchical_df),
      by = c("comparison_name",'annot.type',"direction"),
      quiet = FALSE) %>% as.data.frame()

    dmr_totals <- dmr %>% as.data.frame() %>% group_by(comparison_name,direction) %>% summarize(total=n()) %>% as.data.frame()
    dmr_outside_annotation<-lapply(1:nrow(dmr_totals), function(x) {
      sum_in_annot<-sum(dm_catsum[dm_catsum$comparison_name == dmr_totals$comparison_name[x] & dm_catsum$direction == dmr_totals$direction[x],]$n)
      return(data.frame(comparison_name=dmr_totals[x,]$comparison_name,
      direction=dmr_totals[x,]$direction,
      dmr_outside_annotation=dmr_totals[x,]$total-sum_in_annot))
      })
    
    dmr_outside_annotation<-do.call("rbind",dmr_outside_annotation)
    
    print("DMR sites that fall outside of annotation (sanity check)...")
    print(dmr_outside_annotation$dmr_outside_annotation)

    dm_catsum$celltype<-unlist(lapply(strsplit(dm_catsum$comparison_name,"_"),"[",1))
    dm_catsum$celltype<-factor(dm_catsum$celltype,levels=names(celltype_col))
    dm_catsum<-dm_catsum %>% dplyr::group_by(comparison_name,direction) %>% mutate(percent=round((n/sum(n))*100,2)) %>% as.data.frame() 

    print("Plotting randomized set for comparison (and statistical test)...")
    # Annotate the random regions using the same annotations as above
    dm_random_annotated = annotate_regions(
        regions = dm_random_regions,
        annotations = annotations,
        ignore.strand = TRUE,
        quiet = FALSE)

    print("Perform same hierarchical assignment for random features...")
    dm_hierarchical_df_random <- dm_random_annotated %>% as.data.frame() %>%
        group_by(seqnames, start, end) %>%
        mutate(annot.type = factor(annot.type, levels = feature_priority)) %>%
        arrange(annot.type, .by_group = TRUE) %>%
        slice(1) %>%
        ungroup()

   dm_catsum_rando = summarize_categorical(
      annotated_regions = GRanges(dm_hierarchical_df_random),
      by = c('annot.type'),
      quiet = FALSE) %>% as.data.frame()
    dm_catsum_rando$percent<-round(dm_catsum_rando$n/sum(dm_catsum_rando$n),2)

    plt<-ggplot(dm_catsum,aes(x=celltype,y=n,fill=annot.type))+
    geom_bar(stat="identity")+#geom_text(aes(label = percent),  size = 5)+
    facet_wrap(celltype~direction,scale="free",ncol=2)
    
    plt2<-ggplot(dm_catsum_rando,aes(x="random",y=n,fill=annot.type))+
    geom_bar(stat="identity")#+geom_text(aes(label = percent),size = 5)

    ggsave(plt|plt2,file=paste0(output_folder,"/",prefix,".","dmr.cgisland.annotation.pdf"),width=20,height=20)
    saveRDS(dm_catsum_rando,file=paste0(output_folder,"/",prefix,".","dmr.gene_feature.annotation.randomset.rds"))
    saveRDS(dm_catsum,file=paste0(output_folder,"/",prefix,".","dmr.gene_feature.annotation.counts.rds"))

}



```

# 1. Celltype level comparison (one v rest)

### Methylation DMR summary stats
```R
comparison_set="dmr_across_celltype"
#read in DMRs (cell type comparisons)
output_folder<-paste(project_data_directory,"04_dmr",comparison_set,sep="/")
input_dmr_folder<-paste(project_data_directory,"04_dmr",comparison_set,"dmr_out",sep="/")

output_plots_folder<-paste0(output_folder,"/plots")
dir.create(output_plots_folder)

dmr=do.call("rbind",lapply(list.files(input_dmr_folder,full.names=T,pattern=".collapse.rds"),function(i){readRDS(i)}))
dmr$comparison_name=dmr$type

#filter to cell type comparisons (one v rest)
dmr<- dmr %>% filter(endsWith(comparison_name,suffix="_v_rest"))

seqinfo(dmr) <- Seqinfo(
  seqnames = levels(seqnames(dmr)),
  seqlengths = hg38_lengths[levels(seqnames(dmr))],
  isCircular = rep(FALSE,length(levels(seqnames(dmr)))),
  genome = "hg38")

output_dmr_tracks_folder<-paste("/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/04_dmr/bigwig_output_dmr_celltype_group_500bp",comparison_set,"met_dmr",sep="/")
dir.create(output_dmr_tracks_folder,recursive=T)

write_dmr_bed(dmr=dmr,prefix="celltype_v_rest",output_plots_folder=output_dmr_tracks_folder)

#dmr counts per celltype comparison
dmr_counts<-dmr %>% 
            filter(endsWith(type,suffix="_v_rest")) %>% 
            group_by(type,direction) %>% 
            summarize(count=n()) %>% 
            as.data.frame()
dmr_counts$type<-gsub(dmr_counts$type,pattern="_v_rest",repl="")
dmr_counts$type<-factor(dmr_counts$type,levels=names(celltype_col))

#plot number of dmrs per cell type (hyper and hypo)
plt<-ggplot(dmr_counts,aes(x=type,y=count,fill=type))+
      geom_bar(stat="identity")+scale_fill_manual(values=celltype_col)+
      facet_wrap(~direction)+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(plt,file=paste0(output_plots_folder,"/",comparison_set,".dmr_counts.pdf"))

```

## DMR annotation
Compare DMR features over genes and CGIsland annotations.

```
annotate_dmrs(dmr=dmr,prefix=comparison_set,output_folder=output_plots_folder)
```

### Run GSEA on methylation dmrs

```R
output_gsea_folder<-paste0(output_folder,"/gsea")
dir.create(output_gsea_folder)

#run GSEA on the DMR sites
gsea_across_sets(obj, 
                  dmr,
                  sample_name="celltype", 
                  prefix=comparison_set,
                  output_gsea_directory=output_gsea_folder)

#plot GSEA results
#transcription factor targets GSEA
gsea<-readRDS(paste0(output_gsea_folder,"/",comparison_set,".GSEA_enrichment.celltype.TFT.rds"))
gsea$names<-gsub(gsea$pathway,pattern="_TARGET_GENES",repl="") #make names a bit more readable
plot_gsea(gsea=gsea,out_setname="TFT",prefix=comparison_set,output_directory=output_plots_folder)

#cancer hallmark
gsea<-readRDS(paste0(output_gsea_folder,"/",comparison_set,".GSEA_enrichment.celltype.hallmark.rds"))
#get top hits per group
gsea$names<-gsea$pathway
plot_gsea(gsea=gsea,out_setname="hallmark",prefix=comparison_set,output_directory=output_plots_folder)

#cancer 3C
gsea<-readRDS(paste0(output_gsea_folder,"/",comparison_set,".GSEA_enrichment.celltype.3CA.rds"))
#get top hits per group
gsea$names<-gsea$pathway
plot_gsea(gsea=gsea,out_setname="3CA",prefix=comparison_set,output_directory=output_plots_folder)

#position
gsea<-readRDS(paste0(output_gsea_folder,"/",comparison_set,".GSEA_enrichment.celltype.position.rds"))
#get top hits per group
gsea$names<-gsea$pathway
plot_gsea(gsea=gsea,out_setname="position",prefix=comparison_set,output_directory=output_plots_folder)
```

## Run RNA DE
```R

output_de_folder<-paste(project_data_directory,"04_dmr","dmr_across_celltype","rna_de",sep="/")
dir.create(output_de_folder)

#calculate rna differences per cell type (one v rest)
table(Idents(rna))
markers<-FindAllMarkers(rna)
markers$comparison_name<-paste0(markers$cluster,"_v_rest")
saveRDS(markers,paste0(output_de_folder,"/",comparison_set,".rna.DE.rds"))

markers<-readRDS(paste0(output_de_folder,"/",comparison_set,".rna.DE.rds"))

markers<-markers %>% filter(p_val_adj<0.05) %>% filter(avg_log2FC>1)

output_de_tracks_folder<-paste("/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/04_dmr/bigwig_output_dmr_celltype_group_500bp",comparison_set,"rna_de",sep="/")
dir.create(output_de_tracks_folder)

write_de_bed(markers=markers,prefix="celltype_v_rest",output_folder=output_de_tracks_folder)

marker_counts<-markers %>% 
            filter(p_val_adj<0.05) %>% 
            mutate(direction=ifelse(avg_log2FC>0,"upregulated","downregulated")) %>%
            group_by(cluster,direction) %>% 
            summarize(count=n()) %>% 
            as.data.frame()

marker_counts$cluster<-factor(marker_counts$cluster,levels=names(celltype_col))
#plot number of dmrs per cell type (hyper and hypo)
plt<-ggplot(marker_counts,aes(x=cluster,y=count,fill=cluster))+
      geom_bar(stat="identity")+scale_fill_manual(values=celltype_col)+
      facet_wrap(~direction)+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(plt,file=paste0(output_plots_folder,"/",comparison_set,".celltype_DE_counts.pdf"))
```

## plot GSEA results on RNA
Using SeuratExtend for enrichment analysis
#https://github.com/huayc09/SeuratExtend/blob/master/vignettes/GSEA.md#conduct-gsea-using-the-go-or-reactome-database

```R

calculate_rna_pathway(rna,out_name="position",cores=50,prefix=comparison_set,
  genesets=SeuratExtendData::Genesets_data$human$GSEA[["positional gene sets"]],
  output_directory=output_plots_folder)
  
calculate_rna_pathway(rna,out_name="hallmark",cores=50,prefix=comparison_set,
  genesets=SeuratExtendData::Genesets_data$human$GSEA[["hallmark gene sets"]],
  output_directory=output_plots_folder)

calculate_rna_pathway(rna,out_name="tft",cores=50,prefix=comparison_set,
  genesets=SeuratExtendData::Genesets_data$human$GSEA[["transcription factor targets"]],
  output_directory=output_plots_folder)

calculate_rna_pathway(rna,out_name="canonical",cores=50,prefix=comparison_set,
  genesets=SeuratExtendData::Genesets_data$human$GSEA[["all canonical pathways"]],
  output_directory=output_plots_folder)
```

## Now plot overlap of RNA and MET DMRs
```R
plot_marker_dmr_overlap(markers=markers,
                        dmrs=dmr,
                        prefix=comparison_set,
                        output_plots_folder=output_plots_folder,
                        order_by=paste0(names(celltype_col),"_v_rest"))

```

## Now plot chromvar motif enrichment over 500bp dmrs

```R
input_dmr_folder<-paste(project_data_directory,"04_dmr",comparison_set,"dmr_out",sep="/")

dmr_500=do.call("rbind",lapply(list.files(input_dmr_folder,full.names=T,pattern=".dmr.filt.rds"),function(i){readRDS(i)}))
dmr$comparison_name=dmr$type

#filter to cell type comparisons (one v rest)
dmr<- dmr %>% filter(endsWith(comparison_name,suffix="_v_rest"))

seqinfo(dmr) <- Seqinfo(
  seqnames = levels(seqnames(dmr)),
  seqlengths = hg38_lengths[levels(seqnames(dmr))],
  isCircular = rep(FALSE,length(levels(seqnames(dmr)))),
  genome = "hg38")

#could use plotting tracks for percentages? then partition to hypomethylated

```
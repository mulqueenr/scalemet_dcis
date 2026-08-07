
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
library(Seurat)
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

#read in DMRs
output_folder<-paste(project_data_directory,"04_dmr","dmr_across_celltype",sep="/")
dmr=do.call("rbind",lapply(list.files(output_folder,full.names=T,pattern=".collapse.rds"),function(i){readRDS(i)}))
dmr$comparison_name=dmr$type

#read GTF for RNA annotation
#get overlap of genes and dmrs per celltype (need to add genomic ranges to RNA markers)
gtf<-import("/data/rmulqueen/projects/scalebio_dcis/ref/gencode.v43.annotation.gtf.gz")
gtf$gene<-gtf$gene_name
```

Functions
```R
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

plot_gsea<-function(gsea=gsea,
                    out_setname="hallmark",
                    prefix=prefix){

  #get top hits per group
  features_to_keep_hypo <- gsea %>% filter(padj<0.05) %>% group_by(test) %>% slice_min(NES,n=10) %>% select(names)
  #features_to_keep_hyper <- gsea %>% filter(padj<0.05) %>% group_by(test) %>% slice_max(NES,n=10) %>% select(names)
  features_to_keep<-unique(unlist(features_to_keep_hypo$names))

  if(length(features_to_keep)>0){
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


  # Simple dot plot of GSEA results
  plt<-ggplot(gsea_plot, aes(x = celltype, y = names)) +
    geom_point(aes(size = -log10(padj), color = NES)) +
    scale_radius(range = c(0.1, 6), name = "-log10(p-value)") +
    scale_color_gradient2(low = "#b84d9c", mid = "white", high = "#666666", midpoint = 0, name = "NES") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.minor = element_blank()
    ) + theme_minimal()

    ggsave(plt,file=paste0(prefix,".",out_setname,".NES.heatmap.pdf"),width=10,height=10)
    

    plt<-Heatmap(t(gsea_nes), #adding transpose for plotting, so column and row swapped
        col = col_fun,rect_gp = gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill) {
            #draw a rectangle at all sites
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
            if(!is.na(gsea_nes[i, j]) & !is.na(gsea_pval[i, j])){
              draw a circle sized by NES and colored by pval if significant
              grid.circle(x = x, y = y, 
                    r = (abs(gsea_pval[i, j])/max_size)/2 * min(unit.c(width, height)), 
                    gp = gpar(fill = col_fun(gsea_nes[i, j]), col = NA))}},
      column_order=1:nrow(gsea_nes),
      row_order=column_order)
    
    } else {

  #cluster by features
  gsea_nes <- gsea %>% filter(names %in% features_to_keep) %>% tidyr::pivot_wider(names_from=test, id_cols=names, values_from=NES)  %>% as.data.frame()
  row.names(gsea_nes)<-gsea_nes$names
  gsea_nes<-gsea_nes[,2:ncol(gsea_nes)]
  clus<-dist(gsea_nes, method = "euclidean") %>%  hclust( method = "complete")
  row_order <- row.names(gsea_nes)[clus$order ]
  gsea_plot<-gsea[gsea$names %in% row_order,]
  gsea_plot$names<-factor(gsea_plot$names,levels=row_order)
  gsea_plot$celltype<-gsub(gsea_plot$test,pattern="_v_rest",repl="")
  gsea_plot$celltype<-factor(gsea_plot$celltype,levels=names(celltype_col))
  
  # Simple dot plot of GSEA results
  plt<-ggplot(gsea_plot, aes(x = celltype, y = names)) +
    geom_point(aes(size = -log10(padj), color = NES)) +
    scale_radius(range = c(0.1, 6), name = "-log10(p-value)") +
    scale_color_gradient2(low = "#b84d9c", mid = "white", high = "#666666", midpoint = 0, name = "NES") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.minor = element_blank()
    ) + theme_minimal()

    ggsave(plt,file=paste0(prefix,".",out_setname,".NES.heatmap.pdf"),width=10,height=10)
    

  } else { print(paste("No significant findings in",out_setname))}

  }
}

```

# Celltype level comparison (one v rest)

## Methylation DMR summary stats
```R
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
ggsave(plt,file="celltype_dmr_counts.pdf")

#filter to cell type comparisons (one v rest)
collapsed_dmrs<- dmr %>% filter(endsWith(comparison_name,suffix="_v_rest"))

#run GSEA on the DMR sites
gsea_across_sets(obj, collapsed_dmrs,sample_name="celltype", prefix="celltype_v_rest")

#plot GSEA results
#transcription factor targets GSEA
gsea<-readRDS("/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/04_dmr/test.GSEA_enrichment.test.TFT.rds")
gsea$names<-gsub(gsea$pathway,pattern="_TARGET_GENES",repl="") #make names a bit more readable
plot_gsea(gsea=gsea,out_setname="TFT",prefix="celltype_v_rest")

#cancer hallmark
gsea<-readRDS("/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/04_dmr/test.GSEA_enrichment.test.hallmark.rds")
#get top hits per group
gsea$names<-gsea$pathway
plot_gsea(gsea=gsea,out_setname="hallmark",prefix="celltype_v_rest")

#cancer 3C
gsea<-readRDS("/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/04_dmr/test.GSEA_enrichment.test.hallmark.rds")
#get top hits per group
gsea$names<-gsea$pathway
plot_gsea(gsea=gsea,out_setname="C3",prefix="celltype_v_rest")

```

## RNA markers calculation and summary stats

```
#compare dmr sites with rna marker genes overlap
table(Idents(rna))
markers<-FindAllMarkers(rna)

output_folder<-paste(project_data_directory,"04_dmr","dmr_across_celltype",sep="/")

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
ggsave(plt,file="celltype_DE_counts.pdf")

gtf<- gtf %>% 
  as.data.frame() %>% 
  filter(type=="gene") %>% 
  filter(gene_name %in% markers$gene) 

ensembl_vector <- mapIds(
  x = org.Hs.eg.db, 
  keys = markers$gene, 
  column = "ENSEMBL", 
  keytype = "SYMBOL", 
  multiVals = "first" # Keeps the first Ensembl ID matched if there's a 1-to-many match
)

markers$gene_id<-ensembl_vector[markers$gene]
markers$direction_de<-ifelse(markers$avg_log2FC>0,"upregulated","downregulated")
markers[is.na(markers$gene_id) & startsWith(markers$gene,prefix="ENSG00"),]$gene_id<-markers[is.na(markers$gene_id) & startsWith(markers$gene,prefix="ENSG00"),]$gene
markers[is.na(markers$gene_id),] #some antisense are still lost
gtf$gene_id<-unlist(lapply(strsplit(gtf$gene_id,"[.]"),'[',1))
markers2<-left_join(markers,gtf,"gene_id",relationship="many-to-many")
markers2<-markers2[!is.na(markers2$seqnames),]

markers_grange<-GRanges(markers2)
dmr_grange<-GRanges(dmr)

#get count of dmr overlaps per marker genes
#add 5kb upstream of promoter
dmr_overlap<-lapply(as.character(unique(markers$cluster)),function(celltype){
  #prepare RNA differences
  markers_grange_celltype<-markers_grange %>% as.data.frame() %>% filter(cluster==celltype) %>% GRanges()
  promoter_ranges <- promoters(markers_grange_celltype, upstream = 5000, downstream = 2000)
  
  #prepare dmr
  dmr_grange_celltype<-dmr_grange %>% as.data.frame() %>% filter(dmr_grange$type==paste0(celltype,"_v_rest")) %>% GRanges()
  dmr_grange_celltype$dmr_name<-paste(seqnames(dmr_grange_celltype),start(dmr_grange_celltype),end(dmr_grange_celltype),sep="_")
  dmr_grange_celltype <- split(dmr_grange_celltype,~direction)
  
  overlaps_genebody_hypo <- findOverlaps(query = dmr_grange_celltype$hypo, subject = markers_grange_celltype)
  overlaps_genebody_hyper <- findOverlaps(query = dmr_grange_celltype$hyper, subject = markers_grange_celltype)

  overlaps_promoter_hypo <- findOverlaps(query = dmr_grange_celltype$hypo, subject = promoter_ranges)
  overlaps_promoter_hyper <- findOverlaps(query = dmr_grange_celltype$hyper, subject = promoter_ranges)

  #add collapsed list of DMRs as metadata column into markers_grange_celltype
  markers_grange_celltype$gene_body_dmrs_hypo<-NA
  markers_grange_celltype[overlaps_genebody_hypo@to,]$gene_body_dmrs_hypo<-paste(dmr_grange_celltype$hypo[overlaps_genebody_hypo@from,]$dmr_name,collapse=",")

  markers_grange_celltype$gene_body_dmrs_hyper<-NA
  markers_grange_celltype[overlaps_genebody_hyper@to,]$gene_body_dmrs_hyper<-paste(dmr_grange_celltype$hyper[overlaps_genebody_hyper@from,]$dmr_name,collapse=",")

  #and promoters
  markers_grange_celltype$gene_promoter_dmrs_hypo<-NA
  markers_grange_celltype[overlaps_promoter_hypo@to,]$gene_promoter_dmrs_hypo<-paste(dmr_grange_celltype$hypo[overlaps_promoter_hypo@from,]$dmr_name,collapse=",")

  markers_grange_celltype$gene_promoter_dmrs_hyper<-NA
  markers_grange_celltype[overlaps_promoter_hyper@to,]$gene_promoter_dmrs_hyper<-paste(dmr_grange_celltype$hyper[overlaps_promoter_hyper@from,]$dmr_name,collapse=",")

  return(markers_grange_celltype)
})


dmr_gene_overlap<-do.call("c",dmr_overlap)

saveRDS(dmr_gene_overlap,"04_DE_DMR_overlaps.celltype_v_rest.rds")

overlap_counts<-dmr_gene_overlap %>% 
  as.data.frame() %>% 
  filter(type=="gene") %>% 
  group_by(cluster,direction_de) %>% 
  summarize(
    gene_count=n(),
    hypo_genebody_overlap=sum(!is.na(gene_body_dmrs_hypo)),
    hyper_genebody_overlap=sum(!is.na(gene_body_dmrs_hyper)),
    hypo_promoter_overlap=sum(!is.na(gene_promoter_dmrs_hypo)),
    hyper_promoter_overlap=sum(!is.na(gene_promoter_dmrs_hyper)),
    hypo_promoter_or_genebody_overlap=sum(!is.na(gene_promoter_dmrs_hypo) | !is.na(gene_body_dmrs_hypo))) %>%
    as.data.frame()


overlap_counts$cluster<-factor(overlap_counts$cluster,levels=names(celltype_col))

library(dplyr)

#can make this a function

overlap_counts_promoter<-overlap_counts %>% filter(direction_de=="upregulated") %>% dplyr::select(cluster,gene_count,hypo_promoter_overlap)
overlap_counts_promoter$hypo_promoter_nonoverlap<-overlap_counts_promoter$gene_count - overlap_counts_promoter$hypo_promoter_overlap

overlap_percent_promoter<-overlap_counts_promoter %>% mutate(perc_cov=(hypo_promoter_overlap/gene_count)*100)
overlap_percent_promoter<-setNames(overlap_percent_promoter$cluster,round(overlap_percent_promoter$perc_cov,2))

overlap_counts_promoter<-pivot_longer(overlap_counts_promoter,cols=c("hypo_promoter_nonoverlap","hypo_promoter_overlap"))
plt1<-ggplot(overlap_counts_promoter,aes(x=cluster,y=value,fill=name))+
      geom_col()+
      ylab("Promoter hypo DMR overlap for DE upregulated genes")
plt1_1<-ggplot()+geom_text(aes(x=names(overlap_percent_promoter),y=0,label=overlap_percent_promoter,size=0.5))

overlap_counts_genebody<-overlap_counts %>% filter(direction_de=="upregulated") %>% dplyr::select(cluster,gene_count,hypo_genebody_overlap)
overlap_counts_genebody$hypo_genebody_nonoverlap<-overlap_counts_genebody$gene_count - overlap_counts_genebody$hypo_genebody_overlap

overlap_percent_genebody<-overlap_counts_genebody %>% mutate(perc_cov=(hypo_genebody_overlap/gene_count)*100)
overlap_percent_genebody<-setNames(overlap_percent_genebody$cluster,round(overlap_percent_genebody$perc_cov,2))

overlap_counts_genebody<-pivot_longer(overlap_counts_genebody,cols=c("hypo_genebody_nonoverlap","hypo_genebody_overlap"))
plt2<-ggplot(overlap_counts_genebody,aes(x=cluster,y=value,fill=name))+
      geom_col()+
      ylab("Gene body hypo DMR overlap for DE upregulated genes")
plt2_2<-ggplot()+geom_text(aes(x=names(overlap_percent_genebody),y=0,label=overlap_percent_genebody,size=0.5))

overlap_counts_genebody_or_promoter<-overlap_counts %>% filter(direction_de=="upregulated") %>% dplyr::select(cluster,gene_count,hypo_promoter_or_genebody_overlap)
overlap_counts_genebody_or_promoter$hypo_genebody_or_promoter_nonoverlap<-overlap_counts_genebody_or_promoter$gene_count - overlap_counts_genebody_or_promoter$hypo_promoter_or_genebody_overlap

overlap_percent_genebody_or_promoter<-overlap_counts_genebody_or_promoter %>% mutate(perc_cov=(hypo_promoter_or_genebody_overlap/gene_count)*100)
overlap_percent_genebody_or_promoter<-setNames(overlap_percent_genebody_or_promoter$cluster,round(overlap_percent_genebody_or_promoter$perc_cov,2))

overlap_counts_genebody_or_promoter<-pivot_longer(overlap_counts_genebody_or_promoter,cols=c("hypo_genebody_or_promoter_nonoverlap","hypo_promoter_or_genebody_overlap")) %>% as.data.frame()
plt3<-ggplot(overlap_counts_genebody,aes(x=cluster,y=value,fill=name))+
      geom_col()+ 
      ylab("Gene body or promoter hypo DMR overlap for DE upregulated genes")
plt3_2<-ggplot()+geom_text(aes(x=names(overlap_percent_genebody_or_promoter),y=0,label=overlap_percent_genebody_or_promoter,size=0.5))

ggsave(plt1/plt1_1/plt2/plt2_2/plt3/plt3_2,file="04_de_gene_dmr_hypo_overlap.pdf")

```

One v one comparisons (to copy DMR analysis) within cell types and across groups.

```R

#save markers (by celltype)
#run markers (celltype x group)
#do similar DMR overlap comparison
rna$celltype_by_group<-paste(rna$coarse_celltype,rna$Group,sep="_")
Idents(rna)<-rna$celltype_by_group

celltype_by_group_markers<-do.call("rbind",lapply(unique(rna$coarse_celltype),function(celltype){
  print(paste("Running celltype x group comparisons for",celltype))
  print(paste("Running one v rest comparisons by group"))

  rna_sub<-subset(rna,cells=row.names(rna@meta.data[rna@meta.data$coarse_celltype==celltype,]))
  #one v rest within cell type
  markers_one_v_rest<-FindAllMarkers(rna_sub)
  markers_one_v_rest$comparison_name<-paste0(markers_one_v_rest$cluster,"_v_rest_within_celltype")
  markers_one_v_rest$ident.1=markers_one_v_rest$cluster
  markers_one_v_rest$ident.2=paste0(celltype,"_rest")

  #one v one within celltype
  markers_one_v_one<-do.call("rbind",
  lapply(unique(Idents(rna_sub)),function(i){
    do.call("rbind",lapply(unique(Idents(rna_sub)),function(j){
        if(i!=j){
          print(paste("Running DE on",paste0(i,"_v_",j)))
          markers_one_v_one<-FindMarkers(rna_sub,
            ident.1=i,
            ident.2=j)
          markers_one_v_one$comparison_name<-paste0(i,"_v_",j)
          markers_one_v_one$ident.1=i
          markers_one_v_one$ident.2=j
          markers_one_v_one$gene=row.names(markers_one_v_one)
          markers_one_v_one$cluster<-markers_one_v_one$comparison_name
          markers_one_v_one<-markers_one_v_one[,colnames(markers_one_v_rest)]
          return(markers_one_v_one)
        }}))
  }))
  out<-rbind(markers_one_v_rest,markers_one_v_one)
  return(out)
}))

#run find all markers per celltype group and do same dmr overlap

```
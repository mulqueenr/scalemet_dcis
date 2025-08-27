
library(fgsea)
library(msigdbr)
library(ggplot2)
library(patchwork)
library(optparse)

gsea_enrichment<-function(annot,
                          species="mouse",
                          category="M3",
                          subcategory="TFT:GTRD",
                          out_setname="TFT",
                          outname=outname,
                          de_features_set,
                          col,group1,group2,
                          assay){
  pathwaysDF <- msigdbr(species=species, 
                        category=category, 
                        subcategory = subcategory)

  #limit pathways to genes in our data
  pathwaysDF<-pathwaysDF[pathwaysDF$ensembl_gene %in% unique(annot[annot$gene_biotype=="protein_coding",]$gene_id),]
  
  pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

  group1_features<-de_features_set %>%
    dplyr::filter(group == group1) %>%
    dplyr::filter(padj<0.05) %>% 
    dplyr::arrange(desc(logFC)) %>%
    dplyr::select(feature, logFC)

  ranks<-setNames(nm=group1_features$feature,group1_features$logFC)

  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = ranks,
                    minSize  = 10,
                    nproc = 1)

  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

  plt1<-plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam=0.5)+
        theme(axis.text.y = element_text( size = rel(0.2)),
        axis.text.x = element_text( size = rel(0.2)))

  # only plot the top 20 pathways NES scores
  nes_plt_dat<-rbind(
    fgseaRes  %>% arrange(desc(NES)) %>% head(n= 10),
    fgseaRes  %>% arrange(desc(NES)) %>% tail(n= 10))
  
  nes_plt_dat$col<-"#808080"
  nes_plt_dat[nes_plt_dat$NES>0 & nes_plt_dat$pval<0.05,]$col<-col[group1]
  nes_plt_dat[nes_plt_dat$NES<0 & nes_plt_dat$pval<0.05,]$col<-col[group2]

  plt2<-ggplot(nes_plt_dat, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill= col)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
        title="Hallmark pathways NES from GSEA") + 
    theme_minimal()+scale_fill_identity()+ggtitle(out_setname)
  return(patchwork::wrap_plots(list(plt1,plt2),ncol=2))
}

plot_gsea<-function(obj,annot,dmrs,
                    outname=outname,
                    assay=assay,col=col,group1,group2){

  #run gsea enrichment on different sets https://www.gsea-msigdb.org/gsea/msigdb/mouse/genesets.jsp?collection=CP:REACTOME
  #REACTOME
  react_plt<-gsea_enrichment(species="mouse",
              category="M2",
              subcategory="CP:REACTOME",
              out_setname="REACTOME",outname=outname,
              de_features_set=dmrs,
              col=col,group1=group1,group2=group2,assay=assay,annot=annot)

  position_plt<-gsea_enrichment(species="mouse",
              category="M1",
              subcategory=NULL,
              out_setname="position",outname=outname,
              de_features_set=dmrs,
              col=col,group1=group1,group2=group2,assay=assay,annot=annot)

  hallmark_plt<-gsea_enrichment(species="mouse",
              category="MH",
              subcategory=NULL,
              out_setname="hallmark",outname=outname,
              de_features_set=dmrs,
              col=col,group1=group1,group2=group2,assay=assay,annot=annot)

  tft_plt<-gsea_enrichment(species="mouse",
              category="M3",
              subcategory="GTRD",
              out_setname="tftargets",outname=outname,
              de_features_set=dmrs,
              col=col,group1=group1,group2=group2,assay=assay,annot=annot)

  plt<-patchwork::wrap_plots(list(react_plt,position_plt,hallmark_plt,tft_plt),nrow=4)
  ggsave(plt,file=paste0("pairwise.",outname,".",assay,".gsea.NES.pdf"),width=20,height=10)
}


plot_gsea(obj=obj,
            dmrs=markers_peaks,
            outname="gsea_peaks",
            assay=assay,
            group1=group1,
            group2=group2,
            col=col,
            annot=annot)
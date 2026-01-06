#```bash
#singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
#source activate
#update data.table
#```


```R
#install.packages("data.table")
#ensure 0.9 version of amethyst
#system('wget https://github.com/lrylaarsdam/amethyst/releases/download/v0.0.0.9000/amethyst_0.0.0.9000.tar.gz')
#system('tar -xf amethyst_0.0.0.9000.tar.gz')
#devtools::install("./amethyst")
#source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
#running locally

set.seed(111)
library(amethyst)
library(data.table)
library(dplyr)
library(msigdbr)
library(fgsea)
library(ggplot2)

options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)
obj<-readRDS(file="06_scaledcis.cnv_clones.amethyst.rds")

#update h5 files for cells
#h5_directory=paste0(project_data_directory,"/h5_files")
#system(paste("mkdir",h5_directory))
#for(i in unique(obj@h5paths$path)){
#   system(paste0("facet convert ",paste0(h5_directory,"/",basename(i))," ",i))
#}


```


Running through all the 500bp window generation split by:
- Celltype (fine_celltype)
- CNV based clones (cnv_clonenames)
- CNV based clones against all lumhr (cnv_clonename_alllumhr)
- Celltype + Diagnosis (Group + fine_celltype)

Saving to folders for further processing.
For CNV based clones:
1. loop through cnv based clones samples to do DMR per sample per clone.
2. group together clones that have shared features (e.g. 1q amp) and compare to others
3. Plot DMR regions across genome per clone (see if clustered around major CNV changes or in trans)

For celltype by diagnosis
1. DMRs for TME cells per diagnosis (especially SYNCH vs non-SYNCH DCIS, and normal vs DCIS)

```R
#############################################
#Make 500bp windows once and slice for comparisons
#############################################
dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")
system(paste("mkdir -p", dmr_outdir))

#clones (limited to lumhr and cancer)
dmr_clones_outdir=paste(sep="/",dmr_outdir,"cnv_clones")
system(paste("mkdir -p", dmr_clones_outdir))
obj_lumhr<-subsetObject(obj,cells=row.names(obj@metadata)[(obj@metadata$fine_celltype %in% c("lumhr","cancer"))])

clone500bpwindows <- amethyst::calcSmoothedWindows(obj_lumhr, 
                                        type = "CG", 
                                        threads = 20,
                                        step = 500,
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "cnv_clonename",
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)
saveRDS(clone500bpwindows,file=paste0(dmr_clones_outdir,"/","dmr_analysis.cnv_clones.500bp_windows.rds"))

#clones (all lumhr treated as same group)
dmr_clones_outdir=paste(sep="/",dmr_outdir,"cnv_clones_alllumhr")
system(paste("mkdir -p", dmr_clones_outdir))
obj_lumhr_all<-subsetObject(obj,cells=row.names(obj@metadata)[(obj@metadata$fine_celltype %in% c("lumhr","cancer"))])
obj_lumhr_all@metadata$cnv_clonename_alllumhr<-obj_lumhr_all@metadata$cnv_clonename
obj_lumhr_all@metadata[obj_lumhr_all@metadata$fine_celltype=="lumhr",]$cnv_clonename_alllumhr<-"lumhr"
clone500bpwindows <- amethyst::calcSmoothedWindows(obj_lumhr_all, 
                                        type = "CG", 
                                        threads = 20,
                                        step = 500,
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "cnv_clonename_alllumhr",
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)
saveRDS(clone500bpwindows,file=paste0(dmr_clones_outdir,"/","dmr_analysis.cnv_clones_alllumhr.500bp_windows.rds"))

#celltype
dmr_celltype_outdir=paste(sep="/",dmr_outdir,"celltype")
system(paste("mkdir -p", dmr_celltype_outdir))

celltype500bpwindows <- calcSmoothedWindows(obj, 
                                        type = "CG", 
                                        threads = 1,
                                        step = 500, # change to 500 for real data unless you have really low coverage
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "fine_celltype",
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)
saveRDS(celltype500bpwindows,file=paste0(dmr_celltype_outdir,"/","dmr_analysis.celltype.500bp_windows.rds"))


#celltype by diagnosis
dmr_celltype_by_diag_outdir=paste(sep="/",dmr_outdir,"celltype_by_diagnosis")
system(paste("mkdir -p", dmr_celltype_by_diag_outdir))
obj@metadata$celltype_diag<-paste(sep="_",obj@metadata$Group,obj@metadata$fine_celltype)
celltypediag500bpwindows <- calcSmoothedWindows(obj, 
                                        type = "CG", 
                                        threads = 1,
                                        step = 500, # change to 500 for real data unless you have really low coverage
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "fine_celltype",
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)
saveRDS(celltypediag500bpwindows,file=paste0(dmr_celltype_by_diag_outdir,"/","dmr_analysis.celltype_by_diag.500bp_windows.rds"))

```


```R
#############################################
###### DMR Per CNV Clone Per Sample    ######
#############################################

#add plotting functions for DMR 
#https://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/pbmc_vignette/pbmc_vignette.html
#add dmr locations per clone over heatmap cnv plot

#output to DMR analysis folder
dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")

#run 500bp smoothed windows by clone split once, save and subset

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

dmr_clone_by_clone<-function(obj,sample_name){
    #set dir
    clone_outdir=paste0(dmr_clones_outdir,"/",sample_name)
    system(paste("mkdir -p",clone_outdir))
    
    print(paste("Subset object for sample:",sample_name))
    
    #subset object (sample name match AND fine cell type is lumhr or cancer (by original 500bp summaries above))
    cells_in=row.names(obj@metadata)[(obj@metadata$Sample %in% c(sample_name))]
    obj_sample<-subsetObject(obj,cells=cells_in)
    #subset cell types to just epithelial/cancer

    #read 500bp windows
    clone500bpwindows<-readRDS(paste0(dmr_clones_outdir,"/","dmr_analysis.cnv_clones.500bp_windows.rds"))
    
    #subset to just sample clones
    print(paste("Subset 500bp windows for sample:",sample_name))

    pct_columns_to_keep=c(1,2,3,which(grepl(pattern=sample_name,unlist(lapply(strsplit(colnames(clone500bpwindows[["pct_matrix"]]),"_"),"[",1)))))
    sum_columns_to_keep=c(1,2,3,which(grepl(pattern=sample_name,unlist(lapply(strsplit(colnames(clone500bpwindows[["sum_matrix"]]),"_"),"[",1)))))
    
    pct_mat<-clone500bpwindows[["pct_matrix"]] %>% select(all_of(pct_columns_to_keep))
    sum_mat<-clone500bpwindows[["sum_matrix"]] %>% select(all_of(sum_columns_to_keep))

    print(paste("Saving amethyst object for clones:",sample_name))
    #save clone object for future genome track plotting
    obj_sample@genomeMatrices[["cg_clone_tracks"]] <- pct_mat #load it into amethyst object for plotting
    saveRDS(obj_sample,file=paste0(clone_outdir,"/",sample_name,".clones_amethyst.rds"))

    print(paste("Running DMRs across clones for:",sample_name))
    dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
            eachVsAll = TRUE, # If TRUE, each group found in the sumMatrix will be tested against all others
            nminTotal = 3, # Min number observations across all groups to include the region in calculations
            nminGroup = 3) # Min number observations across either members or nonmembers to include the region
    saveRDS(dmrs,file=paste0(clone_outdir,"/",sample_name,".clones_dmr.rds"))
   
    print(paste("Filtering DMRs across clones for:",sample_name))
    dmrs <- filterDMR(dmrs, 
                method = "bonferroni", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
                filter = FALSE, # If TRUE, removes insignificant results
                pThreshold = 0.05, # Maxmimum adjusted p value to allow if filter = TRUE
                logThreshold = 1.5) # Minimum absolute value of the log2FC to allow if filter = TRUE
    collapsed_dmrs <- collapseDMR(obj_sample, 
                            dmrs, 
                            maxDist = 2000, # Max allowable overlap between DMRs to be considered adjacent
                            minLength = 2000, # Min length of collapsed DMR window to include in the output
                            reduce = T, # Reduce results to unique observations (recommended)
                            annotate = T) # Add column with overlapping gene names
    saveRDS(collapsed_dmrs,file=paste0(clone_outdir,"/",sample_name,".clones_dmr_filt_collapse.rds"))

    rename_dmr_output<-setNames(nm=1:length(unique(collapsed_dmrs$test)),gsub(colnames(sum_mat)[grepl(colnames(sum_mat),pattern="_t$")],pattern="_t",replacement=""))
    collapsed_dmrs$type <- rename_dmr_output[collapsed_dmrs$test]
    saveRDS(collapsed_dmrs,file=paste0(clone_outdir,"/",sample_name,".clones_dmr_filt_collapse.rds"))
    
    #plot dmr counts per clone
    pal <- colorRampPalette(colors = c("#F9AB60", "#E7576E", "#630661", "#B5DCA5"))
    COLS <- pal(length(unique(collapsed_dmrs$type)))

    plt<-ggplot(collapsed_dmrs |> dplyr::group_by(type, direction) |> dplyr::summarise(n = n()), 
        aes(y = type, x = n, fill = type)) + 
        geom_col() + 
        facet_grid(vars(direction), scales = "free_y") + 
        scale_fill_manual(values = COLS) + 
        theme_classic()
    ggsave(plt,file=paste0(clone_outdir,"/",sample_name,".clones_dmr_counts.pdf"))
}

lapply(clone_v_clone_analysis,function(x) dmr_clone_by_clone(obj,sample_name=x))
```

```R
#add plotting functions for DMR 
#https://htmlpreview.github.io/?https://github.com/lrylaarsdam/amethyst/blob/main/vignettes/pbmc_vignette/pbmc_vignette.html
#add dmr locations per clone over heatmap cnv plot

#output to DMR analysis folder
dmr_outdir=paste(sep="/",project_data_directory,"DMR_analysis")

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

dmr_clone_by_lumhr<-function(obj,sample_name){
    #set dir
    clone_outdir=paste0(dmr_clones_outdir,"/",sample_name)
    system(paste("mkdir -p",clone_outdir))
    
    print(paste("Subset object for sample:",sample_name))
    
    #subset object to all lumhr OR sample specific clones
    cells=row.names(obj@metadata)[obj@metadata$Sample %in% c(sample_name) & obj@metadata$cnv_clonename != "NA"]
    lumhr=row.names(obj@metadata)[obj@metadata$fine_celltype %in% c("lumhr")]
    cells<-unique(c(cells,lumhr))
    obj_sample<-subsetObject(obj,cells=cells)

    #read 500bp windows
    clone500bpwindows<-readRDS(paste0(dmr_clones_outdir,"/","dmr_analysis.cnv_clones_alllumhr.500bp_windows.rds"))
    
    #subset to just sample clones
    print(paste("Subset 500bp windows for sample:",sample_name))
    pct_columns_to_keep=c(1,2,3,which(grepl(pattern=sample_name,unlist(lapply(strsplit(colnames(clone500bpwindows[["pct_matrix"]]),"_"),"[",1)))))
    sum_columns_to_keep=c(1,2,3,4,5,which(grepl(pattern=sample_name,unlist(lapply(strsplit(colnames(clone500bpwindows[["sum_matrix"]]),"_"),"[",1)))))
    
    pct_mat<-clone500bpwindows[["pct_matrix"]] %>% select(all_of(pct_columns_to_keep))
    sum_mat<-clone500bpwindows[["sum_matrix"]] %>% select(all_of(sum_columns_to_keep))

    print(paste("Saving amethyst object for clones:",sample_name))
    #save clone object for future genome track plotting
    obj_sample@genomeMatrices[["cg_clone_tracks"]] <- pct_mat #load it into amethyst object for plotting
    saveRDS(obj_sample,file=paste0(clone_outdir,"/",sample_name,".clones_lumhrall_amethyst.rds"))

    #run tests of clones against lumhr, and all clones v lumhr
    comparisons=colnames(pct_mat)[5:ncol(pct_mat)]
    comparisons<-comparisons[!grepl(comparisons,pattern="_diploid")]
    comparisons<-comparisons[!is.na(comparisons)]

    compare_tests=data.frame(
      name=c(paste0(comparisons,"_lumhr"),paste0(sample_name,"_allclones_lumhr")),
      A=c(comparisons,paste(comparisons,collapse=",")),
      B=rep("lumhr",length(comparisons)+1))

    print(paste("Running DMRs across clones for:",sample_name))
    dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
            comparisons=compare_tests,
            nminTotal = 3, # Min number observations across all groups to include the region in calculations
            nminGroup = 3) # Min number observations across either members or nonmembers to include the region
    saveRDS(dmrs,file=paste0(clone_outdir,"/",sample_name,".clones_lumhrall_dmr.rds"))
   
    print(paste("Filtering DMRs across clones for:",sample_name))
    dmrs <- filterDMR(dmrs, 
                method = "bonferroni", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
                filter = FALSE, # If TRUE, removes insignificant results
                pThreshold = 0.05, # Maxmimum adjusted p value to allow if filter = TRUE
                logThreshold = 1.5) # Minimum absolute value of the log2FC to allow if filter = TRUE
    collapsed_dmrs <- collapseDMR(obj_sample, 
                            dmrs, 
                            maxDist = 2000, # Max allowable overlap between DMRs to be considered adjacent
                            minLength = 2000, # Min length of collapsed DMR window to include in the output
                            reduce = T, # Reduce results to unique observations (recommended)
                            annotate = T) # Add column with overlapping gene names
    saveRDS(collapsed_dmrs,file=paste0(clone_outdir,"/",sample_name,".clones_lumhrall_dmr_filt_collapse.rds"))

    rename_dmr_output<-setNames(nm=1:length(unique(collapsed_dmrs$test)),compare_tests$name)
    collapsed_dmrs$type <- rename_dmr_output[collapsed_dmrs$test]
    saveRDS(collapsed_dmrs,file=paste0(clone_outdir,"/",sample_name,".clones_lumhrall_dmr_filt_collapse.rds"))
    
    #plot dmr counts per clone
    pal <- colorRampPalette(colors = c("#F9AB60", "#E7576E", "#630661", "#B5DCA5"))
    COLS <- pal(length(unique(collapsed_dmrs$type)))

    plt<-ggplot(collapsed_dmrs |> dplyr::group_by(type, direction) |> dplyr::summarise(n = n()), 
        aes(y = type, x = n, fill = type)) + 
        geom_col() + 
        facet_grid(vars(direction), scales = "free_y") + 
        scale_fill_manual(values = COLS) + 
        theme_classic()
    ggsave(plt,file=paste0(clone_outdir,"/",sample_name,".clones_lumhrall_dmr_counts.pdf"))
}

lapply(clone_v_lumhr_analysis,function(x) dmr_clone_by_lumhr(obj,sample_name=x))



```

Run DMR results through gene ontology for interpretation (to be coded)


```R


#gsea of top DE genes
gsea_enrichment<-function(annot,species="human",
                          category="C3",
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
    #dplyr::filter(padj<0.05) %>% 
    dplyr::arrange(logFC) %>%
    dplyr::select(feature, logFC)

  ranks<-setNames(nm=group1_features$feature,group1_features$logFC)

  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = ranks,
                    minSize  = 10,
                    nproc = 1)

  topPathwaysUp <- fgseaRes %>% filter(ES > 0) %>% slice_max(NES,n=10) %>% dplyr::select(pathway)
  topPathwaysDown <- fgseaRes %>% filter(ES < 0) %>% slice_max(abs(NES),n=10) %>% dplyr::select(pathway)
  topPathways <- unlist(c(topPathwaysUp, rev(topPathwaysDown)))

  plt1<-plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam=0.5)+
        theme(axis.text.y = element_text( size = rel(0.2)),
        axis.text.x = element_text( size = rel(0.2)))
  #plot gsea ranking of genes by top pathways
  #pdf(paste0("pairwise.",outname,".",assay,".",out_setname,".gsea.pdf"),width=20,height=10)
  #print(plt)
  #dev.off()

  # only plot the top 20 pathways NES scores
  nes_plt_dat<-rbind(
    fgseaRes  %>% slice_max(NES,n= 10),
    fgseaRes  %>% slice_min(NES,n= 10))
  
  nes_plt_dat$col<-"#808080"
  nes_plt_dat[nes_plt_dat$NES>0 & nes_plt_dat$pval<0.05,]$col<-col[group1]
  nes_plt_dat[nes_plt_dat$NES<0 & nes_plt_dat$pval<0.05,]$col<-col[group2]

  plt2<-ggplot(nes_plt_dat, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill= col)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
        title="Hallmark pathways NES from GSEA") + 
    theme_minimal()+scale_fill_identity()+ggtitle(out_setname)+ylim(c(-4,4))
  return(patchwork::wrap_plots(list(plt1,plt2),ncol=2))
}

plot_gsea<-function(obj,annot,dmrs,
                    outname=outname,
                    assay=assay,col=col,group1,group2,outdir){

  #run gsea enrichment on different sets
  tft_plt<-gsea_enrichment(species="human",
              category="C3",
              subcategory="TFT:GTRD",
              out_setname="TFT",outname=outname,
              de_features_set=dmrs,
              col=col,group1=group1,group2=group2,assay=assay,annot=annot)

  position_plt<-gsea_enrichment(species="human",
              category="C1",
              subcategory=NULL,
              out_setname="position",outname=outname,
              de_features_set=dmrs,
              col=col,group1=group1,group2=group2,assay=assay,annot=annot)

  hallmark_plt<-gsea_enrichment(species="human",
              category="H",
              subcategory=NULL,
              out_setname="hallmark",outname=outname,
              de_features_set=dmrs,
              col=col,group1=group1,group2=group2,assay=assay,annot=annot)
  plt<-patchwork::wrap_plots(list(tft_plt,position_plt,hallmark_plt),nrow=3,axes="collect_x")
  ggsave(plt,
          file=paste("pairwise",outname,"GSEA",assay,"pdf",sep="."),
          path=outdir,
          width=20,height=10)


}
```



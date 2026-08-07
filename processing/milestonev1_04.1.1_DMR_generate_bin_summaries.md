
Running differential methylation across cell types (one v all strategy)

Running DMR comparisons across different sets.

Because power is different (due to different cell counts per group) should limit to FC with sufficient power.

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

## Run 500 or 250bp windows per celltype/group and per clone and lumhr

Then split in different ways to calculate DMR

So the idea is to run these comparisons first on RNA then recapitulate on methylation.

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
library(sparseMatrixStats)
library(irlba)
library(rtracklayer)

project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"

#read in object from directory
task_cpus=300
processing_folder="04_dmr"
wd=paste(sep="/",project_data_directory,processing_folder)
system(paste0("mkdir -p ",wd))
setwd(wd)

```
Functions
```R
#generate DMR input
write_binned_bigwigs<-function(celltype_tracks=celltype_tracks,
                                outdir=celltype_outdir, step=step,
                                i){
    #split bw into 4 files per cell type by methylation/average methylation
    out_dat<-celltype_tracks %>% distinct(chr,start,end, .keep_all=TRUE) %>% select(chr,start,end,i) 

    hg38_seq_info<-Seqinfo(genome="hg38")
    out_dat<-GRanges(out_dat[complete.cases(out_dat),]) #filter NA
    out_dat<-out_dat[out_dat@seqnames %in% hg38_seq_info@seqnames,] #filter chr
    out_dat<-resize(out_dat,width=step)
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
    export.bw(out_dat_hypermet,con=paste0(outdir,"/",i,".hypermet.bw"))
    export.bw(out_dat_met_mid,con=paste0(outdir,"/",i,".midmet.bw"))
    export.bw(out_dat_met_low,con=paste0(outdir,"/",i,".lowmet.bw"))
    export.bw(out_dat_met_hypomet,con=paste0(outdir,"/",i,".hypomet.bw"))
  #make a ucsc/igv trackhub
  col_out=c("hypermet"="#333333","midmet"="#444444","lowmet"="#bcbcbc","hypomet"="#FF13F0")
  
  #write out multiwig trackhub
  multiwig<-c(
      paste(c("track", paste0(i)),collapse=" "),
      paste(c("container", "multiWig"),collapse=" "),
      paste(c("shortLabel",  paste0(i)),collapse=" "),
      paste(c("longLabel",  paste0(i)),collapse=" "),
      paste(c("visibility", "full"),collapse=" "),
      paste(c("aggregate", "transparentOverlay"),collapse=" "),
      paste(c("showSubtrackColorOnUi", "on"),collapse=" "))
  track_lines<-lapply(names(col_out),function(met){
    return(c(
      paste(c("track", paste0(i,".",met,".bw")),collapse=" "),
      paste(c("shortLabel", paste0(i,".",met)),collapse=" "),
      paste(c("longLabel",  paste0(i,".",met)),collapse=" "),
      paste(c("parent",  paste0(i)),collapse=" "),
      paste(c("type", "bigWig"),collapse=" "),
      paste(c("group", i),collapse=" "),
      paste(c("autoscale", "group"),collapse=" "),
      paste(c("graphTypeDefault", "bar"),collapse=" "),
      paste(c("color", paste(col2rgb(col_out[met])[,1],collapse=","),collapse=" "))))
  })
  return(c(list(multiwig),track_lines))
  

}

generate_bigwig<-function(obj=immune,
                      suffix="immune",
                      groupBy="fine_cluster_phenograph",
                      outdir="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/03_fine_celltyping/",
                      step=500){
  bigwig_output_dir=paste0(outdir,"/bigwig_output_",suffix)
  system(paste("mkdir -p ", bigwig_output_dir))
  system(paste("mkdir -p ", paste0(bigwig_output_dir,"/","hg38")))

  obj@h5paths$path<-obj@h5paths$paths
  obj@h5paths$barcode<-row.names(obj@h5paths)

  #update to new smoothed windows
  celltype500bpwindows <- calcSmoothedWindows(obj, 
                                          type = "CG", 
                                          threads = 300,
                                          step = step, 
                                          smooth = 1,
                                          genome = "hg38",
                                          index = "chr_cg",
                                          groupBy = groupBy,
                                          returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                          returnPctMatrix = TRUE)
  saveRDS(celltype500bpwindows,file=paste0(bigwig_output_dir,"/","03.1.VMR_umap.",suffix,".fine_cluster.500bp_windows.rds"))
  obj@genomeMatrices[[paste0("cg_",suffix,"_cells_perc")]] <- celltype500bpwindows[["pct_matrix"]]
  #output tracks as bigwig
  groups<-colnames(celltype500bpwindows[["pct_matrix"]])[4:ncol(celltype500bpwindows[["pct_matrix"]])]
  groups<-groups[groups!="NA"]
  
  #make track hub.txt
  writeLines(c(
    paste(c("hub", groupBy),collapse=" "),
    paste(c("shortLabel", groupBy),collapse=" "),
    paste(c("longLabel", groupBy),collapse=" "),
    paste(c("genomesFile", "genomesFile.txt"),collapse=" "),
    paste(c("email ryan"),collapse=" ")),
  paste0(bigwig_output_dir,"/","hub.txt"))

  #make track genomesFile.txt
  writeLines(c(
    paste(c("genome", "hg38"),collapse=" "),
    paste(c("trackDb", "hg38/trackDb.txt"),collapse=" ")),
  paste0(bigwig_output_dir,"/","genomesFile.txt"))
  
  #write out bigwigs and trackhub data
  trackhub_dat<-lapply(groups,function(i) {
    write_binned_bigwigs(celltype_tracks=celltype500bpwindows[["pct_matrix"]],
                          outdir=paste0(bigwig_output_dir,"/","hg38"),
                          i=i,
                          step=step)})
  file_conn <- file(paste0(bigwig_output_dir,"/hg38/","trackDb.txt"), open = "a")
  for(i in trackhub_dat){writeLines(unlist(i), con = file_conn)}
  close(file_conn)

  return(obj)
}

```

## Run 500 or 250bp windows per celltype/group and per clone and lumhr

### Run per celltype x group
Running at 500bp and 250bp window sizes

```R
output_directory=paste0(wd,"/","dmr_celltype_group")
system(paste0("mkdir -p ",output_directory))

obj<-readRDS(file=paste(project_data_directory,"03_fine_celltyping","03_scaledcis.final_celltypes.amethyst.rds",sep="/"))

obj@metadata$celltype_group<-paste(obj@metadata$celltype,obj@metadata$Group,sep="_")

obj<-generate_bigwig(obj=obj,
                        suffix="dmr_celltype_group_500bp",
                        groupBy="celltype_group",
                        step=500,
                        outdir=getwd())

obj<-generate_bigwig(obj=obj,
                        suffix="dmr_celltype_group_250bp",
                        groupBy="celltype_group",
                        step=250,
                        outdir=getwd())                      
```


## Run 500 or 250bp windows per celltype/group and per clone and lumhr

### Run per clone (within sample only)
Running at 500bp and 250bp window sizes
Only run on cancer and lumhr

#Rerun
```R

obj<-readRDS(file=paste(project_data_directory,"03_fine_celltyping","03_scaledcis.final_celltypes.amethyst.rds",sep="/"))

dat<-subsetObject(obj,cells=row.names(obj@metadata[obj@metadata$celltype %in% c("cancer","lumhr"),]))

#add in lumhr as a class to also act as baseline
dat@metadata[(dat@metadata$Group=="HBCA") & (dat@metadata$celltype=="lumhr"),]$cnv_clonename_500kb<-"HBCA_lumhr"

table(dat@metadata$cnv_clonename_500kb,dat@metadata$celltype)

#note that "diploid" clones will only be lumhr/cancer cells because of initial cell type filters
dat<-subsetObject(dat,cells=row.names(dat@metadata[!(dat@metadata$cnv_clonename_500kb %in% c("NA")),]))

#REMOVED CLONE FILTER, BECAUSE THEY CAN STILL BE USEFUL AS MERGED SETS
#CAN FILTER FOR SMALLER COMPARISONS AS NEEDED
#require at least 50 cells for clones passing filter
#clones_passing_filter <-table(dat@metadata$cnv_clonename_500kb)[table#(dat@metadata$cnv_clonename_500kb)>50]
#note that bigwig output in the dropbox will still be 50 cell requirement (since its hard to tell otherwise)
table(dat@metadata$cnv_clonename_500kb)

#dat<-subsetObject(dat,cells=row.names(dat@metadata[dat@metadata$cnv_clonename_500kb %in% names#(clones_passing_filter),]))

dat<-generate_bigwig(obj=dat,
                        suffix="dmr_clones_500bp",
                        groupBy="cnv_clonename_500kb",
                        step=500,
                        outdir=getwd())

dat<-generate_bigwig(obj=dat,
                        suffix="dmr_clones_250bp",
                        groupBy="cnv_clonename_500kb",
                        step=250,
                        outdir=getwd()) 
```



<!--

### Run per clone (using HBCA lumHR as baseline)
Running at 500bp and 250bp window sizes

- set lumhr as the diploid baseline,
- remove other "diploid" clones for comparison (these include other celltypes)
#Rerun
```R

obj<-readRDS(file=paste(project_data_directory,"03_fine_celltyping","03_scaledcis.final_celltypes.amethyst.rds",sep="/"))
dat<-subsetObject(obj,cells=row.names(obj@metadata[obj@metadata$celltype %in% c("cancer","lumhr"),]))
dat<-subsetObject(dat,cells=row.names(dat@metadata[!(dat@metadata$cnv_clonename_500kb %in% c("NA")),]))

dat@metadata[(dat@metadata$Group=="HBCA") & (dat@metadata$celltype=="lumhr"),]$cnv_clonename_500kb<-"HBCA_lumhr"
clones_passing_filter<-table(dat@metadata$cnv_clonename_500kb)[table(dat@metadata$cnv_clonename_500kb)>50]
dat<-subsetObject(dat,cells=row.names(dat@metadata[dat@metadata$cnv_clonename_500kb %in% names(clones_passing_filter),]))

dat<-generate_bigwig(obj=dat,
                        suffix="dmr_clones_HBCALumHR_500bp",
                        groupBy="cnv_clonename_500kb",
                        step=500,
                        outdir=getwd())

dat<-generate_bigwig(obj=dat,
                        suffix="dmr_clones_HBCALumHR_250bp",
                        groupBy="cnv_clonename_500kb",
                        step=250,
                        outdir=getwd())     
```













<!--
## Now running groupwise comparisons across cell types

First running on 500bp windows
```R


obj<-readRDS(file=paste(project_data_directory,"03_fine_celltyping","03_scaledcis.final_celltypes.amethyst.rds",sep="/"))
dat<-amethyst::subsetObject(obj,cells=row.names(obj@metadata[!(obj@metadata$celltype %in% c("stromal_unknown")) & !is.na(obj@metadata$celltype),]))

dat@metadata$celltype_group<-paste(dat@metadata$celltype,dat@metadata$Group,sep="_")

input_folder=paste0(wd,"/","bigwig_output_dmr_celltype_group_500bp")
suffix="dmr_celltype_group_500bp"

output_directory=paste0(input_folder,"/","dmr_celltype_group")
system(paste0("mkdir -p ",output_directory))

celltype500bpwindows<-readRDS(file=paste0(input_folder,"/","03.1.VMR_umap.",suffix,".fine_cluster.500bp_windows.rds"))

pct_mat<-celltype500bpwindows[["pct_matrix"]] 
sum_mat<-celltype500bpwindows[["sum_matrix"]] 


#do 1v1 and 1vrest for all celltype comparisons possible
#1vrest
comparisons_1vrest<-do.call("rbind",lapply(unique(dat@metadata$celltype),function(celltype){
  comparisons_set<-colnames(pct_mat)[grepl(colnames(pct_mat),pattern=celltype)]
  comparisons_set<-do.call("rbind",lapply(comparisons_set,function(i){
    data.frame(
    name=paste0(i,"_v_",celltype,"_rest"),
    A=i,
    B=paste(comparisons_set[grep(comparisons_set,invert=T,pattern=i)],collapse=","))
  }))
  }))

#1v1
comparisons_1v1<-do.call("rbind",lapply(unique(dat@metadata$celltype),function(celltype){
  comparisons_set<-colnames(pct_mat)[grepl(colnames(pct_mat),pattern=celltype)]
  comparisons_set<-do.call("rbind",lapply(comparisons_set,function(i){
    j_list=comparisons_set[grep(comparisons_set,invert=T,pattern=i)]
    do.call("rbind",lapply(j_list,function(j){
      data.frame(
      name=paste0(i,"_v_",j),
      A=i,
      B=j
      )}))
}))}))

comparisons<-rbind(comparisons_1vrest,comparisons_1v1)
row.names(comparisons)<-comparisons$name

#calculate DMRs for group-wise comparisons across cell types
#run DMR analysis per row in comparisons
dmr_out<-mclapply(row.names(comparisons), function(i) {
  celltype <- unlist(lapply(strsplit(i,"_"),"[",1))[1]
  print(paste("Running DMRs across cell types for:",i))
  comparisons_output_directory<-paste0(output_directory,"/",celltype,"/",i)
  dmr_bed_output_directory<-paste0(output_directory,"/",celltype,"/dmr_bed")

  system(paste0("mkdir -p ",paste0(output_directory,"/",celltype)))
  system(paste0("mkdir -p ",comparisons_output_directory))
  system(paste0("mkdir -p ",dmr_bed_output_directory))
  dmrs <- testDMR(sum_mat, # Sum of c and t observations in each genomic window per group
          comparisons=comparisons[i,],
          nminTotal = 3, # Min number observations across all groups to include the region in calculations
          nminGroup = 3) # Min number observations across either members or nonmembers to include the region
  dmrs$comparison_name <- i
  dmrs$celltype <- celltype
  dmrs$groupA <-comparisons[i,]$A 
  dmrs$groupB <- comparisons[i,]$B
  saveRDS(dmrs,paste0(comparisons_output_directory,"/","04.1.",suffix,".",i,".dmr.rds"))

  print(paste("Filtering DMRs across clones for:",i))
  dmrs <- filterDMR(dmrs, 
              method = "BH", #bonferroni is too strict, lose a lot of significant results
              filter = FALSE, # If TRUE, removes insignificant results
              pThreshold = 0.05, # Maxmimum adjusted p value to allow if filter = TRUE
              logThreshold=0) # Minimum absolute value of the log2FC to allow if filter = TRUE
  dmrs<-dmrs %>% filter(padj<0.05)
  dmrs$comparison_name <- i
  dmrs <- collapseDMR(obj, 
                          dmrs, 
                          maxDist = 0, # Max allowable overlap between DMRs to be considered adjacent
                          minLength = 100, # Min length of collapsed DMR window to include in the output
                          reduce = T, # Reduce results to unique observations (recommended)
                          annotate = T) # Add column with overlapping gene names
  dmrs$comparison_name <- i
  dmrs$celltype <- unlist(lapply(strsplit(i,"_"),"[",1))[1]
  dmrs$groupA <-comparisons[i,]$A 
  dmrs$groupB <- comparisons[i,]$B

  saveRDS(dmrs,paste0(comparisons_output_directory,"/","04.1.",suffix,".",i,".dmr.filt.rds"))

  #output DMRs as bed files
  hypo_dmrs<-dmrs %>% filter(direction=="hypo") %>% select(chr,dmr_start,dmr_end) %>% as.data.frame() %>% dplyr::rename(dmr_start=start,dmr_end=end) %>% GRanges()
  export(hypo_dmrs, con = paste0(dmr_bed_output_directory,"/","04.1.",suffix,".",i,".dmr.filt.hypo.bed"), format = "bed")

  hypo_dmrs<-dmrs %>% filter(direction=="hyper") %>% select(chr,dmr_start,dmr_end) %>% as.data.frame() %>% dplyr::rename(dmr_start=start,dmr_end=end) %>% GRanges()
  export(hyper_dmrs, con = paste0(dmr_bed_output_directory,"/","04.1.",suffix,".",i,".dmr.filt.hyper.bed"), format = "bed")

  return(dmrs)
  },mc.cores=100)


dmr_out<-do.call("rbind",dmr_out)


```

Performing same comparisons using RNA data
```R

#compare dmr sites with rna marker genes overlap
library(Seurat)
library(GenomicRanges)


rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.pf.rds")
rna$celltype_group<-paste(rna$coarse_celltype,rna$Group,sep="_")
rna<-JoinLayers(rna)
input_folder=paste0(wd,"/","bigwig_output_dmr_celltype_group_500bp")
output_directory=paste0(input_folder,"/","dmr_celltype_group")

de_out <- lapply(row.names(comparisons), function(i) {
  celltype <- unlist(lapply(strsplit(i,"_"),"[",1))[1]
  print(paste("Running DE genes across cell types for:",i))
  comparisons_output_directory<-paste0(output_directory,"/",celltype,"/",i)

  #run paired comparisons that DMR underwent
  cell_A<-row.names(rna@meta.data)[rna@meta.data$celltype_group %in% unlist(strsplit(comparisons[i,]$A,","))]
  cell_B<-row.names(rna@meta.data)[rna@meta.data$celltype_group %in% unlist(strsplit(comparisons[i,]$B,","))]
  if(length(cell_A)>50 & length(cell_B)>50){
    print(table(rna@meta.data[cell_A,]$celltype_group))
    print(table(rna@meta.data[cell_B,]$celltype_group))
    rna_markers<-FindMarkers(rna@assays$RNA,cells.1=cell_A,cells.2=cell_B,only.pos=TRUE)
    rna_markers$comparison_name<-comparisons[i,]$name
    rna_markers$celltype<-celltype
    rna_markers$groupA<-comparisons[i,]$A
    rna_markers$groupB<-comparisons[i,]$B
    rna_markers$gene<-row.names(rna_markers)
    saveRDS(rna_markers,file=paste0(comparisons_output_directory,"/","04.1.",suffix,".",i,".de.RNAmarkers.rds"))
    return(rna_markers)
  }
})

de_out<-do.call("rbind",de_out)
#also output bed file of significant RNA for plotting

```

Plotting of DMR count per comparison

```R

dmr_out<-do.call("rbind",lapply(list.files(path=output_directory,recursive=T,full.names=T,pattern="*.dmr.filt.rds"),readRDS))
de_out<-do.call("rbind",lapply(list.files(path=output_directory,recursive=T,full.names=T,pattern="*.de.RNAmarkers.rds"),readRDS))

#color assignment is fluor as cancer associated cell type, rest muted versions
celltype_col=c(
"pericyte"="#FF9900",
"fibroblast"="#FF0000",
"endothelial"="#FFFF66",
"myeloid"="#99FFFF",
"bcell"="#0099CC",
"tcell"="#99FF99",
"basal"="#990099",
"lumsec"="#CC0066",
"lumhr"="#FF00CC",
"cancer"="#00FF99")

plot_list<-mclapply(c("basal","fibroblast","endothelial","tcell","myeloid"),function(i){
  print(paste("Plotting...",i))
  celltype_to_plot=i
  comparison_order_for_plotting=c(
    paste(celltype_to_plot,"HBCA","v",celltype_to_plot,"DCIS",sep="_"),
    paste(celltype_to_plot,"HBCA","v",celltype_to_plot,"Synchronous",sep="_"),
    paste(celltype_to_plot,"HBCA","v",celltype_to_plot,"IDC",sep="_"),
    paste(celltype_to_plot,"DCIS","v",celltype_to_plot,"Synchronous",sep="_"),
    paste(celltype_to_plot,"DCIS","v",celltype_to_plot,"IDC",sep="_"),
    paste(celltype_to_plot,"Synchronous","v",celltype_to_plot,"IDC",sep="_"),

    paste(celltype_to_plot,"HBCA","v",celltype_to_plot,"rest",sep="_"),
    paste(celltype_to_plot,"DCIS","v",celltype_to_plot,"rest",sep="_"),
    paste(celltype_to_plot,"Synchronous","v",celltype_to_plot,"rest",sep="_"),
    paste(celltype_to_plot,"IDC","v",celltype_to_plot,"rest",sep="_")
  )

  test_dat<-dmr_out %>% 
            filter(celltype==celltype_to_plot) %>% 
            filter(dmr_padj<0.05) %>% 
            filter(comparison_name %in% comparison_order_for_plotting) %>% 
            group_by(comparison_name,direction,.drop=FALSE) %>% 
            summarize(count=n(),dmr_kbp=sum(dmr_length)/1000)

  test_dat_rna<-de_out %>% 
            filter(celltype==celltype_to_plot) %>% 
            filter(p_val_adj<0.05) %>% filter(avg_log2FC>0) %>%
            filter(comparison_name %in% comparison_order_for_plotting) %>% 
            group_by(comparison_name,.drop=FALSE) %>% 
            summarize(de_count=n())


  #mutate to one data frame
  test_dat$comparison_name<-factor(test_dat$comparison_name,levels=comparison_order_for_plotting)
  test_dat$direction<-factor(test_dat$direction,levels=c("hyper","hypo"))
  test_dat_rna$comparison_name<-factor(test_dat_rna$comparison_name,levels=comparison_order_for_plotting)

  test_dat_rna<-tidyr::complete(test_dat_rna,comparison_name,fill = list(de_count = 0))
  #test_dat<-tidyr::complete(test_dat,comparison_name,direction,fill = list(dmr_kbp = 0,count=0))

  plot_dat<-data.frame(comparison_name=comparison_order_for_plotting,
  hypo_dmr_count=test_dat %>% filter(direction=="hypo") %>% arrange(comparison_name) %>% pull(count),
  hypo_dmr_length=test_dat %>% filter(direction=="hypo") %>% arrange(comparison_name) %>% pull(dmr_kbp),
  hyper_dmr_count=test_dat %>% filter(direction=="hyper") %>% arrange(comparison_name) %>% pull(count),
  hyper_dmr_length=test_dat %>% filter(direction=="hyper") %>% arrange(comparison_name) %>% pull(dmr_kbp),
  de_count=test_dat_rna %>% arrange(comparison_name) %>% pull(de_count))
  
  plot_dat$comparison_name <- factor(plot_dat$comparison_name ,levels=rev(comparison_order_for_plotting))

  #plot count and size of DMRs and DE
  #if comparison_name is 0 make sure it is still plotted

  plt1<-ggplot(plot_dat ,aes(x=comparison_name,y=hyper_dmr_count,fill=I(unname(celltype_col[celltype_to_plot]))))+
        geom_bar(stat="identity") + 
        geom_text(aes(label=paste0(hyper_dmr_length,"kbp")), vjust=0,color="yellow") +
        theme_minimal() + 
        coord_flip()

  plt2<-ggplot(plot_dat,aes(x=comparison_name,y=hypo_dmr_count,color=I(unname(celltype_col[celltype_to_plot]))))+
        geom_bar(stat="identity",fill="white") + 
        geom_text(aes(label=paste0(hypo_dmr_length,"kbp")), vjust=0.1,color="red") +
        theme_minimal() + 
        coord_flip()

  plt3<-ggplot(plot_dat,aes(x=comparison_name,y=de_count,color=I(unname(celltype_col[celltype_to_plot]))))+
        geom_bar(stat="identity") + 
        geom_text(aes(label=de_count), vjust=0.1,color="red") +
        theme_minimal() + 
        coord_flip()
  system(paste0("mkdir -p ",paste0(output_directory,"/","dmr_plots/")))

  #calculate average distance from DMR site to DE gene
  
  gene_locations<-obj@ref %>% filter(type=="gene") %>% filter(!duplicated(gene_name))

  met_overlap<-do.call("rbind",lapply(comparison_order_for_plotting,function(j){
    print(paste("Overlapping...",j))

    met_markers_hypo<-dmr_out %>% 
              filter(celltype==celltype_to_plot) %>% 
              filter(dmr_padj<0.05) %>% 
              filter(direction=="hypo") %>%
              filter(comparison_name == j) %>% GRanges()

    met_markers_hyper<-dmr_out %>% 
              filter(celltype==celltype_to_plot) %>% 
              filter(dmr_padj<0.05) %>% 
              filter(direction=="hyper") %>%
              filter(comparison_name == j) %>% GRanges()

    rna_markers<- de_out %>% 
              filter(celltype==celltype_to_plot) %>% 
              filter(p_val_adj<0.05) %>% filter(avg_log2FC>0) %>%
              filter(comparison_name == j) %>% mutate(gene_name=gene)
    
    rna_markers<-left_join(rna_markers, gene_locations, by = "gene_name") %>% filter(!is.na(seqid)) %>% GRanges()

    #count DMR distance to each gene marker
    met_distance_hypo<-distanceToNearest(met_markers_hypo,rna_markers)
    met_distance_hyper<-distanceToNearest(met_markers_hyper,rna_markers)

    summary_frame<-data.frame(
        celltype=celltype_to_plot,
        comparison=j,
        hypo_dmr_count=length(met_markers_hypo),
        hyper_dmr_count=length(met_markers_hyper),
        de_count=length(rna_markers),
        hypo_less_than_5kbp_from_marker_gene=met_distance_hypo %>% as.data.frame() %>% filter(distance<5000) %>% pull(subjectHits) %>% unique() %>% length(),
        hypo_overlapping_marker_gene=met_distance_hypo %>% as.data.frame() %>% filter(distance==0) %>% pull(subjectHits) %>% unique() %>% length(),
        hyper_less_than_5kbp_from_marker_gene=met_distance_hyper %>% as.data.frame() %>% filter(distance<5000) %>% pull(subjectHits) %>% unique() %>% length(),
        hyper_overlapping_marker_gene=met_distance_hyper %>% as.data.frame() %>% filter(distance==0) %>% pull(subjectHits) %>% unique() %>% length())
    return(summary_frame)}))

  met_overlap$comparison<-factor(met_overlap$comparison,levels=rev(comparison_order_for_plotting))
  met_overlap$unexplained_by_met<-met_overlap$de_count-(met_overlap$hypo_less_than_5kbp_from_marker_gene+met_overlap$hyper_less_than_5kbp_from_marker_gene)
  
  met_overlap <- met_overlap %>% 
                  select(comparison,hypo_less_than_5kbp_from_marker_gene,hyper_less_than_5kbp_from_marker_gene,unexplained_by_met) %>%
                  tidyr::pivot_longer(cols=c(hypo_less_than_5kbp_from_marker_gene, hyper_less_than_5kbp_from_marker_gene,unexplained_by_met),names_to="group",values_to="count")
  met_overlap$group<-factor(met_overlap$group,levels=c("hypo_less_than_5kbp_from_marker_gene","hyper_less_than_5kbp_from_marker_gene","unexplained_by_met"))
  plt<-ggplot(met_overlap,aes(x=comparison,y=count,fill=group))+
  geom_bar(position="stack", stat="identity")+
  theme_minimal()+
  coord_flip() 
  
  plot<-plt1|plt2|plt
  return(plot)
  #run GSEA on cell type
  #system(paste0("mkdir -p ",output_directory,"/gsea_data"))
  #for(j in comparison_order_for_plotting){
  #  gsea_across_sets(obj, 
  #                  collapsed_dmrs=dmr_out %>% filter(comparison_name == j),
  #                  sample_name=j, 
  #                  prefix=paste0(output_directory,"/gsea_data/","04.1.",suffix,".",celltype_to_plot))
  #}

  #tft_gsea<-do.call("rbind",lapply(
  #  list.files(path=paste0(output_directory,"/gsea_data/"),pattern=paste0("04.1.",suffix,".",celltype_to_plot,".*","TFT.rds"),full.names=T),
  #  readRDS))
  #plot_gsea(gsea=tft_gsea,out_setname="TFT",prefix=paste0(output_directory,"/dmr_plots/","04.1.",suffix,".",celltype_to_plot))

  #hallmark_gsea<-do.call("rbind",lapply(
  #  list.files(path=paste0(output_directory,"/gsea_data/"),pattern=paste0("04.1.",suffix,".",celltype_to_plot,".*","hallmark.rds"),full.names=T),
  #  readRDS))
  #plot_gsea(gsea=hallmark_gsea,out_setname="hallmark",prefix=paste0(output_directory,"/dmr_plots/","04.1.",suffix,".",celltype_to_plot))

  #position_gsea<-do.call("rbind",lapply(
  #  list.files(path=paste0(output_directory,"/gsea_data/"),pattern=paste0("04.1.",suffix,".",celltype_to_plot,".*","position.rds"),full.names=T),
  #  readRDS))
  #plot_gsea(gsea=position_gsea,out_setname="position",prefix=paste0(output_directory,"/dmr_plots/","04.1.",suffix,".",celltype_to_plot))

  #canceratlas_gsea<-do.call("rbind",lapply(
  #  list.files(path=paste0(output_directory,"/gsea_data/"),pattern=paste0("04.1.",suffix,".",celltype_to_plot,".*","3CA.rds"),full.names=T),
  #  readRDS))
  #plot_gsea(gsea=canceratlas_gsea,out_setname="3CA",prefix=paste0(output_directory,"/dmr_plots/","04.1.",suffix,".",celltype_to_plot))
  },mc.cores=50)


plt<-wrap_plots(plot_list,ncol=1) + plot_layout(axes = "collect")
  ggsave(plt,file=paste0(output_directory,"/","dmr_plots/","04.1.",suffix,".allcells.dmr.counts.pdf"),height=20,width=20,units="in")

#calculate DMR proximity to DE genes (or overlap with DE genes)
#DMR overlap with breast cancer gene list
#GSEA enrichment
#overlap with bulk methylation progression markers

#Run chromvar on hypomethylated 500bp windows
#with and without excluding CGI
```

Overlap previous markers of progression with DMRs

```R
library(stringr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(readxl)

input_folder=paste0(wd,"/","bigwig_output_dmr_celltype_group_500bp")
output_directory=paste0(input_folder,"/","dmr_celltype_group")

dmr_out<-do.call("rbind",lapply(list.files(path=output_directory,recursive=T,full.names=T,pattern="*.dmr.filt.rds"),readRDS))
cg_v1<-read.csv(file="/data/rmulqueen/projects/scalebio_dcis/ref/infinium-methylationepic-v-1-0-b5-manifest-file.csv",skip=7,sep=",")

fleischer_2014<-c("cg05809947",
"cg12219311",
"cg26466505",
"cg20691428",
"cg20869305",
"cg22174844",
"cg26225829",
"cg04947065",
"cg16575694",
"cg16559598",
"cg08729004",
"cg13635578",
"cg13744452",
"cg04817034",
"cg07130508",
"cg00226265",
"cg13749939",
"cg25817165")

genes<-cg_v1 %>% filter(IlmnID %in% fleischer_2014) %>% pull("GencodeBasicV12_NAME")
cg_loci<-cg_v1 %>% filter(IlmnID %in% fleischer_2014) %>% select(CHR_hg38,Start_hg38,End_hg38,IlmnID) %>% mutate(chr=CHR_hg38,start=Start_hg38,end=End_hg38) %>% GRanges()

dmr_overlap<-findOverlaps(dmr_out %>% GRanges(),cg_loci)

dmr_out_filtered<-dmr_out[dmr_overlap@from,]
dmr_out_filtered$fleischer_2014_cgprobe<-cg_loci[dmr_overlap@to,]$IlmnID

dmr_fleischer_2014_probe_overlap<-table(dmr_out_filtered$fleischer_2014_cgprobe,dmr_out_filtered$celltype) %>% as.data.frame()
dmr_fleischer_2014_probe_overlap$Freq<-as.integer(as.logical(dmr_fleischer_2014_probe_overlap$Freq>0))
dmr_fleischer_2014_probe_overlap<-tidyr::pivot_wider(dmr_fleischer_2014_probe_overlap,names_from=Var1,values_from=Freq)%>% as.data.frame()
row.names(dmr_fleischer_2014_probe_overlap)<-dmr_fleischer_2014_probe_overlap[,1]
dmr_fleischer_2014_probe_overlap<-dmr_fleischer_2014_probe_overlap[,2:ncol(dmr_fleischer_2014_probe_overlap)]

plt<-Heatmap(dmr_fleischer_2014_probe_overlap)
pdf(file="test_dmr_fleischer_2014_probe_overlap.pdf")
print(plt)
dev.off()

#now johnson 2015
johnson_2015<-read_excel("/data/rmulqueen/projects/scalebio_dcis/ref/13148_2015_94_MOESM4_ESM.xls",skip=1,sheet=1)
cg_v1<-read.csv(file="/data/rmulqueen/projects/scalebio_dcis/ref/infinium-methylationepic-v-1-0-b5-manifest-file.csv",skip=7,sep=",")
cg_loci<-cg_v1 %>% filter(IlmnID %in% johnson_2015$`Illumina CG ID`) %>% select(CHR_hg38,Start_hg38,End_hg38,IlmnID) %>% mutate(chr=CHR_hg38,start=Start_hg38,end=End_hg38) %>% GRanges()

dmr_overlap<-findOverlaps(dmr_out %>% GRanges(),cg_loci)

dmr_out_filtered<-dmr_out[dmr_overlap@from,]
dmr_out_filtered$johnson_2015_cgprobe<-cg_loci[dmr_overlap@to,]$IlmnID

dmr_johnson_2015_probe_overlap<-table(dmr_out_filtered$johnson_2015_cgprobe,dmr_out_filtered$celltype) %>% as.data.frame()
dmr_johnson_2015_probe_overlap$Freq<-as.integer(as.logical(dmr_johnson_2015_probe_overlap$Freq>0))
dmr_johnson_2015_probe_overlap<-tidyr::pivot_wider(dmr_johnson_2015_probe_overlap,names_from=Var1,values_from=Freq) %>% as.data.frame()
row.names(dmr_johnson_2015_probe_overlap)<-dmr_johnson_2015_probe_overlap[,1]
dmr_johnson_2015_probe_overlap<-dmr_johnson_2015_probe_overlap[,2:ncol(dmr_johnson_2015_probe_overlap)]
plt<-Heatmap(dmr_johnson_2015_probe_overlap,cluster_rows=T,cluster_columns=T)
pdf(file="test_dmr_johnson_2015_probe_overlap.pdf")
print(plt)
dev.off()
```


```R
library(chromVAR)
library(motifmatchr)
library(Matrix)q
library(SummarizedExperiment)
library(BiocParallel)
set.seed(2017)
library(JASPAR2024)
library(TFBSTools)
register(MulticoreParam(50, progressbar = TRUE))
library(BSgenome.Hsapiens.UCSC.hg38)
task_cpus=300

input_folder=paste0(wd,"/","bigwig_output_dmr_celltype_group_500bp")
suffix="dmr_celltype_group_500bp"
output_directory=paste0(input_folder,"/","dmr_celltype_group")
system(paste0("mkdir -p ",output_directory))

celltype500bpwindows<-readRDS(file=paste0(input_folder,"/","03.1.VMR_umap.",suffix,".fine_cluster.500bp_windows.rds"))

pct_mat<-celltype500bpwindows[["pct_matrix"]] 
pct_mat<-pct_mat[!duplicated(paste(pct_mat$chr,pct_mat$start,pct_mat$end,sep="_")),]
#set up binary counts matrix per window
binary_met<-pct_mat[,4:ncol(pct_mat)]
binary_met<-(binary_met>20)+0 #hypomethylated regions (<20% met) are 1, rest 0 or NA

windows_passing_filter<-which(rowSums(is.na(binary_met),na.rm=T)>1) #require at least 1 groups to have a hypomethylated region
binary_met<-binary_met[windows_passing_filter,]
ranges<-GRanges(pct_mat[windows_passing_filter,1:3])

binary_met[which(is.na(binary_met),arr.ind=T)]<-0 #set remaining NA to 0

#set up hypomet matrix
binary_met <- SummarizedExperiment(assays = list(counts = binary_met),rowRanges = sort(ranges))
ranges<-sort(ranges)

binary_met <- addGCBias(binary_met, genome = BSgenome.Hsapiens.UCSC.hg38)


#compute deviations for peaks
binary_met <- filterPeaks(binary_met,non_overlapping=FALSE)
bg <- getBackgroundPeaks(object = binary_met)


#set up motifs location calls across ranges
motifs <- getMatrixSet(x = JASPAR2024()@db, opts = list(collection = "CORE", species = "Homo sapiens"))

ranges_list<-split(binary_met@rowRanges, ntile(seq_along(binary_met@rowRanges),task_cpus))
motif_ix <- mclapply(1:length(ranges_list),function(i){
                        matchMotifs(motifs, ranges_list[[i]], genome = BSgenome.Hsapiens.UCSC.hg38)
                        },mc.cores=task_cpus)
motif_ix<-do.call("rbind",motif_ix)

#run chromvar
dev <- computeDeviations(object = binary_met, annotations = motif_ix, background_peaks = bg)
variability <- computeVariability(dev)

pdf("test.chromvar.variability.pdf")
plotVariability(variability, use_plotly = FALSE) 
dev.off()



```
BCMDCIS RNA difference

```r
rna@meta.data[rna@meta.data$celltype %in% c("lumhr_HBCA"),]$rna_clonename<-"lumhr_HBCA"
Idents(rna)<-rna$rna_clonename

marker_plotter<-function(ident.1="lumhr_HBCA",ident.2="BCMDCIS41T_c1"){
  markers<-FindMarkers(rna,ident.1=ident.1,ident.2=ident.2)

  # Define significance thresholds
  log2FC_threshold <- 1
  pval_threshold <- 0.05

  # Categorize genes into Significant vs Not Significant
  de_markers <- markers %>%
    mutate(Gene = rownames(markers)) %>%
    mutate(Significance = case_when(
      avg_log2FC > log2FC_threshold & p_val_adj < pval_threshold ~ "Up-regulated",
      avg_log2FC < -log2FC_threshold & p_val_adj < pval_threshold ~ "Down-regulated",
      TRUE ~ "Not Significant"
    ))

  # Generate the Volcano Plot
  plt<-ggplot(de_markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up-regulated" = "ident.color1", 
                                  "Down-regulated" = "blue", 
                                  "Not Significant" = "grey")) +
    theme_minimal() +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    labs(title = "Volcano Plot: Condition A vs Condition B",
        x = "Average Log2 Fold Change",
        y = "-Log10 Adjusted P-value")
  return(plt)
}

plt1<-marker_plotter(ident.1="lumhr_HBCA",ident.2="BCMDCIS41T_c3")

plt2<-marker_plotter(ident.1="BCMDCIS41T_c3",ident.2="BCMDCIS41T_c2")

plt3<-marker_plotter(ident.1="BCMDCIS41T_c2",ident.2="BCMDCIS41T_c1")

ggsave(wrap_plots(plt1,plt2,plt3,ncol=3),file="test.rna_dcis41t_volcano.pdf",width=30)

```
```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
```


```R
source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
library(Seurat)
library(ComplexHeatmap)

task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)

system(paste0("mkdir -p ",project_data_directory,"/fine_celltyping"))

```

#Subcluster to assess more specific cell types

```R
obj<-readRDS(file="04_scaledcis.amethyst.broad_celltype.rds")

#prepare cgi
cgisland="/data/rmulqueen/projects/scalebio_dcis/ref/cpgIslandExt.bed"
cgi<-rtracklayer::import(cgisland)
cgi<-as.data.frame(cgi)
colnames(cgi)<-c("chr","start","end","strand")

plot_histogram_page<-function(celltype){
    plt<-histograModified(obj, 
        baseline="mean",
        genes = unlist(cell_markers[celltype]),
        colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
        matrix = "cg_cluster_id_tracks", arrowScale = .03, trackScale = .5,
        legend = F, cgisland=cgi)
    return(plt)
}

fine_celltyping_cluster<-function(dat=obj,broad_celltype="stromal",d=20,neigh=50,q=0.01){
    output_directory=paste0(project_data_directory,"/fine_celltyping/",broad_celltype)
    system(paste0("mkdir -p ",output_directory))

    print(paste("Subsetting object based on broad celltype:",broad_celltype))
    dat<-subsetObject(dat,cells=row.names(dat@metadata[dat@metadata$broad_celltype==broad_celltype,]))
    #overcluster by existing DMR
    print(paste("Reclustering on only",broad_celltype,"cells"))
    dat<-cluster_by_windows(
        dat, window_name="dmr_sites",
        outname=paste(broad_celltype,"broad_dmr_sites_clustering",sep="."),
        threads=task_cpus, est_dim=d, neighbors=neigh, dist=q,
        k_pheno=neigh)
    dat@metadata[[paste(broad_celltype,"clusters",sep=".")]]<-dat@metadata$cluster_id  #set celltype clusters
    table(dat@metadata[[paste(broad_celltype,"clusters",sep=".")]])
    #find broad celltype specific DMRs
    print(paste("Identifying DMRs on",broad_celltype,"cells"))
    dat<-dmr_and_1kb_window_gen(dat,
        prefix=broad_celltype,
        groupBy=paste(broad_celltype,"clusters",sep="."),
        threads=10,step=1000)
    collapsed_dmrs<-readRDS(file=paste0(broad_celltype,".dmr.",paste(broad_celltype,"clusters",sep="."),".collapsed.rds"))
    #recluster on broad_celltype DMRs
    print(paste("Recluster",broad_celltype,"cells on cell type DMRs"))
    dmr_bed<-as.data.frame(cbind(collapsed_dmrs$chr,
                            collapsed_dmrs$dmr_start,
                            collapsed_dmrs$dmr_end))
    colnames(dmr_bed)<-c("chr","start","end")
    dmr_bed$start<-as.numeric(dmr_bed$start)
    dmr_bed$end<-as.numeric(dmr_bed$end)

    dat@genomeMatrices[["fine_dmr_sites"]] <- makeWindows(dat, 
                                                     bed = dmr_bed,
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 20, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 
    #cluster on defined DMRs
    dat<-cluster_by_windows(
        dat, window_name="fine_dmr_sites",
        outname=paste(broad_celltype,"broad_dmr_sites_clustering",sep="."),
        threads=task_cpus, est_dim=d, neighbors=neigh, dist=q,
        k_pheno=neigh)
    #overwrite temporary clusters to celltype-dmr defined clusters
    dat@metadata[[paste(broad_celltype,"clusters",sep=".")]]<-dat@metadata$cluster_id  #set celltype clusters

    print(paste("Make 500bp windows on",broad_celltype,"clusters"))
    #summarize celltype-dmr clusters over windows
    dat<-dmr_and_1kb_window_gen(dat,
        prefix=broad_celltype,
        groupBy=paste(broad_celltype,"clusters",sep="."),
        threads=10,step=500)

    print(paste("Saving final amethyst object of",broad_celltype))
    saveRDS(dat,paste0(output_directory,"/",paste("04_scaledcis",broad_celltype,"fine_celltype.amethyst.rds",sep=".")))
}

plot_dmr_RNAmarker_overlap<-function(dmr,rna_markers,celltype){
    #flow through:
    # - select DMRs for plotting that match with RNA overlap
    # - filter to hypo dmr (best to quickly visualize)
    # - filter to sites with RNA marker overlap
    # - table dmr to marker count
    # - plot gene body methylation hypoM plots

    #filter to hypo DMRs (more likely to increase expression)
    dmr<-dmr %>% group_by(test) %>% filter(direction=="hypo") 

    #dmr gene names is a nested list (dmr associated with multiple genes, so unlist them) and remove white space
    dmr_gene_names<-unlist(
        lapply(
            unlist(
                strsplit(unlist(dmr$gene_names),",")),
                function(x) gsub("[[:space:]]", "", x)
            )
        )

    #define ideal markers for methylation from the scRNA data
    markers$gene_name<-row.names(markers)
    markers<-markers %>% filter(gene_name %in% dmr_gene_names)

    #find which hypoM overlaps markers to make counts table
    dmr_marker_counts<-lapply(unique(dmr$test), 
        function(x){
        dmr_subset<-dmr %>% filter(test==x)

        dmr_subset_gene_names<-unlist(
        lapply(
            unlist(
                strsplit(unlist(dmr_subset$gene_names),",")),
                function(x) gsub("[[:space:]]", "", x)
            )
        )

        dmr_marker_overlap_count<-markers %>% 
            group_by(cluster,.drop=FALSE) %>% 
            filter(gene_name %in% dmr_subset_gene_names) %>% 
            tally() %>% 
            column_to_rownames(var="cluster") %>% 
            select(n)
        return(dmr_marker_overlap_count)
    })
    marker_overlap<-do.call("cbind",dmr_marker_counts)
    colnames(marker_overlap)<-unique(dmr$test)

    #normalize by dmrs per cluster
    mat_percent <- sweep(marker_overlap, 2, colSums(marker_overlap), FUN = "/")

    #normalize by markers by celltype
    marker_count<-unlist(markers %>% group_by(cluster,.drop=FALSE) %>% tally() %>% select(n))
    mat_percent <- sweep(mat_percent, 1, marker_count, FUN = "/")
    mat_percent<-mat_percent[complete.cases(mat_percent),]

    #make percentage by normalized sums
    mat_percent <- sweep(mat_percent, 2, colSums(mat_percent), FUN = "/") * 100

    print("Plotting percentage overlap of DMRs...")
    pdf(paste(celltype,"dmr_marker_overlap.heatmap.pdf",sep="."))
    print(draw(Heatmap(mat_percent)))
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
    filter(gene_length<25000)

    print(paste(nrow(rna_markers_merge),"RNA markers remaining..."))
    #so now we have a list of RNA markers that are set to size thresholds
    #pick based on significant DMRs (across all clusters) per fine_celltype rna marker

    dmr<-dmr %>% filter(direction=="hypo" & dmr_logFC< -2 & dmr_length > 2000)
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
    
    plotting_list
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



```

Stromal
```R



celltype="stromal"
fine_celltyping_cluster(dat=obj,
                        broad_celltype=celltype,
                        d=20,
                        neigh=50,
                        q=0.01,
                        cell_markers)

#read in subset amethyst object and dmrs
dat_fine<-readRDS(
    paste(project_data_directory,"fine_celltyping",celltype,
    paste("04_scaledcis",celltype,"fine_celltype.amethyst.rds",sep="."),
    sep="/"))

#using apriori list
cell_markers<-list()
cell_markers[["fibro"]]<-c("DCN","APOD","LUM","COL1A2","COL1A1","FAP")
cell_markers[["endo"]]=c("CCL21","TFF3","MMRN1","CLDN5","AL357507.1","PKHD1L1","KLHL4","LINC02147","RHOJ","ST6GALNAC3","MMRN1","MECOM","BTNL9","MCTP1","PTPRB","VWF","ADGRL4","LDB2")
cell_markers[["perivasc"]]=c("RGS6","KCNAB1","COL25A1","ADGRL3","PRKG1","MYL9","ADIRF","NR2F2-AS1","AC012409.2")
cell_markers[["adipo"]]=c("PDE3B","ACACB","WDPCP","PCDH9","CLSTN2","ADIPOQ","TRHDE")

plt_list<-lapply(names(cell_markers),
    function(celltype){
            genes<-unlist(cell_markers[celltype])
            print(paste("Plotting:",celltype))
            print(paste("Genes to plot:",genes))
            plt<-histograModified(obj=dat_fine, 
                baseline="mean",
                genes = unlist(cell_markers[celltype]),
                order=order,
                matrix = paste0("cg_",broad_celltype,".clusters_tracks"), arrowScale = .03, trackScale = .5,
                legend = F, cgisland=cgi) + ggtitle(celltype)
            ggsave(plt,
            file=paste(broad_celltype,"clusters",celltype,"fine_celltype.markers.pdf",sep="."),
            width=length(genes)*5,
            height=length(met_clusters)*20,
            limitsize=F)
        })

#dmr set
dmr<-readRDS(
        paste(project_data_directory,"fine_celltyping",celltype,
        paste(celltype,"dmr",celltype,"clusters.collapsed.rds",sep="."),
        sep="/"))

#rna markers defined by paired 10x RNA analysis
rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.rds")
table(rna$broad_celltype)
rna<-subset(rna,broad_celltype %in% c("fibro","lymphatic","perivasc","vascular"))
rna<-subset(rna,fine_celltype %in% c("suspected_doublet"),invert=TRUE)
table(rna$fine_celltype)

#trying to match the correct level of cell typing
rna@meta.data[rna@meta.data$fine_celltype %in% c("VSMC","peri","periVSMC_unknown"),]$broad_celltype<-"peri"
rna@meta.data[rna@meta.data$fine_celltype %in% c("endo_artery","endo_capillary","endo_lymphatic","endo_unknown","endo_vein"),]$broad_celltype<-"endo"
rna@meta.data[rna@meta.data$fine_celltype %in% c("fibro_SFRP4","fibro_major","fibro_matrix","fibro_prematrix","endo_vein"),]$broad_celltype<-"fibro"
rna@meta.data[rna@meta.data$fine_celltype %in% c("fibro_CAF"),]$broad_celltype<-"CAF"
rna@meta.data[rna@meta.data$fine_celltype %in% c("endo_TEC"),]$broad_celltype<-"TEC"

rna<-JoinLayers(rna)
markers<-FindAllMarkers(rna,group.by="broad_celltype")

#plot overlap of RNA markers for fine cell types and our cluster hypo DMRs
plot_dmr_RNAmarker_overlap(dmr=dmr,
                            rna_markers=markers,
                            celltype=celltype)

#plot RNA defined markers that overlap with large DMRs
plot_dmr_RNAMarker_hypoM(dat=dat_fine,
                            dmr=dmr,
                            rna_markers=markers,
                            broad_celltype=celltype,
                            order=c("1","4","2","7","5","9","3","8","6")) #loose ordering by umap grouping

#ASSIGN CELL TYPES BY PLOT
dat_fine@metadata$fine_celltypes<-"unknown"
dat_fine@metadata[dat_fine@metadata[[paste(celltype,"clusters",sep=".")]] %in% c("0"),]$fine_celltypes<-"unknown"

#SAVE FINE CELL TYPES
saveRDS(dat_fine,paste("04_scaledcis",celltype,"fine_celltype.amethyst.rds",sep="."))

```


Immune
Myeloid/B Cells/T cells

```R

celltype="immune"
fine_celltyping_cluster(dat=obj,broad_celltype=celltype,d=20,neigh=50,q=0.01,cell_markers)

#read in subset amethyst object and dmrs
dat_fine<-readRDS(
    paste(project_data_directory,"fine_celltyping",celltype,
    paste("04_scaledcis",celltype,"fine_celltype.amethyst.rds",sep="."),
    sep="/"))

cell_markers<-list()
cell_markers[["myeloid"]]<-c("HLA-DRA","HLA-DPA1","CD74")
cell_markers[["tcell"]]<-c("PTPRC","IKZF1","IL7R","GNLY")
cell_markers[["mast"]]<-c("NTM","SYTL3","SLC24A3","TPSB2","HDC")
cell_markers[["bcell"]]=c("CD37","TCL1A","LTB","HLA-DPB1","HLA-DRA","HLA-DPA1")
cell_markers[["plasma"]]=c("IGHA2","IGHA1","JCHAIN","IGHM","IGHG1","IGHG4","IGHG3","IGHG2")

broad_celltype="immune"
met_clusters=unique(dat_fine@metadata$cluster)
plt_list<-lapply(names(cell_markers),
    function(celltype){
            genes<-unlist(cell_markers[celltype])
            print(paste("Plotting:",celltype))
            print(paste("Genes to plot:",genes))
            plt<-histograModified(obj=dat_fine, 
                baseline="mean",
                genes = unlist(cell_markers[celltype]),
                order=order,
                matrix = paste0("cg_",broad_celltype,".clusters_tracks"), arrowScale = .03, trackScale = .5,
                legend = F, cgisland=cgi) + ggtitle(celltype)
            ggsave(plt,
            file=paste(broad_celltype,"clusters",celltype,"fine_celltype.markers.pdf",sep="."),
            width=length(genes)*5,
            height=length(met_clusters)*20,
            limitsize=F)
        })

#dmr set
dmr<-readRDS(
        paste(project_data_directory,"fine_celltyping",celltype,
        paste(celltype,"dmr",celltype,"clusters.collapsed.rds",sep="."),
        sep="/"))

#rna markers defined by paired 10x RNA analysis
rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.rds")
rna<-subset(rna,broad_celltype %in% c("myeloid","bcell","tcell","plasma"))
rna<-subset(rna,fine_celltype %in% c("suspected_doublet"),invert=TRUE)
table(rna$fine_celltype)


#trying to match the correct level of cell typing
rna@meta.data[rna@meta.data$fine_celltype %in% c("bcell"),]$broad_celltype<-"bcell"
rna@meta.data[rna@meta.data$fine_celltype %in% c("myeloid_DC","myeloid_cycling","myeloid_macro","myeloid_mast","myeloid_mono","myeloid_neutrophil"),]$broad_celltype<-"myeloid"
rna@meta.data[rna@meta.data$fine_celltype %in% c("plasma"),]$broad_celltype<-"plasma"
rna@meta.data[rna@meta.data$fine_celltype %in% c("tcell_cd4","tcell_cd8","tcell_nk","tcell_treg"),]$broad_celltype<-"tcell"

rna<-JoinLayers(rna)

#by broad cell types first
markers<-FindAllMarkers(rna,group.by="broad_celltype")

#plot overlap of RNA markers for fine cell types and our cluster hypo DMRs
plot_dmr_RNAmarker_overlap(dmr=dmr,
                            rna_markers=markers,
                            celltype=celltype)

#plot RNA defined markers that overlap with large DMRs
plot_dmr_RNAMarker_hypoM(dat=dat_fine,
                            dmr=dmr,
                            rna_markers=markers,
                            order=c("7","6","8","3","10","1","4","2","5"),
                            broad_celltype=celltype)

```

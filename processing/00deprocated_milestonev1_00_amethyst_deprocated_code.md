
/**
#set celltype plotting order
obj@metadata$celltype<-factor(obj@metadata$celltype,levels=c(
"basal",
"lumhr",
"lumsec",
"fibro",
"endothelial",
"myeloid",
"tcell"))

plt<-ggplot2::ggplot(obj@metadata, aes(x = sample)) +
    ggplot2::geom_bar(aes(fill = celltype), position = "fill") +
    ggplot2::scale_fill_manual(values = cell_colors ) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "sample identity", y = "percentage", title = "sample composition") +
    ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(plt,file="scaledcis.amethyst.celltype_comp.pdf")

hbca_samp<-c("BCMHBCA03R",
"BCMHBCA04R",
"BCMHBCA10L-3h",
"BCMHBCA12R-3h",
"BCMHBCA15R-3h",
"BCMHBCA22R-4h",
"BCMHBCA26L-24hTis-4h",
"BCMHBCA29L-2h",
"BCMHBCA38L-3h",
"BCMHBCA85L-3h",
"HBCA-83L",
"HBCA-17T",
"HBCA-19T",
"HBCA-16R")

puredcis_samp<-c("BCMDCIS07T",
"BCMDCIS102T-4h",
"BCMDCIS41T",
"BCMDCIS66T",
"BCMDCIS79T-24hTis_DCIS",
"BCMDCIS80T",
"BCMDCIS94T-24hTis",
"BCMDCIS99T",
"ECIS25T",
"ECIS26T",
"ECIS36T",
"ECIS48T",
"DCIS-66T",
"DCIS-79T")

obj@metadata$clin<-"SYNCH"
obj@metadata[obj@metadata$sample %in% puredcis_samp,]$clin<-"DCIS"
obj@metadata[obj@metadata$sample %in% hbca_samp,]$clin<-"HBCA"

celltype_counts<-obj@metadata %>% group_by(clin,sample) %>% count(celltype)

plt<-ggplot2::ggplot(celltype_counts, aes(x = sample,)) +
    ggplot2::geom_bar(aes(fill = celltype,y=n), position = "fill",stat="identity") +
    ggplot2::scale_fill_manual(values = cell_colors ) +
    ggplot2::theme_classic() +
    facet_wrap(~clin,scales="free")+
    ggplot2::labs(x = "sample identity", y = "percentage", title = "sample composition") +
    ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

ggsave(plt,file="scaledcis.amethyst.celltype_comp.pdf")

celltype_counts<-obj@metadata %>% group_by(clin) %>% count(celltype)
plt<-ggplot2::ggplot(celltype_counts, aes(x = clin)) +
    ggplot2::geom_bar(aes(fill = celltype,y=n), position = "fill",stat="identity") +
    ggplot2::scale_fill_manual(values = cell_colors ) +
    ggplot2::theme_classic() +
    facet_wrap(~clin,scales="free")+
    ggplot2::labs(x = "sample identity", y = "percentage", title = "sample composition") +
    ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

ggsave(plt,file="scaledcis.amethyst.celltype_comp.clin.pdf")

sample_counts<-obj@metadata %>% group_by(clin,sample) %>% count(sample)

plt<-ggplot(sample_counts,aes(x=sample,y=n))+geom_bar(stat = "identity")+theme_minimal()+theme(axis.text.x = element_text(angle=90))+facet_wrap(~clin,scales="free_x")
ggsave(plt,file="scaledcis.cellcounts.pdf")
#summarize over stages

#run subclustering per cell type

#add in scbam locations to each cell file for CNV calling
sc_bams=list.files(path=project_data_directory,pattern="*bam",recursive=TRUE,full.names=TRUE)
sc_bams<-sc_bams[grep(sc_bams,pattern="cnv/sc_bam/")]
sample<-basename(unlist(lapply(strsplit(sc_bams,"[.]"),"[",1)))
tn5_well<-unlist(lapply(strsplit(sc_bams,"[.]"),"[",2))
barc<-unlist(lapply(strsplit(sc_bams,"[.]"),"[",3))
run<-unlist(lapply(strsplit(sc_bams,"/"),"[",7))
sc_bam_df<-data.frame(
    files=sc_bams,run=run,sample=sample,
    tn5_well=tn5_well,barc=barc)
sc_bam_df$cellid=paste(sep="+",sc_bam_df$barc,sc_bam_df$run)
sc_bam_df<-sc_bam_df[!duplicated(sc_bam_df$cellid),]
row.names(sc_bam_df)<-sc_bam_df$cellid

obj@metadata$sc_bam<-sc_bam_df[row.names(obj@metadata),]$files

write.table(as.data.frame(obj@metadata),sep="\t",
    col.names=T,row.names=T,file="scaledcis.metadata.tsv")
    

```
Cluster by cell type
```R
#subset
#cluster on all DMRs
#call DMRs
#cluster on cell type DMRs



############################
######CELL SUBTYPING##########
############################
#
celltype_clustering<-function(celltype,d=12,neigh=25,q=0.000001){
    objsub<-subsetObject(obj, cells = row.names(obj@metadata[obj@metadata$celltype %in% celltype,]))
    #all cells DMRs for clustering
    objsub<-cluster_by_windows(objsub,window_name="dmr_sites",
            outname=paste("subcelltype",celltype,d,neigh,q,sep="_"),
            threads.=task_cpus,est_dim=d,neighbors.=neigh,dist.=q,
            k_pheno=15)
    #in clusters define cell subtype DMRs
    out<-dmr_and_1kb_window_gen(objsub,
        prefix=paste(celltype,"dmr_sites",sep="."),
        threads.=50)
    dmr_out<-out[["dmr"]]
    objsub<-out[["object"]]

    dmr_out<-dmr_out |> 
        dplyr::filter(dmr_padj<0.01) |> 
        dplyr::filter(dmr_logFC < -1 | dmr_logFC > 1) |> 
        select(chr,dmr_start,dmr_end) |> 
        distinct(chr,dmr_start,dmr_end) |> 
        makeGRangesFromDataFrame() |> 
        reduce()
    dmr_bed<-data.frame(seqnames=seqnames(dmr_out),start=start(dmr_out),end=end(dmr_out))
    #cluster on cell subtype DMRs
    objsub@genomeMatrices[["dmr_sites_celltype"]] <- makeWindows(objsub, 
                                                        bed = dmr_bed,
                                                        type = "CG", 
                                                        metric = "score", 
                                                        threads = 20, 
                                                        index = "chr_cg", 
                                                        nmin = 2) 
    objsub<-cluster_by_windows(objsub,window_name="dmr_sites_celltype",
            outname=paste("subcelltype",celltype,d,neigh,q,sep="_"),
            threads.=task_cpus,est_dim=d,neighbors.=neigh,dist.=q,
            k_pheno=15)
    #make tracks from cell subtype subclustering
    out<-dmr_and_1kb_window_gen(objsub,
        prefix=paste(celltype,"dmr_sites_celltype",sep="."),
        threads.=50)
    dmr_out<-out[["dmr"]]
    objsub<-out[["object"]]
    saveRDS(objsub,file=paste0("./merged_data/scaledcis.",celltype,".amethyst.pf.celltypedmr.rds"))
    saveRDS(dmr_out,file=paste0("./merged_data/scaledcis.",celltype,".nakshatri_cluster.celltypedmr.rds"))
}

#rerunning function with additional processing
for(i in list.files(path="./merged_data/",pattern="scaledcis.*.amethyst.pf.nakshatri_clus.rds")){
    if(i!="scaledcis.amethyst.pf.nakshatri_clus.rds"){
        print(i)
        celltype<-gsub(pattern="scaledcis.",replacement="",basename(i))
        celltype<-gsub(pattern=".amethyst.pf.nakshatri_clus.rds",replacement="",celltype)
        print(celltype)
        objsub<-readRDS(paste0("./merged_data/",i))
        print("read in")
        objsub<-cluster_by_windows(objsub,window_name="dmr_sites",
                outname=paste("subcelltype_dmr",celltype,d,neigh,q,sep="_"),
                threads.=task_cpus,est_dim=d,neighbors.=neigh,dist.=q,
                k_pheno=15)
        out<-dmr_and_1kb_window_gen(objsub,
            prefix=paste(celltype,"dmr_sites_celltype",sep="."),
            threads.=50)
        dmr_out<-out[["dmr"]]
        objsub<-out[["object"]]
        print("saving")
        saveRDS(objsub,file=paste0("./merged_data/scaledcis.",celltype,".amethyst.pf.celltypedmr.rds"))
        saveRDS(dmr_out,file=paste0("./merged_data/scaledcis.",celltype,".nakshatri_cluster.celltypedmr.rds"))
}} 

lapply(unique(obj@metadata$celltype),celltype_clustering)
lapply(c("lumhr","myeloid"),celltype_clustering)

cell_subtype<-list()
celltype="lumsec"      
objsub<-readRDS(file=paste0("./merged_data/scaledcis.",celltype,".amethyst.pf.celltypedmr.rds"))
plt_list<-lapply(names(cell_colors),
    function(celltype){
        genes<-unlist(cell_markers[celltype])
        genes<-genes[genes %in% objsub@ref$gene_name]
        if(length(genes)>1){
        print(paste("Plotting:",celltype))
        print(paste("Genes to plot:",genes))
        plt<-histograModified(objsub, 
            baseline="mean",
            genes = unlist(cell_markers[celltype]),
            colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
            matrix = "cg_cluster_id_tracks", arrowScale = .03, trackScale = .5,
            legend = F, cgisland=cgi) + ggtitle(celltype)
        return(plt)
        }else{return(ggplot()) }})
plt_out<-wrap_plots(plt_list,ncol=1) 
ggsave(plt_out ,file=paste0(celltype,".celltype.marker.pdf"),width=30,height=length(cell_colors)*30,limitsize=F)
objsub@metadata$cell_subtype<-"lumsec"
cell_subtype[[celltype]]<-setNames(nm=row.names(objsub@metadata),objsub@metadata$cell_subtype)

celltype="lumhr"   
objsub<-readRDS(file=paste0("./merged_data/scaledcis.",celltype,".amethyst.pf.celltypedmr.rds"))
plt_list<-lapply(names(cell_colors),
    function(celltype){
        genes<-unlist(cell_markers[celltype])
        genes<-genes[genes %in% objsub@ref$gene_name]
        if(length(genes)>1){
        print(paste("Plotting:",celltype))
        print(paste("Genes to plot:",genes))
        plt<-histograModified(objsub, 
            baseline="mean",
            genes = unlist(cell_markers[celltype]),
            colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
            matrix = "cg_cluster_id_tracks", arrowScale = .03, trackScale = .5,
            legend = F, cgisland=cgi) + ggtitle(celltype)
        return(plt)
        }else{return(ggplot()) }})
plt_out<-wrap_plots(plt_list,ncol=1) 
ggsave(plt_out ,file=paste0(celltype,".celltype.marker.pdf"),width=30,height=length(cell_colors)*30,limitsize=F)
objsub@metadata$cell_subtype<-"lumhr"
cell_subtype[[celltype]]<-setNames(nm=row.names(objsub@metadata),objsub@metadata$cell_subtype)

celltype="fibro"  
objsub<-readRDS(file=paste0("./merged_data/scaledcis.",celltype,".amethyst.pf.celltypedmr.rds"))
plt_list<-lapply(names(cell_colors),
    function(celltype){
        genes<-unlist(cell_markers[celltype])
        genes<-genes[genes %in% objsub@ref$gene_name]
        if(length(genes)>1){
        print(paste("Plotting:",celltype))
        print(paste("Genes to plot:",genes))
        plt<-histograModified(objsub, 
            baseline="mean",
            genes = unlist(cell_markers[celltype]),
            colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
            matrix = "cg_cluster_id_tracks", arrowScale = .03, trackScale = .5,
            legend = F, cgisland=cgi) + ggtitle(celltype)
        return(plt)
        }else{return(ggplot()) }})
plt_out<-wrap_plots(plt_list,ncol=1) 
ggsave(plt_out ,file=paste0(celltype,".celltype.marker.pdf"),width=30,height=length(cell_colors)*30,limitsize=F)
objsub@metadata$cell_subtype<-"fibro"
cell_subtype[[celltype]]<-setNames(nm=row.names(objsub@metadata),objsub@metadata$cell_subtype)

celltype="endothelial" 
objsub<-readRDS(file=paste0("./merged_data/scaledcis.",celltype,".amethyst.pf.celltypedmr.rds"))
plt_list<-lapply(names(cell_colors),
    function(celltype){
        genes<-unlist(cell_markers[celltype])
        genes<-genes[genes %in% objsub@ref$gene_name]
        if(length(genes)>1){
        print(paste("Plotting:",celltype))
        print(paste("Genes to plot:",genes))
        plt<-histograModified(objsub, 
            baseline="mean",
            genes = unlist(cell_markers[celltype]),
            colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
            matrix = "cg_cluster_id_tracks", arrowScale = .03, trackScale = .5,
            legend = F, cgisland=cgi) + ggtitle(celltype)
        return(plt)
        }else{return(ggplot()) }})
plt_out<-wrap_plots(plt_list,ncol=1) 
ggsave(plt_out ,file=paste0(celltype,".celltype.marker.pdf"),width=30,height=length(cell_colors)*30,limitsize=F)
objsub@metadata$cell_subtype<-"endothelial"
cell_subtype[[celltype]]<-setNames(nm=row.names(objsub@metadata),objsub@metadata$cell_subtype)

celltype="myeloid"
objsub<-readRDS(file=paste0("./merged_data/scaledcis.",celltype,".amethyst.pf.celltypedmr.rds"))
plt_list<-lapply(names(cell_colors),
    function(celltype){
        genes<-unlist(cell_markers[celltype])
        genes<-genes[genes %in% objsub@ref$gene_name]
        if(length(genes)>1){
        print(paste("Plotting:",celltype))
        print(paste("Genes to plot:",genes))
        plt<-histograModified(objsub, 
            baseline="mean",
            genes = unlist(cell_markers[celltype]),
            colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
            matrix = "cg_cluster_id_tracks", arrowScale = .03, trackScale = .5,
            legend = F, cgisland=cgi) + ggtitle(celltype)
        return(plt)
        }else{return(ggplot()) }})
plt_out<-wrap_plots(plt_list,ncol=1) 
ggsave(plt_out ,file=paste0(celltype,".celltype.marker.pdf"),width=30,height=length(cell_colors)*30,limitsize=F)
objsub@metadata$cell_subtype<-"bcell"
objsub@metadata[objsub@metadata$cluster_id %in% c("11","10","9","1"),]$cell_subtype<-"myeloid"
cell_subtype[[celltype]]<-setNames(nm=row.names(objsub@metadata),objsub@metadata$cell_subtype)

celltype="tcell"  
objsub<-readRDS(file=paste0("./merged_data/scaledcis.",celltype,".amethyst.pf.celltypedmr.rds"))
plt_list<-lapply(names(cell_colors),
    function(celltype){
        genes<-unlist(cell_markers[celltype])
        genes<-genes[genes %in% objsub@ref$gene_name]
        if(length(genes)>1){
        print(paste("Plotting:",celltype))
        print(paste("Genes to plot:",genes))
        plt<-histograModified(objsub, 
            baseline="mean",
            genes = unlist(cell_markers[celltype]),
            colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
            matrix = "cg_cluster_id_tracks", arrowScale = .03, trackScale = .5,
            legend = F, cgisland=cgi) + ggtitle(celltype)
        return(plt)
        }else{return(ggplot()) }})
plt_out<-wrap_plots(plt_list,ncol=1) 
ggsave(plt_out ,file=paste0(celltype,".celltype.marker.pdf"),width=30,height=length(cell_colors)*30,limitsize=F)
objsub@metadata$cell_subtype<-"tcell"
objsub@metadata[objsub@metadata$cluster_id %in% c("5","4"),]$cell_subtype<-"myeloid2"
cell_subtype[[celltype]]<-setNames(nm=row.names(objsub@metadata),objsub@metadata$cell_subtype)

celltype="basal" 
objsub<-readRDS(file=paste0("./merged_data/scaledcis.",celltype,".amethyst.pf.celltypedmr.rds"))
plt_list<-lapply(names(cell_colors),
    function(celltype){
        genes<-unlist(cell_markers[celltype])
        genes<-genes[genes %in% objsub@ref$gene_name]
        if(length(genes)>1){
        print(paste("Plotting:",celltype))
        print(paste("Genes to plot:",genes))
        plt<-histograModified(objsub, 
            baseline="mean",
            genes = unlist(cell_markers[celltype]),
            colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
            matrix = "cg_cluster_id_tracks", arrowScale = .03, trackScale = .5,
            legend = F, cgisland=cgi) + ggtitle(celltype)
        return(plt)
        }else{return(ggplot()) }})
plt_out<-wrap_plots(plt_list,ncol=1) 
ggsave(plt_out ,file=paste0(celltype,".celltype.marker.pdf"),width=30,height=length(cell_colors)*30,limitsize=F)
objsub@metadata$cell_subtype<-"basal"
cell_subtype[[celltype]]<-setNames(nm=row.names(objsub@metadata),objsub@metadata$cell_subtype)

p1 <- dimFeature(obj, colorBy = celltype, reduction = "umap")
ggsave(p1,file="celltype.umap.pdf")

obj<-readRDS(file="scaledcis.amethyst.pf.celltype.rds")
obj@metadata$cell_subtype<-NA
obj@metadata[unlist(lapply(cell_subtype,names)),]$cell_subtype<-unname(unlist(cell_subtype))
#hotfix!
obj@metadata[obj@metadata$cluster_id=="19",]$cell_subtype<-"lumhr"
obj@metadata[obj@metadata$cluster_id=="9",]$cell_subtype<-"myeloid"
saveRDS(obj,file="scaledcis.amethyst.pf.celltype.rds")

p1 <- dimFeature(obj, colorBy = cell_subtype, reduction = "umap")
ggsave(p1,file="subcelltype.umap.pdf")


#set celltype plotting order
obj@metadata$cell_subtype<-factor(obj@metadata$cell_subtype,
levels=c(
"basal",
"lumsec",
"lumhr",
"fibro",
"endothelial",
"myeloid",
"myeloid2",
"bcell",
"tcell"))

cell_colors=c(
"basal"="#844c9d",
"lumsec"="#f68d3e",
"lumhr"="#e23e96",
"fibro"="#68b1af",
"endothelial"="#b5d564",
"myeloid"="#8088c2",
"myeloid2"="#7ecdc2",
"bcell"="#65cbe4",
"tcell"="#1d87c8")

p1 <- dimFeature(obj, colorBy = cell_subtype, reduction = "umap",colors=cell_colors)
ggsave(p1,file="subcelltype.test.umap.pdf")

celltype_counts<-obj@metadata %>% group_by(clin,sample) %>% count(cell_subtype)
plt<-ggplot2::ggplot(celltype_counts, aes(x = sample,)) +
    ggplot2::geom_bar(aes(fill = cell_subtype,y=n), position = "fill",stat="identity") +
    ggplot2::scale_fill_manual(values = cell_colors ) +
    ggplot2::theme_classic() +
    facet_wrap(~clin,scales="free")+
    ggplot2::labs(x = "sample identity", y = "percentage", title = "sample composition") +
    ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

ggsave(plt,file="scaledcis.amethyst.subcelltype_comp.pdf")

celltype_counts<-obj@metadata %>% group_by(clin) %>% count(cell_subtype)
plt<-ggplot2::ggplot(celltype_counts, aes(x = clin)) +
    ggplot2::geom_bar(aes(fill = cell_subtype,y=n), position = "fill",stat="identity") +
    ggplot2::scale_fill_manual(values = cell_colors ) +
    ggplot2::theme_classic() +
    facet_wrap(~clin,scales="free")+
    ggplot2::labs(x = "sample identity", y = "percentage", title = "sample composition") +
    ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

ggsave(plt,file="scaledcis.amethyst.celltype_comp.clin.pdf")

met<-obj@metadata %>% group_by(sample,clin,cell_subtype) %>% summarize(count=n()) %>% mutate(prcnt_of_total = count/sum(count)) %>% as.data.frame()

plt<-ggplot(met,aes(y=prcnt_of_total,x=clin,fill=cell_subtype))+geom_boxplot()+theme_minimal()+facet_wrap(.~cell_subtype,nrow=1)
ggsave(plt,file="scaledcis.amethyst.celltype_barplot.pdf")

sample_counts<-obj@metadata %>% group_by(clin,sample) %>% count(sample)

plt<-ggplot(sample_counts,aes(x=sample,y=n))+geom_bar(stat = "identity")+theme_minimal()+theme(axis.text.x = element_text(angle=90))+facet_wrap(~clin,scales="free_x")
ggsave(plt,file="scaledcis.cellcounts.pdf")
#summarize over stages


out<-dmr_and_1kb_window_gen(obj,step.=500,groupBy="cell_subtype",prefix="dmr_sites_subcelltype",threads.=50)
dmr_out<-out[["dmr"]]
obj<-out[["object"]]
saveRDS(dmr_out,file="scaledcis.amethyst.pf.subcelltype.dmr.rds")
saveRDS(obj,file="scaledcis.amethyst.pf.subcelltype.rds")


objsub<-subsetObject(obj, cells = row.names(obj@metadata[obj@metadata$celltype %in% "lumhr",]))

#in clusters define cell subtype DMRs
outsub<-dmr_and_1kb_window_gen(objsub,groupBy="clin",
    prefix=paste("clin","dmr_sites",sep="."),
    threads.=50)
dmr_out<-outsub[["dmr"]]
objsub<-outsub[["object"]]
saveRDS(dmr_out,file="scaledcis.amethyst.pf.lumhr_clin.dmr.rds")
saveRDS(objsub,file="scaledcis.amethyst.pf.lumhr_clin.rds")


    dmr_out<-dmr_out |> 
        dplyr::filter(dmr_padj<0.01) |> 
        dplyr::filter(dmr_logFC < -1 | dmr_logFC > 1) |> 
        select(chr,dmr_start,dmr_end) |> 
        distinct(chr,dmr_start,dmr_end) |> 
        makeGRangesFromDataFrame() |> 
        reduce()
    dmr_bed<-data.frame(seqnames=seqnames(dmr_out),start=start(dmr_out),end=end(dmr_out))
    #cluster on cell subtype DMRs
    objsub@genomeMatrices[["dmr_sites_celltype"]] <- makeWindows(objsub, 
                                                        bed = dmr_bed,
                                                        type = "CG", 
                                                        metric = "score", 
                                                        threads = 20, 
                                                        index = "chr_cg", 
                                                        nmin = 2) 
    objsub<-cluster_by_windows(objsub,window_name="dmr_sites_celltype",
            outname=paste("subcelltype",celltype,d,neigh,q,sep="_"),
            threads.=task_cpus,est_dim=d,neighbors.=neigh,dist.=q,
            k_pheno=15)
out<-dmr_and_1kb_window_gen(obj,step.=500,groupBy="cell_subtype",prefix="dmr_sites_subcelltype",threads.=50)
dmr_out<-out[["dmr"]]
obj<-out[["object"]]
saveRDS(dmr_out,file="scaledcis.amethyst.pf.subcelltype.dmr.rds")
saveRDS(obj,file="scaledcis.amethyst.pf.subcelltype.rds")


obj<-readRDS(file="scaledcis.amethyst.pf.celltype.rds")
met<-read.csv("/data/rmulqueen/projects/scalebio_dcis/sample_selection/simplified_patient_metadata.csv",header=T)
row.names(met)<-met$Sample
met$Age<-as.numeric(met$Age)
sample_heatmap<-met[c("Diagnosis","Age","Ethnicity","ER","PR","Her2","Grade")]
sample_heatmap$Diagnosis<-factor(sample_heatmap$Diagnosis,levels=c("HBCA","DCIS","IDC"))
sample_heatmap<-sample_heatmap[order(sample_heatmap$Diagnosis),]

library(circlize)
library(ComplexHeatmap,lib.loc="/home/rmulqueen/R/x86_64-pc-linux-gnu-library/4.3")
age_col=colorRamp2(breaks=c(min(sample_heatmap$Age,na.rm=T),max(sample_heatmap$Age,na.rm=T)),c("#d789d7","#2a3d66"))
#plot metadata
ha = rowAnnotation(age=sample_heatmap$Age,
                      ethnicity=sample_heatmap$Ethnicity,
                      grade=sample_heatmap$Grade,
                      col = list(age=age_col),
                      show_legend=TRUE)

pdf("sample_metadata.heatmap.pdf")
plt<-Heatmap(sample_heatmap[c("Diagnosis","ER","PR","Her2")],
 cluster_columns=F,cluster_rows=F,
 left_annotation=ha)
draw(plt)
dev.off()

sample_counts<-obj@metadata %>% group_by(sample) %>% count(sample)
plt<-ggplot(sample_counts,aes(x=sample,y=n))+geom_bar(stat = "identity")+theme_minimal()+theme(axis.text.x = element_text(angle=90))
ggsave(plt,file="scaledcis.cellcounts.pdf")

obj@metadata$run<-unlist(lapply(strsplit(row.names(obj@metadata),"[+]"),"[",4))
obj@metadata$method<-"homebrew"
scale_oligoes<-read.csv("/data/rmulqueen/projects/scalebio_dcis/tools/ScaleMethyl/references/i7.txt",sep="\t",header=F)
i7_idx<-unlist(lapply(strsplit(row.names(obj@metadata),"[+]"),"[",1))
obj@metadata[i7_idx %in% scale_oligoes$V1,]$method<-"scale"

dat<- obj@metadata %>% filter(run=="250329_RM_scalebio_batch1_initseq") 

plt<-ggplot(dat,aes(x=log10(cg_cov),y=method))+geom_jitter()+geom_boxplot()+theme_minimal()+xlim(c(0,7))
ggsave(plt,file="method.cg_cov.pdf")

plt<-ggplot(dat,aes(x=log10(unique_reads),y=method))+geom_jitter()+geom_boxplot()+theme_minimal()+xlim(c(0,7))
ggsave(plt,file="method.unique_reads.pdf")

plt<-ggplot(dat,aes(x=mcg_pct,y=method))+geom_jitter()+geom_boxplot()+theme_minimal()+xlim(c(0,100))
ggsave(plt,file="method.mcg_pct.pdf")

plt<-ggplot(dat,aes(x=tss_enrich,y=method))+geom_jitter()+geom_boxplot()+theme_minimal()
ggsave(plt,file="method.tss_enrich.pdf")

plt<-ggplot(dat,aes(x=mito_reads/unique_reads,y=method))+geom_jitter()+geom_boxplot()+theme_minimal()+xlim(c(0,100))
ggsave(plt,file="method.mito_reads.pdf")

plt<-ggplot(dat,aes(x=unique_reads,y=cg_cov,color=method))+geom_point()+theme_minimal()
ggsave(plt,file="method.read_per_cg.pdf")

dat%>%group_by(method) %>% summarize(count=n(),mean_reads=mean(unique_reads),mito=mean(mito_reads)/mean(unique_reads),cg_per_read=mean(cg_cov/unique_reads))


#plot umap by methods
p1 <- dimFeature(obj, colorBy = method, reduction = "umap")
ggsave(p1,file="method.umap.pdf")

p1 <- dimFeature(obj, colorBy = clin, reduction = "umap")
ggsave(p1,file="clin.umap.pdf")

#plot umap by celltype, 


cell_markers<-list()
cell_markers[["basal"]]<-c("CARMN")
cell_markers[["lumsec"]]<-c("KRT15")
cell_markers[["lumhr"]]<-c("ANKRD30A")
cell_markers[["fibro"]]<-c("COL1A1")
cell_markers[["endothelial"]]=c("CLDN5")
cell_markers[["myeloid"]]<-c("HLA-DRA","HLA-DPA1","CD74")
cell_markers[["tcell"]]<-c("PTPRC")
cell_markers[["mast"]]<-c("NTM","SYTL3","SLC24A3","TPSB2","HDC")
cell_markers[["bcell"]]=c("CD37","TCL1A","LTB","HLA-DPB1","HLA-DRA","HLA-DPA1")
cell_markers[["plasma"]]=c("IGHA2","IGHA1","JCHAIN","IGHM","IGHG1","IGHG4","IGHG3","IGHG2")
cell_markers[["adipo"]]=c("PDE3B","ACACB","WDPCP","PCDH9","CLSTN2","ADIPOQ","TRHDE")
cell_markers[["dcis_diff"]]=c("HOBX13","EN1","TBX15","DLX4")

cell_colors=c(
"basal"="#844c9d",
"lumsec"="#f68d3e",
"lumhr"="#e23e96",
"fibro"="#68b1af",
"endothelial"="#b5d564",
"myeloid"="#8088c2",
"myeloid2"="#7ecdc2",
"bcell"="#65cbe4",
"tcell"="#1d87c8",
"dcis_diff"="black")

plt_list<-lapply(names(cell_colors),
    function(celltype){
        genes<-unlist(cell_markers[celltype])
        genes<-genes[genes %in% obj@ref$gene_name]
        if(length(genes)>0){
        print(paste("Plotting:",celltype))
        print(paste("Genes to plot:",genes))
        plt<-histograModified(obj, 
            baseline="mean",
            genes = unlist(cell_markers[celltype]),
            colors= c(cell_colors[celltype], "#dbdbdb","#cccccc", "#999999"),
            matrix = "cg_cell_subtype_tracks", arrowScale = .03, trackScale = .5,
            legend = F, cgisland=cgi) + ggtitle(celltype)
        return(plt)
        }else{return(ggplot()) }})
plt_out<-wrap_plots(plt_list,ncol=1) 
ggsave(plt_out ,file=paste0("allcelltypes",".celltype.marker.pdf"),width=10,height=length(cell_colors)*50,limitsize=F)


library(dplyr)
library(stringr)
library(msigdbr,lib.loc="/home/users/rmulqueen/R/x86_64-conda-linux-gnu-library/4.4") #add this
library(fgsea,lib.loc="/home/users/rmulqueen/R/x86_64-conda-linux-gnu-library/4.4") #add this
dmrs <- dmr_out |> dplyr::filter(direction=="hypo") |> dplyr::filter(gene_names != "NA") 

#set up gene universe
gene_universe<-dmrs |> pull(gene_names)
gene_universe<-str_replace_all(unlist(lapply(strsplit(gene_universe, ","),unlist)), " ", "") #flatten and remove 
gene_universe<-gene_universe[!duplicated(gene_universe)]


#H is hallmark cancer #C1 is cytoband #C3 is TF targets
run_gsea_enrichment <- function(obj,sample="IDC_12",assay="RNA",category="H") {
    gsea<-lapply(unique(de_genes$cluster), 
    function(x){
        program_genes=unlist(de_genes[de_genes$cluster==x,]$gene)
        out<-runGSEA(program_genes, universe=row.names(obj@assays[[assay]]), category = category) 
        out<-as.data.frame(out)
        out$cluster<-x
        return(out)})
    
    gsea_out<-do.call("rbind",gsea)
    pltdat<-gsea_out %>% group_by(cluster) %>% slice_max(order_by = -padj, n = 5)
    plt<-ggplot(pltdat,aes(x=cluster,y=pathway))+
    geom_point(aes(size = -log10(padj), fill = overlap), shape=21)+
    theme_minimal() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    ggsave(plt,file=paste(sample,"clone_cnv",category,"RNA.pdf",sep="_"),width=20)
}

  #find cancer hallmark signatures
    top_p_H <- do.call("rbind",
    lapply(unique(peaks_out$metaprogram_cluster), 
    function(i) {
    program_name=i
    program_genes=unlist(peaks_out[peaks_out$metaprogram_cluster==program_name,]$SYMBOL)
    out<-runGSEA(program_genes, universe=row.names(dat@assays$RNA), category = "H")
    out$program<-paste0(program_name)
    return(out)
    }
    ))
  pltdat<-top_p_H %>% group_by(program) %>% slice_max(order_by = -padj, n = 5)
  plt<-ggplot(pltdat,aes(x=program,y=pathway))+
  geom_point(aes(size = -log10(padj), fill = overlap), shape=21)+
  theme_minimal() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(plt,file=paste0(prefix,"_metaprograms","_H","_","cistopic",".pdf"),width=20)

   #run gsea enrichment at gene level on different sets
   gsea_enrichment(prefix=prefix,dmrs=dmrs,gene_universe=gene_universe,
   category="C3",subcategory="TFT:GTRD",out_setname="TFT") #find enrichment in tft (transcription factor targets)

   gsea_enrichment(prefix=prefix,dmrs=dmrs,gene_universe=gene_universe,
   category="C1",subcategory=NULL,out_setname="position") #find enrichment in c1 signatures (positional)

   gsea_enrichment(prefix=prefix,dmrs=dmrs,gene_universe=gene_universe,
   category="H",subcategory=NULL,out_setname="hallmark") #find cancer hallmark signatures

```

```


#dmr around promoter regions for cell type calling too (add to plot list)
#methyltree output
#scage output with celltypes



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











#promoter plus into gene body
    promoter_into_genebody_bed <- dat@ref |> dplyr::filter( type == "gene") |>
            dplyr::select(seqid, start, end, gene_name, strand) |>
            dplyr::distinct(gene_name, .keep_all = TRUE) |>
            dplyr::rename("chr" = "seqid", "gene" = "gene_name") |>
            data.table::data.table()
    
    #adding another 5kb into gene body from promoter regions
    promoter_into_genebody_bed[, `:=`(p_start = ifelse(strand == "+", start - 1500, end - 1500),
                    p_end = ifelse(strand == "+", start + 1500 + 5000, end + 1500 + 5000))] # make bed file for promoter cg
        
    promoter_into_genebody_bed[, `:=`(start = NULL, end = NULL, strand = NULL)]
    setnames(promoter_into_genebody_bed, c("chr", "gene", "start", "end"))
    row.names(promoter_into_genebody_bed)<-promoter_into_genebody_bed$gene
    options(future.globals.maxSize= 50*1000*1024^2) #50gb
    dat@genomeMatrices[["promoter_5kb_genebody"]] <- makeWindows(dat, 
                                                        bed = promoter_into_genebody_bed[,c(1,3,4)],
                                                        type = "CG", 
                                                        metric = "score", 
                                                        threads = 5, 
                                                        index = "chr_cg", 
                                                        nmin = 2) 
histograM(dat, 
          genes = c("ELANE", "MPEG1", "SPI1", "CD2", "CD3D"), 
          matrix = "cg_cluster_tracks",
          legend = F)
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



#add additional methylation windows

dat@genomeMatrices[["cg_promoter"]] <- makeWindows(dat, 
                                                     genes = protein_coding,
                                                     promoter = TRUE, 
                                                     type = "CG", 
                                                     metric = "percent", 
                                                     threads = 300, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 



dat@genomeMatrices[["cg_genebody"]] <- makeWindows(dat, 
                                                     genes = protein_coding,
                                                     promoter = FALSE, 
                                                     type = "CG", 
                                                     metric = "percent", 
                                                     threads = 300, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

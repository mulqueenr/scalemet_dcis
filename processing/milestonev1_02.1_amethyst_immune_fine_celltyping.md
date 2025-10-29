```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
```


```R
source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
library(Seurat)
library(ComplexHeatmap)
library(GeneOverlap)

task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)
setwd(wd)

system(paste0("mkdir -p ",project_data_directory,"/fine_celltyping"))

obj<-readRDS(file="04_scaledcis.broad_celltype.amethyst.rds")

#prepare cgi
cgisland="/data/rmulqueen/projects/scalebio_dcis/ref/cpgIslandExt.bed"
cgi<-rtracklayer::import(cgisland)
cgi<-as.data.frame(cgi)
colnames(cgi)<-c("chr","start","end","strand")

broad_celltype="immune"
dat<-obj

output_directory=paste0(project_data_directory,"/fine_celltyping/",broad_celltype)
system(paste0("mkdir -p ",output_directory))


print(paste("Subsetting object based on broad celltype:",broad_celltype))
dat<-subsetObject(dat,cells=row.names(dat@metadata[dat@metadata$broad_celltype %in% c("myeloid","lymphoid"),]))

#overcluster by existing DMR
print(paste("Reclustering on only",broad_celltype,"cells"))
dat<-cluster_by_windows(
    dat, 
    window_name="nakshatri_dmr_sites",
    outname=paste(output_directory,
            paste(broad_celltype,"broad_dmr_sites_clustering",sep="."),
            sep="/"),
    threads=task_cpus, 
    est_dim=10, 
    neighbors=8, 
    dist=1E-8,
    k_pheno=50)

dat@metadata[[paste(broad_celltype,"clusters",sep=".")]]<-dat@metadata$cluster_id  #set celltype clusters
table(dat@metadata[[paste(broad_celltype,"clusters",sep=".")]])

#find broad celltype specific DMRs
print(paste("Identifying DMRs on",broad_celltype,"cells"))
dat<-dmr_and_1kb_window_gen(dat,
    prefix=paste(output_directory,broad_celltype,sep="/"),
    groupBy=paste(broad_celltype,"clusters",sep="."),
    threads=10,step=1000)

collapsed_dmrs<-readRDS(file=paste0(output_directory,"/",
    broad_celltype,".dmr.",paste(broad_celltype,"clusters",sep="."),".collapsed.rds"))

table(collapsed_dmrs$celltype)

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

#cluster on defined DMRs within broad celltype
dat<-cluster_by_windows(
    dat, window_name="fine_dmr_sites",
    outname=paste0(output_directory,"/",paste(broad_celltype,"broad_dmr_sites_clustering",sep=".")),
    threads=10, 
    est_dim=15, 
    neighbors=20, 
    dist=1E-8,
    k_pheno=100)

dat@metadata[[paste(broad_celltype,"clusters",sep=".")]]<-dat@metadata$cluster_id  #set celltype clusters
table(dat@metadata[[paste(broad_celltype,"clusters",sep=".")]])


print(paste("Make 500bp windows on",broad_celltype,"clusters"))

#summarize celltype-dmr clusters over windows
dat<-dmr_and_1kb_window_gen(dat,
    prefix=paste(output_directory,broad_celltype,sep="/"),
    groupBy=paste(broad_celltype,"clusters",sep="."),
    threads=10,step=500)
    
print(paste("Saving final amethyst object of",broad_celltype))
saveRDS(dat,paste0(output_directory,"/",paste("04_scaledcis",broad_celltype,"fine_celltype.amethyst.rds",sep=".")))


#read in subset amethyst object and dmrs
dat_fine<-readRDS(
    paste(project_data_directory,"fine_celltyping",broad_celltype,
    paste("04_scaledcis",broad_celltype,"fine_celltype.amethyst.rds",sep="."),
    sep="/"))


#UMAP PLOTTING
p1 <- dimFeature(dat_fine, colorBy = cluster_id, reduction = "umap") + ggtitle(broad_celltype)
p2 <- dimFeature(dat_fine, colorBy = mcg_pct, reduction = "umap") + ggtitle("%metCG")+scale_color_gradientn(colors = c("black", "turquoise", "gold", "red"),guide="colourbar")
p3 <- dimFeature(dat_fine, colorBy = log10(cov), reduction = "umap") + ggtitle("cov")+ scale_color_gradientn(colors = c("black", "turquoise", "gold", "red"),guide="colourbar")
p4 <- dimFeature(dat_fine, colorBy = Group, reduction = "umap") + ggtitle("Group")
p5 <- dimFeature(dat_fine, colorBy = Race_Ethnicity, reduction = "umap") + ggtitle("Race/Ethnicity")
p6 <- dimFeature(dat_fine, colorBy = Age, reduction = "umap") + ggtitle("Age")
ggsave((p1|p2)/(p3|p4)/(p5|p6),file=paste0(output_directory,"/",paste("04_scaledcis",broad_celltype,"fine_celltype.umap.pdf",sep=".")),width=20,height=20)  



#STACKED BARPLOT PER CLUSTER PLOTTING
p1<-ggplot(dat_fine@metadata, aes(fill=Group,x=cluster_id)) + 
    geom_bar(position="fill", stat="count")+
    scale_x_discrete(limits=factor(order))
p2<-ggplot(dat_fine@metadata, aes(fill=batch,x=cluster_id)) + 
    geom_bar(position="fill", stat="count")+
    scale_x_discrete(limits=factor(order))
p3<-ggplot(dat_fine@metadata, aes(fill=method,x=cluster_id)) + 
    geom_bar(position="fill", stat="count")+
    scale_x_discrete(limits=factor(order))
p4<-ggplot(dat_fine@metadata, aes(fill=Sample,x=cluster_id)) + 
    geom_bar(position="fill", stat="count")+
    scale_x_discrete(limits=factor(order))
p5<-ggplot(dat_fine@metadata, aes(fill=Menopause,x=cluster_id)) + 
    geom_bar(position="fill", stat="count")+
    scale_x_discrete(limits=factor(order))
p6<-ggplot(dat_fine@metadata, aes(color=cluster_id,y=mcg_pct,x=cluster_id)) + 
    geom_jitter()+geom_violin()+ylim(c(0,100))+
    scale_x_discrete(limits=factor(order))
p7<-ggplot(dat_fine@metadata, aes(color=cluster_id,y=cov,x=cluster_id)) + 
    geom_jitter()+geom_violin()+
    scale_x_discrete(limits=factor(order))
ggsave(p1/p2/p3/p4/p5/p6/p7,file=paste0(output_directory,"/",paste("04_scaledcis",broad_celltype,"fine_celltype.stacked_barplots.pdf",sep=".")),width=10,height=40)  


#more broad cell markers
cell_markers<-list()
cell_markers[["myeloid"]]<-c("HLA-DRA","HLA-DPA1","CD74")
cell_markers[["tcell"]]<-c("PTPRC","IKZF1","IL7R","GNLY")
cell_markers[["bcell"]]=c("CD37","TCL1A","LTB","HLA-DPB1")
cell_markers[["plasma"]]=c("IGHA2","IGHA1","JCHAIN","IGHM","IGHG1","IGHG4","IGHG3","IGHG2")
cell_markers[["adipo"]]=c("PDE3B","ACACB","WDPCP","PCDH9","CLSTN2","ADIPOQ","TRHDE")

#set order by umap groupings
order<-c("2","4","8","6","5","1","7","3","9")
cell_markers<-cell_markers
plt_list<-lapply(names(cell_markers),
    function(celltype){
            genes<-unlist(cell_markers[celltype])
            print(paste("Plotting:",celltype))
            print(paste("Genes to plot:",genes))
            plt<-histograModified(obj=dat_fine, 
                baseline="mean",
                genes = unlist(cell_markers[celltype]),
                order=order, promoter_focus=TRUE,trackOverhang=10000,
                matrix = paste0("cg_",broad_celltype,".clusters_tracks"), arrowScale = .03, trackScale = .5,
                legend = F, cgisland=cgi) + ggtitle(celltype)
            ggsave(plt,
            file=paste(project_data_directory,"fine_celltyping",broad_celltype,
            paste("04_scaledcis",broad_celltype,celltype,"broad_celltype.markers.pdf",sep="."),
            sep="/"),
            width=length(genes)*2,
            height=length(order)*2,
            limitsize=F)
        })


immune_markers<-list()
immune_markers[["myeloid_mono_classical"]]<-c('SERPINB2', 'AQP9', 'FCN1', 'VCAN', 'AC245128.3', 'EREG', 'TNIP3', 'CD300E', 'S100A8', 'S100A9', 'THBS1', 'CYP1B1', 'SLC39A8', 'MAP3K20', 'EHD1', 'LUCAT1', 'FCAR', 'S100A12', 'PID1', 'ACOD1')
immune_markers[["myeloid_mono_nonclassical"]]<-c('CCDC68', 'CLEC12A', 'CD48', 'OLIG1', 'CDKN1C', 'S100A4', 'KRTAP2-3', 'CD52', 'WARS', 'LINC02432', 'CFP', 'UNC45B', 'APOBEC3A', 'TKT', 'RIPOR2', 'CDA', 'PLAC8', 'FFAR2', 'PAPSS2', 'GK5')
immune_markers[["myeloid_macro_C3"]]<-c('C3', 'GPR158', 'C1QC', 'C1QB', 'AC079760.2', 'FCGR3A', 'SDS', 'A2M', 'AL590705.1', 'MT3', 'TREM2', 'HLA-DQB2', 'OLR1', 'ZBED9', 'IL4I1', 'DOK5', 'IGLC7', 'PMEPA1', 'IBSP', 'CD9')
immune_markers[["myeloid_macro_LYVE1"]]<-c('PDK4', 'SELENOP', 'FOLR2', 'LYVE1', 'LILRB5', 'MAF', 'KLF4', 'RNASE1', 'SLC40A1', 'RHOB', 'F13A1', 'MRC1', 'TSC22D3', 'PLTP', 'STAB1', 'NR4A2', 'MS4A4A', 'DAB2', 'CD163', 'CTSC')
immune_markers[["myeloid_macro_CCL4"]]<-c('CCL4', 'SELENOP', 'FOLR2', 'RNASE1', 'MRC1', 'F13A1', 'LYVE1', 'CD163', 'STAB1', 'CXCL3', 'HMOX1', 'CCL3', 'CXCL2', 'AC010980.2', 'GCLM', 'CD93', 'CTSL', 'PLTP', 'MAFB', 'THBD')
immune_markers[["myeloid_macro_APOC1"]]<-c('ACP5', 'CYP27A1', 'LIPA', 'PLD3', 'TREM2', 'PLA2G7', 'CD36', 'FABP4', 'GPNMB', 'CCL18', 'NCEH1', 'RARRES1', 'FUCA1', 'APOC1', 'LPL', 'OTOA', 'BLVRB', 'HS3ST2', 'NR1H3', 'CHIT1')
immune_markers[["myeloid_macro_interferon"]]<-c( 'EPSTI1', 'CXCL10', 'GBP1', 'STAT1', 'CXCL11', 'CXCL9', 'GBP4', 'IFI44L', 'FCGR1A', 'SERPING1', 'IFI6', 'IGSF6', 'IFIH1', 'IDO1', 'RASSF4', 'XAF1', 'APOL6', 'ISG15', 'FYB1', 'PSME2')
immune_markers[["myeloid_macro_stress"]]<-c( 'HSPA1B', 'BAG3', 'DNAJB1', 'HSPA6', 'MRPL18', 'ZFAND2A', 'CHORDC1', 'SDS', 'CACYBP', 'FKBP4', 'DNAJA4', 'HSPE1', 'RGS1', 'HSPH1', 'HSPA1A', 'OLR1', 'HSPB1', 'CLEC2B', 'CKS2', 'PLEK')
immune_markers[["myeloid_macro_CTSK"]]<-c( 'AP000904.1', 'CTSK', 'ACP5', 'SLC9B2', 'MMP9', 'DGKI', 'SPP1', 'TCIRG1', 'ATP6V0D2', 'CKB', 'CST3', 'CD109', 'TIMP2', 'JDP2', 'SPARC', 'C9orf16', 'GRN', 'SIGLEC15', 'AK5', 'SNX10')
immune_markers[["myeloid_cycling"]]<-c( 'MKI67', 'PCLAF', 'CKS1B', 'DIAPH3', 'TYMS', 'TOP2A', 'HIST1H1B', 'ASPM', 'TK1', 'DHFR', 'NUSAP1', 'NCAPG', 'DLGAP5', 'CENPM', 'TPX2', 'RRM2', 'CEP55', 'MYBL2', 'BIRC5', 'CENPF')
immune_markers[["myeloid_neutrophil"]]<-c('CRISP3', 'ANKRD34B', 'ABCA13', 'CD177', 'ORM1', 'SERPINB4', 'DEFA3', 'DPYS', 'C7', 'AZU1', 'CLC', 'FOXQ1', 'RAB6B', 'SFRP1', 'TWIST1', 'MAP1B', 'SLPI', 'ADH1B', 'FOLR3', 'BTNL9')
immune_markers[["myeloid_mast"]]<-c('CPA3', 'TPSAB1', 'MS4A2', 'TPSB2', 'HDC', 'IL1RL1', 'GATA2', 'HPGDS', 'GCSAML', 'HPGD', 'ADCYAP1', 'KIT', 'KRT1', 'PRG2', 'CTSG', 'CLU', 'TPSD1', 'CALB2', 'RAB27B', 'SLC18A2')
immune_markers[["myeloid_cDC1"]]<-c('CLEC9A', 'XCR1', 'DNASE1L3', 'IDO1', 'LGALS2', 'CPNE3', 'WDFY4', 'C1orf54', 'DAPP1', 'RAB11FIP1', 'PPA1', 'CPVL', 'C9orf135', 'FOXB1', 'AC016717.2', 'CST3', 'PLEKHA5', 'GPR157', 'FOXN2', 'NCOA7')
immune_markers[["myeloid_cDC2"]]<-c('FCER1A', 'IL1R2', 'CLEC10A', 'CD1C', 'CST7', 'GPAT3', 'DAPP1', 'CFP', 'EREG', 'CCL22', 'LGALS2', 'JAML', 'PID1', 'AREG', 'IL7R', 'AC020656.1', 'IL1R1', 'MARCKSL1', 'ADAM19', 'SLC7A11')
immune_markers[["myeloid_mDC"]]<-c('LAMP3', 'BIRC3', 'LACRT', 'NUB1', 'CCR7', 'MARCKSL1', 'IDO1', 'DAPP1', 'POGLUT1', 'LINC01539', 'GPR157', 'IL12B', 'LAD1', 'KIF2A', 'FSCN1', 'IL7R', 'TXN', 'DUSP5', 'FOXD4L1', 'CD200', 'RAB9A')
immune_markers[["myeloid_pDC"]]<-c('CD5', 'LYPD8', 'PRL', 'AC136475.3', 'LINC01087', 'AC006058.1', 'BDKRB2', 'POTEI', 'DSP', 'AC026369.3', 'GZMB', 'TNFSF4', 'CD2', 'TIGIT', 'LTB', 'TMEM45A', 'PGR', 'CD3G', 'AC015936.1', 'ACOT7')

immune_markers[["tcell_cd4_memory"]]<-c('ADAM23', 'NEFL', 'LINC02273', 'ANTXR2', 'MFHAS1', 'AP3M2', 'GPR183', 'PASK', 'S1PR1', 'FTH1', 'GLIPR1', 'SESN3', 'IL7R', 'CCR7', 'SLC2A3', 'DDIT4', 'CD28', 'ANXA1', 'TRAT1', 'KLF2')
immune_markers[["tcell_cd4_stress"]]<-c('MYADM', 'LMNA', 'ANXA1', 'KLF6', 'HSP90AA1', 'RGCC', 'FAM107B', 'AHNAK', 'HSPH1', 'VIM', 'HSPA8', 'GPR183', 'HSP90AB1', 'S100A10', 'S100A11', 'EZR', 'HSPD1', 'IL7R', 'ANKRD12', 'TUBB4B')
immune_markers[["tcell_cd4_naive"]]<-c('CA6', 'ADTRP', 'GTSCR1', 'LINC00402', 'LINC01481', 'LEF1', 'MAL', 'SELL', 'AL391097.1', 'ANKRD55', 'NEXMIF', 'CCR7', 'LINC00861', 'TESPA1', 'ANK1', 'LTB', 'SESN3', 'TSHZ2', 'PASK', 'LDHB')
immune_markers[["tcell_treg"]]<-c('CD177', 'FANK1', 'LINC02099', 'IL1R2', 'CCR8', 'FOXP3', 'PTGIR', 'LAYN', 'IKZF4', 'F5', 'RTKN2', 'CARD16', 'MAGEH1', 'IL2RA', 'LINC01943', 'BATF', 'IKZF2', 'HTATIP2', 'GK', 'CTLA4')
immune_markers[["tcell_cd4_Tfh"]]<-c('IGFL2', 'CXCL13', 'CPM', 'CD200', 'KSR2', 'NMB', 'HMSD', 'TUBA3D', 'PTPN13', 'GNG4', 'SERPINE2', 'AC008011.2', 'HMHB1', 'LINC01229', 'TSHZ2', 'GK', 'PGM2L1', 'SESN3', 'MFHAS1', 'FKBP5')
immune_markers[["tcell_cd8_TEM"]]<-c('DKK3', 'ENC1', 'GZMK', 'EOMES', 'YBX3', 'GGA2', 'CMC1', 'CD8A', 'CST7', 'CD8B', 'SH2D1A', 'HLA-DPB1', 'CRTAM', 'DUSP2', 'LYST', 'COTL1', 'ZEB2', 'PIK3R1', 'RNF19A', 'TRAT1')
immune_markers[["tcell_cd8_TRMZNF683"]]<-c('CSN2', 'AL590434.1', 'CCL4', 'OR5AU1', 'FGF10-AS1', 'CCL4L2', 'TNFSF9', 'AL133163.2', 'DUSP6', 'FOSB', 'AC020916.1', 'IER5L', 'MZF1-AS1', 'NR4A2', 'IFNG', 'RASGEF1B', 'PLCH1', 'POLQ', 'GZMA', 'GZMK')
immune_markers[["tcell_cd8_TRMAUTS"]]<-c('AUTS2', 'PRSS3', 'GFPT2', 'DAPK2', 'PAGE5', 'ITGA1', 'KLRC1', 'IGFBP3', 'AC020571.1', 'HPSE2', 'CTH', 'COL25A1', 'LINC01871', 'FGF9', 'PARD6G', 'SPRY1', 'SNTG2', 'CD8A', 'ADAMTS6', 'TAC1')
immune_markers[["tcell_cd8_stress"]]<-c( 'HSPA1A', 'DNAJB1', 'HSPH1', 'HSP90AA1', 'HSPB1', 'HSPE1', 'DNAJA1', 'TUBA4A', 'HSP90AB1', 'HSPD1', 'CDKN1C', 'B3GNT5', 'CSF1', 'HSPA8', 'KLF6', 'EGR1', 'POU3F1', 'ANKRD37', 'ZNF80', 'MIR222HG')
immune_markers[["tcell_gdT"]]<-c('KIR2DL3', 'KLRC2', 'KIR2DL1', 'KIR3DL1', 'CHPT1', 'SMC4', 'TRGC2', 'KLRD1', 'XCL1', 'CEMIP2', 'CD7', 'METRNL', 'IKZF2', 'AHI1', 'CTSW', 'KIR2DL4', 'CAST', 'LINC01871', 'ID3', 'PRKX')
immune_markers[["tcell_cycling"]]<-c('DLGAP5', 'E2F8', 'HJURP', 'SPC25', 'RRM2', 'GTSE1', 'AC007240.1', 'BIRC5', 'MYBL2', 'MKI67', 'MCM10', 'KIF4A', 'CENPA', 'SPC24', 'CDC45', 'PCLAF', 'ESCO2', 'KIF20A', 'UBE2C', 'DIAPH3')
immune_markers[["tcell_cd8_interferon"]]<-c( 'IFIT3', 'IFIT1', 'IFI44L', 'HERC5', 'MX2', 'MX1', 'IFI6', 'XAF1', 'ISG15', 'STAT1', 'EIF2AK2', 'LY6E', 'EPSTI1', 'IFI44', 'PARP14', 'OASL', 'OAS2', 'RSAD2', 'SAMD9L', 'MT2A')
immune_markers[["tcell_cd8_exhausted"]]<-c( 'REG4', 'RYR2', 'HAVCR2', 'NPW', 'VCAM1', 'ETV1', 'PNOC', 'KRT86', 'TNS3', 'MYO7A', 'PHEX', 'LRRN3', 'AC243829.4', 'FNDC9', 'BUB1B', 'BEND4', 'AKAP5', 'ASB2', 'NHS', 'RGS13', 'ENTPD1' )
immune_markers[["tcell_nk_FCGR3A"]]<-c('MGAM', 'MYOM2', 'FCGR3A', 'KLRF1', 'PRSS57', 'MLC1', 'AKR1C3', 'FGFBP2', 'SPON2', 'RBPMS2', 'FCER1G', 'KIR3DL1', 'CX3CR1', 'KIR2DL1', 'SH2D1B', 'S1PR5', 'CHST2', 'CCL3', 'CFAP97D2', 'KIR2DL3')
immune_markers[["tcell_nk_cytotoxic"]]<-c('LILRB1', 'CES1', 'PROK2', 'FGFBP2', 'ASCL2', 'CX3CR1', 'FCRL6', 'LINC02384', 'GZMH', 'CFAP97D2', 'LINC01936', 'FCGR3A', 'S1PR5', 'AL590434.1', 'PLEK', 'NKG7', 'ADGRG1', 'AC011933.3', 'PRSS23', 'KIR2DL1')
immune_markers[["tcell_ilc1"]]<-c('ADGRG3', 'CCNJL', 'LINC00996', 'KIR2DL4', 'CSF2', 'B3GNT7', 'FAM166B', 'FES', 'SCT', 'SH2D1B', 'CFAP46', 'AL157895.1', 'KRT86', 'ATP8B4', 'FCER1G', 'GNLY', 'ZNF660', 'TYROBP', 'GOLIM4', 'HOXA6')
immune_markers[["tcell_ilc23"]]<-c('CD300LF', 'KIT', 'SOST', 'LINC00299', 'HOXA6', 'PLA2G4A', 'NPTX2', 'VSTM2L', 'FAM30A', 'SCT', 'XCL1', 'SOX4', 'XCL2', 'AREG', 'TNFRSF18', 'SCN1B', 'TNFSF4', 'MAP3K8', 'SSBP2', 'PLCH1')
immune_markers[["tcell_MAIT"]]<-c('KLRB1', 'RORA', 'KISS1', 'CEBPD', 'TRDV2', 'AC022392.1', 'AC026461.3', 'FEZ1', 'REL', 'GPR65', 'MYBL1', 'LTB', 'MT3', 'SLC4A10', 'NFKB1', 'IGHV3-23', 'CCL20', 'RORC', 'NCR3', 'DPP4')
immune_markers[["tcell_GZMK"]]<-c('ZMAT4', 'KRT13', 'COL24A1', 'CD69', 'COX6A2', 'AL157895.1', 'CCNJL', 'SOST', 'MLC1', 'ADGRG3', 'SMIM24', 'HOXA6', 'SH2D1B', 'TYROBP', 'XCL2', 'GNLY', 'TAL1', 'XCL1', 'AREG', 'FCER1G')

immune_markers[["bcell_naive"]]<-c('IGHD', 'YBX3', 'TCL1A', 'SPRY1', 'CD83', 'BCL7A', 'TXNIP', 'FCER2', 'CLEC2D', 'HLA-DQA2', 'MYB', 'KLF2', 'ESCO2', 'STAG3', 'HLA-DRB5', 'GPR18', 'CMSS1', 'ACSM3', 'GBP4')
immune_markers[["bcell_memory"]]<-c('REL', 'CD83', 'LTB', 'CD69', 'GPR183', 'CCR7', 'CD48', 'ID3', 'BCL2A1', 'KYNU', 'HLA-DRB5', 'PIKFYVE', 'SAMSN1', 'IER5', 'CD80', 'CHORDC1', 'SCIMP', 'DUSP2', 'KLF6', 'NFKB1')
immune_markers[["bcell_interferon"]]<-c('CCL17', 'STAT1', 'IFI44L', 'GBP4', 'ISG15', 'SERPINA9', 'MX1', 'AICDA', 'DUSP4', 'XAF1', 'GBP1', 'ZBED2', 'GBP5', 'GBP2', 'AIM2', 'SAMD9L', 'IFNG', 'MX2', 'ITGAX', 'PLIN2')
immune_markers[["bcell_stress"]]<-c('HSPA1B', 'RRBP1', 'DNAJB1', 'HSPA1A', 'FOS', 'ANKRD28', 'SAMD7', 'H1FX', 'IGKV1-37', 'RGS1', 'HSPB1', 'AC026369.3', 'AMPD1', 'IGKV2-26', 'FAM92B', 'AC021074.3', 'IGHV3-41', 'LINC01405', 'CPEB4', 'IGKV2-4')

immune_markers[["plasma_IGA"]]<-c('JCHAIN', 'IGHA2', 'IGHA1', 'LGALSL', 'NCAM1', 'CPVL', 'CCR10', 'AC106028.4', 'ZNF215', 'IGHV7-81', 'IGKV2-28', 'CADPS2', 'LHX8', 'IGLV2-18', 'IGHV3OR16-8', 'OTP', 'CTSW', 'AC136428.1', 'IGHV3-72')
immune_markers[["plasma_IGG"]]<-c('MZB1', 'IGHG3', 'DERL3', 'FKBP11', 'IGHG4', 'ITM2C', 'IGHG1', 'PRDX4', 'SEC11C', 'XBP1', 'SSR4', 'SSR3', 'JSRP1', 'ERLEC1', 'CYTOR', 'IGKV3D-11', 'SELENOS', 'DPEP1', 'CD38', 'MYDGF')


#set order by umap groupings
cell_markers<-immune_markers
plt_list<-lapply(names(cell_markers),
    function(celltype){
            genes<-unlist(cell_markers[celltype])
            print(paste("Plotting:",celltype))
            print(paste("Genes to plot:",genes))
            plt<-histograModified(obj=dat_fine, 
                baseline="mean",
                genes = unlist(cell_markers[celltype]),
                order=order, promoter_focus=TRUE,trackOverhang=10000,
                matrix = paste0("cg_",broad_celltype,".clusters_tracks"), arrowScale = .03, trackScale = .5,
                legend = F, cgisland=cgi) + ggtitle(celltype)
            ggsave(plt,
            file=paste(project_data_directory,"fine_celltyping",broad_celltype,
            paste("04_scaledcis",broad_celltype,celltype,"fine_celltype.markers.pdf",sep="."),
            sep="/"),
            width=length(genes)*2,
            height=length(order)*2,
            limitsize=F)
        })

#dmr set
dmr<-readRDS(
        paste(project_data_directory,"fine_celltyping",broad_celltype,
        paste(broad_celltype,"dmr",broad_celltype,"clusters.collapsed.rds",sep="."),
        sep="/"))

#rna markers defined by paired 10x RNA analysis
rna<-readRDS("/data/rmulqueen/projects/scalebio_dcis/rna/tenx_dcis.rds")
table(rna$broad_celltype)
rna<-subset(rna,broad_celltype %in% c("tcell","bcell","plasma","myeloid")) #using all stromal types
rna<-subset(rna,fine_celltype %in% c("suspected_doublet"),invert=TRUE)
table(rna$fine_celltype)

rna<-JoinLayers(rna)
markers<-FindAllMarkers(rna,group.by="fine_celltype")

#plot overlap of overexpressed RNA markers for fine cell types and our cluster hypo DMRs
plot_dmr_RNAmarker_overlap(dmr=dmr,
                            rna_markers=markers,
                            celltype=broad_celltype,
                            clus_order=order,
                            plot_value="pval")



#grab cluster dmr hypo markers
clus_markers<-dmr %>% filter(celltype=="7") %>% filter(direction=="hypo") %>% slice_min(dmr_logFC,n=100)
clus_markers<-unlist(
        lapply(
            unlist(
                strsplit(
                    unlist(clus_markers$gene_names),",")),
                function(x) gsub("[[:space:]]", "", x)
            )
        )
#overlap with RNA markers
markers %>% filter(gene %in% clus_markers) %>% filter(avg_log2FC>2) %>% select(cluster) %>% table()



######FINAL IMMUNE CELL TYPING###########
dat_fine@metadata$fine_celltypes<-"tcell"
dat_fine@metadata[dat_fine@metadata[[paste(broad_celltype,"clusters",sep=".")]] %in% c("9"),]$fine_celltypes<-"bcell"
dat_fine@metadata[dat_fine@metadata[[paste(broad_celltype,"clusters",sep=".")]] %in% c("3","7"),]$fine_celltypes<-"myeloid1"
dat_fine@metadata[dat_fine@metadata[[paste(broad_celltype,"clusters",sep=".")]] %in% c("1","8"),]$fine_celltypes<-"myeloid2"


#summarize celltype-dmr clusters over windows
dat_fine<-dmr_and_1kb_window_gen(dat_fine,
    prefix=paste(output_directory,broad_celltype,sep="/"),
    groupBy="fine_celltypes",
    threads=10,step=500)

print(paste("Saving final amethyst object of",broad_celltype))
saveRDS(dat_fine,paste0(output_directory,"/",paste("04_scaledcis",broad_celltype,"fine_celltype.amethyst.rds",sep=".")))



#fix rna marker overlap (i think heatmap doesnt work, but pulling with dplyr does?)
#fix atac overlap of markers





# #plot with atac peak overlaps
# atac_markers<-read.table("/data/rmulqueen/projects/scalebio_dcis/ref/CEDAR_multiome.celltypes.markers.peaks.df.tsv",header=T)
# atac_markers$chr<-unlist(lapply(strsplit(atac_markers$feature,"-"),"[",1))
# atac_markers$start<-unlist(lapply(strsplit(atac_markers$feature,"-"),"[",2))
# atac_markers$end<-unlist(lapply(strsplit(atac_markers$feature,"-"),"[",3))
# clus_order<-order

# plot_dmr_ATACmarker_overlap<-function(dmr,atac_markers,celltype,clus_order,plot_value="Jaccard"){
#     #plot_value can be one of pval, odds.ratio, intersection, union, Jaccard
#     #flow through:
#     # - filter DMRs to HYPO
#     # - filter RNA markers to overexpressed
#     # - table dmr to marker count, and plot odds ratio
#     # - plot gene body methylation hypoM plots
#     output_directory=paste0(project_data_directory,"/fine_celltyping/",celltype)

#     #PLOT HYPO
#     #filter to hypo DMRs (more likely to increase expression)
#     dmr_hypo<- dmr %>% group_by(test) %>% filter(direction=="hypo") 
#     dmr_granges<-makeGRangesFromDataFrame(dmr_hypo,keep.extra.columns=TRUE)
#     dmr_granges<-split(dmr_granges,dmr_granges$test)

#     atac_markers <- atac_markers %>% filter(logFC>0) #look only for increased expression markers
#     atac_markers_granges<-makeGRangesFromDataFrame(atac_markers,keep.extra.columns=TRUE)
#     atac_markers_granges<-split(atac_markers_granges,atac_markers_granges$group)

#     overlap<-as.data.frame(do.call("rbind",lapply(1:length(dmr_granges),calc_overlap_stats)))
#     colnames(overlap)<-c("cluster","celltype","odds.ratio","pval","jaccard")

#     dmr_count<-dmr_hypo %>% group_by(test,.drop=FALSE) %>% tally() 
#     dmr_count<-setNames(dmr_count$n,nm=dmr_count$test)

#     atac_count<-atac_markers %>% group_by(group,.drop=FALSE) %>% tally() 
#     atac_count<-setNames(atac_count$n,nm=atac_count$group)

#     or_mat<-as.data.frame(data.table::dcast(data=as.data.table(overlap),formula=cluster~celltype,value.var=plot_value))
#     row.names(or_mat)<-or_mat$cluster
#     or_mat<-or_mat[,2:ncol(or_mat)]
#     or_mat <- or_mat %>%
#                 mutate_all(~as.numeric(as.character(.)))    
    
#     print("Plotting hypo overlap of DMRs...")
#     pdf(paste0(output_directory,"/","04_scaledcis.",celltype,".dmr_ATAC_overlap.heatmap.pdf"))
#     print(
#         Heatmap(or_mat,
#             bottom_annotation=HeatmapAnnotation(bar=anno_barplot(atac_count)),
#             row_order=clus_order,cluster_rows=FALSE,
#             show_row_names = TRUE, show_column_names = TRUE, row_names_side="right",column_names_side="bottom",
#             column_title = "Hypomethylated and ATAC Peaks",
#             left_annotation=rowAnnotation(bar=anno_barplot(dmr_count)),
#             name=plot_value,
#             cell_fun = function(j, i, x, y, width, height, fill) {
#                                 grid.text(format(round(or_mat,digits=2),nsmall=2)[i, j], x, y, gp = gpar(fontsize = 10))
#                                 }
#             ))
#     dev.off()
# }


# plot_dmr_ATACmarker_overlap(dmr=dmr,
#                             rna_markers=markers,
#                             celltype=broad_celltype,
#                             plot_value="pval")

# #plot percentage cells in fine celltype per patient(rna), and percentage cells in cluster per patient (met)



# # #plot RNA defined markers that overlap with large DMRs
# # plot_dmr_RNAMarker_hypoM(dat=dat_fine,
# #                             dmr=dmr,
# #                             rna_markers=markers,
# #                             broad_celltype=celltype,
# #                             order=c("1","4","2","7","5","9","3","8","6")) #loose ordering by umap grouping

# # #ASSIGN CELL TYPES BY PLOT
# # dat_fine@metadata$fine_celltypes<-"unknown"
# # dat_fine@metadata[dat_fine@metadata[[paste(celltype,"clusters",sep=".")]] %in% c("0"),]$fine_celltypes<-"unknown"

# # #SAVE FINE CELL TYPES
# # saveRDS(dat_fine,paste("04_scaledcis",celltype,"fine_celltype.amethyst.rds",sep="."))

# ```

# #run for lymphoid, myeloid

# ```R

# ######RUNNING###########

# #read in subset amethyst object and dmrs
# dat_fine<-readRDS(
#     paste(project_data_directory,"fine_celltyping",broad_celltype,
#     paste("04_scaledcis",broad_celltype,"fine_celltype.amethyst.rds",sep="."),
#     sep="/"))


# ```


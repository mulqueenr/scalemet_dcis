
#set colors
celltype_col=c(
    "pericyte"="#FF6600",
    "fibroblast"="#FF0066",
    "endothelial"="#FFCC00",
    "unknown"="#666666",
    "suspected_doublet"="#666666",

    "myeloid"="#00FFFF",
    "bcell"="#0099FF",
    "tcell"="#0033FF",

    "basal"="#6600FF",
    "lumsec"="#CC00FF",
    "lumhr"="#FF00CC",
    "cancer"="#00FF99")


race_ethnicity_col=c(
    "African_American"="#0F6FC6",
    "Asian"="#10CF9B",           
    "White"="#A5C249",  
    "Hispanic_or_Latino"="#DBEFF9")


age_col=colorRamp2(breaks=c(min(meta_cat$Age,na.rm=T),
                                median(meta_cat$Age,na.rm=T),
                                max(meta_cat$Age,na.rm=T)),
                                c("#53FFFF","#90A2FF","#FF7BFF"))

class_col=c("+"="black",
            "-"="grey",
            "N/A"="white")

grade_col=c("N/A"="white",
            "G1"="#CCDF92",
            "G2"="#8A9A5B",
            "G3"="#45503B")

differentiation_col=c("N/A"="white",
            "well"="#B7957C",
            "moderate"="#734939",
            "poor"="#A6432D")

menopause_col=c("Hysterectomy (perimenopausal)"="#E2BD6B",
                "Perimenopausal"="#E2BD6B",
                "Post-menopausal(Hysterectomy)"="#4D067B",
                "Post-menopausal"="#4D067B",
                "Pre-menopausal"="#B984DB",
                "Unknown"="white",
                "s/p hysterectomy"="#E2BD6B")

group_col=c("HBCA"="#C7DBD8FF",
            "DCIS"="#C81B74FF",
            "Synchronous"="#993382FF",
            "IDC"="#3D227FFF")

ploidy_col=c("diploid"="#CCCCCC",
            "aneuploid"="#00FF99",
            "NA"="#666666")

#met_col=scale_color_gradient2(low="#ff70ff",mid="#CCCCCC",high="#000000",midpoint=median(epithelial@metadata$mcg_pct))

#coverage_col=scale_color_gradient2(low="#CCCCFF",mid="#3333FF",high="#000066",midpoint=median(log10(epithelial@metadata$unique_reads)))


#Sample colors, too many to make sense of, so using random assignment from high contrast interpolation of Rcolorbrewer
cols <- brewer.pal(n = 12, name = "Paired")
col_func<-colorRampPalette(cols)
sample_col <- col_func(length(unique(obj@metadata$Sample))) 
set.seed(123)
sample_col<-setNames(nm=unique(obj@metadata$Sample),sample(sample_col))



for(i in 1:length(sample_col)){
    print(paste0("'",names(sample_col)[i],"'='",sample_col[[i]],"',"))
}

sample_col=c(
'BCMDCIS05T'='#B294C7',
'BCMDCIS07T'='#DD9A88',
'BCMDCIS102T_24hTis'='#A69C6A',
'BCMDCIS124T'='#5B9EC9',
'BCMDCIS22T'='#F0EB99',
'BCMDCIS28T'='#DBB466',
'BCMDCIS32T'='#FE982C',
'BCMDCIS35T'='#FE870D',
'BCMDCIS41T'='#F78620',
'BCMDCIS49T'='#9774B6',
'BCMDCIS52T'='#2D82AF',
'BCMDCIS65T'='#E42022',
'BCMDCIS66T'='#7D54A5',
'BCMDCIS70T'='#9E8099',
'BCMDCIS74T'='#E8945A',
'BCMDCIS79T_24hTis_DCIS'='#98D277',
'BCMDCIS79T_24hTis_IDC'='#F1E185',
'BCMDCIS80T_24hTis'='#A6D78D',
'BCMDCIS82T_24hTis'='#7EBA98',
'BCMDCIS92T_24hTis'='#75C15D',
'BCMDCIS94T_24hTis'='#754B99',
'BCMDCIS97T'='#CBB0D0',
'BCMDCIS99T'='#3687BC',
'BCMHBCA03R'='#C7B699',
'BCMHBCA04R'='#F16667',
'BCMHBCA09R-3h'='#52AF43',
'BCMHBCA12R-3h'='#FDBB69',
'BCMHBCA16R-3h'='#B15928',
'BCMHBCA17R-3h'='#389F2E',
'BCMHBCA19R-4h'='#F7995C',
'BCMHBCA22R-4h'='#C68647',
'BCMHBCA26L-24hTis-4h'='#F06C45',
'BCMHBCA29L-2h'='#6F9E4C',
'BCMHBCA38L-3h'='#80B6D6',
'BCMHBCA83L-3h'='#F88A89',
'BCMHBCA85L-3h'='#D9A295',
'ECIS25T'='#E83F2E',
'ECIS26T'='#A6CEE3',
'ECIS36T'='#569EA4',
'ECIS48T'='#EA4344',
'ECIS57T'='#FDAA4A')
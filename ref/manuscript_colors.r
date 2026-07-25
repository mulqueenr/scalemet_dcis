
#set colors
celltype_col=c(
    "pericyte"="#FF6600",
    "fibroblast"="#FF0066",
    "endothelial"="#FFCC00",
    "unknown"="#666666",

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

set_pokemon<-list("HBCA"="lapras",
                "DCIS"="abra",
                "Synchronous"="ditto",
                "IDC"="cloyster")

sample_col<-unlist(lapply(1:length(set_pokemon),function(x){
    samples<-unique(obj@metadata[obj@metadata$Group==names(set_pokemon)[x],]$Sample)
    cols <- colorRampPalette(palettetown::pokepal(pokemon = set_pokemon[[x]]))
    sample_col<-setNames(nm=samples,cols(length(samples)))
    return(sample_col)
    }))

for(i in 1:length(sample_colors)){
    print(paste0("'",names(sample_colors)[i],"'='",sample_colors[[i]],"'"))
}

sample_col=c(
'BCMHBCA03R'='#C0E0E8',
'BCMHBCA04R'='#91C5CD',
'BCMHBCA09R-3h'='#527A82',
'BCMHBCA12R-3h'='#436F8B',
'BCMHBCA16R-3h'='#87ADCD',
'BCMHBCA17R-3h'='#D1DEEC',
'BCMHBCA19R-4h'='#1060B0',
'BCMHBCA22R-4h'='#B6DEEB',
'BCMHBCA26L-24hTis-4h'='#58789D',
'BCMHBCA29L-2h'='#2F5BAF',
'BCMHBCA38L-3h'='#6F60A5',
'BCMHBCA83L-3h'='#C62616',
'BCMHBCA85L-3h'='#E84838',
'BCMDCIS05T'='#88C0C8',
'BCMDCIS07T'='#87B3D0',
'BCMDCIS41T'='#385860',
'BCMDCIS66T'='#7898B4',
'BCMDCIS82T_24hTis'='#4880F0',
'BCMDCIS99T'='#83ABD3',
'ECIS25T'='#C02010',
'ECIS26T'='#D3847C',
'ECIS36T'='#F89888',
'BCMDCIS32T'='#68B0E0',
'BCMDCIS35T'='#4C6960',
'BCMDCIS79T_24hTis_DCIS'='#4B5074',
'BCMDCIS80T_24hTis'='#896E39',
'BCMDCIS92T_24hTis'='#8E8C96',
'ECIS48T'='#D0D0D0',
'BCMDCIS102T_24hTis'='#785080',
'BCMDCIS124T'='#986AA5',
'BCMDCIS22T'='#885A95',
'BCMDCIS28T'='#482050',
'BCMDCIS49T'='#A77FB5',
'BCMDCIS52T'='#A58AB0',
'BCMDCIS65T'='#404040',
'BCMDCIS70T'='#757575',
'BCMDCIS74T'='#B2B2B2',
'BCMDCIS79T_24hTis_IDC'='#F8F8F8',
'BCMDCIS94T_24hTis'='#D2D2D2',
'BCMDCIS97T'='#A0A0A0',
'ECIS57T'='#606060')
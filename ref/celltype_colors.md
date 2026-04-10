
#set colors
celltype_col=c(
"perivascular"="#FF9900",
"fibroblast"="#FF0000",
"endothelial"="#FFFF66",
"unknown"="#FF6699",

"monocyte"="#99FFFF",
"macrophage"="#0066FF",
"bcell"="#0099CC",
"tcell_treg"="#99FF99",
"tcell_cd4"="#009966",
"tcell_cd8"="#66FF00",

"basal"="#990099",
"lumsec"="#CC0066",
"lumhr"="#FF00CC",
"cancer"="#00FF99")


race_ethnicity_col=c(
    "African_American"="#0F6FC6",
    "Asian"="#10CF9B",           
    "White"="#A5C249",  
    "Hispanic_or_Latino"="#DBEFF9")


age_col=colorRamp2(breaks=c(min(meta_cat$Age,na.rm=T),median(meta_cat$Age,na.rm=T),max(meta_cat$Age,na.rm=T)),c("#53FFFF","#90A2FF","#FF7BFF"))

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

group_col=c("DCIS"="#278192",
            "HBCA"="#20223E",
            "Synchronous"="#00B089",
            "IDC"="#8FF7BD")
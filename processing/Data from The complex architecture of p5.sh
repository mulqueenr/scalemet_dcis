Data from The complex architecture of p53 binding sites
Alon Senitzki 1, Jessy Safieh 1, Vasundhara Sharma 2, Dmitrij Golovenko 1, Yael Danin-Poleg 1, Alberto Inga 2, Tali E Haran 1

They did not release a motif or PWM formatted file so recreated using their strategy (input into weblogo using their sup table to check). Splitting LHS and RHS of p53.

>Downloaded Supp Table 1, with list of 250 sequences for both LHS and RHS. 
>Input into excel to format as fasta.
>Converted fasta to transfac http://phisite.org/cgi-bin/phisite/pssm-convert.pl
>Input transfac into MEME suite transfac2meme command (CLI)
> Input LHS and RHS merged together through FIMO

fimo --oc . --verbosity 1 --bgfile db/UCSCMammal/mm39.fna.bfile --thresh 0.01 p53_0gap.meme db/UCSCMammal/mm39.fna
fimo --oc . --verbosity 1 --bgfile db/UCSCMammal/mm39.fna.bfile --thresh 0.01 p53_18gap.meme db/UCSCMammal/mm39.fna

used the methods from the paper to generate LHS and RHS p53 motifs
merged together for 0-gap version
added 0.25 per base for each line to extend the gap for 1-18bp gap versions
ran on FIMO, default settings for mm39
convert p53 gapped motifs from results file to bed file

mkdir -p /home/rmulqueen/p53_motif_scan_for_mv/FIMO_bed

for i in $(ls /home/rmulqueen/p53_motif_scan_for_mv/FIMO_Results/fimo_p53_*tsv); 
    do 
    motif_name=$(basename $i | sed s/".tsv"// | sed s/fimo_//)
    awk  -v var="$motif_name" 'NR > 1{ print $2,$3,$4,$5,$6,$7,$8,$9,var}' $i | grep "^chr" > /home/rmulqueen/p53_motif_scan_for_mv/FIMO_bed/${motif_name}.bed; done

#merge them all into a bed file
cat /home/rmulqueen/p53_motif_scan_for_mv/FIMO_bed/*bed | awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' | sort -k1,1 -k2,2n > /home/rmulqueen/p53_motif_scan_for_mv/p53_merged.bed

#define a function to take the gene name from a text file list, expand 5kb up and downstream from gene location, and report overlapping p53 motifs
overlap_motifs(){
    gene_name=$1 #$1 is argument input
    p53_bed="/home/rmulqueen/p53_motif_scan_for_mv/p53_merged.bed" #merged gapped motif scanning 
    mm39_ref="/home/rmulqueen/ref/refdata-gex-GRCm39-2024-A/" #from 10x ref data download
    zcat $mm39_ref/genes/genes.gtf.gz \
    | grep "$gene_name" \
    | awk 'OFS="\t" {print $1,$4,$5,$7}' \
    | sort -k1,1n -k2,2n \
    | bedtools merge | head -n 1 \
    | bedtools slop -i stdin -g $mm39_ref/star/chrNameLength.txt  -l 5000 -r 5000 \
    | bedtools intersect -a stdin -b $p53_bed -wb | awk -v gene=$gene_name 'OFS="\t" {print $0,gene}'
}

export -f overlap_motifs

parallel -j 50 overlap_motifs < p53_targets.txt > p53_gapped_motif_target_overlap.tsv #using parallel to run with 50 jobs at once

#output file is

#Plot motif distributions to define cutoffs

```R
setwd("/home/rmulqueen/p53_motif_scan_for_mv")
library(ggplot2)
library(patchwork)
library(dplyr)

dat_genomewide<-read.csv("/home/rmulqueen/p53_motif_scan_for_mv/p53_merged.bed",sep="\t",col.names=c("motif_chr","motif_start",
 "motif_end","motif_strand","motif_score","motif_pval","motif_qval","seq","gap_size"))

row.names(dat_genomewide)<-paste(dat_genomewide$motif_chr, dat_genomewide$motif_start, dat_genomewide$motif_end, dat_genomewide$motif_strand, dat_genomewide$gap_size,sep="_")

dat_targets<-read.csv("p53_gapped_motif_target_overlap.tsv",sep="\t",col.names=c("gene_chr","gene_start_5kb","gene_end_5kb","motif_chr","motif_start",
 "motif_end","motif_strand","motif_score","motif_pval","motif_qval","seq","gap_size","gene_name"))

dat_targets$name<-paste(dat_targets$motif_chr, dat_targets$motif_start, dat_targets$motif_end, dat_targets$motif_strand, dat_targets$gap_size,sep="_")
dat_targets<-dat_targets[!duplicated(dat_targets$name),]
row.names(dat_targets)<-dat_targets$name

dat_genomewide$in_target<-ifelse(row.names(dat_genomewide) %in% row.names(dat_targets),yes="target",no="not_target")

dat_genomewide$gap_size<-factor(dat_genomewide$gap_size,levels=paste0("p53_",seq(0,18,1),"gap"))
plt<-ggplot(dat_genomewide,aes(x=gap_size,y=motif_score,fill=gap_size))+geom_violin(aes(fill=in_target))+theme_minimal()
ggsave(plt,file="motif_score_estimates_genomewide.pdf",width=20)


dat_targets  %>% filter(motif_qval<=0.05) %>% 
    summarize(mean_qval=mean(motif_qval), mean_pval=mean(motif_pval),mean_motif_score=mean(motif_score), 
                            sd_motif_score=sd(motif_score), median_motif_score=median(motif_score),nrow=n(), genes_with_motifs=n_distinct(gene_name))

group_by(gap_size,in_target) %>% 


```
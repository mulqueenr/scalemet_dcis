[See notebook for wet lab processing.](https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio%20sciMETv2.one%7CD3F046A2-B151-0443-938E-82A415D420EB%2F240109%20Scale%20Met%20Alpha%20Test%20Kit%7C2C45CEED-7824-7C4B-8CFC-697EE8D6A947%2F%29)


```bash
#export environment variables for working/scratch directories
export SCRATCH="/data/rmulqueen/projects/scalebio_dcis/scratch/scalemet_work"
export TMPDIR="/data/rmulqueen/projects/scalebio_dcis/scratch"
export NXF_SINGULARITY_CACHEDIR="/data/rmulqueen/projects/scalebio_dcis/singularity"
export SINGULARITY_BINDPATH="/data/rmulqueen/projects/scalebio_dcis/tools/ScaleMethyl/bin" 

#set up directories and variables
projDir="/data/rmulqueen/projects/scalebio_dcis"
scalebio_nf="/data/rmulqueen/projects/scalebio_dcis/tools/ScaleMethyl" 
params="/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/dcis_runParams.yml"
runDir="${projDir}/data/240202_prelim1"
bclDir="/data/rmulqueen/projects/scalebio_dcis/seq/PM2517/"

```

SAMPLE SPECIFIC TN5 SHEET

```bash

mkdir -p ${runDir}
mkdir -p ${runDir}/samplesheets
cd ${runDir}

echo """sample,barcodes,libName
MCF10A,1A01-1B12,ScaleMethyl
MCF7,1C01-1E12,ScaleMethyl
MDA-MB-231,1F01-1H12,ScaleMethyl
HBCA-16R,2A01-2D12,ScaleMethyl
HBCA-83L,2E01-3H12,ScaleMethyl""" > ${runDir}/samplesheets/samples.csv

```

MAKE SAMPLE SHEETS FOR PCR

```bash

#generate sample sheet using a modified version of bcl_convert_sheet.py to allow for pcr plate specifications.
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.csv \
${projDir}/tools/ScaleMethyl/references/lib.json \
${bclDir}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/scalebio_prelim1_samplesheet.csv

#make nf-core input sheet
cd $runDir
echo """id,samplesheet,flowcell
prelim1,${runDir}/samplesheets/scalebio_prelim1_samplesheet.csv,${bclDir}""" > pipeline_samplesheet.csv

```

BCL to FASTQ

```bash
mkdir -p $SCRATCH/scalemet_prelim1
cd $runDir
#use nf-core to split out fastqs
nextflow run nf-core/demultiplex \
    --input pipeline_samplesheet.csv \
    --outdir fastq \
    --trim_fastq false \
     --remove_adapter false \
    --skip_tools fastp,fastqc,kraken,multiqc,checkqc,falco,md5sum,samshee \
    -profile singularity \
    -w $SCRATCH/scalemet_prelim1

```

RUN SCALE_METHYL PIPELINE

```bash
#remove undetermined ones and empty files
rm -rf ${runDir}/fastq/Undetermined*fastq.gz
find ${runDir}/fastq/ -type f -size 1M -exec rm {} \;

cd ${runDir}
nextflow run ${scalebio_nf} \
--fastqDir ${runDir}/fastq/prelim1 \
--samples ${runDir}/samplesheets/samples.csv \
--outDir ${runDir}/scale_dat \
--maxMemory 300.GB \
--bamOut true \
--trimOut true \
--maxCpus 200 \
-profile singularity \
-params-file ${params} \
-w $SCRATCH/scalemet_prelim1 \
-resume

```

Correct all symlink with actual links and split out merged bam files into single-cell bam files.

```bash

#also copy all symlink files from scalebio nextflow by following symlinks (so we don't need work dir maintained)
find ${runDir} -maxdepth 7 -type l -exec bash -c 'cp -L -R "$(readlink -m "$0")" "$0".dereferenced' {} \; #copy files
find ${runDir} -maxdepth 7 -name "*.dereferenced" -type f -exec bash -c 'mv $0 $(echo $0 | sed -e 's/".dereferenced"//g' -)' {} \; #move to old file names
```

Run initial postprocessing to fit into amethyst object for future merging.

```bash
cd ${runDir}

nextflow run ${projDir}/tools/scalemet_dcis/src/scaleDCIS_postprocessing.nf.groovy \
--runDir ${runDir}/scale_dat \
--maxMemory 300.GB \
--maxForks 100 \
--maxCpus 200 \
--outputPrefix scale \
-w $SCRATCH/scalemet_prelim1 \
-resume
```

```bash
#make splitting barcode list from whitelist
outdir="/data/rmulqueen/projects/scalebio_dcis/data/240202_prelim1/scale_dat/sc_bams_nodedup"
mkdir -p $outdir
indir="/data/rmulqueen/projects/scalebio_dcis/data/240202_prelim1/scale_dat/alignments/bam"

function filter_cells_for_splitting(){
    file_in="$1"
    samtools view $file_in \
    | awk 'OFS="\t" {split($1,a,":");print a[8]}' \
    | sort \
    | uniq -c \
    | awk -v file_in=$file_in ' $1>50000 {print file_in"\t"$2}'
}

export -f filter_cells_for_splitting
parallel -j 100 filter_cells_for_splitting ::: $(find $indir -maxdepth 2 -name "MDA-MB-231*.bam") > $outdir/cell_barcodes.tsv

cd $outdir

function split_bam(){
    file_in=$(echo $1 | awk '{print $1}')
    cellid=$(echo $1 | awk '{print $2}')
    echo $file_in
    echo $cellid
    ((samtools view -H $file_in) & (samtools view $file_in | grep "${cellid}")) | samtools view -b -o MDA231_${cellid}_sc.bam
}
export -f split_bam
parallel -j 100 -a $outdir/cell_barcodes.tsv split_bam

#generate library complexity based on 10% downsample rates
#count unique chr:start sites
function proj_complexity() {
cellid="${1::-4}"
for i in $(seq 0.1 0.1 1.0); do
uniq_count=$(samtools view -F 3332 -s $i $1 \
| awk 'OFS="\t"{print $3,$4}' \
| sort \
| uniq -c \
| wc -l)
total_count=$(samtools view -F 3332 -s $i $1 | wc -l)
echo "${cellid},${i},${total_count},${uniq_count}"; done > ${cellid}.projected_metrics.txt
}

export -f proj_complexity
parallel -j 100 proj_complexity ::: $(ls *.bam)

cat *projected_metrics.txt > ../proj_read_counts.txt
#excluding reads that meet any below conditions:
#read unmapped (0x4)
#not primary alignment (0x100)
#read is PCR or optical duplicate (0x400)
#supplementary alignment (0x800)




```
```R
library(ggplot2)
library(patchwork)
library(drc)
library(parallel)
setwd("/volumes/USR2/Ryan")

scale<-read.table("scalebio_proj_read_counts.txt",sep=",")
colnames(scale)<-c("cellid","downsamp_perc","total_reads","uniq_reads")
scale$sample<-unlist(lapply(strsplit(scale$cellid,"[.]"),"[",1))
scale$method<-"scale"
scale<-scale[!(scale$cellid=="MDA231__sc"),]
kismet<-read.table("kismet_proj_read_counts.txt",sep=",")
colnames(kismet)<-c("cellid","downsamp_perc","total_reads","uniq_reads")
kismet$sample<-"MDA-MB-231"
kismet$method<-"kismet"

compl<-rbind(scale,kismet)

michaelis_menten_fit<-function(x){
    mm<-compl[compl$cellid==x,]
    colnames(mm)<-c("cellid","downsample_perc","S","v","sample","method")
    model.drm <- drm(v ~ S, data = mm, fct = MM.2())
    km_uniq <- data.frame(S = coef(model.drm)[2])
    km_uniq$v <- predict(model.drm, newdata = km_uniq)
    vmax<-as.numeric(coef(model.drm)[1])
    km<-as.numeric(coef(model.drm)[2])
    current_total_reads<-as.numeric(mm[mm$downsample_perc==1.0,]$S)
    current_uniq_reads<-as.numeric(mm[mm$downsample_perc==1.0,]$v)
    method<-mm[mm$downsample_perc==1.0,]$method
    return(c(x,vmax,km,km_uniq$v,current_total_reads,current_uniq_reads,method))
}

projdat<-as.data.frame(do.call("rbind",mclapply(mc.cores=100,unique(compl$cellid),michaelis_menten_fit)))


colnames(projdat)<-c("sample",
"projected_total_fragments",
"projected_optimal_seq_effort",
"projected_reads_at_optimal_effort",
"current_total_reads",
"current_uniq_reads",
"method")

#filter cells with less than 10000 current reads
projdat$current_uniq_reads<-as.numeric(projdat$current_uniq_reads)

projdat<-projdat[projdat$current_uniq_reads>50000,]
projdat<-projdat[!is.na(projdat$current_uniq_reads),]

#these are already sequence exhausted, so don't need projections (also DRC isnt installed)
#dat<-dat[dat$downsample_perc==0.50,]
#dat$perc_uniq<-dat$uniq_frag/dat$total_frag

plt1<-ggplot(projdat,aes(x=method,y=log10(as.numeric(projected_reads_at_optimal_effort)),alpha=0.5,color=method))+
geom_jitter()+
geom_boxplot(outlier.shape=NA)+
theme_minimal()+
ylim(c(0,7))
ggsave(plt1,file="uniq_reads_per_method.pdf")

library(dplyr)
projdat %>% group_by(method) %>% summarize(count=n(),reads=mean(as.numeric(projected_reads_at_optimal_effort),na.rm=T))
```
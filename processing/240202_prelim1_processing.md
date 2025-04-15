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
--maxCpus 200 \
-profile singularity \
-params-file ${params} \
-w $SCRATCH/scalemet_prelim1


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
--maxCpus 200 \
--outputPrefix scale \
-w $SCRATCH/scalemet_prelim1 \
-resume
```







Split bams to single-cells and run copykit

```bash
#set up functions
#count reads export
count_reads() { 
        samtools view $1 | awk -v b=$1 '{split($1,a,":"); print a[8],b}' | sort | uniq -c | sort -k1,1n
}

#split bams export
split_bams() { 
        test=$1
        idx=$(echo $test | cut -d ' ' -f 2 )
        bam=$(echo $test | cut -d ' ' -f 3)
        outprefix=$(echo $bam | awk -F/ '{print $NF}')
        outprefix=$(echo $outprefix | sed -e 's/.dedup.bam//g' -)
        echo ./sc_bams/${outprefix}.${idx}.bam
        ((samtools view -H $bam) && (samtools view $bam | awk -v i=$idx '{split($1,a,":"); if(a[8]==i); print $0}')) | samtools view -bS > ./sc_bams/${outprefix}.${idx}.bam
}

export -f count_reads
export -f split_bams

#filter to bam files with >100000 unique reads
cd ${runDir}/scale_dat
mkdir -p ${runDir}/scale_dat/cnv
mkdir -p ${runDir}/scale_dat/sc_bams

parallel -j 100 count_reads ::: $(find ${runDir}/scale_dat/alignments -maxdepth 5 -name '*bam') | sort -k1,1n > ${runDir}/scale_dat/cnv/scale_unique_read_counts.tsv

awk '$1>100000 {print $0}' ${runDir}/scale_dat/cnv/scale_unique_read_counts.tsv > ${runDir}/scale_dat/cnv/scale_cells.pf.txt

#split bam files to scbams
parallel -j 300 -a ${runDir}/scale_dat/cnv/scale_cells.pf.txt split_bams
```


Run CNV calling per sample

```bash
singularity exec \
--bind /data/rmulqueen/projects/scalebio_dcis/ \
~/singularity/copykit.sif \
Rscript /data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/copykit_cnvcalling.R \
--input_dir ${runDir}/scale_dat/sc_bams \
--output_dir ${runDir}/scale_dat/cnv \
--output_prefix scale \
--task_cpus 125 &

```

Run AMETHYST object initiation per sample
Generates METHYLTREE input for processing as well

```bash
singularity \
exec \
--bind /data/rmulqueen/projects/scalebio_dcis  \
~/singularity/amethyst.sif \
Rscript /data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_initial_processing.R \
--input_dir ${runDir}/scale_dat \
--task_cpus 150

```

Run METHYLTREE per sample, and on all merged

```bash
singularity \
exec \
--bind /data/rmulqueen/projects/scalebio_dcis  \
~/singularity/methyltree.sif \
Rscript /data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_initial_processing.R \
--input_dir ${runDir}/scale_dat/methyltree \
--task_cpus 150
```
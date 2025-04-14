[See notebook for wet lab processing.](https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio%20sciMETv2.one%7CD3F046A2-B151-0443-938E-82A415D420EB%2F240910%20ScaleBio%20DCIS%20Samples%7C30534461-040E-C54F-BB40-7D53F8115495%2F%29)

[See updated experiment notes for homebrew comparison.](https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio%20sciMETv2.one%7CD3F046A2-B151-0443-938E-82A415D420EB%2F250214%20ScaleBio%20Homebrew%20Sorting%20Extra%20Plates%7CA5A091ED-32C0-D24E-81B9-1EFF5EA1252B%2F%29)

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
runDir="${projDir}/data/241007_prelim3"
bclDir="/data/rmulqueen/projects/scalebio_dcis/seq/241004_A01819_0637_BHY5MJDMXY"


```

SAMPLE SPECIFIC TN5 SHEET

```bash

mkdir -p ${runDir}
mkdir -p ${runDir}/samplesheets
cd ${runDir}

echo """sample,barcodes,libName
DCIS-92T,1A01-1D12,ScaleMethyl
DCIS-66T,1E01-1H12,ScaleMethyl
DCIS-79T,2A01-2D12,ScaleMethyl
IDC-79T,2E01-2H12,ScaleMethyl
HBCA-19T,3A01-3D12,ScaleMethyl
HBCA-17T,3E01-3H12,ScaleMethyl""" > ${runDir}/samplesheets/samples.csv

```

MAKE SAMPLE SHEETS FOR PCR

```bash

#generate sample sheet using a modified version of bcl_convert_sheet.py to allow for pcr plate specifications.
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.csv \
${projDir}/tools/ScaleMethyl/references/lib.json \
${bclDir}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/scalebio_prelim3_samplesheet.csv

#make nf-core input sheet
cd $runDir
echo """id,samplesheet,flowcell
prelim3,${runDir}/samplesheets/scalebio_prelim3_samplesheet.csv,${bclDir}""" > pipeline_samplesheet.csv

```

BCL to FASTQ

```bash
mkdir -p $SCRATCH/scalemet_prelim3
cd $runDir
#use nf-core to split out fastqs
nextflow run nf-core/demultiplex \
    --input pipeline_samplesheet.csv \
    --outdir fastq \
    --trim_fastq false \
     --remove_adapter false \
    --skip_tools fastp,fastqc,kraken,multiqc,checkqc,falco,md5sum,samshee \
    -profile singularity \
    -w $SCRATCH/scalemet_prelim3

```

RUN SCALE_METHYL PIPELINE

```bash
#remove undetermined ones and empty files
rm -rf ${runDir}/fastq/prelim3/Undetermined*fastq.gz
find ${runDir}/fastq/prelim3 -type f -size 1M -exec rm {} \;

cd ${runDir}
nextflow run ${scalebio_nf} \
--fastqDir ${runDir}/fastq/prelim3 \
--samples ${runDir}/samplesheets/samples.csv \
--outDir ${runDir}/scale_dat \
--maxMemory 300.GB \
--maxCpus 200 \
-profile singularity \
-params-file ${params} \
-w $SCRATCH/scalemet_prelim3 \
-resume

```

Correct all symlink with actual links and split out merged bam files into single-cell bam files.

```bash

#also copy all symlink files from scalebio nextflow by following symlinks (so we don't need work dir maintained)
find ${runDir} -maxdepth 7 -type l -exec bash -c 'cp -L -R "$(readlink -m "$0")" "$0".dereferenced' {} \; #copy files
find ${runDir} -maxdepth 7 -name "*.dereferenced" -type f -exec bash -c 'mv $0 $(echo $0 | sed -e 's/".dereferenced"//g' -)' {} \; #move to old file names

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
parallel -j 200 -a ${runDir}/scale_dat/cnv/scale_cells.pf.txt split_bams
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
--task_cpus 125

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
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

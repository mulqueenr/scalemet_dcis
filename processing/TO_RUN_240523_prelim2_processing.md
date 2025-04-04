[See notebook for wet lab processing.](https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio%20sciMETv2.one%7CD3F046A2-B151-0443-938E-82A415D420EB%2F240422%20Scale%20Met%20Alpha%20Test%20Kit%202%7CC8153809-4616-E444-A791-2800DD23B717%2F%29)


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
runDir="${projDir}/data/240523_prelim2"
bclDir="/data/rmulqueen/projects/scalebio_dcis/seq/240523_VH00219_594_AAFLYGNM5"


```

SAMPLE SPECIFIC TN5 SHEET

```bash

mkdir -p ${runDir}
mkdir -p ${runDir}/samplesheets
cd ${runDir}

echo """sample,barcodes,libName
BCMDCIS41T,3A01-3A08,ScaleMethyl
BCMDCIS41T,3B01-3B08,ScaleMethyl
BCMDCIS41T,3C01-3C08,ScaleMethyl
BCMDCIS41T,3D01-3D08,ScaleMethyl
BCMDCIS41T,3E01-3E08,ScaleMethyl
BCMDCIS41T,3F01-3F08,ScaleMethyl
BCMDCIS41T,3G01-3G08,ScaleMethyl
BCMDCIS41T,3H01-3H08,ScaleMethyl
BCMDCIS66T,3A09-3A11,ScaleMethyl
BCMDCIS66T,3B09-3B11,ScaleMethyl
BCMDCIS66T,3C09-3C11,ScaleMethyl
BCMDCIS66T,3D09-3D11,ScaleMethyl
BCMDCIS66T,3E09-3E11,ScaleMethyl
BCMDCIS66T,3F09-3F11,ScaleMethyl
BCMDCIS66T,3G09-3G11,ScaleMethyl
BCMDCIS66T,3H09-3H11,ScaleMethyl
BCMDCIS81T,3A12-3A12,ScaleMethyl
BCMDCIS81T,3B12-3B12,ScaleMethyl
BCMDCIS81T,3C12-3C12,ScaleMethyl
BCMDCIS81T,3D12-3D12,ScaleMethyl
BCMDCIS81T,3E12-3E12,ScaleMethyl
BCMDCIS81T,3F12-3F12,ScaleMethyl
BCMDCIS81T,3G12-3G12,ScaleMethyl
BCMDCIS81T,3H12-3H12,ScaleMethyl""" > ${runDir}/samplesheets/samples.csv

```

MAKE SAMPLE SHEETS FOR PCR

```bash

#generate sample sheet using a modified version of bcl_convert_sheet.py to allow for pcr plate specifications.
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.csv \
${projDir}/tools/ScaleMethyl/references/lib.json \
${bclDir}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/scalebio_prelim2_samplesheet.csv

#make nf-core input sheet
cd $runDir
echo """id,samplesheet,flowcell
prelim2,${runDir}/samplesheets/scalebio_prelim2_samplesheet.csv,${bclDir}""" > pipeline_samplesheet.csv

```

BCL to FASTQ

```bash
mkdir -p $SCRATCH/scalemet_prelim2
cd $runDir
#use nf-core to split out fastqs
nextflow run nf-core/demultiplex \
    --input pipeline_samplesheet.csv \
    --outdir fastq \
    --trim_fastq false \
     --remove_adapter false \
    --skip_tools fastp,fastqc,kraken,multiqc,checkqc,falco,md5sum,samshee \
    -profile singularity \
    -w $SCRATCH/scalemet_prelim2

```

RUN SCALE_METHYL PIPELINE

```bash
#remove undetermined ones and empty files
rm -rf ${runDir}/fastq/Undetermined*fastq.gz
find ${runDir}/fastq/ -type f -size 1M -exec rm {} \;

cd ${runDir}
nextflow run ${scalebio_nf} \
--fastqDir ${runDir}/fastq/prelim2 \
--samples ${runDir}/samplesheets/samples.csv \
--outDir ${runDir}/scale_dat \
--maxMemory 300.GB \
--maxCpus 200 \
-profile singularity \
-params-file ${params} \
-w $SCRATCH/scalemet_prelim2


```

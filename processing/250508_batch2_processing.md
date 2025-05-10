[See notebook for wet lab processing.](https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/Documents/metACT/scalebio_dcis_processing.one#250501%20Batch%202%20Plates%201,2&section-id={A0170AC0-CBA6-B848-B58D-A9AF996E46F1}&page-id={B1598264-521A-9243-9A8E-C85848922058}&end)


[Sample selection and plate processing notes:](/volumes/USR2/Ryan/projects/scalebio_dcis/sample_selection/dcis_sample_selection.xlsx)

Set env variables.

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
runDir="${projDir}/data/250508_RM_scalebio_batch2_initseq"
bclDir="/data/rmulqueen/projects/scalebio_dcis/seq/250506_VH01788_104_222CF7LNX"

```

SAMPLE SPECIFIC TN5 SHEET

```bash
mkdir -p ${runDir}
mkdir -p ${runDir}/samplesheets

cd ${runDir}
samples="${runDir}/samplesheets/samples.csv"

echo """sample,barcodes,libName
BCMHBCA12R-3h,1A01-1A06;2A07-2A12;3A01-3A06,ScaleMethyl
BCMDCIS65T,1B01-1B06;2B07-2B12;3B01-3B06;1C07-1C12;2B01-2B06;3C07-3C12,ScaleMethyl
BCMDCIS05T,1C01-1C06;2C07-2C12;3C01-3C06;3D11-3D12,ScaleMethyl
BCMDCIS102T-4h,1D01-1D06;2D07-2D12;3D01-3D06,ScaleMethyl
BCMHBCA04R,1E01-1E06;2E07-2E12;3E01-3E06,ScaleMethyl
ECIS25T,1F01-1F06;2F07-2F12;3F01-3F06,ScaleMethyl
BCMDCIS35T-3h,1G01-1G06;2G07-2G12;3G01-3G06,ScaleMethyl
BCMDCIS97T,1H01-1H06;2H07-2H12;3H01-3H06,ScaleMethyl
BCMHBCA29L-2h,1A07-1A12;2A01-2A06;3A07-3A12,ScaleMethyl
BCMHBCA85L-3h,1B07-1B12;2B01-2B06;3B07-3B12,ScaleMethyl
BCMDCIS49T-24hTis,1D07-1D12;2D01-2D06;3D07-3D10,ScaleMethyl
ECIS48T,1E07-1E12;2E01-2E06;3E07-3E12,ScaleMethyl
BCMDCIS70T,1F07-1F12;2F01-2F06;3F07-3F12,ScaleMethyl
BCMHBCA22R-4h,1G07-1G12;2G01-2G06;3G07-3G12,ScaleMethyl
ECIS26T,1H07-1H12;2H01-2H06;3H07-3H12,ScaleMethyl""" > ${samples}
```

MAKE SAMPLE SHEETS FOR PCR

```bash
#generate sample sheet using a modified version of bcl_convert_sheet.py to allow for pcr plate specifications.
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.csv \
${projDir}/tools/ScaleMethyl/references/lib.json \
${bclDir}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/scalebio_batch2_plate3-4_samplesheet.csv

#sample sheets for homebrew plates, run one plate at a time and merge
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.csv \
${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
${bclDir}/RunInfo.xml \
--splitFastq \
--i7Set D \
--i5Set 1 > ${runDir}/samplesheets/homebrew_batch2_plate1_samplesheet.csv

python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.csv \
${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
${bclDir}/RunInfo.xml \
--splitFastq \
--i7Set D \
--i5Set 2 > ${runDir}/samplesheets/homebrew_batch2_plate2_samplesheet.csv

cat ${runDir}/samplesheets/homebrew_batch2_plate1_samplesheet.csv > ${runDir}/samplesheets/homebrew_batch2_plate1-2_samplesheet.csv
grep "^ScaleMethyl_" ${runDir}/samplesheets/homebrew_batch2_plate2_samplesheet.csv >> ${runDir}/samplesheets/homebrew_batch2_plate1-2_samplesheet.csv

#make nf-core input sheet
cd $runDir
echo """id,samplesheet,flowcell
homebrew_batch2_plate1-2,${runDir}/samplesheets/homebrew_batch2_plate1-2_samplesheet.csv,${bclDir}
scale_batch2_plate3-4,${runDir}/samplesheets/scalebio_batch2_plate3-4_samplesheet.csv,${bclDir}""" > pipeline_samplesheet.csv

```

BCL to FASTQ

```bash
mkdir -p $SCRATCH/scalemet_work

#use nf-core to split out fastqs
nextflow run nf-core/demultiplex \
    --input pipeline_samplesheet.csv \
    --outdir fastq \
    --trim_fastq false \
     --remove_adapter false \
    --skip_tools fastp,fastqc,kraken,multiqc,checkqc,falco,md5sum,samshee \
    -profile singularity \
    -w ${SCRATCH}/scalemet_work

```

RUN SCALE_METHYL PIPELINE

```bash
mkdir -p ${SCRATCH}/scalemet_native_work

cd ${runDir}
nextflow run ${scalebio_nf} \
--fastqDir ${runDir}/fastq/homebrew_batch2_plate1-2 \
--samples ${runDir}/samplesheets/samples.csv \
--outDir ${runDir}/homebrew_dat \
--libStructure ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
--maxMemory 300.GB \
--maxCpus 200 \
-profile singularity \
-params-file ${params} \
-w ${SCRATCH}/scalemet_homebrew_work


cd ${runDir}
nextflow run ${scalebio_nf} \
--fastqDir ${runDir}/fastq/scale_batch2_plate3-4 \
--samples ${runDir}/samplesheets/samples.csv \
--outDir ${runDir}/scale_dat \
--maxMemory 300.GB \
--maxCpus 200 \
-profile singularity \
-params-file ${params} \
-w ${SCRATCH}/scalemet_native_work
```


Correct all symlink with actual links and split out merged bam files into single-cell bam files.

```bash

#also copy all symlink files from scalebio nextflow by following symlinks (so we don't need work dir maintained)
find ${runDir} -maxdepth 7 -type l -exec bash -c 'cp -L -R "$(readlink -m "$0")" "$0".dereferenced' {} \; #copy files
find ${runDir} -maxdepth 7 -name "*.dereferenced" -type f -exec bash -c 'mv $0 $(echo $0 | sed -e 's/".dereferenced"//g' -)' {} \; #move to old file names
```


```bash
cd ${runDir}

nextflow run ${projDir}/tools/scalemet_dcis/src/scaleDCIS_postprocessing.nf.groovy \
--runDir ${runDir}/homebrew_dat \
--maxMemory 300.GB \
--maxForks 100 \
--maxCpus 200 \
--outputPrefix homebrew \
-w $SCRATCH/scalemet_batch2_homebrew

nextflow run ${projDir}/tools/scalemet_dcis/src/scaleDCIS_postprocessing.nf.groovy \
--runDir ${runDir}/scale_dat \
--maxMemory 300.GB \
--maxForks 100 \
--maxCpus 200 \
--outputPrefix scale \
-w $SCRATCH/scalemet_batch2_scale \
-resume

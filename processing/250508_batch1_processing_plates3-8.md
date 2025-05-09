[See notebook for wet lab processing.](https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio_dcis_processing.one%7CA0170AC0-CBA6-B848-B58D-A9AF996E46F1%2F%29)

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
runDir="${projDir}/data/250508_RM_scalebio_batch1_plates3_8"
bclDir="/data/rmulqueen/projects/scalebio_dcis/seq/PM2563"

```

SAMPLE SPECIFIC TN5 SHEET

```bash
mkdir -p ${runDir}
mkdir -p ${runDir}/samplesheets

cd ${runDir}
samples="${runDir}/samplesheets/samples.csv"

echo """sample,barcodes,libName
BCMDCIS32T,1A01-1A06;2A07-2A12;3A01-3A03,ScaleMethyl
BCMHBCA03R,1B01-1B06;2B07-2B12;3B01-3B03;3E04-3E06;3E10-3E12,ScaleMethyl
BCMDCIS94T-24hTis,1C01-1C06;2C07-2C12;3C01-3C03,ScaleMethyl
BCMHBCA26L-24hTis-4h,1D01-1D06;2D07-2D12;3D01-3D03;3G01-3G06;3G10-3G12;1G01-1G06;2G07-2G12,ScaleMethyl
ECIS36T,1E01-1E06;2E07-2E12;3E01-3E03;3F04-3F06;3F10-3F12,ScaleMethyl
BCMDCIS74T,1F01-1F06;2F07-2F12;3F01-3F03,ScaleMethyl
BCMDCIS22T,1H01-1H06;2H07-2H12;3H01-3H03,ScaleMethyl
BCMDCIS66T,1A07-1A12;2A01-2A06;3A04-3A12,ScaleMethyl
BCMDCIS28T,1B07-1B12;2B01-2B06;3B04-3B12,ScaleMethyl
BCMDCIS99T,1C07-1C12;2C01-2C06;3C07-3C09,ScaleMethyl
BCMDCIS52T,1D07-1D12;1F07-1F12;2D01-2D06;2F01-2F06;3H04-3H06;3D07-3D09;3F07-3F09;3H10-3H12,ScaleMethyl
BCMDCIS80T,1E07-1E12;2E01-2E06;3D04-3D06;3E07-3E09;3D10-3D12,ScaleMethyl
BCMDCIS07T,1G07-1G12;2G01-2G06;3C04-3C06;3G07-3G09;3C10-3C12,ScaleMethyl
BCMHBCA38L-3h,1H07-1H12;2H01-2H06;3H07-3H09,ScaleMethyl""" > ${samples}
```

MAKE SAMPLE SHEETS FOR PCR

```bash
#generate sample sheet using a modified version of bcl_convert_sheet.py to allow for pcr plate specifications.
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.csv \
${projDir}/tools/ScaleMethyl/references/lib.json \
${bclDir}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/scalebio_batch1_plate3-8_samplesheet.csv

#sample sheets for homebrew plates, run one plate at a time and merge
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.csv \
${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
${bclDir}/RunInfo.xml \
--splitFastq \
--i7Set A \
--i5Set 3 > ${runDir}/samplesheets/homebrew_batch1_plate3_samplesheet.csv

python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.csv \
${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
${bclDir}/RunInfo.xml \
--splitFastq \
--i7Set C \
--i5Set 3 > ${runDir}/samplesheets/homebrew_batch1_plate8_samplesheet.csv

cat ${runDir}/samplesheets/homebrew_batch1_plate3_samplesheet.csv > ${runDir}/samplesheets/homebrew_batch1_plate3-8_samplesheet.csv
grep "^ScaleMethyl_" ${runDir}/samplesheets/homebrew_batch1_plate8_samplesheet.csv >> ${runDir}/samplesheets/homebrew_batch1_plate3-8_samplesheet.csv

#replace override cycles (just using sed) and also only run on lane 5 of run (from pipeline samplesheet)
sed -i 's/OverrideCycles,Y148;I34;I10;Y146/OverrideCycles,Y148;I10N24;I10;Y146/g' \
${runDir}/samplesheets/homebrew_batch1_plate3-8_samplesheet.csv

#make nf-core input sheet
cd $runDir
echo """id,samplesheet,lane,flowcell
homebrew_batch1_plate3-8,${runDir}/samplesheets/homebrew_batch1_plate3-8_samplesheet.csv,5,${bclDir}""" > pipeline_samplesheet.csv

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

nextflow run ${scalebio_nf} \
--fastqDir ${runDir}/fastq/homebrew_batch1_plate3-8/L005 \
--samples ${runDir}/samplesheets/samples.csv \
--libStructure ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
--outDir ${runDir}/homebrew_dat \
--maxMemory 300.GB \
--maxCpus 200 \
-profile singularity \
-params-file ${params} \
-w ${SCRATCH}/scalemet_homebrew_work

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
-w $SCRATCH/scalemet_batch1_homebrew \
-resume

nextflow run ${projDir}/tools/scalemet_dcis/src/scaleDCIS_postprocessing.nf.groovy \
--runDir ${runDir}/scale_dat \
--maxMemory 300.GB \
--maxForks 100 \
--maxCpus 200 \
--outputPrefix scale \
-w $SCRATCH/scalemet_batch1_scale \
-resume

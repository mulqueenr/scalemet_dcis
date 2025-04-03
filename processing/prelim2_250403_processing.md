[See notebook for wet lab processing.](https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio%20sciMETv2.one%7CD3F046A2-B151-0443-938E-82A415D420EB%2F240422%20Scale%20Met%20Alpha%20Test%20Kit%202%7CC8153809-4616-E444-A791-2800DD23B717%2F%29)


```bash
#export environment variables for working/scratch directories
export SCRATCH="/data/rmulqueen/projects/scalebio_dcis/scratch/scalemet_work"
export TMPDIR="/data/rmulqueen/projects/scalebio_dcis/scratch"
export NXF_SINGULARITY_CACHEDIR="/data/rmulqueen/projects/scalebio_dcis/singularity"
export SINGULARITY_BINDPATH="/data/rmulqueen/projects/scalebio_dcis/tools/ScaleMethyl/bin" 
mkdir -p $SCRATCH

#set up directories and variables
projDir="/data/rmulqueen/projects/scalebio_dcis"
scalebio_nf="/data/rmulqueen/projects/scalebio_dcis/tools/ScaleMethyl" 
params="/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/dcis_runParams.yml"
runDir="${projDir}/data/241007_RM_scalebio_dcis2"
bclDir="/data/rmulqueen/projects/scalebio_dcis/seq/240523_VH00219_594_AAFLYGNM5"

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

#generate sample sheet using a modified version of bcl_convert_sheet.py to allow for pcr plate specifications.
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.csv \
${projDir}/tools/ScaleMethyl/references/lib.json \
${bclDir}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/scalebio_batchprelim_plate1-2_samplesheet.csv

#run bcl-convert within amethyst sif
cd ${runDir}
singularity shell \
--bind ${bclDir} \
--bind ${runDir} \
--bind ${bclDir}:/var/log/bcl-convert \
~/singularity/amethyst.sif

projDir="/data/rmulqueen/projects/scalebio_dcis"
bclDir="/data/rmulqueen/projects/scalebio_dcis/seq/241004_A01819_0637_BHY5MJDMXY"
runDir="${projDir}/data/241007_RM_scalebio_dcis2"
cd ${runDir}

bcl-convert \
--bcl-input-directory ${bclDir} \
--output-directory ${runDir}/fastq \
--no-lane-splitting true \
--bcl-num-parallel-tiles 1 \
--bcl-num-conversion-threads 10 \
--bcl-num-compression-threads 10 \
--bcl-num-decompression-threads 10 \
--sample-sheet ${runDir}/samplesheets/scalebio_batchprelim_plate1-2_samplesheet.csv \
--force

#remove undetermined ones and empty files
rm -rf ${runDir}/fastq/Undetermined*fastq.gz
find ${runDir}/fastq/ -type f -size 1M -exec rm {} \;

cd $runDir
nextflow run ${scalebio_nf} \
--fastqDir ${runDir}/fastq \
--samples ${runDir}/samplesheets/samples.csv \
--outDir ${runDir}/dat \
--maxMemory 200.GB \
--maxCpus 100 \
-profile singularity \
-params-file ${params} \
-w ${SCRATCH}/scalemet_work \
-resume


```

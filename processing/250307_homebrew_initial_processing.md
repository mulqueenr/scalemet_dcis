[See notebook for wet lab processing.](https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio%20sciMETv2.one%7CD3F046A2-B151-0443-938E-82A415D420EB%2F240910%20ScaleBio%20DCIS%20Samples%7C30534461-040E-C54F-BB40-7D53F8115495%2F%29)

[See updated experiment notes for homebrew comparison.](https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio%20sciMETv2.one%7CD3F046A2-B151-0443-938E-82A415D420EB%2F250214%20ScaleBio%20Homebrew%20Sorting%20Extra%20Plates%7CA5A091ED-32C0-D24E-81B9-1EFF5EA1252B%2F%29)

```bash
#Should still only need to identify samples at the tagmentation level, and expanding the i5.txt and i7.txt should take care of itself.
#export environment variables for working/scratch directories
export SCRATCH="/home/rmulqueen/projects/scalebio_dcis/scratch/scalemet_work"
export TMPDIR="/home/rmulqueen/projects/scalebio_dcis/scratch"
export NXF_SINGULARITY_CACHEDIR="/home/rmulqueen/projects/scalebio_dcis/singularity"
export SINGULARITY_BINDPATH="/home/rmulqueen/projects/scalebio_dcis/tools/ScaleMethyl/bin" 
mkdir -p $SCRATCH

#set up directories and variables
projDir="/home/rmulqueen/projects/scalebio_dcis"
scalebio_nf="/home/rmulqueen/projects/scalebio_dcis/tools/ScaleMethyl" 
params="/home/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/dcis_runParams.yml"
runDir="${projDir}/data/250307_RM_scalebio_dcis2_homebrew"
bclDir="/home/rmulqueen/projects/scalebio_dcis/seq/250227_RM_CurioWGS_scalemet"

mkdir -p ${runDir}
cd ${runDir}
samples="${runDir}/samples.csv"

echo """sample,barcodes,libName
DCIS-92T,1A01-1D12,ScaleMethyl
DCIS-66T,1E01-1H12,ScaleMethyl
DCIS-79T,2A01-2D12,ScaleMethyl
IDC-79T,2E01-2H12,ScaleMethyl
HBCA-19T,3A01-3D12,ScaleMethyl
HBCA-17T,3E01-3H12,ScaleMethyl""" > ${samples}

#generate sample sheet using a modified version of bcl_convert_sheet.py to allow for pcr plate specifications.
source activate 
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
samples.csv \
${projDir}/tools/ScaleMethyl/references/lib.json \
${bclDir}/RunInfo.xml --splitFastq \
--i7Set C \
--i5Set 3 > samplesheet.csv

#update Override Cycles in samplesheet to fit sequencing run
#OverrideCycles,Y50;I10;I10N14;Y47

bcl-convert \
--bcl-input-directory ${bclDir} \
--output-directory ${runDir}/fastq \
--bcl-num-conversion-threads 50 \
--bcl-num-compression-threads 50 \
--bcl-num-decompression-threads 50 \
--no-lane-splitting true \
--sample-sheet samplesheet.csv \
--force

#remove empty ones
rm -rf ${runDir}/fastq/Undetermined*fastq.gz

source activate conda #(to use more recent java version)
cd $runDir

nextflow run ${scalebio_nf} \
--fastqDir ${runDir}/fastq \
--samples ${runDir}/samples.csv \
--outDir ${runDir} \
--maxMemory 300.GB \
--maxCpus 200 \
-profile singularity \
-params-file ${params} \
-w ${SCRATCH}/scalemet_work \
-resume



```

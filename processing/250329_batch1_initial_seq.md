[See notebook for wet lab processing.](https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio_dcis_processing.one%7CA0170AC0-CBA6-B848-B58D-A9AF996E46F1%2F%29)

[Sample selection and plate processing notes:](/volumes/USR2/Ryan/projects/scalebio_dcis/sample_selection/dcis_sample_selection.xlsx)

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
runDir="${projDir}/data/250329_RM_scalebio_batch1_initseq"
bclDir="/home/rmulqueen/projects/scalebio_dcis/seq/250325_VH01788_97_2227YC2NX"

mkdir -p ${runDir}
cd ${runDir}
samples="${runDir}/samples.csv"

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

#generate sample sheet using a modified version of bcl_convert_sheet.py to allow for pcr plate specifications.
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
samples.csv \
${projDir}/tools/ScaleMethyl/references/lib.json \
${bclDir}/RunInfo.xml \
--splitFastq \
--i7Set A,B \
--i5Set 1 > scalebio_kit_samplesheet.csv

#generate sample sheet using a modified version of bcl_convert_sheet.py to allow for pcr plate specifications.
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
samples.csv \
${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
${bclDir}/RunInfo.xml \
--splitFastq \
--i7Set A,B \
--i5Set 1,2 > homebrew_samplesheet.csv

#run bcl-convert within amethyst sif
cd ${runDir}
singularity shell \
--bind ${bclDir} \
--bind ${runDir} \
--bind ${bclDir}:/var/log/bcl-convert \
~/singularity/amethyst.sif

#scalebio kit fq split
bcl-convert \
--bcl-input-directory ${bclDir} \
--output-directory ${runDir}/scale_fastq \
--bcl-num-conversion-threads 30 \
--bcl-num-compression-threads 30 \
--bcl-num-decompression-threads 30 \
--shared-thread-odirect-output true \
--no-lane-splitting true \
--sample-sheet scalebio_kit_samplesheet.csv \
--force

nextflow run ${scalebio_nf} \
--fastqDir ${runDir}/scale_fastq \
--samples ${runDir}/samples.csv \
--outDir ${runDir} \
--maxMemory 300.GB \
--maxCpus 200 \
-profile singularity \
-params-file ${params} \
-w ${SCRATCH}/scalemet_work \
-resume


#homebrew fq split
bcl-convert \
--bcl-input-directory ${bclDir} \
--output-directory ${runDir}/homebrew_fastq \
--bcl-num-conversion-threads 10 \
--bcl-num-compression-threads 10 \
--bcl-num-decompression-threads 10 \
--bcl-num-parallel-tiles 1 \
--shared-thread-odirect-output true \
--no-lane-splitting true \
--sample-sheet homebrew_samplesheet.csv \
--force




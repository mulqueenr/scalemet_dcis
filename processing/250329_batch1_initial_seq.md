[See notebook for wet lab processing.](https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio_dcis_processing.one%7CA0170AC0-CBA6-B848-B58D-A9AF996E46F1%2F%29)

[Sample selection and plate processing notes:](/volumes/USR2/Ryan/projects/scalebio_dcis/sample_selection/dcis_sample_selection.xlsx)

Note: this should be rerun, the run transfer has some missing filter files which are costing us some tiles.

Set env variables.

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
```

Sample sheet setup for Tn5

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

Scalebio kit processing plates 1 and 6
```bash
#generate sample sheet using a modified version of bcl_convert_sheet.py to allow for pcr plate specifications.
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.csv \
${projDir}/tools/ScaleMethyl/references/lib.json \
${bclDir}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/scalebio_batch1_plate1-6_samplesheet.csv

#run bcl-convert within amethyst sif
cd ${runDir}
singularity shell \
--bind ${bclDir} \
--bind ${runDir} \
--bind ${bclDir}:/var/log/bcl-convert \
~/singularity/amethyst.sif

projDir="/home/rmulqueen/projects/scalebio_dcis"
scalebio_nf="/home/rmulqueen/projects/scalebio_dcis/tools/ScaleMethyl" 
params="/home/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/dcis_runParams.yml"
runDir="${projDir}/data/250329_RM_scalebio_batch1_initseq"
bclDir="/home/rmulqueen/projects/scalebio_dcis/seq/250325_VH01788_97_2227YC2NX"

#scalebio kit fq split
bcl-convert \
--bcl-input-directory ${bclDir} \
--output-directory ${runDir}/scale_fastq \
--shared-thread-odirect-output true \
--no-lane-splitting true \
--bcl-num-parallel-tiles 1 \
--bcl-num-conversion-threads 10 \
--bcl-num-compression-threads 10 \
--bcl-num-decompression-threads 10 \
--sample-sheet ${runDir}/samplesheets/scalebio_batch1_plate1-6_samplesheet.csv \
--force

mkdir -p ${SCRATCH}/scalemet_native_work
nextflow run ${scalebio_nf} \
--fastqDir ${runDir}/scale_fastq \
--samples ${runDir}/samplesheets/samples.csv \
--outDir ${runDir}/scale_dat \
--maxMemory 300.GB \
--maxCpus 100 \
-profile singularity \
-params-file ${params} \
-w ${SCRATCH}/scalemet_native_work \
-resume

```

Homebrew processing plates 2 and 7

```bash
#maybe rewrite the python script to just be pairwise combo of plates rather than all combos or just run for each combo and append together
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
samples.csv \
${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
${bclDir}/RunInfo.xml \
--splitFastq \
--i7Set A \
--i5Set 1 > ${runDir}/samplesheets/homebrew_batch1_plate2_samplesheet.csv

python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
samples.csv \
${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
${bclDir}/RunInfo.xml \
--splitFastq \
--i7Set B \
--i5Set 2 > ${runDir}/samplesheets/homebrew_batch1_plate7_samplesheet.csv

cat ${runDir}/samplesheets/homebrew_batch1_plate2_samplesheet.csv > ${runDir}/samplesheets/homebrew_batch1_plate2-7_samplesheet.csv
grep "^ScaleMethyl_" ${runDir}/samplesheets/homebrew_batch1_plate7_samplesheet.csv >> ${runDir}/samplesheets/homebrew_batch1_plate2-7_samplesheet.csv

#set up directories and variables
projDir="/home/rmulqueen/projects/scalebio_dcis"
runDir="${projDir}/data/250329_RM_scalebio_batch1_initseq"
scalebio_nf="/home/rmulqueen/projects/scalebio_dcis/tools/ScaleMethyl" 
params="/home/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/dcis_runParams.yml"
bclDir="/home/rmulqueen/projects/scalebio_dcis/seq/250325_VH01788_97_2227YC2NX"

#run bcl-convert within amethyst sif
cd ${runDir}
singularity shell \
--bind ${bclDir} \
--bind ${runDir} \
--bind ${bclDir}:/var/log/bcl-convert \
~/singularity/amethyst.sif

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
--sample-sheet ${runDir}/samplesheets/homebrew_batch1_plate2-7_samplesheet.csv \
--force


nextflow run ${scalebio_nf} \
--fastqDir ${runDir}/homebrew_fastq \
--samples ${runDir}/samplesheets/samples.csv \
--outDir ${runDir}/homebrew_dat \
--libStructure ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
--maxMemory 300.GB \
--maxCpus 100 \
-profile singularity \
-params-file ${params} \
-w ${SCRATCH}/scalemet_work \
-resume

```

Correct all symlink with actual links and split out merged bam files into single-cell bam files.

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
        outprefix=$(echo $bam | cut -d '/' -f 4)
        outprefix=$(echo $outprefix | sed -e 's/.dedup.bam//g' -)
        echo ./sc_bams/${outprefix}.${idx}.bam
        ((samtools view -H $bam) && (samtools view $bam | awk -v i=$idx '{split($1,a,":"); if(a[8]==i); print $0}')) | samtools view -bS > ./sc_bams/${outprefix}.${idx}.bam
}

export -f count_reads
export -f split_bams

```
Run on scalebio plates

```bash
scale_runDir="/home/rmulqueen/projects/scalebio_dcis/data/250329_RM_scalebio_batch1_initseq/scale_dat"

#filter to bam files with >100000 unique reads
cd $scale_runDir
parallel -j 100 count_reads ::: $(find ./ -maxdepth 4 -name '*bam' -type l ) | sort -k1,1n > scale_unique_read_counts.tsv
awk '$1>100000 {print $0}' scale_unique_read_counts.tsv > scale_cells.pf.txt
mkdir -p sc_bams
parallel -j 100 -a scale_cells.pf.txt split_bams

#also copy all symlink files from scalebio nextflow by following symlinks (so we don't need work dir maintained)
find ${scale_runDir} -maxdepth 5 -type l -exec bash -c 'cp -L -R "$(readlink -m "$0")" "$0".dereferenced' {} \; #copy files
find ${scale_runDir} -maxdepth 5 -name "*.dereferenced" -type f -exec bash -c 'mv $0 $(echo $0 | sed -e 's/".dereferenced"//g' -)' {} \; #move to old file names

#clear space from scratch
rm -rf ${SCRATCH}/scalemet_native_work

mkdir -p ${scale_runDir}/copykit_cnv
singularity exec \
--bind /home/rmulqueen/projects/scalebio_dcis/ \
~/singularity/copykit.sif \
Rscript ~/projects/scalebio_dcis/tools/scalemet_dcis/src/copykit_cnvcalling.R \
--input_dir ${scale_runDir}/sc_bams \
--output_prefix ${scale_runDir}/copykit_cnv/scale \
--task_cpus 150

```

Run on homebrew plates

```bash
homebrew_runDir="/home/rmulqueen/projects/scalebio_dcis/data/250329_RM_scalebio_batch1_initseq/homebrew_dat"

#filter to bam files with >100000 unique reads
cd $homebrew_runDir
parallel -j 100 count_reads ::: $(find ./ -maxdepth 4 -name '*bam' -type l ) | sort -k1,1n > homebrew_unique_read_counts.tsv
awk '$1>100000 {print $0}' homebrew_unique_read_counts.tsv > homebrew_cells.pf.txt
mkdir -p sc_bams
parallel -j 100 -a homebrew_cells.pf.txt split_bams

#also copy all symlink files from scalebio nextflow by following symlinks (so we don't need work dir maintained)
find ${homebrew_runDir} -maxdepth 5 -type l -exec bash -c 'cp -L -R "$(readlink -m "$0")" "$0".dereferenced' {} \; #copy files
find ${homebrew_runDir} -maxdepth 5 -name "*.dereferenced" -type f -exec bash -c 'mv $0 $(echo $0 | sed -e 's/".dereferenced"//g' -)' {} \; #move to old file names

rm -rf ${SCRATCH}/scalemet_work

mkdir -p ${homebrew_runDir}/copykit_cnv
singularity exec \
--bind /home/rmulqueen/projects/scalebio_dcis/ \
~/singularity/copykit.sif \
Rscript ~/projects/scalebio_dcis/tools/scalemet_dcis/src/copykit_cnvcalling.R \
--input_dir ${homebrew_runDir}/sc_bams \
--output_dir ${homebrew_runDir}/copykit_cnv/ \
--output_prefix homemade \
--task_cpus 150

```

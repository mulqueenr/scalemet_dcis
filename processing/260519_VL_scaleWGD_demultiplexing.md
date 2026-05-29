# Example processing of library pools 

- using custom python script to make demultiplexing barcodes
- using nf-core to demultiplex
- using scalebio pipeline for methylation analysis


All processed plates
| Plate| Chemistry |Index i5 |Index I7|
|---|---|---|---|
|Plate1|Scalebio Chemistry|Scale Bio Plate 1|Scale Bio Plate 1|
|Plate2|Scalebio Chemistry|Scale Bio Plate 2|Scale Bio Plate 2|
|Plate3|Homebrew|Homebrew 1|Homebrew A|
|Plate4|Homebrew|Homebrew 2|Homebrew B|
|Plate5|Homebrew|NOT RUN|NOT RUN|
|Plate6|Homebrew|NOT RUN|NOT RUN|
|Plate11|Homebrew|Homebrew 1|Homebrew B|
|Plate12|Homebrew|Homebrew 2|Homebrew A|
|Plate13|Homebrew|Homebrew 3|Homebrew B|
|Plate14|Homebrew|Homebrew 3|Homebrew C|
|Plate15|Homebrew|NOT RUN|NOT RUN|
|Plate16|Homebrew|NOT RUN|NOT RUN|

Plates 5,6 15 and 16 stored in case we need more methylomes.


Setting up environment variables on SBIO server for project. 
Using scalebio_dcis project folder for already installed tools.

```bash
#export environment variables for working/scratch directories
export SCRATCH="/data/rmulqueen/projects/subclonal_wgd/scratch/scalemet_work"
export TMPDIR="/data/rmulqueen/projects/subclonal_wgd/scratch"
export NXF_SINGULARITY_CACHEDIR="/data/rmulqueen/projects/subclonal_wgd/singularity"
export SINGULARITY_BINDPATH="/data/rmulqueen/projects/scalebio_dcis/tools/ScaleMethyl/bin" #using preinstalled ScaleMethyl pipeline here from scale git repo https://github.com/ScaleBio/ScaleMethyl

#set up directories and variables
projDir="/data/rmulqueen/projects/subclonal_wgd"
scalebio_nf="/data/rmulqueen/projects/scalebio_dcis/tools/ScaleMethyl" #using preinstalled ScaleMethyl pipeline here from scale git repo
params="/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/dcis_runParams.yml"
runDir="${projDir}/data/260519_initial_run"

bclDir="/data/instruments/nextseq2000/NavinLab/Ryan/260518_VH01788_170_AAJ5K5CM5" #preliminary run on plate 12 only to check QC
```

SAMPLE SPECIFIC TN5 SHEETS
Make a sheet per plate processed so it can be run lane-specific on future novaseq runs.

|Processing ID |Sample	ID|	Alternative ID|
|---|---|---|
|1|GEM2.1_LR_1|2A6|
|2|GEM2.1_LR_2|2A7|
|3|GEM2.8_PT_2|2B8|
|4|GEM2.8_PT_3|2B9|
|5|GEM2.12_PT_1|2D6|
|6|GEM2.12_PT_2|2D7|

Samples 1+2 were combined (due to low cell counts)
Will be referred to as GEM2.1_LR_1-2.

*All samples were FACS sorted into D (diploid) and A (aneuloid/WGD gates, and were tagmented with that a priori knowledge)*

Generating a csv of tagmentation wells following scalebio input requirements. Samples loaded by column.
Name should only contain [a-z],[A-Z],[0-9], dash (-)
[GEM2]-[.##]-[LR|PT]-[##]-[A|D]

```bash
mkdir -p ${runDir}
mkdir -p ${runDir}/samplesheets

samples="${runDir}/samplesheets/samples.csv"

echo """sample,barcodes,libName
GEM2-1-LR-12-A,1A07;1B07;1C07;1D07;1E07;1F07;1G07;1H07;1A12;1B12;1C12;1D12;1E12;1F12;1G12;1H12;2A01;2B01;2C01;2D01;2E01;2F01;2G01;2H01;3F07;3G07;3H07;3A09;3B09;3C09;3D09;3E09;3F09;3G09;3H09,ScaleMethyl
GEM2-1-LR-12-D,1A04;1B04;1C04;1D04;1E04;1F04;1G04;1H04;2A12;2B12;2C12;2D12;2E12;2F12;2G12;2H12,ScaleMethyl
GEM2-8-PT-2-A,1A10;1B10;1C10;1D10;1E10;1F10;1G10;1H10;2A04;2B04;2C04;2D04;2E04;2F04;2G04;2H04;2A06;2B06;2C06;2D06;2E06;2F06;2G06;2H06;3A07;3B07,ScaleMethyl
GEM2-8-PT-2-D,1A06;1B06;1C06;1D06;1E06;1F06;1G06;1H06;2A03;2B03;2C03;2D03;2E03;2F03;2G03;2H03;3A02;3B02;3C02;3D02;3E02;3F02;3G02;3H02;3A06;3B06;3C06;3D06;3E06;3F06;3G06;3H06;3C07;3D07;3E07,ScaleMethyl
GEM2-8-PT-3-A,1A05;1B05;1C05;1D05;1E05;1F05;1G05;1H05;2A08;2B08;2C08;2D08;2E08;2F08;2G08;2H08,ScaleMethyl
GEM2-8-PT-3-D,1A02;1B02;1C02;1D02;1E02;1F02;1G02;1H02;2A07;2B07;2C07;2D07;2E07;2F07;2G07;2H07;3A01;3B01;3C01;3D01;3E01;3F01;3G01;3H01;3A11;3B11;3C11;3D11;3E11;3F11;3G11;3H11,ScaleMethyl
GEM2-12-PT-1-A,1A03;1B03;1C03;1D03;1E03;1F03;1G03;1H03;1A09;1B09;1C09;1D09;1E09;1F09;1G09;1H09;2A05;2B05;2C05;2D05;2E05;2F05;2G05;2H05;3A10;3B10;3C10;3D10;3E10;3F10;3G10;3H10,ScaleMethyl
GEM2-12-PT-1-D,1A11;1B11;1C11;1D11;1E11;1F11;1G11;1H11;2A02;2B02;2C02;2D02;2E02;2F02;2G02;2H02;2A10;2B10;2C10;2D10;2E10;2F10;2G10;2H10;3A05;3B05;3C05;3D05;3E05;3F05;3G05;3H05,ScaleMethyl
GEM2-12-PT-2-A,1A08;1B08;1C08;1D08;1E08;1F08;1G08;1H08;3A03;3B03;3C03;3D03;3E03;3F03;3G03;3H03;3A04;3B04;3C04;3D04;3E04;3F04;3G04;3H04;3A12;3B12;3C12;3D12;3E12;3F12;3G12;3H12,ScaleMethyl
GEM2-12-PT-2-D,1A01;1B01;1C01;1D01;1E01;1F01;1G01;1H01;2A09;2B09;2C09;2D09;2E09;2F09;2G09;2H09;2A11;2B11;2C11;2D11;2E11;2F11;2G11;2H11;3A08;3B08;3C08;3D08;3E08;3F08;3G08;3H08,ScaleMethyl""" > ${samples}

```

Pull git repo with code and homebrew indexes

```bash
mkdir -p ${projDir}/tools/
cd ${projDir}/tools/
git clone https://github.com/mulqueenr/scalemet_dcis.git ${projDir}/tools/scalemet_dcis #clone dcis project repo (has homebrew indexes)

```

## For scalebio indexes, use ScaleMethyl bcl_convert_sheet.py

```bash
#example
#NOT RUN, since i tested seq on a homebrew plate
#python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
#${runDir}/samplesheets/tn5.csv \
#${projDir}/tools/ScaleMethyl/references/lib.json \
#${bclDir}/RunInfo.xml \
#--splitFastq > ${runDir}/samplesheets/plate1-2_samplesheet.csv
```

## For homebrew indexes, use custom script.
Can be pulled from here:

```bash
#sample sheets for homebrew plates, run one plate at a time and merge
#we ran plate 12 which is 
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.csv \
${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
${bclDir}/RunInfo.xml \
--splitFastq \
--i7Set A \
--i5Set 2 > ${runDir}/samplesheets/plate12_samplesheet.csv

```

Make samplesheet input for nf-core demultiplexing.

```bash
cd $runDir

#make nf-core input sheet
echo """id,samplesheet,lane,flowcell
plate12,${runDir}/samplesheets/plate12_samplesheet.csv,,${bclDir}
""" > pipeline_samplesheet1.csv 
```

# BCL to FASTQ
Using nf-core demultiplex

```bash
mkdir -p $SCRATCH
cd $runDir

#use nf-core to split out fastqs
nextflow run nf-core/demultiplex \
    --input pipeline_samplesheet1.csv \
    --outdir fastq \
    --trim_fastq false \
    --remove_adapter false \
    --skip_tools fastp,fastqc,kraken,multiqc,checkqc,falco,md5sum,samshee \
    -profile singularity \
    -w ${SCRATCH}

```

## Clean up FASTQ Files
Removing undetermined reads, so they aren't included in scalebio processing.

```bash
#clean up undetermined reads for each plate
parallel -j 100 rm {} ::: $(find . -type f -name "Undetermined*.fastq.gz")
```

# Run SCALEMETHYL Pipeline
Splitting plates by 
- homebrew and scale

Homebrew plates require a custom lib.json structure that points to the i7 and i5 custom PCR indexes.

For homebrew running, need to provide the homebrew indexes for splitting:
Library structure with homebrew indexes:
https://github.com/mulqueenr/scalemet_dcis/blob/main/ref/homebrew.lib.json
https://github.com/mulqueenr/scalemet_dcis/blob/main/ref/homebrew_i5.txt
https://github.com/mulqueenr/scalemet_dcis/blob/main/ref/homebrew_i7.txt
https://github.com/mulqueenr/scalemet_dcis/blob/main/ref/tgmt.txt


Example params:
https://github.com/mulqueenr/scalemet_dcis/blob/main/src/dcis_runParams.yml

```bash
mkdir -p ${SCRATCH}
mkdir -p ${runDir}/scalemethyl_pipeline_out
cd ${runDir}

plate="plate12"
mkdir -p ${runDir}/scalemethyl_pipeline_out/$plate

#note for lane-specified splitting fastqdir is one deeper (have to specify lane)
nice nextflow run ${scalebio_nf} \
--fastqDir ${runDir}/fastq/$plate/ \
--samples ${runDir}/samplesheets/samples.csv \
--outDir ${runDir}/scalemethyl_pipeline_out/$plate \
--libStructure ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
--maxMemory 500.GB \
--maxCpus 200 \
-profile singularity \
-params-file ${params} \
-w ${SCRATCH}/scalemet_milestone_work;

```

## Example run of scale standard processed plates 1+2

using the defaults lib structure

```bash
#NOT RUN
#plate="plate1-2"
#mkdir -p ${runDir}/scalemethyl_pipeline_out/$plate

#nice nextflow run ${scalebio_nf} \
#--fastqDir ${runDir}/fastq/$plate/ \
#--samples ${runDir}/samplesheets/samples.csv \
#--outDir ${runDir}/scalemethyl_pipeline_out/$plate \
#--maxMemory 500.GB \
#--maxCpus 200 \
#-profile singularity \
#-params-file ${params} \
#-w ${SCRATCH}/scalemet_milestone_work;

```

Correct all symlink with actual files so scratch can be cleared.

```bash

#also copy all symlink files from scalebio nextflow by following symlinks (so we don't need work dir maintained)
find ${runDir}/scalemethyl_pipeline_out/ -maxdepth 7 -type l -exec bash -c 'cp -L -R "$(readlink -m "$0")" "$0".dereferenced' {} \; #copy files
find ${runDir}/scalemethyl_pipeline_out/ -maxdepth 7 -name "*.dereferenced" -type f -exec bash -c 'mv $0 $(echo $0 | sed -e 's/".dereferenced"//g' -)' {} \; #move to old file names

```

#could potentially estimate library complexity via picard tools, but for now just passing back to Lesky.
```bash

#make function
#go to folder, count indexes per bam, take indexes over 10k reads, split out bam with that index into sc_bam folder, run picard tools to estimate complexity
mkdir -p $runDir/sc_bams

dir_list=$(find ${runDir}/scalemethyl_pipeline_out/plate12/alignments/dedup -maxdepth 7 -type d)

#goes through each tagmentation well directory and splits out cells from bams with >10k reads

for dir_in in $dir_list; do
echo $dir_in
bam_in=$(ls ${dir_in}/*.dedup.bam)
echo $bam_in
index_filt=$(samtools view $bam_in | awk '{split($1,a,":"); print a[8]}' | sort | uniq -c | sort -k1,1n | awk '$1>10000 {print $2}')
for idx in $index_filt; do
((samtools view -H $bam_in) && (samtools view $bam_in | awk -v idx=$idx '{split($1,a,":"); if(a[8]==idx){print $0}}')) | samtools view -O BAM - > $runDir/sc_bams/$idx.dedup.bam
done
done



#projecting complexity
cd $runDir/sc_bams

function proj_complexity() {
bam=$1
cellid="${bam::-10}"
#project complexity bam, just mark duplicates
samtools sort -m 10G -n $bam | \
samtools fixmate -p -m - - | \
samtools sort -m 10G | \
samtools markdup --mode t -s - ${cellid}.mkdup.bam

#${cellid}.mkdup.bam

java -jar /home/rmulqueen/tools/picard.jar \
EstimateLibraryComplexity \
MAX_OPTICAL_DUPLICATE_SET_SIZE=-1 \
TMP_DIR="." \
READ_NAME_REGEX="^[^:]+:[^:]+:[^:]+:[^:]+:([0-9]+):([0-9]+):([0-9]+):.+$" \
I=${cellid}.mkdup.bam \
VALIDATION_STRINGENCY=LENIENT \
O=${cellid}.complex_metrics.txt

#format a bit
grep "^Unknown" ${cellid}.complex_metrics.txt | \
awk -v cellid=${cellid} 'OFS="," {print cellid,$3,$9,$10}' > ${cellid}.projected_metrics.picard.txt
}

export -f proj_complexity
parallel -j 300 proj_complexity ::: $(ls *.bam)

cat *.projected_metrics.picard.txt > plate12.projected_metrics.txt

```

```R
library(dplyr)
dat<-read.table("plate12.projected_metrics.txt",sep=",")
colnames(dat)<-c("idx","read_count","percent_duplication","estimated_size")
dat %>% filter(read_count>50000)
```
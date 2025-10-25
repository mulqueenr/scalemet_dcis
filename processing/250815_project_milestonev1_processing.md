# Milestone v1 Preprocessing

Combining all sequencing runs into the same demultiplexing script, to merge cells at the fastq level prior to running scalemethyl pipeline per plate. See individual flowcell processing for one note link on library generation.

All processed plates

| Batch | Plate | Processing | i5 Index | i7 Index | NovaSeq Lane |
|---|---|---|---|---|---|
| 1 | 1 | scale kit | scale plate 1 | scale plate 1 | 1 |
| 1 | 6 | scale kit | scale plate 2 | scale plate 2 | 1 |
| 1 | 12 | homebrew | D (4) | B | 1 |
| 3 | 1 | homebrew | A (1) | E | 2 |
| 2 | 3 | scale kit | scale plate 1 | scale plate 1 | 2 |
| 2 | 4 | scale kit | scale plate 2 | scale plate 2 | 2 |
| 3 | 2 | scale kit | scale plate 1 | scale plate 1 | 3 |
| 3 | 5 | scale kit | scale plate 2 | scale plate 2 | 3 |
| 3 | 11 | homebrew | E (5) | A | 3 |
| prelim1-2 | 1 | scale kit | scale plate 1 | scale plate 1 | 4 |
| prelim1-2 | 2 | scale kit | scale plate 2 | scale plate 2 | 4 |
| 1 | 3 | homebrew | C (3) | A | 4 |
| prelim3 | 1 | scale kit | scale plate 1 | scale plate 1 | 5 |
| prelim3 | 2 | scale kit | scale plate 2 | scale plate 2 | 5 |
| 3 | 13 | homebrew | E (5) | B | 5 |
| 2 | 2 | homebrew | B (2) | D | 6 |
| 3 | 15 | homebrew | D (4) | C | 6 |
| 3 | 16 | homebrew | D (4) | D | 6 |
| 3 | 17 | homebrew | E (5) | D | 6 |
| 2 | 1 | homebrew | A (1) | D | 7 |
| 3 | 6 | homebrew | B (2) | E | 7 |
| 1 | 8 | homebrew | C (3) | C | 7 |
| 1 | 11 | homebrew | D (4) | A | 7 |
| 1 | 2 | homebrew | A (1) | A | 8 |
| 1 | 4 | homebrew | C (3) | B | 8 |
| 1 | 7 | homebrew | B (2) | B | 8 |
| 1 | 9 | homebrew | C (3) | D | 8 |

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
runDir="${projDir}/data/250815_milestone_v1"

bclDir_240202="/data/rmulqueen/projects/scalebio_dcis/seq/PM2517/" #prelim1
bclDir_240523="/data/rmulqueen/projects/scalebio_dcis/seq/240523_VH00219_594_AAFLYGNM5" #prelim2
bclDir_241004="/data/rmulqueen/projects/scalebio_dcis/seq/241004_A01819_0637_BHY5MJDMXY" #prelim3
bclDir_250325="/data/rmulqueen/projects/scalebio_dcis/seq/250325_VH01788_97_2227YC2NX" #batch1 plates 1-2-6-7
bclDir_250424="/data/rmulqueen/projects/scalebio_dcis/seq/PM2563" #batch1 plates 3-8
bclDir_250506="/data/rmulqueen/projects/scalebio_dcis/seq/250506_VH01788_104_222CF7LNX" #batch 2
bclDir_250808="/data/rmulqueen/projects/scalebio_dcis/seq/20250808_LH00503_0129_A22WHTCLT4" #novaseq batch 1-2-3-prelim1-2-3

#modified /data/rmulqueen/projects/scalebio_dcis/tools/ScaleMethyl/modules/input_reads.nf to decrease alignment thread usage
```


SAMPLE SPECIFIC TN5 SHEETS

```bash
mkdir -p ${runDir}
mkdir -p ${runDir}/samplesheets

################PRELIM 1-2################
samples="${runDir}/samplesheets/samples.prelim1-2.csv"
echo """sample,barcodes,libName
MCF10A,1A01-1B12,ScaleMethyl
MCF7,1C01-1E12,ScaleMethyl
MDA-MB-231,1F01-1H12,ScaleMethyl
HBCA-16R,2A01-2D12,ScaleMethyl
HBCA-83L,2E01-3H12,ScaleMethyl
BCMDCIS41T,3A01-3A08;3B01-3B08;3C01-3C08;3D01-3D08;3E01-3E08;3F01-3F08;3G01-3G08;3H01-3H08,ScaleMethyl
BCMDCIS66T,3A09-3A11;3B09-3B11;3C09-3C11;3D09-3D11;3E09-3E11;3F09-3F11;3G09-3G11;3H09-3H11,ScaleMethyl""" > ${samples}


################PRELIM 3################
samples="${runDir}/samplesheets/samples.prelim3.csv"

echo """sample,barcodes,libName
DCIS-92T,1A01-1D12,ScaleMethyl
DCIS-66T,1E01-1H12,ScaleMethyl
DCIS-79T,2A01-2D12,ScaleMethyl
IDC-79T,2E01-2H12,ScaleMethyl
HBCA-19T,3A01-3D12,ScaleMethyl
HBCA-17T,3E01-3H12,ScaleMethyl""" > ${samples}

################BATCH 1################
samples="${runDir}/samplesheets/samples.batch1.csv"

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

################BATCH 2################
samples="${runDir}/samplesheets/samples.batch2.csv"

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

################BATCH 3################
samples="${runDir}/samplesheets/samples.batch3.csv"

echo """sample,barcodes,libName
BCMHBCA22R-4h,1A01-1A06;2A07-2A12;3A01-3A05,ScaleMethyl
ECIS57T,1B01-1B06;2B07-2B12;3B01-3B06,ScaleMethyl
BCMDCIS82T-24hTis,1C01-1C06;2C07-2C12;3C01-3C06,ScaleMethyl
ECIS26T,1D01-1D06;2D07-2D12;3D01-3D06;3A06,ScaleMethyl
BCMDCIS97T,1E01-1E06;2E07-2E12;3E01-3E06;1C07-1C12;2C01-2C06;3C07-3C12,ScaleMethyl
BCMDCIS102T-4h,1F01-1F12;2F01-2F12;3F01-3F11,ScaleMethyl
BCMHBCA29L-2h,1G01-1G06;1D07-1D12;2D01-2D06;2G07-2G12;3G01-3G06;3D07-3D12,ScaleMethyl
BCMHBCA09R-3h,1H01-1H06;2E01-2E06;3H01-3H06;1E07-1E12;2H07-2H12;3E07-3E12,ScaleMethyl
BCMDCIS35T-3h,1A07-1A12;1G07-1G12;2A01-2A06;2G01-2G06;3A07-3A12;3G07-3G12,ScaleMethyl
BCMDCIS124T,1B07-1B12;1H07-1H12;2B01-2B06;2H01-2H06;3B07-3B12;3H07-3H12,ScaleMethyl""" > ${samples}

```

Pull git repo with code and homebrew indexes

```bash
cd ${projDir}/tools/
git clone https://github.com/mulqueenr/scalemet_dcis.git /data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis #this repo
```

## Nextseq 240202
```bash
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.prelim1-2.csv \
${projDir}/tools/ScaleMethyl/references/lib.json \
${bclDir_240202}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/prelim1-2_scalebio_plate1_samplesheet.csv
```

## Nextseq 240523
```bash
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.prelim1-2.csv \
${projDir}/tools/ScaleMethyl/references/lib.json \
${bclDir_240523}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/prelim1-2_scalebio_plate2_samplesheet.csv
```

## Novaseq 241004
```bash
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.prelim3.csv \
${projDir}/tools/ScaleMethyl/references/lib.json \
${bclDir_241004}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/prelim3_scalebio_plate1_samplesheet.csv
```

## Nextseq 250325
```bash

python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.batch1.csv \
${projDir}/tools/ScaleMethyl/references/lib.json \
${bclDir_250325}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/batch1_scalebio_plate1-6_samplesheet.csv

#sample sheets for homebrew plates, run one plate at a time and merge
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv \
${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
${bclDir_250325}/RunInfo.xml \
--splitFastq \
--i7Set A \
--i5Set 1 > ${runDir}/samplesheets/batch1_homebrew_plate2_samplesheet.csv

python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv \
${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
${bclDir_250325}/RunInfo.xml \
--splitFastq \
--i7Set B \
--i5Set 2 > ${runDir}/samplesheets/batch1_homebrew_plate7_samplesheet.csv
```

## Novaseq 250424
Lane 5 only

```bash
#sample sheets for homebrew plates, run one plate at a time and merge
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv \
${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
${bclDir_250424}/RunInfo.xml \
--splitFastq \
--i7Set A \
--i5Set 3 > ${runDir}/samplesheets/batch1_homebrew_plate3_samplesheet.csv

python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv \
${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
${bclDir_250424}/RunInfo.xml \
--splitFastq \
--i7Set C \
--i5Set 3 > ${runDir}/samplesheets/batch1_homebrew_plate8_samplesheet.csv

#replace override cycles (just using sed) and also only run on lane 5 of run (from pipeline samplesheet)
sed -i 's/OverrideCycles,Y148;I34;I10;Y146/OverrideCycles,Y148;I10N24;I10;Y146/g' \
${runDir}/samplesheets/batch1_homebrew_plate3_samplesheet.csv

#replace override cycles (just using sed) and also only run on lane 5 of run (from pipeline samplesheet)
sed -i 's/OverrideCycles,Y148;I34;I10;Y146/OverrideCycles,Y148;I10N24;I10;Y146/g' \
${runDir}/samplesheets/batch1_homebrew_plate8_samplesheet.csv
```

## NextSeq 250506
```bash
#generate sample sheet using a modified version of bcl_convert_sheet.py to allow for pcr plate specifications.
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.batch2.csv \
${projDir}/tools/ScaleMethyl/references/lib.json \
${bclDir_250506}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/batch2_scalebio_plate3-4_samplesheet.csv

#sample sheets for homebrew plates, run one plate at a time and merge
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch2.csv \
${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
${bclDir_250506}/RunInfo.xml \
--splitFastq \
--i7Set D \
--i5Set 1 > ${runDir}/samplesheets/batch2_homebrew_plate1_samplesheet.csv

python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch2.csv \
${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
${bclDir_250506}/RunInfo.xml \
--splitFastq \
--i7Set D \
--i5Set 2 > ${runDir}/samplesheets/batch2_homebrew_plate2_samplesheet.csv
```

## NovaSeq 250808
Lane 1
```bash
#generate sample sheet using a modified version of bcl_convert_sheet.py to allow for pcr plate specifications.
#scalebio indexes batch 1
cd $runDir/samplesheets
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/ScaleMethyl/references/lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/novaseq_lane1_batch1_scalebio_plate1-6.csv

#sample sheets for homebrew plates, run one plate at a time and merge
#batch 1 plate 12
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq --i7Set B --i5Set 4 > ${runDir}/samplesheets/novaseq_lane1_batch1_homebrew_plate12.csv

```
Lane 2
```bash
#scalebio indexes batch 2
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.batch2.csv ${projDir}/tools/ScaleMethyl/references/lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/novaseq_lane2_batch2_scalebio_plate3-4.csv

#batch 3 plate 1
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch3.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq --i7Set E --i5Set 1 > ${runDir}/samplesheets/novaseq_lane2_batch3_homebrew_plate1.csv
```

Lane 3
```bash
#scalebio indexes batch 3
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.batch3.csv ${projDir}/tools/ScaleMethyl/references/lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/novaseq_lane3_batch3_scalebio_plate2-5.csv

#batch 3 plate 11
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch3.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq --i7Set A --i5Set 5 > ${runDir}/samplesheets/novaseq_lane3_batch3_homebrew_plate11.csv
```

Lane 4
```bash
#scalebio indexes batch prelim1 and prelim2
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.prelim1-2.csv ${projDir}/tools/ScaleMethyl/references/lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/novaseq_lane4_prelim1-2_scalebio_plate1-2.csv

#batch 1 plate 3
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch3.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq --i7Set A --i5Set 3 > ${runDir}/samplesheets/novaseq_lane4_batch1_homebrew_plate3.csv
```

Lane 5
```bash
#scalebio indexes batch prelim3
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.prelim3.csv ${projDir}/tools/ScaleMethyl/references/lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/novaseq_lane5_prelim3_scalebio_plate1-2.csv

#batch 3 plate 13
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch3.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq --i7Set B --i5Set 5 > ${runDir}/samplesheets/novaseq_lane5_batch3_homebrew_plate13.csv
```

Lane 6
```bash
#batch 2 plate 2
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch2.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq --i7Set D --i5Set 2 > ${runDir}/samplesheets/novaseq_lane6_batch2_homebrew_plate2.csv

#batch 3 plate 15
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq --i7Set C --i5Set 4 > ${runDir}/samplesheets/novaseq_lane6_batch3_homebrew_plate15.csv

#batch 3 plate 16
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq --i7Set D --i5Set 4 > ${runDir}/samplesheets/novaseq_lane6_batch3_homebrew_plate16.csv

#batch 3 plate 17
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq --i7Set D --i5Set 5 > ${runDir}/samplesheets/novaseq_lane6_batch3_homebrew_plate17.csv
```

Lane 7
```bash
#batch 2 plate 1
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq --i7Set D --i5Set 1 > ${runDir}/samplesheets/novaseq_lane7_batch2_homebrew_plate1.csv

#batch 3 plate 6
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq --i7Set E --i5Set 2 > ${runDir}/samplesheets/novaseq_lane7_batch3_homebrew_plate6.csv

#batch 1 plate 8
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq --i7Set C --i5Set 3 > ${runDir}/samplesheets/novaseq_lane7_batch1_homebrew_plate8.csv

#batch 1 plate 11
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq --i7Set A --i5Set 4 > ${runDir}/samplesheets/novaseq_lane7_batch1_homebrew_plate11.csv
```

Lane 8
```bash
#batch 1 plate 2
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch2.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq --i7Set A --i5Set 1 > ${runDir}/samplesheets/novaseq_lane8_batch1_homebrew_plate2.csv

#batch 1 plate 4
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch3.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq --i7Set B --i5Set 3 > ${runDir}/samplesheets/novaseq_lane8_batch1_homebrew_plate4.csv

#batch 1 plate 7
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch3.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq --i7Set B --i5Set 2 > ${runDir}/samplesheets/novaseq_lane8_batch1_homebrew_plate7.csv

#batch 1 plate 9
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir_250808}/RunInfo.xml \
--splitFastq --i7Set D --i5Set 3 > ${runDir}/samplesheets/novaseq_lane8_batch1_homebrew_plate9.csv
```

Novaseq 250808 is the only instance samples are sequenced again, so samplesheets for novaseq have a "novaseq" prefix (to account for proper index complementing etc.)

```bash
cd $runDir

#make nf-core input sheet
echo """id,samplesheet,lane,flowcell
prelim1-2_scalebio_plate1,${runDir}/samplesheets/prelim1-2_scalebio_plate1_samplesheet.csv,,${bclDir_240202}
prelim1-2_scalebio_plate2,${runDir}/samplesheets/prelim1-2_scalebio_plate2_samplesheet.csv,,${bclDir_240523}
batch1_scalebio_plate1-6,${runDir}/samplesheets/batch1_scalebio_plate1-6_samplesheet.csv,,${bclDir_250325}
batch1_homebrew_plate2,${runDir}/samplesheets/batch1_homebrew_plate2_samplesheet.csv,,${bclDir_250325}
batch1_homebrew_plate7,${runDir}/samplesheets/batch1_homebrew_plate7_samplesheet.csv,,${bclDir_250325}
batch1_homebrew_plate3,${runDir}/samplesheets/batch1_homebrew_plate3_samplesheet.csv,5,${bclDir_250424}
batch1_homebrew_plate8,${runDir}/samplesheets/batch1_homebrew_plate8_samplesheet.csv,5,${bclDir_250424}
batch2_scalebio_plate3-4,${runDir}/samplesheets/batch2_scalebio_plate3-4_samplesheet.csv,,${bclDir_250506}
batch2_homebrew_plate1,${runDir}/samplesheets/batch2_homebrew_plate1_samplesheet.csv,,${bclDir_250506}
batch2_homebrew_plate2,${runDir}/samplesheets/batch2_homebrew_plate2_samplesheet.csv,,${bclDir_250506}
""" > pipeline_samplesheet1.csv 

#make nf-core input sheet
echo """id,samplesheet,lane,flowcell
batch1_scalebio_plate1-6,${runDir}/samplesheets/novaseq_lane1_batch1_scalebio_plate1-6.csv,1,${bclDir_250808}
batch1_homebrew_plate12,${runDir}/samplesheets/novaseq_lane1_batch1_homebrew_plate12.csv,1,${bclDir_250808}
batch2_scalebio_plate3-4,${runDir}/samplesheets/novaseq_lane2_batch2_scalebio_plate3-4.csv,2,${bclDir_250808}
batch3_homebrew_plate1,${runDir}/samplesheets/novaseq_lane2_batch3_homebrew_plate1.csv,2,${bclDir_250808}
batch3_scalebio_plate2-5,${runDir}/samplesheets/novaseq_lane3_batch3_scalebio_plate2-5.csv,3,${bclDir_250808}
batch3_homebrew_plate11,${runDir}/samplesheets/novaseq_lane3_batch3_homebrew_plate11.csv,3,${bclDir_250808}
prelim1-2_scalebio_plate1-1,${runDir}/samplesheets/novaseq_lane4_prelim1-2_scalebio_plate1-1.csv,4,${bclDir_250808}
batch1_homebrew_plate3,${runDir}/samplesheets/novaseq_lane4_batch1_homebrew_plate3.csv,4,${bclDir_250808}
prelim3_scalebio_plate1-2,${runDir}/samplesheets/novaseq_lane5_prelim3_scalebio_plate1-2.csv,5,${bclDir_250808}
batch3_homebrew_plate13,${runDir}/samplesheets/novaseq_lane5_batch3_homebrew_plate13.csv,5,${bclDir_250808}
batch2_homebrew_plate2,${runDir}/samplesheets/novaseq_lane6_batch2_homebrew_plate2.csv,6,${bclDir_250808}
batch3_homebrew_plate15,${runDir}/samplesheets/novaseq_lane6_batch3_homebrew_plate15.csv,6,${bclDir_250808}
batch3_homebrew_plate16,${runDir}/samplesheets/novaseq_lane6_batch3_homebrew_plate16.csv,6,${bclDir_250808}
batch3_homebrew_plate17,${runDir}/samplesheets/novaseq_lane6_batch3_homebrew_plate17.csv,6,${bclDir_250808}
""" > pipeline_samplesheet2.csv

echo """id,samplesheet,lane,flowcell
batch2_homebrew_plate1,${runDir}/samplesheets/novaseq_lane7_batch2_homebrew_plate1.csv,7,${bclDir_250808}
batch3_homebrew_plate6,${runDir}/samplesheets/novaseq_lane7_batch3_homebrew_plate6.csv,7,${bclDir_250808}
batch1_homebrew_plate8,${runDir}/samplesheets/novaseq_lane7_batch1_homebrew_plate8.csv,7,${bclDir_250808}
batch1_homebrew_plate11,${runDir}/samplesheets/novaseq_lane7_batch1_homebrew_plate11.csv,7,${bclDir_250808}
batch1_homebrew_plate2,${runDir}/samplesheets/novaseq_lane8_batch1_homebrew_plate2.csv,8,${bclDir_250808}
batch1_homebrew_plate4,${runDir}/samplesheets/novaseq_lane8_batch1_homebrew_plate4.csv,8,${bclDir_250808}
batch1_homebrew_plate7,${runDir}/samplesheets/novaseq_lane8_batch1_homebrew_plate7.csv,8,${bclDir_250808}
batch1_homebrew_plate9,${runDir}/samplesheets/novaseq_lane8_batch1_homebrew_plate9.csv,8,${bclDir_250808}
""" > pipeline_samplesheet3.csv

echo """id,samplesheet,lane,flowcell
batch1_homebrew_plate9,${runDir}/samplesheets/novaseq_lane8_batch1_homebrew_plate9.csv,8,${bclDir_250808}
batch1_homebrew_plate11,${runDir}/samplesheets/novaseq_lane7_batch1_homebrew_plate11.csv,7,${bclDir_250808}
""" > pipeline_rerun.csv
```

BCL to FASTQ

Running in 3 batches to limit strain of I/O on server.

```bash
mkdir -p $SCRATCH/scalemet_work
cd $runDir

#use nf-core to split out fastqs
nextflow run nf-core/demultiplex \
    --input pipeline_rerun.csv \
    --outdir fastq \
    --trim_fastq false \
    --remove_adapter false \
    --skip_tools fastp,fastqc,kraken,multiqc,checkqc,falco,md5sum,samshee \
    -profile singularity \
    -w ${SCRATCH}/scalemet_work

#use nf-core to split out fastqs
nextflow run nf-core/demultiplex \
    --input pipeline_samplesheet1.csv \
    --outdir fastq \
    --trim_fastq false \
    --remove_adapter false \
    --skip_tools fastp,fastqc,kraken,multiqc,checkqc,falco,md5sum,samshee \
    -profile singularity \
    -w ${SCRATCH}/scalemet_work

#use nf-core to split out fastqs
nextflow run nf-core/demultiplex \
    --input pipeline_samplesheet2.csv \
    --outdir fastq \
    --trim_fastq false \
    --remove_adapter false \
    --skip_tools fastp,fastqc,kraken,multiqc,checkqc,falco,md5sum,samshee \
    -profile singularity \
    -w ${SCRATCH}/scalemet_work

#use nf-core to split out fastqs
nextflow run nf-core/demultiplex \
    --input pipeline_samplesheet3.csv \
    --outdir fastq \
    --trim_fastq false \
    --remove_adapter false \
    --skip_tools fastp,fastqc,kraken,multiqc,checkqc,falco,md5sum,samshee \
    -profile singularity \
    -w ${SCRATCH}/scalemet_work

```

# Clean up FASTQ Files

```bash
#clean up undetermined reads for each plate
parallel -j 100 rm {} ::: $(find . -type f -name "Undetermined*.fastq.gz")

#merge fastq files from previous runs to the novaseq resequenced samples
#all samplesheet1 samples were resequenced on novaseq

#make a new folder per plate of merged_fq
#set Lane to wild card for matching wells
#output concat files with new prefix (ScaleMethylMerged)
#parallelize so it doesnt take forever

concatenate_fastq() {
    well_new=$(echo $2 | sed 's/_L00\*_/_L001_/' |  sed 's/ScaleMethyl/ScaleMethylMerged/')
    find ./fastq/$1/ -name "$2" -type f -exec cat {} + > ./fastq/$1/merged_fq/$well_new; 
}
export -f concatenate_fastq

for plate in prelim1-2_scalebio_plate1-1 prelim3_scalebio_plate1 batch1_scalebio_plate1-6 batch1_homebrew_plate2 batch1_homebrew_plate7 batch1_homebrew_plate3 batch1_homebrew_plate8 batch2_scalebio_plate3-4 batch2_homebrew_plate1 batch2_homebrew_plate2;
    do echo "Combing FASTQ for $plate..."
    mkdir -p ./fastq/$plate/merged_fq
    files_in=$(find ./fastq/$plate/ -name ScaleMethyl_*fastq.gz -exec basename {}  \; | sed 's/_L00[1-8]_/_L00*_/' | sort | uniq)
    parallel -j 200 concatenate_fastq $plate {} ::: $files_in
done

```

# Run SCALEMETHYL Pipeline
Splitting plates by 
1. homebrew vs scale
2. one sequencing run vs multiple for processing

```bash
mkdir -p ${SCRATCH}/scalemet_milestone_work
mkdir -p ${runDir}/scalemethyl_pipeline_out
cd ${runDir}
```

## reseq plate, homebrew
DONE

```bash
for plate in batch1_homebrew_plate2 batch1_homebrew_plate7 batch1_homebrew_plate3 batch1_homebrew_plate8 batch2_homebrew_plate1 batch2_homebrew_plate2;
do 
mkdir -p ${runDir}/scalemethyl_pipeline_out/$plate
batch=$(echo $plate | cut -d'_' -f1)
sed 's/ScaleMethyl/ScaleMethylMerged/g' ${runDir}/samplesheets/samples.$batch.csv > ${runDir}/scalemethyl_pipeline_out/$plate/samples.$plate.csv

nice nextflow run ${scalebio_nf} \
--fastqDir ${runDir}/fastq/$plate/merged_fq \
--samples ${runDir}/scalemethyl_pipeline_out/$plate/samples.$plate.csv \
--outDir ${runDir}/scalemethyl_pipeline_out/$plate \
--libStructure ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
--maxMemory 100.GB \
--maxCpus 20 \
-profile singularity \
-params-file ${params} \
-w ${SCRATCH}/scalemet_milestone_work -resume -qs 50;
done
```

## new plates, homebrew
 
```bash
for plate in batch1_homebrew_plate11 batch1_homebrew_plate4 batch3_homebrew_plate1 batch3_homebrew_plate11 batch3_homebrew_plate13 batch3_homebrew_plate15 batch3_homebrew_plate16 batch3_homebrew_plate17 batch3_homebrew_plate6 batch1_homebrew_plate9 batch1_homebrew_plate12;
do 
mkdir -p ${runDir}/scalemethyl_pipeline_out/$plate
batch=$(echo $plate | cut -d'_' -f1)
lane=$(ls -d ${runDir}/fastq/$plate/L* | awk -F/ '{print $NF}')
cp ${runDir}/samplesheets/samples.$batch.csv ${runDir}/scalemethyl_pipeline_out/$plate/samples.$plate.csv

nice nextflow run ${scalebio_nf} \
--fastqDir ${runDir}/fastq/$plate/$lane \
--samples ${runDir}/scalemethyl_pipeline_out/$plate/samples.$plate.csv \
--outDir ${runDir}/scalemethyl_pipeline_out/$plate \
--libStructure ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
--maxMemory 100.GB \
--maxCpus 1 \
-profile singularity \
-params-file ${params} \
-w ${SCRATCH}/scalemet_milestone_work -resume;
done
```

## new plate, scale

```bash
plate="batch3_scalebio_plate2-5"
batch=$(echo $plate | cut -d'_' -f1)
lane="L003"
mkdir -p ${runDir}/scalemethyl_pipeline_out/$plate
cp ${runDir}/samplesheets/samples.$batch.csv ${runDir}/scalemethyl_pipeline_out/$plate/samples.$plate.csv

nice nextflow run ${scalebio_nf} \
--fastqDir ${runDir}/fastq/$plate/$lane \
--samples ${runDir}/scalemethyl_pipeline_out/$plate/samples.$plate.csv \
--outDir ${runDir}/scalemethyl_pipeline_out/$plate \
--maxMemory 100.GB \
--maxCpus 20 \
-profile singularity \
-params-file ${params} \
-w ${SCRATCH}/scalemet_milestone_work -resume;
```

## reseq plate, scale
```bash
for plate in batch2_scalebio_plate3-4 prelim3_scalebio_plate1-2; 
do batch=$(echo $plate | cut -d'_' -f1)
mkdir -p ${runDir}/scalemethyl_pipeline_out/$plate
sed 's/ScaleMethyl/ScaleMethylMerged/g' ${runDir}/samplesheets/samples.$batch.csv > ${runDir}/scalemethyl_pipeline_out/$plate/samples.$plate.csv

nice nextflow run ${scalebio_nf} \
--fastqDir ${runDir}/fastq/$plate/merged_fq \
--samples ${runDir}/scalemethyl_pipeline_out/$plate/samples.$plate.csv \
--outDir ${runDir}/scalemethyl_pipeline_out/$plate \
-profile singularity \
-params-file ${params} \
-w ${SCRATCH}/scalemet_milestone_work -resume;
done
```

# fix fastq for 
#prelim1-2_scalebio_plate1-2
#prelim3_scalebio_plate1-2
#batch1_scalebio_plate1-6
#batch1_homebrew_plate9
#batch1_homebrew_plate11

```bash

# DONE
plate="prelim3_scalebio_plate1-2"; batch=$(echo $plate | cut -d'_' -f1); lane="L005"
mkdir -p ${runDir}/scalemethyl_pipeline_out/$plate
cp ${runDir}/samplesheets/samples.$batch.csv ${runDir}/scalemethyl_pipeline_out/$plate/samples.$plate.csv
nice nextflow run ${scalebio_nf} --fastqDir ${runDir}/fastq/$plate/$lane --samples ${runDir}/scalemethyl_pipeline_out/$plate/samples.$plate.csv \
--outDir ${runDir}/scalemethyl_pipeline_out/$plate -profile singularity -params-file ${params} -w ${SCRATCH}/scalemet_milestone_work -resume

#DONE
plate="batch1_scalebio_plate1-6"; batch=$(echo $plate | cut -d'_' -f1); lane="L001"
mkdir -p ${runDir}/scalemethyl_pipeline_out/$plate
cp ${runDir}/samplesheets/samples.$batch.csv ${runDir}/scalemethyl_pipeline_out/$plate/samples.$plate.csv
nice nextflow run ${scalebio_nf} --fastqDir ${runDir}/fastq/$plate/$lane --samples ${runDir}/scalemethyl_pipeline_out/$plate/samples.$plate.csv \
--outDir ${runDir}/scalemethyl_pipeline_out/$plate -profile singularity -params-file ${params} -w ${SCRATCH}/scalemet_milestone_work -resume

#DONE
plate="prelim1-2_scalebio_plate1-1"; batch=$(echo $plate | cut -d'_' -f1); lane="L004"
mkdir -p ${runDir}/scalemethyl_pipeline_out/$plate
cp ${runDir}/samplesheets/samples.$batch.csv ${runDir}/scalemethyl_pipeline_out/$plate/samples.$plate.csv
nice nextflow run ${scalebio_nf} --fastqDir ${runDir}/fastq/$plate/$lane --samples ${runDir}/scalemethyl_pipeline_out/$plate/samples.$plate.csv \
--outDir ${runDir}/scalemethyl_pipeline_out/$plate -profile singularity -params-file ${params} -w ${SCRATCH}/scalemet_milestone_work -resume


#expecting 25G of R1 per lane. Getting 3-9G instead
#problem with splitting, not pipeline

for plate in batch1_homebrew_plate9 batch1_homebrew_plate11;
do 
mkdir -p ${runDir}/scalemethyl_pipeline_out/$plate
batch=$(echo $plate | cut -d'_' -f1)
lane=$(ls -d ${runDir}/fastq/$plate/L* | awk -F/ '{print $NF}')
cp ${runDir}/samplesheets/samples.$batch.csv ${runDir}/scalemethyl_pipeline_out/$plate/samples.$plate.csv

nice nextflow run ${scalebio_nf} \
--fastqDir ${runDir}/fastq/$plate/$lane \
--samples ${runDir}/scalemethyl_pipeline_out/$plate/samples.$plate.csv \
--outDir ${runDir}/scalemethyl_pipeline_out/$plate \
--libStructure ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json \
--maxMemory 100.GB \
--maxCpus 1 \
-profile singularity \
-params-file ${params} \
-w ${SCRATCH}/scalemet_milestone_work;
done



```

Correct all symlink with actual links and split out merged bam files into single-cell bam files.

```bash

#also copy all symlink files from scalebio nextflow by following symlinks (so we don't need work dir maintained)
find ${runDir}/scalemethyl_pipeline_out/ -maxdepth 7 -type l -exec bash -c 'cp -L -R "$(readlink -m "$0")" "$0".dereferenced' {} \; #copy files
find ${runDir}/scalemethyl_pipeline_out/ -maxdepth 7 -name "*.dereferenced" -type f -exec bash -c 'mv $0 $(echo $0 | sed -e 's/".dereferenced"//g' -)' {} \; #move to old file names


```

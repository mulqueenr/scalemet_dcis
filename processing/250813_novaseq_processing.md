[See notebook for wet lab processing.](https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/Doc.aspx?sourcedoc={32B0B85E-D559-4659-AF43-88CF3D29307C}&wd=target%28scalebio_dcis_processing.one%7C52B52D1F-07A0-48D1-A89A-CBFF2E600AC2%2FSequencing%20Preparation%7C13D73B59-E7F1-5E40-849F-97E71808BF91%2F%29&wdpartid={CAA2D319-2B15-8942-852A-064E78A416C3}{1}&wdsectionfileid={52B52D1F-07A0-48D1-A89A-CBFF2E600AC2})


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
runDir="${projDir}/data/250813_RM_scalebio_novaseq"
bclDir="/data/rmulqueen/projects/scalebio_dcis/seq/20250808_LH00503_0129_A22WHTCLT4"

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

MAKE SAMPLE SHEETS FOR PCR

Per lane (8 total) we need to run sample sheet generation for both:
* scalebio
* homebrew indexes

| NovaSeq Lane | Batch	| Plate | Processing | i5 Index | i7 Index |
|---|---|---|---|---|---|
| 1 | 1 | 1 | scale kit | Scale bio Plate 1 | Scale bio Plate 1 |
| 1 | 1 | 2 | homebrew | A | A |
| 1 | 1 | 6 | scale kit | Scale bio Plate 2 | Scale bio Plate 2 |
| 2 | 2 | 3 | scale kit | Scalebio Plate 1 | Scalebio Plate 1 |
| 2 | 2 | 4 | scale kit | Scalebio Plate 2 | Scalebio Plate 2 |
| 2 | 3 | 1 | homebrew | A | E |
| 3 | 3 | 2 | scale kit | scale plate 1 | scale plate 1 |
| 3 | 3 | 5 | scale kit | scale plate 2 | scale plate 2 |
| 3 | 3 | 6 | homebrew | B | E |
| 4 | 3 | 16 | homebrew | D | D |
| 4 | prelim1 | 1 | scale kit | scale plate 1 | scale plate 1 |
| 4 | prelim2 | 1 | scale kit | scale plate 2 | scale plate 2 |
| 5 | 3 | 17 | homebrew | E | D |
| 5 | prelim3 | 1 | scale kit | scale plate 1 | scale plate 1 |
| 5 | prelim3 | 2 | scale kit | scale plate 2 | scale plate 2 |
| 6 | 1 | 3 | homebrew | C | A |
| 6 | 1 | 4 | homebrew | C | B |
| 6 | 1 | 7 | homebrew | B | B |
| 6 | 1 | 8 | homebrew | C | C |
| 7 | 1 | 9 | homebrew | C | D |
| 7 | 1 | 11 | homebrew | D | A |
| 7 | 1 | 12 | homebrew | D | B |
| 7 | 2 | 1 | homebrew | A | D |
| 8 | 2 | 2 | homebrew | B | D |
| 8 | 3 | 11 | homebrew | E | A |
| 8 | 3 | 13 | homebrew | E | B |
| 8 | 3 | 15 | homebrew | D | C |


########Waiting for data to transfer########

Pull git repo with code and homebrew indexes

```bash
cd ${projDir}/tools/scalemet_dcis/
git clone https://github.com/mulqueenr/scalemet_dcis.git /data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis #this repo
```
Lane 1
```bash
#generate sample sheet using a modified version of bcl_convert_sheet.py to allow for pcr plate specifications.
#scalebio indexes batch 1
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/ScaleMethyl/references/lib.json ${bclDir}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/novaseq_lane1_batch1_scalebio_plate1-6.csv

#sample sheets for homebrew plates, run one plate at a time and merge
#batch 1 plate 2
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir}/RunInfo.xml \
--splitFastq --i7Set A --i5Set 1 > ${runDir}/samplesheets/novaseq_lane1_batch1_homebrew_plate2.csv
```
Lane 2
```bash
#scalebio indexes batch 2
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.batch2.csv ${projDir}/tools/ScaleMethyl/references/lib.json ${bclDir}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/novaseq_lane2_batch2_scalebio_plate3-4.csv

#batch 3 plate 1
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch3.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir}/RunInfo.xml \
--splitFastq --i7Set E --i5Set 1 > ${runDir}/samplesheets/novaseq_lane2_batch3_homebrew_plate1.csv
```

Lane 3
```bash
#scalebio indexes batch 3
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.batch3.csv ${projDir}/tools/ScaleMethyl/references/lib.json ${bclDir}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/novaseq_lane3_batch3_scalebio_plate2-5.csv

#batch 3 plate 6
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch3.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir}/RunInfo.xml \
--splitFastq --i7Set E --i5Set 2 > ${runDir}/samplesheets/novaseq_lane3_batch3_homebrew_plate6.csv
```

Lane 4
```bash
#scalebio indexes batch prelim1 and prelim2
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.prelim1.csv ${projDir}/tools/ScaleMethyl/references/lib.json ${bclDir}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/novaseq_lane4_prelim1-2_scalebio_plate1-1.csv

#batch 3 plate 16
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch3.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir}/RunInfo.xml \
--splitFastq --i7Set D --i5Set 4 > ${runDir}/samplesheets/novaseq_lane4_batch3_homebrew_plate16.csv
```

Lane 5
```bash
#scalebio indexes batch prelim3
python ${projDir}/tools/ScaleMethyl/bin/bcl_convert_sheet.py \
${runDir}/samplesheets/samples.prelim3.csv ${projDir}/tools/ScaleMethyl/references/lib.json ${bclDir}/RunInfo.xml \
--splitFastq > ${runDir}/samplesheets/novaseq_lane5_prelim3_scalebio_plate1-2.csv

#batch 3 plate 17
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch3.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir}/RunInfo.xml \
--splitFastq --i7Set D --i5Set 5 > ${runDir}/samplesheets/novaseq_lane5_batch3_homebrew_plate17.csv
```

Lane 6
```bash
#batch 1 plate 3
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir}/RunInfo.xml \
--splitFastq --i7Set C --i5Set 1 > ${runDir}/samplesheets/novaseq_lane6_batch1_homebrew_plate3.csv

#batch 1 plate 4
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir}/RunInfo.xml \
--splitFastq --i7Set C --i5Set 2 > ${runDir}/samplesheets/novaseq_lane6_batch1_homebrew_plate4.csv

#batch 1 plate 7
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir}/RunInfo.xml \
--splitFastq --i7Set B --i5Set 2 > ${runDir}/samplesheets/novaseq_lane6_batch1_homebrew_plate7.csv

#batch 1 plate 8
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir}/RunInfo.xml \
--splitFastq --i7Set C --i5Set 3 > ${runDir}/samplesheets/novaseq_lane6_batch1_homebrew_plate8.csv

#if multiple homebrew plates, combine with this
cat ${runDir}/samplesheets/novaseq_lane6_batch1_homebrew_plate3.csv > ${runDir}/samplesheets/novaseq_lane6_batch1_homebrew_plate3-4-7-8.csv
grep "^ScaleMethyl_" ${runDir}/samplesheets/novaseq_lane6_batch1_homebrew_plate4.csv >> ${runDir}/samplesheets/novaseq_lane6_batch1_homebrew_plate3-4-7-8.csv
grep "^ScaleMethyl_" ${runDir}/samplesheets/novaseq_lane6_batch1_homebrew_plate7.csv >> ${runDir}/samplesheets/novaseq_lane6_batch1_homebrew_plate3-4-7-8.csv
grep "^ScaleMethyl_" ${runDir}/samplesheets/novaseq_lane6_batch1_homebrew_plate8.csv >> ${runDir}/samplesheets/novaseq_lane6_batch1_homebrew_plate3-4-7-8.csv
```

Lane 7
```bash
#batch 1 plate 9
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir}/RunInfo.xml \
--splitFastq --i7Set D --i5Set 3 > ${runDir}/samplesheets/novaseq_lane7_batch1_homebrew_plate9.csv

#batch 1 plate 11
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir}/RunInfo.xml \
--splitFastq --i7Set A --i5Set 4 > ${runDir}/samplesheets/novaseq_lane7_batch1_homebrew_plate11.csv

#batch 1 plate 12
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch1.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir}/RunInfo.xml \
--splitFastq --i7Set B --i5Set 4 > ${runDir}/samplesheets/novaseq_lane7_batch1_homebrew_plate12.csv

#batch 2 plate 1
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch2.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir}/RunInfo.xml \
--splitFastq --i7Set D --i5Set 1 > ${runDir}/samplesheets/novaseq_lane7_batch2_homebrew_plate1.csv

#if multiple homebrew plates, combine with this
cat ${runDir}/samplesheets/novaseq_lane7_batch1_homebrew_plate9.csv > ${runDir}/samplesheets/novaseq_lane7_batch1-2_homebrew_plate9-11-12-1.csv
grep "^ScaleMethyl_" ${runDir}/samplesheets/novaseq_lane7_batch1_homebrew_plate11.csv >> ${runDir}/samplesheets/novaseq_lane7_batch1-2_homebrew_plate9-11-12-1.csv
grep "^ScaleMethyl_" ${runDir}/samplesheets/novaseq_lane7_batch1_homebrew_plate12.csv >> ${runDir}/samplesheets/novaseq_lane7_batch1-2_homebrew_plate9-11-12-1.csv
grep "^ScaleMethyl_" ${runDir}/samplesheets/novaseq_lane7_batch2_homebrew_plate1.csv >> ${runDir}/samplesheets/novaseq_lane7_batch1-2_homebrew_plate9-11-12-1.csv
```

Lane 8
```bash
#batch 2 plate 2
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch2.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir}/RunInfo.xml \
--splitFastq --i7Set D --i5Set 2 > ${runDir}/samplesheets/novaseq_lane8_batch2_homebrew_plate2.csv

#batch 3 plate 11
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch3.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir}/RunInfo.xml \
--splitFastq --i7Set A --i5Set 5 > ${runDir}/samplesheets/novaseq_lane8_batch3_homebrew_plate11.csv

#batch 3 plate 13
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch3.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir}/RunInfo.xml \
--splitFastq --i7Set B --i5Set 5 > ${runDir}/samplesheets/novaseq_lane8_batch3_homebrew_plate13.csv

#batch 3 plate 15
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
${runDir}/samplesheets/samples.batch3.csv ${projDir}/tools/scalemet_dcis/ref/homebrew.lib.json ${bclDir}/RunInfo.xml \
--splitFastq --i7Set C --i5Set 4 > ${runDir}/samplesheets/novaseq_lane8_batch3_homebrew_plate15.csv

#if multiple homebrew plates, combine with this
cat ${runDir}/samplesheets/novaseq_lane8_batch2_homebrew_plate2.csv > ${runDir}/samplesheets/novaseq_lane8_batch2-3_homebrew_plate2-11-13-15.csv
grep "^ScaleMethyl_" ${runDir}/samplesheets/novaseq_lane8_batch3_homebrew_plate11.csv >> ${runDir}/samplesheets/novaseq_lane8_batch2-3_homebrew_plate2-11-13-15.csv
grep "^ScaleMethyl_" ${runDir}/samplesheets/novaseq_lane8_batch3_homebrew_plate13.csv >> ${runDir}/samplesheets/novaseq_lane8_batch2-3_homebrew_plate2-11-13-15.csv
grep "^ScaleMethyl_" ${runDir}/samplesheets/novaseq_lane8_batch3_homebrew_plate15.csv >> ${runDir}/samplesheets/novaseq_lane8_batch2-3_homebrew_plate2-11-13-15.csv
```


```bash
#make nf-core input sheet
cd $runDir
echo """id,samplesheet,lane,flowcell
lane1_batch1_homebrew_plate2,${runDir}/samplesheets/novaseq_lane1_batch1_homebrew_plate2.csv,1,${bclDir}
lane1_batch1_scalebio_plate1-6,${runDir}/samplesheets/novaseq_lane1_batch1_scalebio_plate1-6.csv,1,${bclDir}
lane2_batch3_homebrew_plate1,${runDir}/samplesheets/novaseq_lane2_batch3_homebrew_plate1.csv,2,${bclDir}
lane2_batch2_scalebio_plate3-4,${runDir}/samplesheets/novaseq_lane2_batch2_scalebio_plate3-4.csv,2,${bclDir}
lane3_batch3_homebrew_plate6,${runDir}/samplesheets/novaseq_lane3_batch3_homebrew_plate6.csv,3,${bclDir}
lane3_batch3_scalebio_plate2-5,${runDir}/samplesheets/novaseq_lane3_batch3_scalebio_plate2-5.csv,3,${bclDir}
lane4_batch3_homebrew_plate16,${runDir}/samplesheets/novaseq_lane4_batch3_homebrew_plate16.csv,4,${bclDir}
lane4_prelim1-2_scalebio_plate1-1,${runDir}/samplesheets/novaseq_lane4_prelim1-2_scalebio_plate1-1.csv,4,${bclDir}
lane5_batch3_homebrew_plate17,${runDir}/samplesheets/novaseq_lane5_batch3_homebrew_plate17.csv,5,${bclDir}
lane5_prelim3_scalebio_plate1-2.,${runDir}/samplesheets/novaseq_lane5_prelim3_scalebio_plate1-2.csv,5,${bclDir}
lane6_batch1_homebrew_plate3-4-7-8,${runDir}/samplesheets/novaseq_lane6_batch1_homebrew_plate3-4-7-8.csv,6,${bclDir}
lane7_batch1-2_homebrew_plate9-11-12-1,${runDir}/samplesheets/novaseq_lane7_batch1-2_homebrew_plate9-11-12-1.csv,7,${bclDir}
lane8_batch2-3_homebrew_plate2-11-13-15,${runDir}/samplesheets/novaseq_lane8_batch2-3_homebrew_plate2-11-13-15.csv,8,${bclDir}
""" > pipeline_samplesheet.csv


echo """id,samplesheet,lane,flowcell
lane1_batch1_homebrew_plate2,${runDir}/samplesheets/novaseq_lane1_batch1_homebrew_plate2.csv,1,${bclDir}
""" > pipeline_samplesheet_homebrewrerun.csv
```

BCL to FASTQ

```bash
mkdir -p $SCRATCH/scalemet_work

#use nf-core to split out fastqs
nextflow run nf-core/demultiplex \
    --input pipeline_samplesheet_homebrewrerun.csv \
    --outdir fastq \
    --trim_fastq false \
    --remove_adapter false \
    --skip_tools fastp,fastqc,kraken,multiqc,checkqc,falco,md5sum,samshee \
    -profile singularity \
    --first-tile-only true \
    -w ${SCRATCH}/scalemet_work


#get percentage of reads unassigned
#actually i think i can just get this from pipeline output
assigned_reads(){
    lane_plate=$(ls $1/L*/)
    assigned=$(zcat $lane_plate/ScaleMethyl*_I1_*fq.gz | grep "^@" | wc -l)
    undetermined=$(zcat $lane_plate/Undetermined*_I1_*fq.gz | grep "^@" | wc -l)
    echo "$(basename $lane_plate),$assigned,undetermined"
    }

export -f lib_complexity
parallel -j 200 lib_complexity ::: $(ls *-1.bam)

ls -d ${runDir}/fastq/lane*/
```

RUN SCALE_METHYL PIPELINE

```bash
mkdir -p ${SCRATCH}/scalemet_native_work

cd ${runDir}

homebrew_lane_plate=$(ls ${runDir}/fastq/lane*homebrew*/L*/
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

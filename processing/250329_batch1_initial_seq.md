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
BCMDCIS32T,1A01-1A06;2A07-2A12;3A01-3A03,BCMDCIS32T
BCMHBCA03R,1B01-1B06;2B07-2B12;3B01-3B03;3E04-3E06;3E10-3E12,BCMHBCA03R
BCMDCIS94T-24hTis,1C01-1C06;2C07-2C12;3C01-3C03,BCMDCIS94T-24hTis
BCMHBCA26L-24hTis_4h,1D01-1D06;2D07-2D12;3D01-3D03;3G01-3G06;3G10-3G12;1G01-1G06;2G07-2G12;,BCMHBCA26L-24hTis-4h
ECIS36T,1E01-1E06;2E07-2E12;3E01-3E03;3F04-3F06;3F10-3F12,ECIS36T
BCMDCIS74T,1F01-1F06;2F07-2F12;3F01-3F03,BCMDCIS74T
BCMDCIS22T,1H01-1H06;2H07-2H12;3H01-3H03,BCMDCIS22T
BCMDCIS66T,1A07-1A12;2A01-2A06;3A04-3A12,BCMDCIS66T
BCMDCIS28T,1B07-1B12;2B01-2B06;3B04-3B12,BCMDCIS28T
BCMDCIS99T,1C07-1C12;2C01-2C06;3C07-3C09,BCMDCIS99T
BCMDCIS52T,1D07-1D12;1F07-1F12;2D01-2D06;2F01-2F06;3H04-3H06;3D07-3D09;3F07-3F09;3H10-3H12,BCMDCIS52T
BCMDCIS80T,1E07-1E12;2E01-2E06;3D04-3D06;3E07-3E09;3D10-3D12,BCMDCIS80T
BCMDCIS07T,1G07-1G12;2G01-2G06;3C04-3C06;3G07-3G09;3C10-3C12,BCMDCIS07T
BCMHBCA38L-3h,1H07-1H12;2H01-2H06;3H07-3H09,BCMHBCA38L-3h""" > ${samples}

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

#i5 is revcomp
#i7 is seq
Sample_ID,index,index2
BCMDCIS32T_A01,GGAGGCCTCC,ATGGAGCTAC
BCMDCIS32T_A01,CAGCAGTATC,ATGGAGCTAC

#update Override Cycles in samplesheet to fit sequencing run
#OverrideCycles,Y50;I10;I10N14;Y47

#run bcl-convert within amethyst sif
cd ${runDir}
singularity shell \
--bind ${bclDir} \
--bind ${runDir} \
--bind ${bclDir}:/var/log/bcl-convert \
~/singularity/amethyst.sif
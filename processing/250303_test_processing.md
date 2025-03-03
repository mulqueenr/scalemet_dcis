[See notebook for wet lab processing.](https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio%20sciMETv2.one%7CD3F046A2-B151-0443-938E-82A415D420EB%2F240910%20ScaleBio%20DCIS%20Samples%7C30534461-040E-C54F-BB40-7D53F8115495%2F%29)

[See updated experiment notes for homebrew comparison.](https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio%20sciMETv2.one%7CD3F046A2-B151-0443-938E-82A415D420EB%2F250214%20ScaleBio%20Homebrew%20Sorting%20Extra%20Plates%7CA5A091ED-32C0-D24E-81B9-1EFF5EA1252B%2F%29)

```bash
#Should still only need to identify samples at the tagmentation level, and expanding the i5.txt and i7.txt should take care of itself.


echo """sample,barcodes,libName
DCIS-92T,1A01-1D12,ScaleMethyl
DCIS-66T,1E01-1H12,ScaleMethyl
DCIS-79T,2A01-2D12,ScaleMethyl
IDC-79T,2E01-2H12,ScaleMethyl
HBCA-19T,3A01-3D12,ScaleMethyl
HBCA-17T,3E01-3H12,ScaleMethyl""" > ${runDir}/samples.csv


#build proper formated singularity container
singularity build ~/singularity/scalemethyl_v1.6.sif ~/singularity/public.ecr.aws-o5l3p3e4-scale-methyl-tools@sha256-6fd63db48e8786ed1cfc17d7e3effd3fd696ccb8e5e54803959e2dcd2f794aec.img

#set up directories and variables
proj_dir="/volumes/USR2/Ryan/projects/metact"
scalebio_nf="${proj_dir}/tools2/ScaleMethyl" #tools is a symlink directory so it wasn't mounting properly for singularity
runDir="${proj_dir}/241007_RM_scalebio_dcis2"
genome="${proj_dir}/ref/reference/genome.json"
fastqDir="${runDir}/241004_A01819_0637_BHY5MJDMXY/241004_A01819_0637_BHY5MJDMXY"
samples="${runDir}/samples.csv"

export SCRATCH="/volumes/USR2/Ryan/scratch/scalemet_work"
export TMPDIR="/volumes/USR2/Ryan/scratch"
export NXF_SINGULARITY_CACHEDIR="/volumes/USR2/Ryan/singularity"
export SINGULARITY_BINDPATH="/volumes/seq/projects/metACT/tools/ScaleMethyl/bin" 

mkdir -p $SCRATCH

source activate conda #(to use more recent java version)
cd $runDir
nextflow run ${scalebio_nf} \
--runFolder ${fastqDir} \
--samples ${runDir}/samples.csv \
--outDir ${runDir} \
--genome ${genome} \
--maxMemory 500.GB \
--maxCpus 300 \
-profile singularity \
-w ${SCRATCH}/scalemet_work \
-resume

```

[See notebook for wet lab processing.](https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio%20sciMETv2.one%7CD3F046A2-B151-0443-938E-82A415D420EB%2F240910%20ScaleBio%20DCIS%20Samples%7C30534461-040E-C54F-BB40-7D53F8115495%2F%29)

[See updated experiment notes for homebrew comparison.](https://mdandersonorg-my.sharepoint.com/personal/rmulqueen_mdanderson_org/_layouts/OneNote.aspx?id=%2Fpersonal%2Frmulqueen_mdanderson_org%2FDocuments%2FmetACT&wd=target%28scalebio%20sciMETv2.one%7CD3F046A2-B151-0443-938E-82A415D420EB%2F250214%20ScaleBio%20Homebrew%20Sorting%20Extra%20Plates%7CA5A091ED-32C0-D24E-81B9-1EFF5EA1252B%2F%29)

```bash
#Should still only need to identify samples at the tagmentation level, and expanding the i5.txt and i7.txt should take care of itself.
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
runDir="${projDir}/data/250307_RM_scalebio_dcis2_homebrew"
bclDir="/data/rmulqueen/projects/scalebio_dcis/seq/250227_RM_CurioWGS_scalemet"

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
python ${projDir}/tools/scalemet_dcis/src/bcl_convert_sheet_pcr.py \
samples.csv \
${projDir}/tools/ScaleMethyl/references/lib.json \
${bclDir}/RunInfo.xml --splitFastq \
--i7Set C \
--i5Set 3 > samplesheet.csv

#update Override Cycles in samplesheet to fit sequencing run
#OverrideCycles,Y50;I10;I10N14;Y47

#run bcl-convert within amethyst sif
cd ${runDir}
singularity shell \
--bind ${bclDir} \
--bind ${runDir} \
--bind ${bclDir}:/var/log/bcl-convert \
~/singularity/amethyst.sif

projDir="/data/rmulqueen/projects/scalebio_dcis"
bclDir="/data/rmulqueen/projects/scalebio_dcis/seq/250227_RM_CurioWGS_scalemet"
runDir="${projDir}/data/250307_RM_scalebio_dcis2_homebrew"
cd ${runDir}
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


Correct all symlink with actual links

```bash

#set up functions
#count reads export
count_reads() { 
        samtools view $1 | awk -v b=$1 '{split($1,a,":"); print a[1],b}' | sort | uniq -c | sort -k1,1n
}

#split bams export
split_bams() { 
        test=$1
        idx=$(echo $test | cut -d ':' -f 8 )
        outidx=$(echo $idx | sed -e 's/+/_/g' -)
        bam=$(echo $test | cut -d ' ' -f 3)
        outprefix=$(echo $bam | cut -d '/' -f 2)
        ((samtools view -H $bam) && (samtools view $bam | awk -v i=$idx '{split($1,a,":"); if(a[1]==i); print $0}')) | samtools view -bS > ./sc_bams/${outprefix}.${idx}.bam
}

export -f count_reads
export -f split_bams

#parallelize it
scale_runDir="/data/rmulqueen/projects/scalebio_dcis/data/250329_RM_scalebio_batch1_initseq/scale_dat"
homebrew_runDir="/data/rmulqueen/projects/scalebio_dcis/data/250329_RM_scalebio_batch1_initseq/homebrew_dat"

cd $scale_runDir
parallel -j 100 count_reads ::: $(find ./ -maxdepth 4 -name '*bam' -type l ) | sort -k1,1n > scale_unique_read_counts.tsv

#filter to bam files with >100000 unique reads
awk '$1>100000 {print $0}' unique_read_counts.tsv > cells_pf.txt
mkdir -p sc_bams
parallel -j 60 -a cells_pf.txt split_bams

scale_runDir="/data/rmulqueen/projects/scalebio_dcis/data/250329_RM_scalebio_batch1_initseq/scale_dat"


find ${scale_runDir} -type l -exec bash -c 'cp -R "$(readlink -m "$0")" ./dedup_bams' {} \; #scale pipeline makes empty files, this throws errors for empty files (but can be ignored)

find $runDir -xtype l -exec bash -c 'target="$(readlink "{}")"; link="{}"; target="$(echo "$target" | sed "s/\/rsrch4\/home\/genetics\/rmulqueen\/projects\/metact/\/volumes\/seq\/projects\/metACT/g")"; ln -Tfs "$target" "$link"' \;
```
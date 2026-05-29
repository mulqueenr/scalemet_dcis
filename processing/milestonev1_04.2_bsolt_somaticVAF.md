cd /data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/scalemethyl_pipeline_out/prelim1-2_scalebio_plate1-1/alignments/dedup/BCMDCIS41T.3B04

#script to make a merged bam per clone, and somatic lineage

bam="BCMDCIS41T.3B04.dedup.bam"
bsbolt BamIndex \
-I $bam

bsbolt CallVariation \
-t 10 \
-I $bam \
-DB /data/rmulqueen/projects/scalebio_dcis/ref/hg38_bsbolt \
-O ${bam::-4}

#make temp bams per clone (and somatic by non-epithelial cells) and call variants
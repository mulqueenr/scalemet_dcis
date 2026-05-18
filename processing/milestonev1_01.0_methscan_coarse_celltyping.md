# Milestone v1 Preprocessing
Processing AllC files into intial filtering and clustering via methscan.
Import cluster info into amethyst for further processing.


Make symlinks of allC files for all samples. Add naming scheme to avoid index collisions across lanes.

```bash
mkdir -p /data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/allc
mkdir -p /data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/allc/CH
mkdir -p /data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/allc/CG

cd /data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/allc

allc_output_link() {
    plate=$(awk -v name="$1" 'BEGIN {split(name,a,"/"); print a[9]}')
    sample=$(awk -v name="$1" 'BEGIN {split(name,a,"/"); print a[13]}')
    met=$(awk -v name="$1" 'BEGIN {split(name,a,"/"); print a[14]}')
    file_name=$(basename $1)
    tbi_file_name="$1.tbi"
    echo "Processing file: $1"
    echo "$plate $sample $met $file_name"
    if [[ "$met" == "CG" ]]; then
        echo "CG"
        ln -s $1 /data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/allc/CG/${plate}.${sample}.${met}.${file_name}
        ln -s $tbi_file_name /data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/allc/CG/${plate}.${sample}.${met}.${file_name}.tbi
    else
        echo "CH"
        ln -s $1 /data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/allc/CH/${plate}.${sample}.${met}.${file_name}
        ln -s $tbi_file_name /data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/allc/CH/${plate}.${sample}.${met}.${file_name}.tbi
    fi
}
export -f allc_output_link
#    

find /data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/scalemethyl_pipeline_out/ -maxdepth 7 -type f -name "*.allc.tsv.gz" | parallel -j 50 allc_output_link {} #symlink of all allc files to directory

```

```bash
mkdir -p /data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/methscan/sc_met

cd /data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/allc/CG/
ulimit -s 65536 #expand stack size for long argument list

methscan prepare \
--round-sites \
--input-format 'allc' \
--chunksize 1000000000 \
*.allc.tsv.gz \
/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/methscan/sc_met

```

Following tutorial for initial processing
https://anders-biostat.github.io/MethSCAn/tutorial.html

```R
library(tidyverse)
setwd("/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/methscan/sc_met")
cell_stats <- read_csv("cell_stats.csv")

head(cell_stats)

plt<-cell_stats %>% 
  ggplot(aes(x = global_meth_frac * 100, y = log10(n_obs))) +
  geom_point(size=0.1,alpha=0.5) + xlim(c(0,100)) +
  labs(x = "global DNA methylation %", y = "# of observed CpG sites") + 
  geom_rect(ymin=log10(20000),ymax=log10(1e7),xmin=20,xmax=90,color="red",fill=NA) +
  theme_minimal()
ggsave(plt,file="methylation_by_obs.unfiltered.pdf")
```

Now filter cells, then smooth to identify VMRs
```bash
methscan filter --min-sites 20000 --min-meth 20 --max-meth 90 sc_met filtered_data #loose filtering of data
methscan smooth filtered_data
methscan scan --threads 100 filtered_data VMRs.bed
#Wrote 68452 VMRs with sequencing coverage in at least 6 cells to VMRs.bed
#The columns in this file correspond to:
#chromosome, VMR start, VMR end, variance, # of CpG sites, # of cells with coverage in the VMR
#Finished scan on Tue May 05 09:40:25 2026. Total runtime was 0:17:19 (hour:min:s).

methscan matrix --threads 200 --sparse VMRs.bed filtered_data VMR_matrix

```

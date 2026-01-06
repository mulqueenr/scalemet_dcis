```bash
mkdir -p /data/rmulqueen/projects/scalebio_dcis/rna
mkdir -p /data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output
cd /data/rmulqueen/projects/scalebio_dcis/rna/cellranger_output

fq_dir=(
"/geo_seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/Fastq_output/BCMHBCA03R"
"/geo_seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/Fastq_output/BCMHBCA04R_2h"
"/geo_seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/Fastq_output/BCMHBCA09R"
"/geo_seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/Fastq_output/BCMHBCA12R_3h"
"/geo_seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/Fastq_output/BCMHBCA16R_3h"
"/geo_seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/Fastq_output/BCMHBCA16R_nu"
"/geo_seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/Fastq_output/BCMHBCA17R"
"/geo_seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/Fastq_output/BCMHBCA17R-NU"
"/geo_seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/Fastq_output/BCMHBCA19R_4h"
"/geo_seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/Fastq_output/BCMHBCA19R-NU"
"/geo_seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/Fastq_output/BCMHBCA22R_4h"
"/geo_seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/Fastq_output/BCMHBCA22R_nu"
"/geo_seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/Fastq_output/BCMHBCA26L_24hTis_4h"
"/geo_seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/Fastq_output/BCMHBCA26L_nu"
"/geo_seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/Fastq_output/BCMHBCA29L_2h"
"/geo_seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/Fastq_output/BCMHBCA38R"
"/geo_seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/Fastq_output/BCMHBCA83L_3h"
"/geo_seq/projects/HBCA/BCMHBCA/BCMHBCA_10x_RNA/Fastq_output/BCMHBCA85L_3h"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS05T"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS07T"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS41T"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS66T"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS80T"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS82T"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS99T"
"/geo_seq/projects/ECIS/ECIS_cellranger_v3_10/fastq/ECIS25T"
"/geo_seq/projects/ECIS/ECIS_cellranger_v3_10/fastq/ECIS36T_24hTis"
"/geo_seq/projects/ECIS/ECIS_cellranger_v3_10/fastq/ECIS26T"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS79T_DCIS"
"/geo_USR3/ctang4/10X_processing_RNA/090924_BCMDCIS_3_Rui/Fastq/outs/fastq_path/2223MNMNX/BCMDCIS124T"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS22T"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS28T"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS32T"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS35T"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS49T_24htis"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS52T"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS52T_novaseq"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS65T"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS70T"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS74T"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS79T_IDC"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS92T_24h"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS94T_24h"
"/geo_seq/projects/ECIS/BCMDCIS_cellranger_v3_10/Fastq/BCMDCIS97T"
"/geo_seq/projects/ECIS/ECIS_cellranger_v3_10/fastq/ECIS48T"
"/geo_seq/projects/ECIS/ECIS_cellranger_v3_10/fastq/ECIS57T")

samp_name=(
'BCMDCIS102T_24hTis'
'BCMDCIS124T'
'BCMDCIS05T'
'BCMDCIS07T'
'BCMDCIS22T'
'BCMDCIS28T'
'BCMDCIS32T'
'BCMDCIS35T'
'BCMDCIS41T'
'BCMDCIS49T'
'BCMDCIS52T'
'BCMDCIS65T'
'BCMDCIS66T'
'BCMDCIS70T'
'BCMDCIS74T'
'BCMDCIS79T_24hTis_DCIS'
'BCMDCIS79T_24hTis_IDC'
'BCMDCIS80T_24hTis'
'BCMDCIS82T_24hTis'
'BCMDCIS92T_24hTis'
'BCMDCIS94T_24hTis'
'BCMDCIS97T'
'BCMDCIS99T'
'ECIS25T'
'ECIS26T'
'ECIS36T'
'ECIS48T'
'ECIS57T'
'BCMHBCA03R'
'BCMHBCA04R'
'BCMHBCA09R-3h'
'BCMHBCA12R-3h'
'BCMHBCA16R-3h'
'BCMHBCA16R-3h-nuc'
'BCMHBCA17R-3h'
'BCMHBCA17R-3h-nuc'
'BCMHBCA19R-4h'
'BCMHBCA19R-4h-nuc'
'BCMHBCA22R-4h'
'BCMHBCA22R-4h-nuc'
'BCMHBCA26L-24hTis-4h'
'BCMHBCA26L-24hTis-4h-nuc'
'BCMHBCA29L-2h'
'BCMHBCA38L-3h'
'BCMHBCA83L-3h'
'BCMHBCA85L-3h'
)

for i in $(seq ${#samp_name[@]}); do
/home/rmulqueen/tools/cellranger-9.0.1/cellranger count \
--fastqs=${fq_dir[i]} \
--id=${samp_name[i]} \
--transcriptome=/home/rmulqueen/ref/refdata-gex-GRCh38-2024-A \
--create-bam=true \
--include-introns=true \
--localcores=200 \
--localmem=100; 
done
```


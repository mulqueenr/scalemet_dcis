transposable element methylation analysis

copy these to ref folder
- check cnv overlap with haploinsuffiency and triplosensitivity markers
- check methylation scores for repeats (expecting lower methylation in cancer for repeats)

## cgi
https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=regulation&hgta_track=cpgIslandExt&hgta_table=cpgIslandExt&hgta_doSchema=describe+table+schema
```bash
cd /data/rmulqueen/projects/scalebio_dcis/ref
wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz"
/data/rmulqueen/projects/scalebio_dcis/ref/cpgIslandExt.bed
```

## haploinsufficiency and triplosensitivity tracks
https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=3596528659_3AWPx4Dyqv9CtPn24flp4abaaA3K&db=hg38&c=chr6&g=dosageSensitivity
```bash
cd /data/rmulqueen/projects/scalebio_dcis/ref
wget "https://zenodo.org/records/6347673/files/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz?download=1"
mv Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz?download=1 Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz
/data/rmulqueen/projects/scalebio_dcis/ref/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz
```

## repeats, lines sines ltr etc
https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_doSchema=describe+table+schema

```bash
cd /data/rmulqueen/projects/scalebio_dcis/ref
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
/data/rmulqueen/projects/scalebio_dcis/ref/rmsk.txt.gz
```
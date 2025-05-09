# scalemet_dcis
 Repo for processing of dcis samples generated through scalebiomethyl kit and homebrew version

**TO DO:**
* Train scAGE on patient age samples https://github.com/alex-trapp/scAge/blob/main/notebooks/construct_reference_notebook.ipynb
* Methyltree https://github.com/ShouWenWang-Lab/MethylTree
* Add Plate metadata column for sorting gate tracking (all plates are from one respective gate)

See processing folder for details on intial processing for each sample batch, and sequencing_runs.md within that for raw data locations and links to notebook pages.
The pipeline is set up as a nextflow module and pulls from multiple singularity SIFs for environments.

## Directory Structure
* ref : Reference files, including i5 and i7 index additions, and oligo design files. As well as markdown for how to prepare bsbolt reference and amethyst processing singularity file.

* processing : Markdown files of how samples were processed on our compute server (SYSBIO). 

* src : custom scripts used in processing, including Rscript of custom amethyst functions, addition of MethylTree etc.

## Set up scalemet methylation pipeline
Adding homebrew indexes for i5 and i7 so we can use sciMETv2 on excess tagmented nuclei from Scale Bio Methyl kits.

```bash
#clone scalemethyl pipeline v 1.2.2
mkdir -p /data/rmulqueen/projects/scalebio_dcis/tools
cd /data/rmulqueen/projects/scalebio_dcis/tools
git clone https://github.com/ScaleBio/ScaleMethyl.git /data/rmulqueen/projects/scalebio_dcis/tools/ScaleMethyl #scalemethyl repo
git clone https://github.com/mulqueenr/scalemet_dcis.git /data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis #this repo

##Add zenodo hosted SIFs here later

```

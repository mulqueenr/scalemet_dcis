# scalemet_dcis
 Repo for processing of dcis samples generated through scalebiomethyl kit and homebrew version

## Directory Structure
* ref : Reference files, including i5 and i7 index additions, and oligo design files. As well as markdown for how to prepare bsbolt reference and amethyst processing singularity file.

* processing : Markdown files of how samples were processed on our compute server GEO. Maybe re-write for seadragon if I'm hogging the cluster too much.

* src : custom scripts used in processing, including Rscript of custom amethyst functions, addition of MethylTree etc.

## Set up scalemet methylation pipeline
Adding homebrew indexes for i5 and i7 so we can use sciMETv2 on excess tagmented nuclei from Scale Bio Methyl kits.

```bash
#clone scalemethyl pipeline v 1.2.2
mkdir -p ~/projects/scalemet_dcis/tools
cd ~/projects/scalemet_dcis/tools
git clone https://github.com/ScaleBio/ScaleMethyl.git ~/projects/scalebio_dcis/tools/ScaleMethyl #scalemethyl repo
git clone https://github.com/mulqueenr/scalemet_dcis.git ~/projects/scalebio_dcis/tools/scalemet_dcis #this repo


#Run ScaleMethyl as normal for processing
#A default custom config and runparameters file for processing.

```
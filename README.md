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
git clone https://github.com/ScaleBio/ScaleMethyl.git #scalemethyl repo
git clone https://github.com/mulqueenr/scalemet_dcis.git #this repo

#modify i5 and i7 to include homebrew indexes as well
#i5 set goes from 1-4, so adding on as 5+
#i7 set goes from A-B, so adding on as C+
#I'm aware there is a collision between homebrew G10 i7 and kit A10 i7, but my plan is to run them on separate flow cells anyway.

cp scalemet_dcis/ref/i5_with_homebrew.txt ScaleMethyl/references/i5.txt
cp scalemet_dcis/ref/i7_with_homebrew.txt ScaleMethyl/references/i7.txt

#Run ScaleMethyl as normal for processing
#A default custom config and runparameters file for processing.

```
# scalemet_dcis
 Repo for processing of dcis samples generated through scalebiomethyl kit and homebrew version

# Set up scalemet methylation pipeline

```bash
#clone scalemethyl pipeline v 1.2.2
mkdir -p ~/projects/scalemet_dcis/tools
cd ~/projects/scalemet_dcis/tools
git clone https://github.com/ScaleBio/ScaleMethyl.git
git clone https://github.com/ScaleBio/ScaleMethyl.git

#modify i5 and i7 to include homebrew indexes as well
#i5 set goes from 1-4, so adding on as 5+
#i7 set goes from A-B, so adding on as C+
/ref/i5_with_homebrew.txt
/ref/i7_with_homebrew.txt

mv 
```
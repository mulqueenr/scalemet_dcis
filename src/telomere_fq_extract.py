import gzip
from Bio import SeqIO
import sys
import os
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from scipy.spatial import distance
import numpy as np
import pandas as pd
from collections import Counter
import argparse
from multiprocessing import Pool

parser = argparse.ArgumentParser()
parser.add_argument('--project_dir',default="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1")
parser.add_argument('--plate',default="batch1_homebrew_plate11")
parser.add_argument('--outdir',default="telomere_fq")
parser.add_argument('--cellline_meta',default="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/merged_data/01_celllines.metadata.tsv")
parser.add_argument('--dcis_meta',default="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/merged_data/05_scaledcis.fine_celltype.metadata.tsv")
parser.add_argument('--cores',default=50)
args = parser.parse_args()

fq_dir=args.project_dir+"/fastq" #fastq directory based on project dir structure
file_list=[os.path.join(dirpath,f) for (dirpath, dirnames, filenames) in os.walk(fq_dir+"/"+args.plate) for f in filenames if "_R1_001.fastq.gz" in f] #find all R1 files in plate directory

output_directory=args.project_dir+"/"+args.outdir+"/"+args.plate
os.system("mkdir -p "+args.project_dir+"/"+args.outdir)   #make output directory
os.system("mkdir -p "+output_directory)  #make plate output directory

df_cellline = pd.read_table(args.cellline_meta)
df_dcis = pd.read_table(args.dcis_meta)
df_idx_whitelist = pd.concat([df_cellline, df_dcis])
df_idx_whitelist = df_idx_whitelist[df_idx_whitelist['plate_info']==args.plate]

i7_whitelist=list(set([x.split("+")[0] for x in df_idx_whitelist.index]))
i5_whitelist=list(set([x.split("+")[1] for x in df_idx_whitelist.index]))
tn5_whitelist=list(set([x.split("+")[2] for x in df_idx_whitelist.index]))

#possible telomere motifs (all aleast a repeat of 2 motifs)
telo_fwd="GGGATTGGGATT"
telo_rev="GGGATTGGGATT"
telo_fwdcomp="CCCTAACCCTAA"
telo_revcomp="AATCCCAATCCC"
telo_fwdcomp_conv="TTTTAATTTTAA"
telo_revcomp_conv="AATTTTAATTTT"
telo_fwd_conv="AAAATTAAAATT"
telo_rev_conv="TTAAAATTAAAA"
telo_list=[telo_fwd,telo_rev,telo_fwdcomp,telo_revcomp,telo_fwdcomp_conv,telo_revcomp_conv,telo_fwd_conv,telo_rev_conv]

#count telomere reads
#write out telomere reads fq to another directory
def count_telo_fq(r1):
    fq1_file=r1
    fq2_file=r1.replace("_R1_001.fastq.gz","_R2_001.fastq.gz")
    read_counter=[]
    telo_counter=[]
    print("Running FastQ split for: " + os.path.basename(r1))
    with gzip.open(fq1_file, "rt") as handle1, \
        gzip.open(fq2_file, "rt") as handle2, \
        open(output_directory+"/"+os.path.basename(fq1_file).replace("R1_001.","R1_001.TELO.")[:-3], "w") as outfile_fq1, \
        open(output_directory+"/"+os.path.basename(fq2_file).replace("R2_001.","R2_001.TELO.")[:-3], "w") as outfile_fq2:
            for (title1, seq1, qual1), (title2, seq2, qual2) in \
            zip(FastqGeneralIterator(handle1), FastqGeneralIterator(handle2)):
                idx1_2=title1.split(":")[-1]
                idx1=idx1_2.split("+")[0]
                idx2=idx1_2.split("+")[1]
                idx3=seq2[0:8]
                idx="%s+%s+%s" % (idx1_2, idx3, args.plate)
                read_counter.append(idx)
                if(any(motif in seq1 for motif in telo_list) or any(motif in seq2 for motif in telo_list)): 
                    if "%s+%s+%s+%s" % (idx1,idx2,idx3,args.plate) in df_idx_whitelist.index:
                        title1="%s+%s" % (title1, idx3)
                        title2="%s+%s" % (title2, idx3)
                        fq1="@%s\n%s\n+\n%s\n" % (title1, seq1, qual1)
                        fq2="@%s\n%s\n+\n%s\n" % (title2, seq2, qual2)
                        telo_counter.append(idx)
                        outfile_fq1.write(fq1_file)
                        outfile_fq2.write(fq2_file)
                    elif min([round(distance.hamming(Seq(idx1), Seq(i))*10) for i in i7_whitelist]) <= 2:
                        i7idx_hamming=[round(distance.hamming(Seq(idx1), Seq(i))*10) for i in i7_whitelist]
                        idx1=i7_whitelist[i7idx_hamming.index(min(i7idx_hamming))]
                        if (min([round(distance.hamming(Seq(idx2), Seq(i))*10) for i in i5_whitelist]) <= 2):
                            i5idx_hamming=[round(distance.hamming(Seq(idx2), Seq(i))*10) for i in i5_whitelist]
                            idx2=i5_whitelist[i5idx_hamming.index(min(i5idx_hamming))]
                            if (min([round(distance.hamming(Seq(idx3), Seq(i))*8) for i in tn5_whitelist]) <= 2):
                                tn5idx_hamming=[round(distance.hamming(Seq(idx3), Seq(i))*8) for i in tn5_whitelist]
                                idx3=tn5_whitelist[tn5idx_hamming.index(min(tn5idx_hamming))]
                                title1=":".join(title1.split(":")[0:len(title1.split(":"))-1])+":"+idx1+"+"+idx2+"+"+idx3
                                title2=":".join(title2.split(":")[0:len(title2.split(":"))-1])+":"+idx1+"+"+idx2+"+"+idx3
                                idx="%s+%s+%s+%s" % (idx1, idx2, idx3, args.plate)
                                fq1="@%s\n%s\n+\n%s\n" % (title1, seq1, qual1)
                                fq2="@%s\n%s\n+\n%s\n" % (title2, seq2, qual2)
                                telo_counter.append(idx)
                                outfile_fq1.write(fq1_file)
                                outfile_fq2.write(fq2_file)
    if len(telo_counter)>0 and len(read_counter)>0:
        print("Counting total reads and telomere reads: " + os.path.basename(r1))
        idx_count, val_count = np.unique(read_counter, return_counts=True)
        read_result=pd.DataFrame(list(zip(idx_count, val_count)))
        read_result.columns = ['cellID', 'totalReads']
        idx_telo, val_telo = np.unique(telo_counter, return_counts=True)
        telo_result=pd.DataFrame(list(zip(idx_telo, val_telo)))
        telo_result.columns = ['cellID', 'telomereReads']
        merged_df = pd.merge(read_result, telo_result, on='cellID', how='inner')
        merged_df.to_csv(output_directory+"/"+os.path.basename(fq1_file).replace("R1_001.fastq.gz","readCounts.tsv"), index=True)
        os.system('gzip -f '+ output_directory+"/"+os.path.basename(fq1_file).replace("R1_001.","R1_001.TELO.")[:-3])
        os.system('gzip -f '+ output_directory+"/"+os.path.basename(fq2_file).replace("R2_001.","R2_001.TELO.")[:-3])

p = Pool(int(args.cores))
with p:
    p.map(count_telo_fq,file_list)
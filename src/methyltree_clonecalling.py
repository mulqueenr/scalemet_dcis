import cospar as cs
cs.settings.set_figure_params()

import sys
sys.path.append('/container_src/MethylTree') #load in methyltree module
import methyltree 
import pandas as pd
import scanpy as sc
import os
import numpy as np
from matplotlib import pyplot as plt
import argparse
import cospar as cs

parser = argparse.ArgumentParser()
parser.add_argument('-i', "--input",default="DCIS-41T")  #description='sample name input (matching prefix to methyltree output from amethyst processing)'
args = parser.parse_args()
in_dir= "/volumes/USR2/Ryan/projects/metact/amethyst_processing"
os.chdir(in_dir)
prefix=args.input

print('Read in H5 Data') # 
dat=pd.read_hdf(prefix+"_methyltree_input.h5", key="data", mode='r')
dat.rename(columns={'cell_id': 'index'}, inplace=True)
df_out=dat.pivot(index='index',columns='genomic_region_id',values='value')
adata=sc.AnnData(df_out) #This adata has not annotation yet. We will use the df_sample from raw_data folder to annotate this object later

print('Read in metadata')
metadat=pd.read_hdf(prefix+"_methyltree_input.h5", key="metadata", mode='r')
metadat.to_csv("sample_sheet.tsv.gz",sep='\t',compression='gzip')

print('Set up output directories')
out_dir=prefix+'_methyltree'+'/data'
save_data_des=prefix
clone_key='large_clone_id'
data_des=prefix
figure_path=prefix+'_methyltree'+'/figure'
data_path=os.getcwd()
os.makedirs(out_dir,exist_ok=True)
os.makedirs(figure_path,exist_ok=True)
df_sample=methyltree.hf.load_sample_info(data_path)

adata_final,stat_out=methyltree.analysis.comprehensive_lineage_analysis(
                                                        out_dir,data_path,save_data_des,clone_key,
                                                        adata_orig=adata,
                                                        clone_color_dict=None,heatmap_additional_key_list=['celltype'],
                                                        compute_similarity=True,
                                                        update_sample_info=True,
                                                        # correct the correlation bias
                                                        similarity_correction=True,
                                                        # remove cell-type specific DNA methylation signals
                                                        remove_celltype_signal=True,cell_type_key='celltype',
                                                        # optimize the lineage tree
                                                        optimize_tree=True,optimize_background_cutoff=0.8,
                                                        # better heatmap visualization
                                                        heatmap_vmax_percentile=99,heatmap_vmin_percentile=70,heatmap_figsize=(10, 9.5),heatmap_show_label=True,
                                                        heatmap_show_legend=True,heatmap_fontsize=2,
                                                        # related to saving figures
                                                        fig_dir=figure_path,data_des=data_des,
                                                        # coarse-graining
                                                        perform_coarse_graining=False,#coarse_grain_figsize=(6, 5),
                                                        # infer the clones based on similarity cutoff
                                                        #perform_clone_inference=True,clone_inference_threshold=0.8,clone_inference_print=True,
                                                        perform_memory_analysis=True, save_adata=True,perform_depth_analysis=True,)

methyltree.plotting.plot_similarity_heatmap_with_multiple_colorbars(adata_final,save_name="test.svg",additional_key_list=['celltype'],heatmap_vmax_percentile=99.9,heatmap_vmin_percentile=80,)
fig = plt1.get_figure()
plt1.savefig("out.png") 

methyltree.lineage.bootstrap_lineage_tree(adata_final,out_dir,
    save_data_des,
    clone_key=clone_key,
    cell_type_key='celltype',remove_celltype_signal=True,
    num_iterations=10,  
    sample_fraction=0.8,  
    similarity_correction=True,
    heatmap_ordering_method='UPGMA',
    recompute=False,
    cores=50)

cs.settings.set_figure_params(fontsize=3)
method='UPGMA'
my_tree_path=f"{out_dir}/lineage_tree_{clone_key}_{save_data_des}_{method}_support.txt"
methyltree.plotting.plot_tree_with_support(my_tree_path,leaf_name_map=None,figsize=(10,20))#leaf_name_map=lambda x: x.split('-')[-1],figsize=(20,50))
plt.savefig(f'{figure_path}/{data_des}_lineage_tree_support.pdf')
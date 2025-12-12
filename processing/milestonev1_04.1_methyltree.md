```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/methyltree.sif
source activate
ipython
#have to run interactive methyltree analysis with ipython
```

```python
#run in ipython
import sys
sys.path.append('/MethylTree') #load in methyltree module
import methyltree 
import pandas as pd
import scanpy as sc
import os
import numpy as np
from matplotlib import pyplot as plt
import argparse
import cospar as cs
import re
from collections import Counter
from ete3 import Tree, PhyloTree
from tqdm import tqdm

#cs.settings.set_figure_params()

parser = argparse.ArgumentParser()
parser.add_argument('-i', "--input",default="BCMDCIS41T")  #description='sample name input (matching prefix to methyltree output from amethyst processing)'
parser.add_argument('-c', "--cpu_cores",default="50")  #description='sample name input (matching prefix to methyltree output from amethyst processing)'

args = parser.parse_args()
in_dir= "/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/methyltree/"
sample_name=args.input
os.chdir(in_dir+"/"+sample_name)
in_dat=in_dir+"/"+sample_name+"/"+"methyltree."+sample_name+".methyltree_input.h5"

print('Read in H5 Data') # 
dat=pd.read_hdf(in_dat, key="data", mode='r')
print('Write out data as methyltree format.')
dat.rename(columns={'cell_id': 'index'}, inplace=True)
df_out=dat.pivot(index='index',columns='genomic_region_id',values='value')
#dat.to_csv("data.tsv.gz",sep='\t',compression='gzip')

print('Read in metadata')
metadat=pd.read_hdf(in_dat, key="metadata", mode='r')
metadat['cellid']=metadat['sample']
metadat['HQ']=metadat['HQ'].astype(bool)
metadat.to_csv("sample_sheet.tsv.gz",sep='\t',compression='gzip')

adata=sc.AnnData(df_out) 
#This adata has not annotation yet. We will use the df_sample from raw_data folder to annotate this object later

print('Set up output directories')
#set directories
out_dir=os.getcwd()+"/methyltree_data"
save_data_des="methyltree"
data_des=sample_name
data_path=os.getcwd()
figure_path=os.getcwd()+"/methyltree_figure"

#make directories
os.makedirs(out_dir,exist_ok=True)
os.makedirs(figure_path,exist_ok=True)

df_sample=methyltree.hf.load_sample_info(data_path)
df_sample.head()
clone_key='cnv_clonename'
adata.obs=df_sample

#compute cell similarities
similarity_method="correlation"
raw_similarity_key = f"X_similarity_{similarity_method}_raw"
X_similarity, shared_site_matrix = methyltree.similarity.compute_similarity_matrix(adata.X,method=similarity_method,)
adata.obsm[raw_similarity_key] = X_similarity
adata.obsm["X_shared_sites"] = shared_site_matrix
adata.obsm["X_similarity"] = adata.obsm[raw_similarity_key].copy()

#removing optimize tree allows to run
adata_final,stat_out=methyltree.analysis.comprehensive_lineage_analysis(
                                                        out_dir,
                                                        data_path,
                                                        save_data_des,
                                                        clone_key,
                                                        adata_orig=adata,
                                                        clone_color_dict=None,
                                                        heatmap_additional_key_list=['celltype','cnv_clonename'],
                                                        #use existing similarity
                                                        compute_similarity=False, 
                                                        update_sample_info=False,
                                                        # correct the correlation bias
                                                        similarity_correction=True,
                                                        # remove cell-type specific DNA methylation signals
                                                        remove_celltype_signal=True,cell_type_key='celltype',
                                                        # optimize the lineage tree
                                                        optimize_tree=False,optimize_background_cutoff=0.8,
                                                        # better heatmap visualization
                                                        heatmap_vmax_percentile=99,heatmap_vmin_percentile=40,heatmap_figsize=(10, 9.5),heatmap_show_label=True,
                                                        heatmap_show_legend=True,heatmap_fontsize=2,
                                                        # related to saving figures
                                                        fig_dir=figure_path,data_des=data_des,
                                                        perform_clone_inference=True,clone_inference_threshold=0.8,clone_inference_print=True,
                                                        perform_memory_analysis=True, save_adata=True,perform_depth_analysis=True,)

#optimize lineage tree
background_cutoff=0.8
method='UPGMA'
my_tree_path=f"{out_dir}/lineage_tree_{clone_key}_{save_data_des}_{method}_support.txt"

with open(my_tree_path, "r") as f:
        my_tree = Tree(f.read(), support=True)


t = PhyloTree(my_tree_path)

from ete3 import ClusterTree


# To obtain all the evolutionary events involving a given leaf node we
# use get_my_evol_events method for a given cell
matches = t.search_nodes(name="119-cancer-BCMDCIS41T_c5")
human_seq = matches[0]
# Obtains its evolutionary events
events = human_seq.get_my_evol_events()

# Alternatively, you can scan the whole tree topology
events = t.get_descendant_evol_events()


#https://etetoolkit.org/docs/latest/tutorial/tutorial_phylogeny.html#detecting-evolutionary-events



```
for leaf in my_tree.iter_leaves():
    leaf.name=leaf.name.split('-')[0]

X_similarity = methyltree.similarity.rescale_similarity_matrix(adata_final.obsm["X_similarity"].copy())
order_x = np.array(my_tree.get_leaf_names()).astype(int)

## determine clone size scale
for i in range(X_similarity.shape[0]):
    X_similarity[i, i] = 0
high_idx = np.nonzero(X_similarity[order_x][:, order_x] > background_cutoff)
temp = abs(high_idx[0] - high_idx[1])
clone_size_scale = np.median(temp)
print(f"Clone size scale: {clone_size_scale}")


## define the weight matrix
for i in range(X_similarity.shape[0]):
    X_similarity[i, i] = 1

factor = 2
size = X_similarity.shape[0]
weight_matrix = np.zeros((size, size))
for i in range(weight_matrix.shape[0]):
    for j in range(weight_matrix.shape[0]):
        if abs(i - j) > factor * clone_size_scale:
            weight_matrix[i, j] = (
                np.exp(-factor) * (factor * clone_size_scale) / abs(i - j)
            )  # np.exp(-factor)*(weight_matrix.shape[0]-abs(i-j)+factor*clone_size_scale)/
        else:
            weight_matrix[i, j] = np.exp(-abs(i - j) / clone_size_scale)
        # weight_matrix[i,j]=np.exp(-abs(i-j)/clone_size_scale)
## perform optimization
current_score = np.sum(X_similarity[order_x][:, order_x] * weight_matrix)
result = []
iterations=10000

for j in tqdm(range(iterations)):
    # my_tree_new=omic_ana.random_subtree_pruning_and_regrafting(my_tree)
    my_tree_new = methyltree.lineage.random_subtree_swapping(my_tree)
    order_x_new = np.array(my_tree_new.get_leaf_names()).astype(int)
    next_score = np.sum(X_similarity[order_x_new][:, order_x_new] * weight_matrix)
    # print(next_score-current_score)
    if next_score > current_score:
        my_tree = my_tree_new
        current_score = next_score
        ordered_clone_array = adata.obs[clone_key].astype(str)[order_x_new]
        unique_clone_set = list(set(ordered_clone_array))
        for x0 in unique_clone_set:
            if np.sum(ordered_clone_array == x0) == 1:
                ordered_clone_array[ordered_clone_array == x0] = "nan"
        df_accuracy = methyltree.metric.lineage_accuracy_from_leaf_order(ordered_clone_array)
        mean_accuracy = df_accuracy[df_accuracy["clone_size"] > 1].iloc[:, -4:].mean() #this is the line with the error, just sliced to last 4 col to fix
        mean_accuracy = mean_accuracy[["accuracy", "entropy", "continuity", "wassertein"]]
        print(f"iter-{j}; score {current_score}; continuity {mean_accuracy[2]}")
        result.append(
            [
                current_score,
                mean_accuracy[0],
                mean_accuracy[1],
                mean_accuracy[2],
                mean_accuracy[3],
            ]
        )

## plot the results
plt.clf()
cs.settings.set_figure_params()
df_tmp = pd.DataFrame(result)
fig, ax = plt.subplots()
ax.scatter(df_tmp[0], df_tmp[1], label="accuracy")
ax.scatter(df_tmp[0], df_tmp[2], label="entropy")
ax.scatter(df_tmp[0], df_tmp[3], label="continuity")
ax.scatter(df_tmp[0], df_tmp[4], label="wassertein")
plt.xlabel("Optimization score")
plt.ylabel("Accuracy")
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.savefig(f"{figure_path}/lineage_optimized.pdf")


vmin = np.percentile(X_similarity.flatten(), 50)
vmax = np.max(np.triu(X_similarity, k=1))
methyltree.plotting.plot_similarity_heatmap_with_clone_label(
    X_similarity[order_x_new][:, order_x_new],
    label=adata.obs[clone_key].astype(str)[order_x_new],
    vmin=vmin,
    vmax=vmax,
    color_bar=False,
    save_name=f"{figure_path}/lineage_optimized.heatmap.pdf",
    cell_index=None,
    show_legend=True,)

#save the optimized tree
my_tree.write(format=1,outfile=f"{out_dir}/lineage_tree_{clone_key}_{save_data_des}_{method}.txt",)

methyltree.lineage.bootstrap_lineage_tree(
    adata_final,
    out_dir,
    save_data_des,
    clone_key,
    cell_type_key='celltype',
    selected_idx=None,
    num_iterations=100,
    sampled_replacement=False,
    sample_fraction=0.8,
    similarity_correction=True,
    similarity_normalize=True,
    remove_celltype_signal=True,
    fig_dir=figure_path,
    recompute=False,
    heatmap_ordering_method="UPGMA",
    cores=64,)


#visualize tree with support values
cs.settings.set_figure_params(fontsize=13)
my_tree_path=f"{out_dir}/lineage_tree_{clone_key}_{save_data_des}_{method}_support.txt"
methyltree.plotting.plot_tree_with_support(my_tree_path,leaf_name_map=lambda x: x.split('-')[-1],figsize=(20,50))
plt.savefig(f'{figure_path}/lineage_tree_support.pdf')


#Predict clones based on support value cutoff
df_predict=methyltree.clone.identify_putative_clones_from_trees_with_support(
    my_tree_path, 
    adata_final.obsm['X_similarity'],
    print_clone=True,
)

#determine clone size
my_tree_path=f"{out_dir}/lineage_tree_{clone_key}_{save_data_des}_{method}_support.txt"
df_state_clone = pd.read_csv(f"{out_dir}/state_clone_info_{clone_key}_{save_data_des}.csv")
rename_leaf_map = dict(zip(df_state_clone["state_info"].astype(str), df_state_clone["clone_id"]))
methyltree.clone.infer_clone_number_across_various_threshold(adata_final, out_dir, clone_key, save_data_des, figure_path)

#Predict clones based on support value cutoff
import seaborn as sns


for i in range(X_similarity.shape[0]):
        X_similarity[i, i] = np.nan  # to avoid counting the diagonal terms

# Obtain randomized similarity
X_similarity_rand = X_similarity.copy()
rand_dist = []
cell_orders = np.arange(X_similarity.shape[0])
for __ in tqdm(range(100)):
    np.random.shuffle(cell_orders)
    methyltree.lineage.traverse_tree_to_update_dist(
        my_tree, X_similarity_rand[cell_orders][:, cell_orders], percentile=50
    )
    for tmp_tree in my_tree.traverse():
        rand_dist.append(tmp_tree.dist)

# Obtain true similarity
methyltree.lineage.traverse_tree_to_update_dist(my_tree, X_similarity, percentile=50)
true_dist = []
for tmp_tree in my_tree.traverse():
    true_dist.append(tmp_tree.dist)

rand_dist = np.array(rand_dist)
rand_dist = rand_dist[~np.isnan(rand_dist)]
true_dist = np.array(true_dist)
true_dist = true_dist[~np.isnan(true_dist)]

plt.subplots(1, 1, figsize=(5.5, 3.5))
color_random = "#a6bddb"
color_data = "#fdbb84"

all_data = list(true_dist) + list(rand_dist)
bins = np.linspace(np.min(all_data), np.max(all_data), 50)

ax = sns.histplot(
    data=true_dist,
    label="Observed",
    bins=bins,
    stat="probability",
    color=color_data,
    alpha=0.5,
)

ax = sns.histplot(
    data=rand_dist,
    label="Random",
    bins=bins,
    stat="probability",
    color=color_random,
)

# ax.legend()
plt.legend(loc=[1.05, 0.4])
# ax.set_xlabel("Sister-cell distance")
plt.xlabel("Tree weight")
plt.ylabel("Normalized frequency")
# plt.yscale('log')
plt.tight_layout()
plt.savefig(f"{figure_path}/lineage_clone_signal.pdf")
true_dist
df_predict=methyltree.clone.identify_putative_clones(
    adata_final.obsm['X_similarity'],
    my_tree, 
    plot_signal_noise=True,
    weight_threshold=np.percentile(true_dist, 10) , #get lower percentile of true distribution
    print_clone=True,
)

Counter(df_predict['predicted_clone'])


# #Predict clones based on support value cutoff
# df_predict2=methyltree.clone.identify_putative_clones_from_trees_with_support(
#     my_tree_path, 
#     adata_final.obsm['X_similarity'],
#     support_threshold=np.percentile(true_dist, 10), 
#     print_clone=True,
#     similarity_percentile=53,
#     method='high', 
# )


df_predict['cell_id_new']=df_predict['cell_id']
df_predict.to_csv(f'{out_dir}/predicted_clone_size.csv')
adata_final.obs["inferred_clone"] = df_predict.sort_values("cell_id_new")["predicted_clone"].to_numpy()

singleton_fraction=np.mean(df_predict['clone_size']==1)
df_freq=methyltree.clone.estimate_clone_number(adata_final, clone_key="inferred_clone")


plt.subplots(figsize=(3.5,3.5))
plt.bar(df_freq['clone_size'],df_freq['frequency'],color='#7fc97f')
plt.title(
    f"Observed cells: {adata_final.shape[0]}"+ "\n"
    +  f"Singleton cell fraction: {singleton_fraction:.2f}"
)
plt.xlabel("Methylation-defined clone size")
#plt.xticks([1,2,3,4,6,8,10,12],[1,2,3,4,6,8,10,12])
plt.ylabel("Count")
plt.tight_layout()
plt.savefig(f'{figure_path}/methylation_clone_size_distribution.pdf')

pd.DataFrame(adata_final.obs['inferred_clone']).reset_index().rename(columns={'index':'sample'}).to_csv(f'{out_dir}/{save_data_des}_predicted_clone_size_with_sample_ID.csv')

size_cutoff=1
true_labels=adata_final.obs[clone_key].to_numpy()[adata_final.uns['order_x']]
labels_pred=adata_final.obs['inferred_clone'].to_numpy()[adata_final.uns['order_x']]

size_cutoff=1
true_labels=adata_final.obs[clone_key].to_numpy()[adata_final.uns['order_x']]
labels_pred=adata_final.obs['inferred_clone'].to_numpy()[adata_final.uns['order_x']]
df=pd.DataFrame({'true_label':true_labels,'predict_label':labels_pred}).reset_index()
valid_true_labels=df.groupby('true_label').agg({'index':'count'}).rename(columns={'index':'count'}).query(f'count>{size_cutoff}').index
valid_pred_labels=df.groupby('predict_label').agg({'index':'count'}).rename(columns={'index':'count'}).query(f'count>{size_cutoff}').index
df_sub=df[df['true_label'].isin(valid_true_labels) & df['predict_label'].isin(valid_pred_labels)]
df_sub_2=df_sub[(df_sub['true_label']!='nan') & (df_sub['predict_label']!='nan')]

frac_predicted_large_clones=len(df[df['predict_label'].isin(valid_pred_labels)])/len(df)
print(f'Fraction of cells in the non-singleton predicted clones: {frac_predicted_large_clones:.2f}')

#update with inferred clones
df_clone_all=pd.read_csv(f'{out_dir}/{save_data_des}_predicted_clone_size_with_sample_ID.csv',index_col=0)
df_clone_all=df_clone_all.set_index('sample')
df_sample=methyltree.metadata.load_sample_info(data_path)
#del df_sample['inferred_clone']
df_sample['inferred_clone']=df_clone_all['inferred_clone']
df_sample=df_sample.reset_index()
methyltree.metadata.backup_and_save_sample_info(df_sample,data_path)
#
df_sample['inferred_clone'].unique()

import sklearn
ARI=sklearn.metrics.adjusted_rand_score(df_sub_2['true_label'], df_sub_2['predict_label'])
print(f'Adjusted rank index: {ARI:.2f}')

df_predict=pd.read_csv(f'{out_dir}/predicted_clone_size.csv',index_col=0)
df_tmp=df_predict.filter(['predicted_clone','clone_size']).drop_duplicates().sort_values('predicted_clone')
#groupby('predicted_clone').agg({'cell_id':'count'}).reset_index()
clone_size_array=df_tmp['clone_size'].to_list()


import random
sample_name='test'
sample_cell_N_dict={'test':adata_final.shape[0]}
singleton_fraction_dict={'test':singleton_fraction} # using the singleton_cell fraction from the previous plot

df_predict=pd.read_csv(f'{out_dir}/predicted_clone_size.csv',index_col=0)
df_tmp=df_predict.filter(['predicted_clone','clone_size']).drop_duplicates().sort_values('predicted_clone')
#groupby('predicted_clone').agg({'cell_id':'count'}).reset_index()
clone_size_array=df_tmp['clone_size'].to_list()


clone_size=10
sampled_cell_N=sample_cell_N_dict[sample_name]
predicted_clone_list=[]
singleton_list=[]

N_c_array=np.arange(20,400,20)

for N_c in N_c_array:
    singleton_tmp_list=[]
    predicted_clone_tmp_list=[]
    for _ in range(100):
        data_list=[]
        clone_size_array_tmp=clone_size*np.random.choice(clone_size_array,N_c,replace=True).astype(int)
        for i in range(0,N_c):
            tmp=np.zeros(clone_size_array_tmp[i]).astype(int)+i
            data_list.append(tmp.astype(str))
        final_data=np.hstack(data_list)
        sampled_data=np.random.choice(final_data, sampled_cell_N, replace=True)
        df=pd.DataFrame({'clone_id':sampled_data}).reset_index().groupby('clone_id').agg({'index':'count'}).reset_index().rename(columns={'index':'clone_size'}).sort_values('clone_size',ascending=False)
        singleton_fraction_tmp=np.sum(df['clone_size']==1)/df['clone_size'].sum()
        predicted_clone_N=sampled_cell_N/(2*(1-singleton_fraction_tmp)**1.5)
        singleton_tmp_list.append(singleton_fraction_tmp)
        predicted_clone_tmp_list.append(predicted_clone_N)
    predicted_clone_list.append(predicted_clone_tmp_list)
    singleton_list.append(singleton_tmp_list)

np.savez(f'{out_dir}/saved_sampling_clone_data_sampled_cell_N{sampled_cell_N}.npz',N_c_array=N_c_array,singleton_list=np.array(singleton_list))

cs.settings.set_figure_params(format='pdf',figsize=[4,3.5],fontsize=14)

sampled_cell_N=sample_cell_N_dict[sample_name]
file=np.load(f'{out_dir}/saved_sampling_clone_data_sampled_cell_N{sampled_cell_N}.npz')
N_c_array=file['N_c_array']
singleton_list=file['singleton_list']

from sklearn.linear_model import LinearRegression

# Create training data
target_singleton=singleton_fraction_dict[sample_name]
df_tmp=pd.DataFrame({'singleton':list(np.mean(singleton_list,1))+[target_singleton],'clone_N':list(N_c_array)+[np.nan]}).sort_values('singleton')
target_id=np.nonzero(np.array(df_tmp['singleton']==target_singleton))[0][0]
if target_id+2 > df_tmp.shape[0]:
    raise ValueError(f"target_id: {target_id}; df_tmp.shape[0]: {df_tmp.shape[0]}. Please Increased upper limit of N_c_array")
df_train=df_tmp.iloc[[target_id-2,target_id-1,target_id+1,target_id+2]]
X_train=[[x] for x in df_train['singleton'].to_numpy()]
y_train=df_train['clone_N'].to_list()

# Create and train the linear regression model
model = LinearRegression()
model.fit(X_train, y_train)

# X_new is the observed singleton fraction
singleton_frac=singleton_fraction_dict[sample_name]
X_new = [[singleton_frac]]

    
# Predict the output value
y_pred = model.predict(X_new)
estimated_clone_N=int(y_pred[0])

print(f"Sample name: {sample_name}; cell number: {sampled_cell_N}; singleton_frac: {X_new[0][0]}; Predicted output value:", y_pred[0])


plt.subplots(figsize=(4,3.5))
#plt.plot(N_c_array,np.array(singleton_list),'.',color='grey')
# rcParams["axes.spines.right"] = False
# rcParams["axes.spines.top"] = False
error_array=np.array([np.std(x) for x in np.array(singleton_list)])
plt.errorbar(N_c_array,np.mean(singleton_list,1),
             yerr=error_array,marker='*',color='r',
            ecolor='grey',elinewidth=1,linestyle ='')

plt.plot(N_c_array,np.mean(singleton_list,1),'-*r')

plt.xlim([0,estimated_clone_N+50])
plt.ylim([0,1])
plt.xlabel('Total clone number')
plt.ylabel('Singleton cell fraction')

plt.plot([0,estimated_clone_N+20],[singleton_frac,singleton_frac],'--b',linewidth=1)
plt.plot([estimated_clone_N,estimated_clone_N],[0,singleton_frac+0.05],'--b',linewidth=1)
plt.title(f'{sample_name}: sample {sampled_cell_N} cells')
plt.tight_layout()
plt.savefig(f'{figure_path}/{sample_name}_clone_number_prediction.pdf')

df_freq=methyltree.clone.estimate_clone_number(adata_final, clone_key="inferred_clones")

methyltree.plotting.plot_similarity_heatmap_with_multiple_colorbars(
    adata_final,additional_key_list=['celltype','cnv_clonename','inferred_clones'],
    save_name=sample_name+".similarity_heatmap.pdf",
    data_des=sample_name,fig_dir=figure_path,
    heatmap_vmax_percentile=99.9,heatmap_vmin_percentile=60,)


methyltree.analysis.coarse_grain_analysis(
    adata_final,
    clone_key='celltype',
    normalized_X_similarity=True)

methyltree.analysis.coarse_grain_analysis(
    adata_final,
    clone_key='cnv_clonename',
    normalized_X_similarity=True)

X_coupling, mega_cluster_list, mega_cluster_list_2= methyltree.analysis.coarse_grained_coupling(
    adata_final,
    key_1='inferred_clones',key_2='cnv_clonename',)

adata_new=methyltree.analysis.embedding_analysis(adata_final,
    clone_key='inferred_clones',
    UMAP_parameters=[10,10,0.3],
    clone_prefix=None,figure_path=figure_path,
    title='',
    fontsize=13,marker_size=50)

cs.settings.set_figure_params(fontsize=3)
method='UPGMA'
my_tree_path=f"{out_dir}/lineage_tree_{clone_key}_{save_data_des}_{method}_support.txt"
methyltree.plotting.plot_tree_with_support(my_tree_path,leaf_name_map=None,figsize=(10,20))#leaf_name_map=lambda x: x.split('-')[-1],figsize=(20,50))
plt.savefig(f'{figure_path}/{data_des}_lineage_tree_support.pdf')

```
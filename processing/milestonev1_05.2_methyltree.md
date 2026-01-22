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
import scipy
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
                                                        compute_similarity=False, #use existing similarity
                                                        update_sample_info=False,
                                                        similarity_correction=True, # correct the correlation bias
                                                        remove_celltype_signal=True,cell_type_key='celltype', # remove cell-type specific DNA methylation signals
                                                        optimize_tree=False,optimize_background_cutoff=0.2, # optimize the lineage tree
                                                        # better heatmap visualization
                                                        heatmap_vmax_percentile=100,heatmap_vmin_percentile=70,heatmap_figsize=(10, 9.5),heatmap_show_label=True,
                                                        heatmap_show_legend=True,heatmap_fontsize=2,
                                                        # related to saving figures
                                                        fig_dir=figure_path,data_des=data_des,
                                                        perform_clone_inference=True,clone_inference_threshold=0.2,clone_inference_print=True,
                                                        perform_memory_analysis=True, save_adata=True,perform_depth_analysis=True,)

#Compute support value of the tree
methyltree.lineage.bootstrap_lineage_tree(adata_final,out_dir,
    save_data_des,
    clone_key=clone_key,
    num_iterations=20,  
    sample_fraction=0.8,  
    similarity_correction=True,
    remove_celltype_signal=True,cell_type_key='celltype',
    heatmap_ordering_method='UPGMA',
    recompute=False,
    cores=100)

#Visualize the lineage tree with support values
cs.settings.set_figure_params(fontsize=13)
method='UPGMA'
my_tree_path=f"{out_dir}/lineage_tree_{clone_key}_{save_data_des}_{method}_support.txt"
methyltree.plotting.plot_tree_with_support(my_tree_path,leaf_name_map=None,figsize=(10,30))#leaf_name_map=lambda x: x.split('-')[-1],figsize=(20,50))
plt.savefig(f'{figure_path}/{sample_name}_lineage_tree_support.pdf')



#take celltype corrected similarity and cluster by scipy

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster, cut_tree
from sklearn.metrics import silhouette_score, silhouette_samples


https://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_silhouette_analysis.html

range_n_clusters = [2, 3, 4, 5, 6]

X=1-pd.DataFrame(adata_final.obsm['X_similarity'],
            index=adata_final.obs.cellid,
            columns=adata_final.obs.cellid)
for n_clusters in range_n_clusters:
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(18, 7)
    ax1.set_xlim([-0.1, 1])
    ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])
    clusterer = KMeans(n_clusters=n_clusters, random_state=10)
    cluster_labels = clusterer.fit_predict(X)
    silhouette_avg = silhouette_score(X, cluster_labels)
    print(
        "For n_clusters =",
        n_clusters,
        "The average silhouette_score is :",
        silhouette_avg,
    )
    sample_silhouette_values = silhouette_samples(X, cluster_labels)

    y_lower = 10
    for i in range(n_clusters):
        ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]
        ith_cluster_silhouette_values.sort()
        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i
        color = cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(
            np.arange(y_lower, y_upper),
            0,
            ith_cluster_silhouette_values,
            facecolor=color,
            edgecolor=color,
            alpha=0.7,
        )
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
        y_lower = y_upper + 10  # 10 for the 0 samples
    ax1.set_title("The silhouette plot for the various clusters.")
    ax1.set_xlabel("The silhouette coefficient values")
    ax1.set_ylabel("Cluster label")
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")
    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
    colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
    ax2.scatter(
        X[:, 0], X[:, 1], marker=".", s=30, lw=0, alpha=0.7, c=colors, edgecolor="k"
    )
    centers = clusterer.cluster_centers_
    ax2.scatter(
        centers[:, 0],
        centers[:, 1],
        marker="o",
        c="white",
        alpha=1,
        s=200,
        edgecolor="k",
    )
    for i, c in enumerate(centers):
        ax2.scatter(c[0], c[1], marker="$%d$" % i, alpha=1, s=50, edgecolor="k")
    ax2.set_title("The visualization of the clustered data.")
    ax2.set_xlabel("Feature space for the 1st feature")
    ax2.set_ylabel("Feature space for the 2nd feature")
    plt.suptitle(
        "Silhouette analysis for KMeans clustering on sample data with n_clusters = %d"
        % n_clusters,
        fontsize=14,
        fontweight="bold",
    )

plt.show()

dist_mat
#dist_mat=scipy.spatial.distance.squareform(dist_mat)
linkage_matrix = linkage(dist_mat, metric="euclidean",optimal_ordering=True)
tree_cut = cut_tree(linkage_matrix,n_clusters=8)
dendrogram(linkage_matrix)
plt.title("test")
nodes = fcluster(linkage_matrix, 4, criterion="maxclust")
silhouette_samples(dist_mat, tree_cut)
silhouette_score(dist_mat + dist_mat.T , nodes, metric='euclidean')
plt.savefig(f'{figure_path}/scipytest_{sample_name}_lineage_tree_support.pdf')


scipy.cluster.hierarchy.cut_tree(my_tree,n_clusters=10)

#convert tree to ClusterNode
from ete3 import ClusterTree
with open(my_tree_path,"r") as f:
    my_tree = Tree(f.read())


mat_file_path=f"{out_dir}/sim_mat_{clone_key}_{save_data_des}_{method}_support.txt"
out_df=pd.DataFrame(adata_final.obsm['X_similarity'],
            index=adata_final.obs.cellid,
            columns=adata_final.obs.cellid)

out_df.to_csv(mat_file_path,sep="\t",header=True,index=True)

pd.read_csv(mat_file_path,sep="\t")

t = ClusterTree(my_tree_path, text_array=out_df)










#Predict clones based on support value cutoff
df_predict=methyltree.clone.identify_putative_clones_from_trees_with_support(
    my_tree_path, 
    adata_final.obsm['X_similarity'],
    support_threshold=0.8, 
    print_clone=True,
    similarity_percentile=50,
    method='low', 
)

df_predict['cell_id_new']=np.array(df_predict['cell_id'].apply(lambda x:x.split('-')[0])).astype(int)
df_predict.to_csv(f'{out_dir}/predicted_clone_size.csv')
adata_final.obs["inferred_clone"] = df_predict.sort_values("cell_id_new")[
    "predicted_clone"
].to_numpy()

singleton_fraction=np.mean(df_predict['clone_size']==1)
print(f'Singleton fraction: {singleton_fraction:.2f}')

df_freq=methyltree.clone.estimate_clone_number(adata_final, clone_key="inferred_clones")

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

clone_key='cnv_clonename'
methyltree.clone.infer_clone_number_across_various_threshold(
    adata_final,
    out_dir,clone_key,
    save_data_des,
    figure_path)

pd.DataFrame(adata_final.obs['inferred_clone']).reset_index().rename(columns={'index':'sample'}).to_csv(f'{out_dir}/{save_data_des}_predicted_clone_size_with_sample_ID.csv')
#methy_ana.accuracy_after_removing_leaves_with_low_support(df_predict,adata,clone_key)


clone_key='inferred_clone'

adata_new=methyltree.analysis.embedding_analysis(adata_final,clone_key,
                            UMAP_parameters=[10,10,0.3],
                            clone_prefix=None,figure_path=figure_path,
                            title='',
                            fontsize=13,marker_size=50)
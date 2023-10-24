# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: python (assess)
#     language: python
#     name: assess
# ---

# +
import os
import re
import sys 
import json
import glob
import scib
import scipy

import inspect

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import homogeneity_score


# -

def get_latent_files(input_path):
    # Load Data
    methods = [x for x in next(os.walk(input_path))[1] if x[0:3] == 'run']
    latent_paths = []
    for method in methods:
        method_path = os.path.join(input_path, method)
        latent_paths.extend([glob.glob(method_path + '/*latent.csv')])
    latent_paths = list(filter(None, latent_paths))
    
    return(latent_paths)


def paired_metrics(
    adata,
    latent_paths,
    output_path,
    metadata_path,
    metadata_index='barcode',
    label = 'cell_type',
    batch = 'batch',
    cluster_method = 'louvain'
):
    metadata = pd.read_table(metadata_path, sep=',', index_col=metadata_index)
    # create metrics DataFrame
    metrics = pd.DataFrame()
    # calculate metrics
    for latent_path in latent_paths:
        print(latent_path)
        # Set Anndata
#         adata = sc.AnnData(latent, dtype='float64')
        adata.obs = metadata.loc[adata.obs_names.to_list()].astype('category')
        # get latent file prefix
        latent = pd.read_csv(latent_path[0], header=0, index_col=0)
        latent_file_name = re.split('/', latent_path[0])[-1]
        latent_prefix = re.split('-latent.csv', latent_file_name)[0]
        latent_prefix
        # add latent in obsm
        ## rearrangement latent index by adata's index
        latent = latent.loc[adata.obs_names.to_list()]
        adata.obsm[latent_prefix] = latent
        
        # Cluster
        if cluster_method == 'louvain':
            # add 2 obsp：细胞间的distances矩阵和connectivities矩阵
            sc.pp.neighbors(adata, use_rep=latent_prefix, n_neighbors=30) # neighborhood graph based on UMAP
            # 添加1个obs：louvain; 添加2个uns：neighbors字典和louvain字典
        #     getNClusters(adata, n_cluster=num_clusters) # 使用Louvain聚类得到要求的cluster数量
            sc.tl.louvain(adata) # Louvain 
        elif cluster_method == 'kmeans':
            kmeans = KMeans(n_clusters=num_clusters, random_state=2022).fit(adata.X) # 添加1个obs：kmeans
            adata.obs[cluster_method] = pd.Series(kmeans.labels_, index=adata.obs.index).astype('category')
        elif cluster_method == 'hc': # hierachical clustering
            hc = AgglomerativeClustering(n_clusters=num_clusters).fit(adata.X) # 添加1个obs：hc
            adata.obs[cluster_method] = pd.Series(hc.labels_,index=adata.obs.index).astype('category')
        
        # 0. latent dimension
        metrics.loc[latent_prefix+'-'+cluster_method, ['nCell']] = adata.obsm[latent_prefix].shape[0]
        metrics.loc[latent_prefix+'-'+cluster_method, ['nDimension']] = adata.obsm[latent_prefix].shape[1]
        # 1. adjusted rank index
        ari = adjusted_rand_score(adata.obs[label], adata.obs[cluster_method])
        metrics.loc[latent_prefix+'-'+cluster_method, ['ARI']] = [ari]
        # 2. adjusted mutual information
        ami = adjusted_mutual_info_score(adata.obs[label], adata.obs[cluster_method],average_method='arithmetic')
        metrics.loc[latent_prefix+'-'+cluster_method, ['AMI']] = [ami]
        # 3. kNN graph_connectivity
        # embedding output
        sc.pp.neighbors(adata, use_rep=latent_prefix)
        graph_connectivity = scib.metrics.graph_connectivity(adata, batch)
        metrics.loc[latent_prefix+'-'+cluster_method, ['graph_connectivity']] = [graph_connectivity]
        # 4. Batch ASW
        Batch_ASW = scib.metrics.silhouette_batch(adata, batch_key=batch, group_key=label, embed=latent_prefix, metric="euclidean", return_all=False, scale=True, verbose=False)
        metrics.loc[latent_prefix+'-'+cluster_method, ['Batch_ASW']] = [Batch_ASW]
        # 5. silhouette
        silhouette = scib.metrics.silhouette(adata, group_key=cluster_method, embed=latent_prefix, metric="euclidean", scale=True)
        metrics.loc[latent_prefix+'-'+cluster_method, ['silhouette']] = [silhouette]
        # 6. graph iLISI
        adata.obsm["X_emb"] = latent
        graph_iLISI = scib.metrics.ilisi_graph(adata, batch_key=batch, type_="embed", k0=90, subsample=None, scale=True, n_cores=1, verbose=False)
        metrics.loc[latent_prefix+'-'+cluster_method, ['graph_iLISI']] = [graph_iLISI]
        # 7. clisi_graph
        graph_clisi = scib.metrics.clisi_graph(adata, batch_key=batch, label_key=label, type_="embed", k0=90, subsample=None, scale=True, n_cores=1, verbose=False)
        metrics.loc[latent_prefix+'-'+cluster_method, ['graph_clisi']] = [graph_clisi]
        # 8. isolated_labels
        isolated_labels = scib.metrics.isolated_labels(adata, batch_key=batch, label_key=label, embed=latent_prefix, cluster=True, iso_threshold=None, return_all=False, verbose=False)
        metrics.loc[latent_prefix+'-'+cluster_method, ['isolated_labels']] = [isolated_labels]
    metrics.to_csv(output_path + '/metrics.csv')
    adata.uns['metrics'] = metrics
        
    return(adata)


# +
input_path = sys.argv[1]
metadata_path = sys.argv[2]
rep_path = sys.argv[3]
output_path = sys.argv[4]

# Load Data
rna = sc.read_h5ad(glob.glob(input_path + '/*-RNA-counts.h5ad')[0])

adt = sc.read_h5ad(glob.glob(input_path + '/*-ADT-counts.h5ad')[0])

adata = ad.AnnData(rna)
adata.obsm['RNA'] = rna.X.copy()
adata.obsm['ADT'] = adt.X.copy()

latent_paths = get_latent_files(rep_path)

adata = paired_metrics(adata, latent_paths, output_path, metadata_path)
del adata.raw
adata.write_h5ad(output_path + "/metrics.h5ad")
adata.uns['metrics'].to_csv(output_path + '/metrics.csv')
# adata.uns['metrics']
# -



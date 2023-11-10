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
import scglue
from scglue.typehint import RandomState
from scglue.utils import get_rs

import inspect

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

import scipy.spatial
from anndata import AnnData
from scipy.sparse.csgraph import connected_components


import sklearn.metrics
import sklearn.neighbors
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


# # FOSCTTM

from typing import Tuple
def foscttm(x: np.ndarray, y: np.ndarray,
            **kwargs) -> Tuple[np.ndarray, np.ndarray]:
    r"""
    Fraction of samples closer than true match (smaller is better)

    Parameters
    ----------
    x
        Coordinates for samples in modality X
    y
        Coordinates for samples in modality y
    **kwargs
        Additional keyword arguments are passed to
        :func:`scipy.spatial.distance_matrix`

    Returns
    -------
    foscttm_x, foscttm_y
        FOSCTTM for samples in modality X and Y, respectively

    Note
    ----
    Samples in modality X and Y should be paired and given in the same order
    """
    if x.shape != y.shape:
        raise ValueError("Shapes do not match!")
    d = scipy.spatial.distance_matrix(x, y, **kwargs)
    foscttm_x = (d < np.expand_dims(np.diag(d), axis=1)).mean(axis=1)
    foscttm_y = (d < np.expand_dims(np.diag(d), axis=0)).mean(axis=0)
    foscttm_mean = np.mean([foscttm_x.mean(), foscttm_y.mean()])
    return foscttm_mean


# # seurat_alignment_score

def seurat_alignment_score(
        x: np.ndarray, y: np.ndarray, neighbor_frac: float = 0.01,
        n_repeats: int = 4, random_state: RandomState = None, **kwargs
) -> float:
    r"""
    Seurat alignment score

    Parameters
    ----------
    x
        Coordinates
    y
        Batch labels
    neighbor_frac
        Nearest neighbor fraction
    n_repeats
        Number of subsampling repeats
    random_state
        Random state
    **kwargs
        Additional keyword arguments are passed to
        :class:`sklearn.neighbors.NearestNeighbors`

    Returns
    -------
    sas
        Seurat alignment score
    """
    rs = get_rs(random_state)
    idx_list = [np.where(y == u)[0] for u in np.unique(y)]
    min_size = min(idx.size for idx in idx_list)
    repeat_scores = []
    for _ in range(n_repeats):
        subsample_idx = np.concatenate([
            rs.choice(idx, min_size, replace=False)
            for idx in idx_list
        ])
        subsample_x = x[subsample_idx]
        subsample_y = y[subsample_idx]
        k = max(round(subsample_idx.size * neighbor_frac), 1)
        nn = sklearn.neighbors.NearestNeighbors(
            n_neighbors=k + 1, **kwargs
        ).fit(subsample_x)
        nni = nn.kneighbors(subsample_x, return_distance=False)
        same_y_hits = (
            subsample_y[nni[:, 1:]] == np.expand_dims(subsample_y, axis=1)
        ).sum(axis=1).mean()
        repeat_score = (k - same_y_hits) * len(idx_list) / (k * (len(idx_list) - 1))
        repeat_scores.append(min(repeat_score, 1))  # score may exceed 1, if same_y_hits is lower than expected by chance
    return np.mean(repeat_scores).item()


# # nearest cell barcode

from typing import Tuple
def nearest_cell_barcode(
    adata,
    rep_1,
    rep_2):
    if adata.obsm[rep_1].shape != adata.obsm[rep_2].shape:
        raise ValueError("Shapes do not match!")
    d = scipy.spatial.distance_matrix(adata.obsm[rep_1], adata.obsm[rep_2])
    dist = pd.DataFrame(d)
    dist.index = adata.obsm[rep_1].index
    dist.columns = adata.obsm[rep_2].index

    # Get Column names of minimum value in every row 
    min_col_in_row = dist.idxmin(axis = 1)
    nearest_atac_of_rna = min_col_in_row.index == min_col_in_row.values
    nearest_atac_of_rna_ratio = sum(nearest_atac_of_rna)/len(nearest_atac_of_rna)

    # Get row index label of minimum value in every column :
    min_row_in_col = dist.idxmin()
    nearest_rna_of_atac = min_row_in_col.index == min_row_in_col.values
    nearest_rna_of_atac_ratio = sum(nearest_rna_of_atac)/len(nearest_rna_of_atac)

    nearest_cell_barcode_mean = (nearest_atac_of_rna_ratio+nearest_rna_of_atac_ratio)/2
    
    nearest_cell_barcode_list = pd.Series([nearest_atac_of_rna_ratio, nearest_rna_of_atac_ratio, nearest_cell_barcode_mean], 
                                        index=['RNA', 'ATAC', 'Mean'])
    
    return nearest_cell_barcode_list


# # nearest cell celltype

from typing import Tuple
def nearest_cell_celltype(
    adata,
    rep_1,
    rep_2,
    label = 'cell_type'):
    if adata.obsm[rep_1].shape != adata.obsm[rep_2].shape:
        raise ValueError("Shapes do not match!")
    d = scipy.spatial.distance_matrix(adata.obsm[rep_1], adata.obsm[rep_2])
    dist = pd.DataFrame(d)
    dist.index = adata.obs[label]
    dist.columns = adata.obs[label]

    # Get Column names of minimum value in every row 
    min_col_in_row = dist.idxmin(axis = 1)
    nearest_atac_of_rna = min_col_in_row.index == min_col_in_row.values
    nearest_atac_of_rna_ratio = sum(nearest_atac_of_rna)/len(nearest_atac_of_rna)

    # Get row index label of minimum value in every column :
    min_row_in_col = dist.idxmin()
    nearest_rna_of_atac = min_row_in_col.index == min_row_in_col.values
    nearest_rna_of_atac_ratio = sum(nearest_rna_of_atac)/len(nearest_rna_of_atac)

    nearest_cell_cell_type_mean = (nearest_atac_of_rna_ratio+nearest_rna_of_atac_ratio)/2
    
    nearest_cell_cell_type_list = pd.Series([nearest_atac_of_rna_ratio, nearest_rna_of_atac_ratio, nearest_cell_cell_type_mean], 
                                            index=['RNA', 'ATAC', 'Mean'])
    return nearest_cell_cell_type_list



def cross_metrics(
    adata,
    latent_paths,
    output_path,
    label = 'cell_type',
    batch = 'batch',
    cluster_method = 'louvain'
):
    # create metrics DataFrame
    metrics = pd.DataFrame()
    
    # calculate metrics
    for latent_path in latent_paths:
        print(latent_path)
        if len(latent_path) != 2:
            raise ValueError("methods should have two latent!")

        rna_latent_path = [x for x in latent_path if 'RNA-latent' in x][0]
        atac_latent_path = [x for x in latent_path if 'ATAC-latent' in x][0]

        # get latent file prefix
        rna_latent = pd.read_csv(rna_latent_path, header=0, index_col=0)
        atac_latent = pd.read_csv(atac_latent_path, header=0, index_col=0)

        rna_latent_file_name = re.split('/', rna_latent_path)[-1]
        atac_latent_file_name = re.split('/', atac_latent_path)[-1]

        latent_prefix = re.split('-RNA-latent.csv', rna_latent_file_name)[0]
        rna_latent_prefix = re.split('-latent.csv', rna_latent_file_name)[0]
        atac_latent_prefix = re.split('-latent.csv', atac_latent_file_name)[0]

        # add latent in obsm
        ## rearrangement latent index by adata's index
        rna_latent = rna_latent.loc[adata.obs_names.to_list()]
        atac_latent = atac_latent.loc[adata.obs_names.to_list()]

        adata.obsm[rna_latent_prefix] = rna_latent
        adata.obsm[atac_latent_prefix] = atac_latent
        
        # 1. FOSCTTM
        metrics.loc[latent_prefix+'-'+cluster_method, ['FOSCTTM']] = [foscttm(adata.obsm[rna_latent_prefix], adata.obsm[atac_latent_prefix])]
        # 2. nearest_cell_barcode
        nearest_cell_barcode_list = nearest_cell_barcode(adata, rep_1=rna_latent_prefix, rep_2=atac_latent_prefix)
        metrics.loc[latent_prefix+'-'+cluster_method, ['nearest_cell_barcode' + '-' + 'RNA']] = [nearest_cell_barcode_list['RNA']]
        metrics.loc[latent_prefix+'-'+cluster_method, ['nearest_cell_barcode' + '-' + 'ATAC']] = [nearest_cell_barcode_list['ATAC']]
        metrics.loc[latent_prefix+'-'+cluster_method, ['nearest_cell_barcode']] = [nearest_cell_barcode_list['Mean']]
        # 11. nearest_cell_celltype
        nearest_cell_celltype_list = nearest_cell_celltype(adata, rep_1=rna_latent_prefix, rep_2=atac_latent_prefix, label='cell_type')
        metrics.loc[latent_prefix+'-'+cluster_method, ['nearest_cell_celltype' + '-' + 'RNA']] = [nearest_cell_celltype_list['RNA']]
        metrics.loc[latent_prefix+'-'+cluster_method, ['nearest_cell_celltype' + '-' + 'ATAC']] = [nearest_cell_celltype_list['ATAC']]
        metrics.loc[latent_prefix+'-'+cluster_method, ['nearest_cell_celltype']] = [nearest_cell_celltype_list['Mean']]
    metrics.to_csv(output_path + '/metrics.csv')
    adata.uns['metrics'] = metrics
        
    return(adata)


# +
input_path = sys.argv[1]
metadata_path = sys.argv[2]

rep_path = sys.argv[3]
output_path = sys.argv[4]

# +
# Load Data
rna = sc.read_h5ad(glob.glob(input_path + '/*-RNA-counts.h5ad')[0])

atac = sc.read_h5ad(glob.glob(input_path + '/*-ATAC-peaks.h5ad')[0])

adata = ad.AnnData(rna)
adata.obsm['RNA'] = rna.X.copy()
adata.obsm['ATAC'] = atac.X.copy()

metadata = pd.read_table(metadata_path, sep=',', index_col='barcode')
adata.obs = metadata.loc[adata.obs_names.to_list()].astype('category')
# -


latent_paths = get_latent_files(rep_path)
latent_paths
adata = cross_metrics(adata, latent_paths, output_path)
# del adata.raw
adata.write_h5ad(output_path + "/metrics.h5ad")
adata.uns['metrics'].to_csv(output_path + '/metrics.csv')
# adata.uns['metrics']

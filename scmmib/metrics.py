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
import gc
import datatable as dt

from typing import Tuple
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.cluster import adjusted_rand_score, adjusted_mutual_info_score, homogeneity_score,\
    silhouette_samples, silhouette_score
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix
from sklearn.metrics import precision_recall_curve, auc, roc_auc_score
from scipy.sparse.csgraph import connected_components

from scmmib.knn_smooth import knn_smoothing


def foscttm(x: np.ndarray, y: np.ndarray,
            **kwargs) -> Tuple[np.ndarray, np.ndarray]:
    if x.shape != y.shape:
        raise ValueError("Shapes do not match!")
    d = scipy.spatial.distance_matrix(x, y, **kwargs)
    foscttm_x = (d < np.expand_dims(np.diag(d), axis=1)).mean(axis=1)
    foscttm_y = (d < np.expand_dims(np.diag(d), axis=0)).mean(axis=0)
    foscttm_mean = np.mean([foscttm_x.mean(), foscttm_y.mean()])
    return foscttm_mean

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

    min_col_in_row = dist.idxmin(axis = 1)
    nearest_mod2_of_rna = min_col_in_row.index == min_col_in_row.values
    nearest_mod2_of_rna_ratio = sum(nearest_mod2_of_rna)/len(nearest_mod2_of_rna)

    min_row_in_col = dist.idxmin()
    nearest_rna_of_mod2 = min_row_in_col.index == min_row_in_col.values
    nearest_rna_of_mod2_ratio = sum(nearest_rna_of_mod2)/len(nearest_rna_of_mod2)

    nearest_cell_barcode_mean = (nearest_mod2_of_rna_ratio+nearest_rna_of_mod2_ratio)/2
    
    nearest_cell_barcode_list = pd.Series([nearest_mod2_of_rna_ratio, nearest_rna_of_mod2_ratio, nearest_cell_barcode_mean], 
                                        index=['RNA', rep_2 , 'Mean'])
    
    return nearest_cell_barcode_list


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

    min_col_in_row = dist.idxmin(axis = 1)
    nearest_mod2_of_rna = min_col_in_row.index == min_col_in_row.values
    nearest_mod2_of_rna_ratio = sum(nearest_mod2_of_rna)/len(nearest_mod2_of_rna)

    min_row_in_col = dist.idxmin()
    nearest_rna_of_mod2 = min_row_in_col.index == min_row_in_col.values
    nearest_rna_of_mod2_ratio = sum(nearest_rna_of_mod2)/len(nearest_rna_of_mod2)

    nearest_cell_cell_type_mean = (nearest_mod2_of_rna_ratio+nearest_rna_of_mod2_ratio)/2
    
    nearest_cell_cell_type_list = pd.Series([nearest_mod2_of_rna_ratio, nearest_rna_of_mod2_ratio, nearest_cell_cell_type_mean], 
                                            index=['RNA', rep_2, 'Mean'])
    return nearest_cell_cell_type_list


def scmmib_gc(adata, label_key):
    """Graph Connectivity

    Quantify the connectivity of the subgraph per cell type label.
    The final score is the average for all cell type labels :math:`C`, according to the equation:

    .. math::

        GC = \\frac {1} {|C|} \\sum_{c \\in C} \\frac {|{LCC(subgraph_c)}|} {|c|}

    where :math:`|LCC(subgraph_c)|` stands for all cells in the largest connected component and :math:`|c|` stands for all cells of
    cell type :math:`c`.

    :param adata: adata with computed neighborhood graph
    :param label_key: name in adata.obs containing the cell identity labels

    Revised from scib.metrics.graph_connectivity, remove KeyError check
    """
#     if "neighbors" not in adata.uns:
#         raise KeyError(
#             "Please compute the neighborhood graph before running this function!"
#         )

    clust_res = []

    for label in adata.obs[label_key].cat.categories:
        if not np.sum(adata.obs[label_key].isin([label])): # avoid 0 cells in certain batch
            continue

        adata_sub = adata[adata.obs[label_key].isin([label])]
        _, labels = connected_components(
            adata_sub.obsp["connectivities"], connection="strong"
        )
        tab = pd.value_counts(labels)
        clust_res.append(tab.max() / sum(tab))

    return np.mean(clust_res)


def cells_correlation(ndarray_1, 
                      ndarray_2, 
                      metadata,
                      celltype_key='cell_type'):
    # corr_results = pd.Series(dtype='float64')
    cell_types_correlation = None
    if celltype_key in metadata.columns:
        cell_types = metadata[celltype_key].unique()
        tmp = []
        for cell_type in cell_types:

            cells_index = metadata[celltype_key] == cell_type
            N = sum(cells_index)

            ndarray_1_mean = np.mean(ndarray_1[cells_index,], axis=0)
            ndarray_2_mean = np.mean(ndarray_2[cells_index,], axis=0)
            tmp.append(np.corrcoef(list(ndarray_1_mean), list(ndarray_2_mean))[0, 1])
        cell_types_correlation = np.mean(tmp)

    corr_results = pd.Series(dtype='float64')
    N = ndarray_1.shape[0]

    combined_matrix = np.concatenate([ndarray_1, ndarray_2], axis=0)
    combined_matrix = np.array(combined_matrix,dtype="float64")
    corr_matrix = np.corrcoef(combined_matrix)
    corr_cells = corr_matrix[:N, N:2*N]
    cells_correlation = np.mean(np.diag(corr_cells))
    # corr_results['cells_correlation'] = corr_cells_mean
    return [cell_types_correlation, cells_correlation]


def cells_auc(gold_standard,  
                imputation,
                metadata,
                celltype_key='cell_type'):
    auprc_cell_results = pd.Series(dtype='float64')
    cell_types_auprc = None
    cell_types_auroc = None
    if celltype_key in metadata.columns:
        cell_types = metadata[celltype_key].unique()

        tmp = []
        tmp2 = []
        for cell_type in cell_types:
            cells_index = metadata[celltype_key] == cell_type

            gold_standard_mean = np.mean(gold_standard[cells_index,], axis=0)
            imputation_mean = np.mean(imputation[cells_index,], axis=0)
            
            raw_binary = np.where(gold_standard_mean != 0, 1, 0)
            if np.mean(raw_binary) > 0.999:
                tmp.append(1)
                tmp2.append(1)
            else:
                # raw_binary = gold_standard[i].astype(bool).astype(int) 
                precision, recall, _ = precision_recall_curve(raw_binary, imputation_mean)
                # auprc_cell_results[cell_type] = auc(recall, precision)
                tmp.append(auc(recall, precision))
                tmp2.append(roc_auc_score(raw_binary, imputation_mean))

            # del raw_binary, precision, recall
            # gc.collect()
        cell_types_auprc = np.mean(tmp)
        cell_types_auroc = np.mean(tmp2)
        # auprc_cell_results['cell_types_auprc'] = auprc_cells_mean
        
    auprc_values = []
    auroc_values = []
    N = gold_standard.shape[0]
    for i in range(N):
        raw_binary = np.where(gold_standard[i] != 0, 1, 0)
        if raw_binary.mean() > 0.99:
            auprc_values.append(1)
            auroc_values.append(1)
        elif raw_binary.mean() == 0:
            pass
        else:
            # raw_binary = gold_standard[i].astype(bool).astype(int) 
            precision, recall, _ = precision_recall_curve(raw_binary, imputation[i])
            auprc_values.append(auc(recall, precision))
            auroc_values.append(roc_auc_score(raw_binary, imputation[i]))

        # del raw_binary, precision, recall
        # gc.collect()

    cells_auprc = np.mean(auprc_values)
    cells_auroc = np.mean(auroc_values)
    return [cell_types_auprc, cell_types_auroc, cells_auprc, cells_auroc]


def cells_nzero_ratio(gold_standard,  
                     metadata,
                     celltype_key='cell_type'):
    cell_types_nzero_ratio = None
    if celltype_key in metadata.columns:
        cell_types = metadata[celltype_key].unique()
        tmp = []
        for cell_type in cell_types:
            cells_index = metadata[celltype_key] == cell_type
            gold_standard_mean = np.mean(gold_standard[cells_index,], axis=0)            
            tmp.append(1-(gold_standard_mean==0).mean())
        cell_types_nzero_ratio = np.mean(tmp)
        
    cells_nzero_ratio = 1-(gold_standard==0).mean()
    return [cell_types_nzero_ratio, cells_nzero_ratio]


# % SCMMIB metric modules
def paired_latent_metrics(adata,
                          method,
                          cluster,
                          batch = 'batch',
                          label = 'cell_type',
                          outf = './text.txt'
                         ):
    # Cluster (deprecated)
    if cluster == 'louvain':
        sc.pp.neighbors(adata, n_neighbors=15, use_rep='X') # add 2 obsp：distances + connectivities # add 1 uns：neighbors
        sc.tl.louvain(adata)
    elif cluster == 'kmeans':
        kmeans = KMeans(n_clusters=num_clusters, random_state=2022).fit(adata.X) # add obs：kmeans
        adata.obs[cluster_method] = pd.Series(kmeans.labels_, index=adata.obs.index).astype('category')
    elif cluster == 'hc': # hierachical clustering
        hc = AgglomerativeClustering(n_clusters=num_clusters).fit(adata.X) # add obs：hc
        adata.obs[cluster_method] = pd.Series(hc.labels_,index=adata.obs.index).astype('category')

    adata.uns['metrics'] = pd.DataFrame()
    # index = method + '-' + cluster 
    index = method
    adata.obsm['latent'] = adata.X.copy()

    # 0. output
    adata.uns['metrics'].loc[index, ['Output']] = 'Embedding'
    adata.uns['metrics'].loc[index, ['nCell']] = adata.X.shape[0]

    # 1. Batch ASW
    if batch in adata.obs.columns:
        Batch_ASW = scib.metrics.silhouette_batch(adata, batch_key=batch, group_key=cluster, 
                                                  embed='latent', metric="euclidean", 
                                                  return_all=False, scale=True, verbose=False)
        adata.uns['metrics'].loc[index, ['Batch_ASW_batch']] = [Batch_ASW]

    if batch == "batch" and 'Site' in adata.obs.columns:  
        Batch_ASW = scib.metrics.silhouette_batch(adata, batch_key='Site', group_key=cluster, 
                                                  embed='latent', metric="euclidean", 
                                                  return_all=False, scale=True, verbose=False)
        adata.uns['metrics'].loc[index, ['Batch_ASW_site']] = [Batch_ASW]

    # 2. kNN graph_connectivity
    if label in adata.obs.columns:
        adata.obs[label]  = adata.obs[label].cat.remove_unused_categories()

        graph_connectivity = scmmib_gc(adata, label)
        adata.uns['metrics'].loc[index, ['graph_connectivity']] = [graph_connectivity]
    if 'cell_type.l1' in adata.obs.columns:
        adata.obs['cell_type.l1']  = adata.obs['cell_type.l1'].cat.remove_unused_categories()
        graph_connectivity = scmmib_gc(adata, 'cell_type.l1')
        adata.uns['metrics'].loc[index, ['graph_connectivity.l1']] = [graph_connectivity]

    # 3. graph iLISI
    if batch in adata.obs.columns:
        graph_iLISI = scib.metrics.ilisi_graph(adata, batch_key=batch, use_rep='latent',
                                               type_=None, k0=150, subsample=None, scale=True, n_cores=1, verbose=False)
        adata.uns['metrics'].loc[index, ['graph_iLISI']] = [graph_iLISI]

    # 4. adjusted rank index
    if label in adata.obs.columns:
        ari = adjusted_rand_score(adata.obs[label], adata.obs[cluster])
        adata.uns['metrics'].loc[index, ['ARI']] = [ari]
    if 'cell_type.l1' in adata.obs.columns:
        ari = adjusted_rand_score(adata.obs['cell_type.l1'], adata.obs[cluster])
        adata.uns['metrics'].loc[index, ['ARI.l1']] = [ari]
        
    # 5. adjusted mutual information
    if label in adata.obs.columns:
        ami = adjusted_mutual_info_score(adata.obs[label], adata.obs[cluster], average_method='arithmetic')
        adata.uns['metrics'].loc[index, ['AMI']] = [ami]
    if 'cell_type.l1' in adata.obs.columns:
        ami = adjusted_mutual_info_score(adata.obs['cell_type.l1'], adata.obs[cluster], average_method='arithmetic')
        adata.uns['metrics'].loc[index, ['AMI.l1']] = [ami]

    # 6. graph cLISI
    if label in adata.obs.columns:
        graph_cLISI = scib.metrics.clisi_graph(adata, label_key=label, use_rep='latent',
                                               type_=None, k0=150, subsample=None, scale=True, n_cores=1, verbose=False)
        adata.uns['metrics'].loc[index, ['graph_cLISI']] = [graph_cLISI]
    if 'cell_type.l1' in adata.obs.columns:
        graph_cLISI = scib.metrics.clisi_graph(adata, label_key='cell_type.l1', use_rep='latent',
                                               type_=None, k0=150, subsample=None, scale=True, n_cores=1, verbose=False)
        adata.uns['metrics'].loc[index, ['graph_cLISI.l1']] = [graph_cLISI]

    # 7. isolated labels ASW
    if batch in adata.obs.columns and label in adata.obs.columns:
        isolated_labels_ASW = scib.metrics.isolated_labels(adata, batch_key=batch, label_key=label, embed='latent',
                                                           cluster=False, iso_threshold=None, return_all=False, verbose=False)
        adata.uns['metrics'].loc[index, ['isolated_labels_ASW']] = [isolated_labels_ASW]
    if batch in adata.obs.columns and 'cell_type.l1' in adata.obs.columns:
        isolated_labels_ASW = scib.metrics.isolated_labels(adata, batch_key=batch, label_key='cell_type.l1', embed='latent',
                                                           cluster=False, iso_threshold=None, return_all=False, verbose=False)
        adata.uns['metrics'].loc[index, ['isolated_labels_ASW.l1']] = [isolated_labels_ASW]
    if outf:
        print("writing to {}".format(outf))
        adata.uns['metrics'].to_csv(outf)
    else:
        return(adata.uns['metrics'])


def paired_graph_metrics(adata,
                         method,
                         cluster,
                         batch = 'batch',
                         label = 'cell_type',
                         outf = None
                        ):
    # Cluster
    if cluster == 'louvain':
        sc.pp.neighbors(adata, n_neighbors=15) # add 2 obsp：distances + connectivities # add 1 uns：neighbors
        adata.obsp['connectivities'] = csr_matrix(adata.obsp[method])
        sc.tl.louvain(adata, obsp=method) # Louvain 
    elif cluster == 'kmeans':
        kmeans = KMeans(n_clusters=num_clusters, random_state=2022).fit(adata.X) # add obs：kmeans
        adata.obs[cluster_method] = pd.Series(kmeans.labels_, index=adata.obs.index).astype('category')
    elif cluster == 'hc': # hierachical clustering
        hc = AgglomerativeClustering(n_clusters=num_clusters).fit(adata.X) # add obs：hc
        adata.obs[cluster_method] = pd.Series(hc.labels_,index=adata.obs.index).astype('category')
    
    if 'metrics' not in adata.uns:
        adata.uns['metrics'] = pd.DataFrame()
    
    adata.uns['metrics'] = pd.DataFrame()
    # index = method + '-' + cluster 
    index = method
    # 1. output
    adata.uns['metrics'].loc[index, ['Output']] = 'Graph'
    adata.uns['metrics'].loc[index, ['nCell']] = adata.X.shape[0]

    # 2. kNN graph_connectivity
    if label in adata.obs.columns:
        # adata.obsp["connectivities"] = adata.obsp[method]
        score = scmmib_gc(adata, label)
        adata.uns['metrics'].loc[index, ['graph_connectivity']] = [score]
    if 'cell_type.l1' in adata.obs.columns:
        # adata.obsp["connectivities"] = adata.obsp[method]
        score = scmmib_gc(adata, 'cell_type.l1')
        adata.uns['metrics'].loc[index, ['graph_connectivity.l1']] = [score]

    # 3. graph iLISI: "Fixing"
    if batch in adata.obs.columns:
#         graph_iLISI = scib.metrics.ilisi_graph(adata, batch_key=batch, type_="knn",
#                                               k0=150, subsample=None, scale=True, n_cores=1, verbose=False)
        graph_iLISI = scib.metrics.ilisi_graph(adata, batch_key=batch, type_="knn", verbose=False)
        adata.uns['metrics'].loc[index, ['graph_iLISI']] = [graph_iLISI]

    # 4. adjusted rank index
    if label in adata.obs.columns:
        ari = adjusted_rand_score(adata.obs[label], adata.obs[cluster])
        adata.uns['metrics'].loc[index, ['ARI']] = [ari]
    if 'cell_type.l1' in adata.obs.columns:
        ari = adjusted_rand_score(adata.obs['cell_type.l1'], adata.obs[cluster])
        adata.uns['metrics'].loc[index, ['ARI.l1']] = [ari]

    # 5. adjusted mutual information
    if label in adata.obs.columns:
        ami = adjusted_mutual_info_score(adata.obs[label], adata.obs[cluster], average_method='arithmetic')
        adata.uns['metrics'].loc[index, ['AMI']] = [ami]
    if 'cell_type.l1' in adata.obs.columns:
        ami = adjusted_mutual_info_score(adata.obs['cell_type.l1'], adata.obs[cluster], average_method='arithmetic')
        adata.uns['metrics'].loc[index, ['AMI.l1']] = [ami]

    # # 6. graph cLISI
    if label in adata.obs.columns:
        graph_cLISI = scib.metrics.clisi_graph(adata, label_key=label,
                                               type_='knn')
        adata.uns['metrics'].loc[index, ['graph_cLISI']] = [graph_cLISI]
    if 'cell_type.l1' in adata.obs.columns:
        graph_cLISI = scib.metrics.clisi_graph(adata, label_key='cell_type.l1',
                                               type_='knn')
        adata.uns['metrics'].loc[index, ['graph_cLISI.l1']] = [graph_cLISI]
    # return(adata)
    if outf:
        adata.uns['metrics'].to_csv(outf)
    else:
        return(adata.uns['metrics'])


def unpaired_latent_metrics(adata,
                            method,
                            cluster,
                            batch = 'batch',
                            label = 'cell_type',
                            mods = ['RNA', 'ATAC'],
                            outf = None,
                            embed_acc = True
                           ):
    adata.uns['metrics'] = pd.DataFrame()
    # index = method +  '-' + cluster
    index = method
    adata.uns['metrics'].loc[index, ['Output']] = 'Embedding'
    adata.uns['metrics'].loc[index, ['nCell']] = adata.X.shape[0]
    # for mod in ['RNA', 'ATAC']:
    if embed_acc:
        for mod in mods:
            # Cluster
            if cluster == 'louvain':
                sc.pp.neighbors(adata, use_rep=mod, n_neighbors=15) # add 2 obsp：distances + connectivities # add 1 uns：neighbors
                sc.tl.louvain(adata)
            elif cluster == 'kmeans':
                kmeans = KMeans(n_clusters=15, random_state=2022).fit(adata.obsm[mod]) # add obs：kmeans
                adata.obs[index] = pd.Series(kmeans.labels_, index=adata.obs.index).astype('category')
            elif cluster == 'hc': # hierachical clustering
                hc = AgglomerativeClustering(n_clusters=15).fit(adata.obsm[mod]) # add obs：hc
                adata.obs[index] = pd.Series(hc.labels_,index=adata.obs.index).astype('category')

            # Metrics
            # 0. output
            

            # 1. Batch ASW
            # if batch in adata.obs.columns:
            #     Batch_ASW = scib.metrics.silhouette_batch(adata, batch_key=batch, group_key=cluster, 
            #                                             embed=mod, metric="euclidean", return_all=False, scale=True, verbose=False)
            #     adata.uns['metrics'].loc[index, ['Batch_ASW_batch' + '-' + mod]] = [Batch_ASW]
            # if 'Site' in adata.obs.columns:
            #     Batch_ASW = scib.metrics.silhouette_batch(adata, batch_key='Site', group_key=cluster, 
            #                                             embed=mod, metric="euclidean", return_all=False, scale=True, verbose=False)
            #     adata.uns['metrics'].loc[index, ['Batch_ASW_site' + '-' + mod]] = [Batch_ASW]

            # # 2. kNN graph_connectivity
            # if label in adata.obs.columns:
            #     graph_connectivity = scmmib_gc(adata, label)
            #     adata.uns['metrics'].loc[index, ['graph_connectivity' + '-' + mod]] = [graph_connectivity]
            # if 'cell_type.l1' in adata.obs.columns:
            #     graph_connectivity = scmmib_gc(adata, 'cell_type.l1')
            #     adata.uns['metrics'].loc[index, ['graph_connectivity.l1' + '-' + mod]] = [graph_connectivity]

            # # 3. graph iLISI
            # if batch in adata.obs.columns:    
            #     adata.obsm["X_emb"] = adata.obsm[mod]
            #     graph_iLISI = scib.metrics.ilisi_graph(adata, batch_key=batch, 
            #                                         type_=None, k0=150, subsample=None, scale=True, n_cores=1, verbose=False)
            #     adata.uns['metrics'].loc[index, ['graph_iLISI' + '-' + mod]] = [graph_iLISI]

            # 4. adjusted rank index
            if label in adata.obs.columns:
                ari = adjusted_rand_score(adata.obs[label], adata.obs[cluster])
                adata.uns['metrics'].loc[index, ['ARI' + '-' + mod]] = [ari]
            if 'cell_type.l1' in adata.obs.columns:
                ari = adjusted_rand_score(adata.obs['cell_type.l1'], adata.obs[cluster])
                adata.uns['metrics'].loc[index, ['ARI.l1' + '-' + mod]] = [ari]

            # 5. adjusted mutual information
            if label in adata.obs.columns:
                ami = adjusted_mutual_info_score(adata.obs[label], adata.obs[cluster], average_method='arithmetic')
                adata.uns['metrics'].loc[index, ['AMI' + '-' + mod]] = [ami]
            if 'cell_type.l1' in adata.obs.columns:
                ami = adjusted_mutual_info_score(adata.obs['cell_type.l1'], adata.obs[cluster], average_method='arithmetic')
                adata.uns['metrics'].loc[index, ['AMI.l1' + '-' + mod]] = [ami]

            # 6. graph cLISI
            # if batch in adata.obs.columns and label in adata.obs.columns:
            adata.obsm["X_emb"] = adata.obsm[mod]
            graph_cLISI = scib.metrics.clisi_graph(adata, label_key=label, 
                                                type_=None, k0=150, subsample=None, scale=True, n_cores=1, verbose=False)
            adata.uns['metrics'].loc[index, ['graph_cLISI' + '-' + mod]] = [graph_cLISI]
            if 'cell_type.l1' in adata.obs.columns:
                adata.obsm["X_emb"] = adata.obsm[mod]
                graph_cLISI = scib.metrics.clisi_graph(adata, label_key='cell_type.l1', 
                                                    type_=None, k0=150, subsample=None, scale=True, n_cores=1, verbose=False)
                adata.uns['metrics'].loc[index, ['graph_cLISI.l1' + '-' + mod]] = [graph_cLISI]

            # 7. isolated_labels
            if batch in adata.obs.columns and label in adata.obs.columns:
                isolated_labels_ASW = scib.metrics.isolated_labels(adata, batch_key=batch, label_key=label, 
                                                                embed=mod, cluster=False, iso_threshold=None, return_all=False, verbose=False)
                adata.uns['metrics'].loc[index, ['isolated_labels_ASW' + '-' + mod]] = [isolated_labels_ASW]
            if batch in adata.obs.columns and 'cell_type.l1' in adata.obs.columns:
                isolated_labels_ASW = scib.metrics.isolated_labels(adata, batch_key=batch, label_key='cell_type.l1', 
                                                                embed=mod, cluster=False, iso_threshold=None, return_all=False, verbose=False)
                adata.uns['metrics'].loc[index, ['isolated_labels_ASW.l1' + '-' + mod]] = [isolated_labels_ASW]
    
    # 8. FOSCTTM
    score = foscttm(adata.obsm['RNA'], adata.obsm[mods[1]])
    adata.uns['metrics'].loc[index, ['FOSCTTM']] = [score]
    
    # 9. nearest_cell_barcode
    score = nearest_cell_barcode(adata, rep_1='RNA', rep_2=mods[1])
    adata.uns['metrics'].loc[index, ['nearest_cell_barcode']] = [score['Mean']]
    adata.uns['metrics'].loc[index, ['nearest_cell_barcode-RNA']] = [score['RNA']]
    adata.uns['metrics'].loc[index, [f'nearest_cell_barcode-{mods[1]}']] = [score[mods[1]]]
    
    # 10. nearest_cell_celltype
    if label in adata.obs.columns:
        adata.obs[label]  = adata.obs[label].cat.remove_unused_categories()
        score = nearest_cell_celltype(adata, rep_1='RNA', rep_2=mods[1], label=label)
        adata.uns['metrics'].loc[index, ['nearest_cell_celltype']] = [score['Mean']]
        adata.uns['metrics'].loc[index, ['nearest_cell_celltype-RNA']] = [score['RNA']]
        adata.uns['metrics'].loc[index, [f'nearest_cell_celltype-{mods[1]}']] = [score[mods[1]]]
    if 'cell_type.l1' in adata.obs.columns:
        adata.obs['cell_type.l1']  = adata.obs['cell_type.l1'].cat.remove_unused_categories()
        score = nearest_cell_celltype(adata, rep_1='RNA', rep_2=mods[1], label='cell_type.l1')
        adata.uns['metrics'].loc[index, ['nearest_cell_celltype.l1']] = [score['Mean']]
        adata.uns['metrics'].loc[index, ['nearest_cell_celltype.l1-RNA']] = [score['RNA']]
        adata.uns['metrics'].loc[index, [f'nearest_cell_celltype.l1-{mods[1]}']] = [score[mods[1]]]
    
    if outf:
        print("writing to {}".format(outf))
        adata.uns['metrics'].to_csv(outf)
    else:
        return(adata.uns['metrics'])



def mosaic_latent_metrics(latents,metadatas,paired="s1d1", unpaired="s3d10", 
                            mod2="atac", batch="batch",label="cell_type",latent_path="", method='sciPENN', writef=True):
    
    if len(latents) != len(metadatas):
        print("metadata do not match input latent files!")
        return
    
    if mod2 != "atac" and mod2!="adt":
        print("Unknown modality! mod2 must be atac or adt!")
        return

    if len(latents) == 3:
        # latent_rna = latent[latent.index.str.endswith('_rna')]
        # latent_mod2 = latent[latent.index.str.endswith(f'_{mod2}')]
        latent_paried = latents[0].copy()
        latent_rna = latents[1].copy()
        latent_mod2 = latents[2].copy()
        latent_rna.columns = latent_paried.columns
        latent_mod2.columns = latent_paried.columns
        latent = pd.concat([latent_paried, latent_rna, latent_mod2], axis=0)
        metadata_all_origin = pd.concat([metadatas[0],metadatas[1],metadatas[2]])

        metadata_paired = metadatas[0].copy()
        metadata_rna = metadatas[1].copy()
        metadata_mod2 = metadatas[2].copy()
        metadata_rna.index = [x + '_rna' for x in metadata_rna.index]
        metadata_rna[batch] = [x + '_rna' for x in metadata_rna[batch]]
        metadata_rna[batch] = metadata_rna[batch].astype('category')

        metadata_mod2.index = [x + f'_{mod2}' for x in metadata_mod2.index]
        metadata_mod2[batch] = [x + f'_{mod2}' for x in metadata_mod2[batch]]
        metadata_mod2[batch] = metadata_mod2[batch].astype('category')

        metadata_all = pd.concat([metadata_paired, metadata_rna, metadata_mod2], axis=0)
        metadata_all[batch] = metadata_all[batch].astype('category')

        latent_reordered = latent.reindex(metadata_all.index)
        adata_all = sc.AnnData(latent_reordered, obs=metadata_all, dtype='float32')
        paired_metrics_path = latent_path.split('latent')[0] + 'paired-metrics.csv'
        if writef:
            paired_latent_metrics(adata_all, method = method, cluster = 'louvain', batch = batch, label = label,outf=paired_metrics_path)
        else:
            adata_all = paired_latent_metrics(adata_all, method = method, cluster = 'louvain', batch = batch, label = label,outf=None)

        metadata_unpaired = metadata_all_origin[metadata_all_origin[batch] == unpaired]
        latent_rna.index = latent_rna.index.str.replace('_rna$', '', regex=True)
        latent_rna_reordered = latent_rna.reindex(metadata_unpaired.index)

        latent_mod2.index = latent_mod2.index.str.replace(f'_{mod2}$', '', regex=True)
        latent_mod2_reordered = latent_mod2.reindex(metadata_unpaired.index)

        adata_unpaired = sc.AnnData(latent_rna_reordered, obs=metadata_unpaired, dtype='float32')
        adata_unpaired.obsm['RNA'] = latent_rna_reordered
        adata_unpaired.obsm['{}'.format(mod2.upper())] = latent_mod2_reordered

        cell_type_counts = adata_unpaired.obs['cell_type'].value_counts()
        valid_cell_types = cell_type_counts[cell_type_counts > 0].index
        adata_unpaired = adata_unpaired[adata_unpaired.obs['cell_type'].isin(valid_cell_types)].copy()
        
        unpaired_metrics_path = latent_path.split('latent')[0] + 'unpaired-metrics.csv'
        if writef:
            unpaired_latent_metrics(adata_unpaired, method = method, cluster = 'louvain', batch = batch, label = label, mods = ["RNA",mod2.upper()], outf=unpaired_metrics_path, embed_acc=False)
        else:
            adata_unpaired2 = unpaired_latent_metrics(adata_unpaired, method = method, cluster = 'louvain', batch = batch, label = label, mods = ["RNA",mod2.upper()],outf=None, embed_acc=False)

    if mod2 == "adt":
        latent_paried = latents[0].copy()
        latent_rna = latents[1].copy()
        latent_rna.columns = latent_paried.columns
        latent = pd.concat([latent_paried, latent_rna], axis=0)
        
        metadata_paired = metadatas[0].copy()
        metadata_rna = metadatas[1].copy()
        metadata_rna.index = [x + '_rna' for x in metadata_rna.index]
        metadata_rna[batch] = [x + '_rna' for x in metadata_rna[batch]]
        # metadata_rna[batch] = metadata_rna[batch].astype('category')

        metadata_noadt = pd.concat([metadata_paired, metadata_rna], axis=0)
        metadata_noadt[batch] = metadata_noadt[batch].astype('category')

        latent_reordered = latent.reindex(metadata_noadt.index)
        adata_noadt = sc.AnnData(latent_reordered, obs=metadata_noadt, dtype='float32')
        
        noADT_paired_metrics_path = latent_path.split('latent')[0] + 'noADT_paired-metrics.csv'
        if writef:
            paired_latent_metrics(adata_noadt, method = method, cluster = 'louvain', batch = batch, label = label,outf=noADT_paired_metrics_path)
        else:
            adata_noadt = paired_latent_metrics(adata_noadt, method = method, cluster = 'louvain', batch = batch, label = label,outf=None)

    if not writef:
        if mod2 == "adt":
            if len(latents) == 3:
                return([adata_all, adata_unpaired2, adata_noadt])
            else:
                return([adata_noadt])
        else:
            return([adata_all, adata_unpaired2])


def mosaic_cnk_latent_metrics(latents,metadatas,paired="s1d1", unpaired="s3d10", 
                            mod2="atac", batch="batch",label="cell_type",latent_path="", method='sciPENN', writef=True):
    """
    Function for cnk paried dataset size robustness evaluation.
    Here we only assess the accuracy in the unpaired modalities.
    """
    if len(latents) != len(metadatas):
        print("metadata do not match input latent files!")
        return
    if mod2 != "atac" and mod2!="adt":
        print("Unknown modality! mod2 must be atac or adt!")
        return

    if len(latents) == 3:
        latent_paried = latents[0].copy()
        latent_rna = latents[1].copy()
        latent_mod2 = latents[2].copy()
        latent_rna.columns = latent_paried.columns
        latent_mod2.columns = latent_paried.columns
        latent = pd.concat([latent_rna, latent_mod2], axis=0)

        metadata_paired = metadatas[0].copy()
        metadata_rna = metadatas[1].copy()
        metadata_mod2 = metadatas[2].copy()
        metadata_all_origin = pd.concat([metadatas[0],metadatas[1],metadatas[2]])
        metadata_rna.index = [x + '_rna' for x in metadata_rna.index]
        metadata_rna[batch] = [x + '_rna' for x in metadata_rna[batch]]
        metadata_rna[batch] = metadata_rna[batch].astype('category')

        metadata_mod2.index = [x + f'_{mod2}' for x in metadata_mod2.index]
        metadata_mod2[batch] = [x + f'_{mod2}' for x in metadata_mod2[batch]]
        metadata_mod2[batch] = metadata_mod2[batch].astype('category')

        metadata_all = pd.concat([metadata_rna, metadata_mod2], axis=0)
        metadata_all[batch] = metadata_all[batch].astype('category')

        latent_reordered = latent.reindex(metadata_all.index)
        adata_all = sc.AnnData(latent_reordered, obs=metadata_all, dtype='float32')
        paired_metrics_path = latent_path.split('latent')[0] + 'new-paired-metrics.csv'
        if writef:
            paired_latent_metrics(adata_all, method = method, cluster = 'louvain', batch = batch, label = label, outf=paired_metrics_path)
        else:
            adata_all = paired_latent_metrics(adata_all, method = method, cluster = 'louvain', batch = batch, label = label, outf=None)

        metadata_unpaired = metadata_unpaired = metadata_all_origin[metadata_all_origin[batch] == unpaired]
        latent_rna.index = latent_rna.index.str.replace('_rna$', '', regex=True)
        latent_rna_reordered = latent_rna.reindex(metadata_unpaired.index)

        latent_mod2.index = latent_mod2.index.str.replace(f'_{mod2}$', '', regex=True)
        latent_mod2_reordered = latent_mod2.reindex(metadata_unpaired.index)

        adata_unpaired = sc.AnnData(latent_rna_reordered, obs=metadata_unpaired, dtype='float32')
        adata_unpaired.obsm['RNA'] = latent_rna_reordered
        adata_unpaired.obsm['{}'.format(mod2.upper())] = latent_mod2_reordered

        cell_type_counts = adata_unpaired.obs['cell_type'].value_counts()
        valid_cell_types = cell_type_counts[cell_type_counts > 0].index
        adata_unpaired = adata_unpaired[adata_unpaired.obs['cell_type'].isin(valid_cell_types)].copy()
        
        unpaired_metrics_path = latent_path.split('latent')[0] + 'new-unpaired-metrics.csv'
        if writef:
            unpaired_latent_metrics(adata_unpaired, method = method, cluster = 'louvain', batch = None, label = label, mods = ["RNA",mod2.upper()], outf=unpaired_metrics_path, embed_acc=False)
        else:
            adata_unpaired = unpaired_latent_metrics(adata_unpaired, method = method, cluster = 'louvain', batch = None, label = label, mods = ["RNA",mod2.upper()], outf=None, embed_acc=False)

    if mod2 == "adt":
        # latent_paried = latents[0]
        latent_rna = latents[1].copy()
        # latent = pd.concat([latent_paried, latent_rna], axis=0)
        latent = latent_rna.copy()
        # metadata_paired = metadatas[0]
        metadata_rna = metadatas[1].copy()
        metadata_rna.index = [x + '_rna' for x in metadata_rna.index]
        metadata_rna[batch] = [x + '_rna' for x in metadata_rna[batch]]
        metadata_rna[batch] = metadata_rna[batch].astype('category')

        # metadata_noadt = pd.concat([metadata_paired, metadata_rna], axis=0)
        metadata_noadt = metadata_rna.copy()
        # metadata_noadt[batch] = metadata_noadt[batch].astype('category')

        latent_reordered = latent.reindex(metadata_noadt.index)
        adata_noadt = sc.AnnData(latent_reordered, obs=metadata_noadt, dtype='float32')
        
        noADT_paired_metrics_path = latent_path.split('latent')[0] + 'new-noADT_paired-metrics.csv'
        if writef:
            paired_latent_metrics(adata_noadt, method = method, cluster = 'louvain', batch = None, label = label,outf=noADT_paired_metrics_path)
        else:
            adata_noadt = paired_latent_metrics(adata_noadt, method = method, cluster = 'louvain', batch = None, label = label,outf=None)
    if not writef:
        if mod2 == "adt":
            if len(latents) == 3:
                return([adata_all, adata_unpaired, adata_noadt])
            else:
                return([adata_noadt])
        else:
            return([adata_all, adata_unpaired])


def imputation_pair_rna_atac(metadata_path,rna_imp_path, rna_path,atac_imp_path, atac_path, method, outf=None ):
    """
    Calculate metrics for single alogirthm results.
    metadata_path: path of cell metadata.
    rna_imp_path: path of rna imputation csv file.
    rna_path: path of rna raw h5ad file.
    atac_imp_path: path of atac imputation csv file.
    atac_path: path of atac raw h5ad file.
    outf: output to given file name or pd.Dataframe() (None). 
    """
    ## prepare input data and parameter
    # method=rna_imp_path.split("-")[-2]
    cluster="louvain"
    # Load metadata
    metadata = pd.read_csv(metadata_path, index_col='barcode', header=0, dtype='category')
    # RNA impute
    # rna_imp = pd.read_csv(rna_imp_path, header=0, index_col=0)
    rna_imp = dt.fread(rna_imp_path)
    rna_imp = rna_imp.to_pandas()
    rna_imp.index = rna_imp.iloc[:,0]
    rna_imp.drop(rna_imp.columns[[0]],axis=1,inplace=True)

    rna_imp_reordered = rna_imp.reindex(metadata.index)
    rna_imp_features = rna_imp_reordered.columns.values
    rna_imputation = rna_imp_reordered.to_numpy()


    ncells = rna_imp.shape[0]
    
    if ncells > 10000:
        sampling = True
    else:
        sampling = False
    if sampling:
        metadata = metadata.iloc[:10000,]
        rna_imputation = rna_imputation[:10000]
    # RNA raw
    rna = sc.read_h5ad(rna_path)
    rna_subset = rna[:, rna.var_names.isin(rna_imp_features)]
    rna_raw_counts = rna_subset.X.toarray() # raw counts   
    print(rna_raw_counts.shape)
    if sampling:
        rna_raw_counts = rna_raw_counts[:10000]
    
    print(rna_raw_counts.shape)
    rna_knn_smooth = knn_smoothing(rna_raw_counts.transpose(), k=20, d=10, dither=0).transpose()

    # ATAC impute
    # atac_imp = pd.read_csv(atac_imp_path, header=0, index_col=0) # too slow
    atac_imp = dt.fread(atac_imp_path)
    atac_imp = atac_imp.to_pandas()
    atac_imp.index = atac_imp.iloc[:,0]
    atac_imp.drop(atac_imp.columns[[0]],axis=1,inplace=True)

    atac_imp_reordered = atac_imp.reindex(metadata.index)
    atac_imp_features = atac_imp_reordered.columns.values
    atac_imputation = atac_imp_reordered.to_numpy()
    if sampling:
        atac_imputation = atac_imputation[:10000]

    # ATAC raw
    atac = sc.read_h5ad(atac_path)
    # atac = sc.read_h5ad(gam_path)
    atac_subset = atac[:, atac.var_names.isin(atac_imp_features)]
    atac_raw_counts = atac_subset.X.toarray() # raw counts
    print(atac_raw_counts.shape)
    if sampling:
        atac_raw_counts = atac_raw_counts[:10000]

    # knn_smooth
    atac_knn_smooth = knn_smoothing(atac_raw_counts.transpose(), k=20, d=10, dither=0).transpose()

    data_res = pd.DataFrame(index=[f"{method}-{cluster}"],columns=["rna_n_features","atac_n_features"],
                            data=[[len(rna_imp_features),len(atac_imp_features)]])

    print("KNN smooth completed.")
    rna_res =  cells_nzero_ratio(rna_raw_counts, metadata, 'cell_type')+cells_nzero_ratio(rna_knn_smooth, metadata, 'cell_type') + \
    cells_correlation(rna_raw_counts,rna_imputation,metadata, 'cell_type') + \
    cells_correlation(rna_knn_smooth,rna_imputation,metadata, 'cell_type' ) 
    
    print("RNA metrics completed")

    rna_res2 = pd.DataFrame(rna_res).T
    rna_res2.columns = ["rna_raw_nzero_ratio_celltype","rna_raw_nzero_ratio_cell",
                            "rna_knn_nzero_ratio_celltype", "rna_knn_nzero_ratio_cell",
                            "rna_raw_celltype_cor", "rna_raw_cell_cor",
                            "rna_knn_celltype_cor", "rna_knn_cell_cor"
                            ]
    rna_res2.index=[f"{method}-{cluster}"]

    atac_res =  cells_nzero_ratio(atac_raw_counts, metadata, 'cell_type')+cells_nzero_ratio(atac_knn_smooth, metadata, 'cell_type') + \
    cells_auc(atac_raw_counts, atac_imputation,metadata, 'cell_type') + \
    cells_auc(atac_knn_smooth, atac_imputation,metadata, 'cell_type' ) 
    
    print("ATAC metrics completed")
    atac_res2 = pd.DataFrame(atac_res).T

    atac_res2.columns = ["atac_raw_nzero_ratio_celltype","atac_raw_nzero_ratio_cell",
                            "atac_knn_nzero_ratio_celltype", "atac_knn_nzero_ratio_cell",
                            "atac_raw_celltype_auprc", "atac_raw_celltype_auroc", "atac_raw_cell_auprc", "atac_raw_cell_auroc",
                            "atac_knn_celltype_auprc", "atac_knn_celltype_auroc", "atac_knn_cell_auprc", "atac_knn_cell_auroc"
                            ]
    atac_res2.index=[f"{method}-{cluster}"]
    out_metrics = pd.concat([data_res, rna_res2,atac_res2],axis=1)
    if outf:
        # if 'imputation_rna' in rna_imp_path:
        #     out_metrics_path = rna_imp_path.split('imputation_rna')[0] + 'imputation-metrics.csv'
        # elif 'imputation-rna' in rna_imp_path:
        #     out_metrics_path = rna_imp_path.split('imputation-rna')[0] + 'imputation-metrics.csv'
        # else:
        #     out_metrics_path = rna_imp_path + 'imputation-metrics.csv'
        print(f"output to {outf}")
        print(out_metrics)
        out_metrics.to_csv(outf)
    else:
        return(out_metrics)


def imputation_mosaic_rna_atac(metadata_path,rna_imp_path, rna_path,atac_imp_path, atac_path,method,paired="s1d1",unpaired="s3d10",batch="batch", outf=None ):
    """
    Calculate metrics for single alogirthm results.
    metadata_path: path of cell metadata.
    rna_imp_path_list: path of rna imputation csv files with same raw data gold standards.
    rna_path: path of rna raw h5ad file.
    atac_imp_path: path of atac imputation csv file.
    atac_path: path of atac raw h5ad file.
    outf: output to given file  name or pd.Dataframe() (None) Only for batch metrics calculation.
    """
  
    ## prepare input data and parameter
    # method=rna_imp_path.split("-")[-2]
    cluster="louvain"
    # Load metadata
    metadata = pd.read_csv(metadata_path, index_col='barcode', header=0, dtype='category')
    metadata_paired = metadata[metadata[batch] == paired]
    metadata_unpaired = metadata[metadata[batch] == unpaired]

    # RNA impute
    # rna_imp = pd.read_csv(rna_imp_path, header=0, index_col=0)
    rna_imp = dt.fread(rna_imp_path)
    rna_imp = rna_imp.to_pandas()
    rna_imp.index = rna_imp.iloc[:,0]
    rna_imp.drop(rna_imp.columns[[0]],axis=1,inplace=True)

    rna_imp.index = rna_imp.index.str.replace('_atac$', '', regex=True)
    rna_imp_reordered = rna_imp.reindex(metadata_unpaired.index)
    rna_imp_features = rna_imp_reordered.columns.values
    # rna_imputation = rna_imp_reordered.to_numpy()

    # RNA raw
    rna = sc.read_h5ad(rna_path)
    rna_unpaired = rna[rna.obs[batch] == unpaired, :]
    rna_unpaired_subset = rna_unpaired[:, rna_unpaired.var_names.isin(rna_imp_features)]
    rna_raw_counts = rna_unpaired_subset.X.toarray() # raw counts        
    
    rna_imputation = rna_imp_reordered.loc[:,rna_unpaired_subset.var_names].to_numpy()

    rna_knn_smooth = knn_smoothing(rna_raw_counts.transpose(), k=20, d=10, dither=0).transpose()


    # ATAC impute
    # atac_imp = pd.read_csv(atac_imp_path, header=0, index_col=0) # too slow
    atac_imp = dt.fread(atac_imp_path)
    atac_imp = atac_imp.to_pandas()
    atac_imp.index = atac_imp.iloc[:,0]
    atac_imp.drop(atac_imp.columns[[0]],axis=1,inplace=True)

    atac_imp.index = atac_imp.index.str.replace('_rna$', '', regex=True)
    atac_imp_reordered = atac_imp.reindex(metadata_unpaired.index)
    atac_imp_features = atac_imp_reordered.columns.values
    # atac_imputation = atac_imp_reordered.to_numpy()

    # ATAC raw
 
    atac = sc.read_h5ad(atac_path)
    atac_unpaired = atac[atac.obs[batch] == unpaired, :]
    atac_unpaired_subset = atac_unpaired[:, atac_unpaired.var_names.isin(atac_imp_features)]
    atac_raw_counts = atac_unpaired_subset.X.toarray() # raw counts
    # atac_subset = atac[:, atac.var_names.isin(atac_imp_features)]
    # atac_raw_counts = atac_subset.X.toarray() # raw counts

    atac_imputation = atac_imp_reordered.loc[:,atac_unpaired_subset.var_names].to_numpy()
    # knn_smooth
    atac_knn_smooth = knn_smoothing(atac_raw_counts.transpose(), k=20, d=10, dither=0).transpose()

    print("KNN smooth complete.")
    data_res = pd.DataFrame(index=[f"{method}-{cluster}"],columns=["rna_n_features","atac_n_features"],
                            data=[[len(rna_imp_features),len(atac_imp_features)]])


    rna_res =  cells_nzero_ratio(rna_raw_counts, metadata_unpaired, 'cell_type')+cells_nzero_ratio(rna_knn_smooth, metadata_unpaired, 'cell_type') + \
        cells_correlation(rna_raw_counts,rna_imputation,metadata_unpaired, 'cell_type') + \
        cells_correlation(rna_knn_smooth,rna_imputation,metadata_unpaired, 'cell_type' ) 
        
    print("RNA metrics complete.")

    rna_res2 = pd.DataFrame(rna_res).T
    rna_res2.columns = ["rna_raw_nzero_ratio_celltype","rna_raw_nzero_ratio_cell",
                            "rna_knn_nzero_ratio_celltype", "rna_knn_nzero_ratio_cell",
                            "rna_raw_celltype_cor", "rna_raw_cell_cor",
                            "rna_knn_celltype_cor", "rna_knn_cell_cor"
                            ]
    rna_res2.index=[f"{method}-{cluster}"]

    atac_res =  cells_nzero_ratio(atac_raw_counts, metadata_unpaired, 'cell_type')+ \
        cells_nzero_ratio(atac_knn_smooth, metadata_unpaired, 'cell_type') + \
        cells_auc(atac_raw_counts, atac_imputation,metadata_unpaired, 'cell_type') + \
        cells_auc(atac_knn_smooth, atac_imputation,metadata_unpaired, 'cell_type' ) 

    print("ATAC metrics complete.")
    atac_res2 = pd.DataFrame(atac_res).T

    atac_res2.columns = ["atac_raw_nzero_ratio_celltype","atac_raw_nzero_ratio_cell",
                            "atac_knn_nzero_ratio_celltype", "atac_knn_nzero_ratio_cell",
                            "atac_raw_celltype_auprc", "atac_raw_celltype_auroc", "atac_raw_cell_auprc", "atac_raw_cell_auroc",
                            "atac_knn_celltype_auprc", "atac_knn_celltype_auroc", "atac_knn_cell_auprc", "atac_knn_cell_auroc"
                            ]
    atac_res2.index=[f"{method}-{cluster}"]
    out_metrics = pd.concat([data_res, rna_res2,atac_res2],axis=1)
    if outf:
        # if 'imputation_rna' in rna_imp_path:
        #     out_metrics_path = rna_imp_path.split('imputation_rna')[0] + 'imputation-metrics.csv'
        # elif 'imputation-rna' in rna_imp_path:
        #     out_metrics_path = rna_imp_path.split('imputation-rna')[0] + 'imputation-metrics.csv'
        # else:
        #     out_metrics_path = rna_imp_path + 'imputation-metrics.csv'
        print(f"output to {outf}")
        out_metrics.to_csv(outf)

    else:
        print("output to stdout")
        return(out_metrics)


def imputation_mosaic_rna_adt(metadata_path,rna_imp_path=None, rna_path=None,adt_imp_path=None, adt_path=None, method="sciPENN",paired="s3d6", unpaired="s2d1", batch="batch", outf=None ):
    """
    Calculate metrics for single alogirthm results.
    metadata_path: path of cell metadata.
    rna_imp_path_list: path of rna imputation csv files with same raw data gold standards.
    rna_path: path of rna raw h5ad file.
    adt_imp_path: path of adt imputation csv file. If None, adt imputation file is omitted.
    adt_path: path of adt raw h5ad file. If None, adt imputation file is omitted.
    outf: output to given filename) or pd.Dataframe() (None).
    """
  
    ## prepare input data and parameter
    # method=rna_imp_path.split("-")[-2]
    cluster="louvain"
    # Load metadata
    metadata = pd.read_csv(metadata_path, index_col='barcode', header=0, dtype='category')
    metadata_paired = metadata[metadata[batch] == paired]
    metadata_unpaired = metadata[metadata[batch] == unpaired]

    # RNA impute
    if rna_imp_path and rna_path:
        rna_imp = pd.read_csv(rna_imp_path, header=0, index_col=0)
        rna_imp.index = rna_imp.index.str.replace('_adt$', '', regex=True)
        rna_imp_reordered = rna_imp.reindex(metadata_unpaired.index)
        rna_imp_features = rna_imp_reordered.columns.values
        # rna_imputation = rna_imp_reordered.to_numpy()

        
        # RNA raw
        rna = sc.read_h5ad(rna_path)
        rna_unpaired = rna[rna.obs[batch] == unpaired, :]
        rna_unpaired_subset = rna_unpaired[:, rna_unpaired.var_names.isin(rna_imp_features)]
        rna_raw_counts = rna_unpaired_subset.X.toarray() # raw counts        

        rna_imputation = rna_imp_reordered.loc[:,rna_unpaired_subset.var_names].to_numpy()
        
        print("rna_imputation shape {}".format(rna_imputation.shape))
        print("rna_raw shape {}".format(rna_raw_counts.shape))
        # rna_subset = rna[:, rna.var_names.isin(rna_imp_features)]
        # rna_raw_counts = rna_subset.X.toarray() # raw counts        
        rna_knn_smooth = knn_smoothing(rna_raw_counts.transpose(), k=20, d=10, dither=0).transpose()


    # ADT impute
    adt_imp = pd.read_csv(adt_imp_path, header=0, index_col=0)    
    adt_imp.index = adt_imp.index.str.replace('_rna$', '', regex=True)
    adt_imp_reordered = adt_imp.reindex(metadata_unpaired.index)
    adt_imp_features = adt_imp_reordered.columns.values
    # adt_imputation = adt_imp_reordered.to_numpy()
    
    # ADT raw
    adt = sc.read_h5ad(adt_path)
    adt_unpaired = adt[adt.obs[batch] == unpaired, :]
    adt_unpaired_subset = adt_unpaired[:, adt_unpaired.var_names.isin(adt_imp_features)]
    adt_raw_counts = adt_unpaired_subset.X.toarray() # raw counts
    # adt_subset = adt[:, adt.var_names.isin(adt_imp_features)]
    # adt_raw_counts = adt_subset.X.toarray() # raw counts

    adt_imputation = adt_imp_reordered.loc[:,adt_unpaired_subset.var_names].to_numpy()

    if rna_imp_path and rna_path:
        data_res = pd.DataFrame(index=[f"{method}-{cluster}"],columns=["rna_n_features","adt_n_features"],
                            data=[[len(rna_imp_features),len(adt_imp_features)]])
    else:
        data_res = pd.DataFrame(index=[f"{method}-{cluster}"],columns=["adt_n_features"],
                            data=[[len(adt_imp_features)]])

    if rna_imp_path and rna_path:
        rna_res =  cells_nzero_ratio(rna_raw_counts, metadata_unpaired, 'cell_type')+cells_nzero_ratio(rna_knn_smooth, metadata_unpaired, 'cell_type') + \
        cells_correlation(rna_raw_counts,rna_imputation,metadata_unpaired, 'cell_type') + \
        cells_correlation(rna_knn_smooth,rna_imputation,metadata_unpaired, 'cell_type' ) 
        
        rna_res2 = pd.DataFrame(rna_res).T
        rna_res2.columns = ["rna_raw_nzero_ratio_celltype","rna_raw_nzero_ratio_cell",
                                "rna_knn_nzero_ratio_celltype", "rna_knn_nzero_ratio_cell",
                                "rna_raw_celltype_cor", "rna_raw_cell_cor",
                                "rna_knn_celltype_cor", "rna_knn_cell_cor"
                                ]
        rna_res2.index=[f"{method}-{cluster}"]
    
    adt_res =  cells_nzero_ratio(adt_raw_counts, metadata_unpaired, 'cell_type')+ \
    cells_correlation(adt_raw_counts, adt_imputation,metadata_unpaired, 'cell_type')
    adt_res2 = pd.DataFrame(adt_res).T

    adt_res2.columns = ["adt_raw_nzero_ratio_celltype","adt_raw_nzero_ratio_cell",                           
                        "adt_raw_celltype_cor", "adt_raw_cell_cor"
                        ]
    adt_res2.index=[f"{method}-{cluster}"]
    
    if rna_imp_path and rna_path:
        out_metrics = pd.concat([data_res,rna_res2,adt_res2],axis=1)
    else:
        out_metrics = pd.concat([data_res,adt_res2],axis=1)
    if outf:
        # if 'imputation_adt' in adt_imp_path:
        #     out_metrics_path = adt_imp_path.split('imputation_adt')[0] + 'imputation-metrics.csv'
        # elif 'imputation-adt' in rna_imp_path:
        #     out_metrics_path = adt_imp_path.split('imputation-adt')[0] + 'imputation-metrics.csv'
        # else:
        #     out_metrics_path = adt_imp_path + 'imputation-metrics.csv'
        out_metrics.to_csv(outf)
    else:
        return(out_metrics)
    

def imputation_rna(metadata_path,rna_imp_path, rna_path, method, outf=None ):
    """
    Calculate metrics for single alogirthm results.
    metadata_path: path of cell metadata.
    rna_imp_path: path of rna imputation csv file.
    rna_path: path of rna raw h5ad file.
    outf: output to give  filename or pd.Dataframe() (None). Only for batch metrics calculation.
    For scMVAE only.
    """
    ## prepare input data and parameter
    # method=rna_imp_path.split("-")[-2]
    cluster="louvain"
    # Load metadata
    metadata = pd.read_csv(metadata_path, index_col='barcode', header=0, dtype='category')

    # RNA impute
    rna_imp = pd.read_csv(rna_imp_path, header=0, index_col=0)
    rna_imp_reordered = rna_imp.reindex(metadata.index)
    rna_imp_features = rna_imp_reordered.columns.values
    rna_imputation = rna_imp_reordered.to_numpy()

    # RNA raw
    rna = sc.read_h5ad(rna_path)
    rna_subset = rna[:, rna.var_names.isin(rna_imp_features)]
    rna_raw_counts = rna_subset.X.toarray() # raw counts        
    rna_knn_smooth = knn_smoothing(rna_raw_counts.transpose(), k=20, d=10, dither=0).transpose()


    data_res = pd.DataFrame(index=[f"{method}-{cluster}"],columns=["rna_n_features"],
                            data=[[len(rna_imp_features)]])

    rna_res =  cells_nzero_ratio(rna_raw_counts, metadata, 'cell_type')+cells_nzero_ratio(rna_knn_smooth, metadata, 'cell_type') + \
    cells_correlation(rna_raw_counts,rna_imputation,metadata, 'cell_type') + \
    cells_correlation(rna_knn_smooth,rna_imputation,metadata, 'cell_type' ) 
    
    rna_res2 = pd.DataFrame(rna_res).T
    rna_res2.columns = ["rna_raw_nzero_ratio_celltype","rna_raw_nzero_ratio_cell",
                            "rna_knn_nzero_ratio_celltype", "rna_knn_nzero_ratio_cell",
                            "rna_raw_celltype_cor", "rna_raw_cell_cor",
                            "rna_knn_celltype_cor", "rna_knn_cell_cor"
                            ]
    rna_res2.index=[f"{method}-{cluster}"]

    out_metrics = pd.concat([data_res, rna_res2],axis=1)
    if outf:
        out_metrics.to_csv(outf)
    else:
        return(out_metrics)
    
def imputation_stabmap(metadata_path,rna_imp_path, rna_path,atac_imp_path, atac_path,method,paired="s1d1",unpaired="s3d10",batch="batch", outf=None ):
    """
    for stabmap knn imputation only.
    did not knn-smoothing for scATAC
    """
    cluster="louvain"
    metadata = pd.read_csv(metadata_path, index_col='barcode', header=0, dtype='category')
    metadata_paired = metadata[metadata[batch] == paired]
    metadata_unpaired = metadata[metadata[batch] == unpaired]

    # RNA impute
    # rna_imp = pd.read_csv(rna_imp_path, header=0, index_col=0)
    rna_imp = dt.fread(rna_imp_path)
    rna_imp = rna_imp.to_pandas()
    rna_imp.index = rna_imp.iloc[:,0]
    rna_imp.drop(rna_imp.columns[[0]],axis=1,inplace=True)

    rna_imp.index = rna_imp.index.str.replace('_atac$', '', regex=True)
    rna_imp_reordered = rna_imp.reindex(metadata_unpaired.index)
    rna_imp_features = rna_imp_reordered.columns.values
    # rna_imputation = rna_imp_reordered.to_numpy()

    # RNA raw
    rna = sc.read_h5ad(rna_path)
    rna_unpaired = rna[rna.obs[batch] == unpaired, :]
    rna_unpaired_subset = rna_unpaired[:, rna_unpaired.var_names.isin(rna_imp_features)]
    rna_raw_counts = rna_unpaired_subset.X.toarray() # raw counts        
    
    rna_imputation = rna_imp_reordered.loc[:,rna_unpaired_subset.var_names].to_numpy()

    rna_knn_smooth = knn_smoothing(rna_raw_counts.transpose(), k=20, d=10, dither=0).transpose()


    # ATAC impute
    # atac_imp = pd.read_csv(atac_imp_path, header=0, index_col=0) # too slow
    atac_imp = dt.fread(atac_imp_path)
    atac_imp = atac_imp.to_pandas()
    atac_imp.index = atac_imp.iloc[:,0]
    atac_imp.drop(atac_imp.columns[[0]],axis=1,inplace=True)

    atac_imp.index = atac_imp.index.str.replace('_rna$', '', regex=True)
    atac_imp_reordered = atac_imp.reindex(metadata_unpaired.index)
    atac_imp_features = atac_imp_reordered.columns.values
    # atac_imputation = atac_imp_reordered.to_numpy()

    # ATAC raw
 
    atac = sc.read_h5ad(atac_path)
    atac_unpaired = atac[atac.obs[batch] == unpaired, :]
    atac_unpaired_subset = atac_unpaired[:, atac_unpaired.var_names.isin(atac_imp_features)]
    atac_raw_counts = atac_unpaired_subset.X.toarray() # raw counts
    # atac_subset = atac[:, atac.var_names.isin(atac_imp_features)]
    # atac_raw_counts = atac_subset.X.toarray() # raw counts

    atac_imputation = atac_imp_reordered.loc[:,atac_unpaired_subset.var_names].to_numpy()
    # knn_smooth
    # atac_knn_smooth = knn_smoothing(atac_raw_counts.transpose(), k=20, d=10, dither=0).transpose()

    # print("KNN smooth complete.")
    data_res = pd.DataFrame(index=[f"{method}-{cluster}"],columns=["rna_n_features","atac_n_features"],
                            data=[[len(rna_imp_features),len(atac_imp_features)]])

    print(rna_raw_counts.shape)
    print(rna_imputation.shape)
    print(rna_knn_smooth.shape)
    print(metadata_unpaired.shape)
    rna_res =  cells_nzero_ratio(rna_raw_counts, metadata_unpaired, 'cell_type')+cells_nzero_ratio(rna_knn_smooth, metadata_unpaired, 'cell_type') + \
        cells_correlation(rna_raw_counts,rna_imputation,metadata_unpaired, 'cell_type') + \
        cells_correlation(rna_knn_smooth,rna_imputation,metadata_unpaired, 'cell_type' ) 
        
    print("RNA metrics complete.")

    rna_res2 = pd.DataFrame(rna_res).T
    rna_res2.columns = ["rna_raw_nzero_ratio_celltype","rna_raw_nzero_ratio_cell",
                            "rna_knn_nzero_ratio_celltype", "rna_knn_nzero_ratio_cell",
                            "rna_raw_celltype_cor", "rna_raw_cell_cor",
                            "rna_knn_celltype_cor", "rna_knn_cell_cor"
                        ]
    rna_res2.index=[f"{method}-{cluster}"]

    atac_res =  cells_nzero_ratio(atac_raw_counts, metadata_unpaired, 'cell_type')+ \
        cells_auc(atac_raw_counts, atac_imputation,metadata_unpaired, None) 

    print("ATAC metrics complete.")
    atac_res2 = pd.DataFrame(atac_res).T

    atac_res2.columns = ["atac_raw_nzero_ratio_celltype","atac_raw_nzero_ratio_cell",
                         "atac_raw_celltype_auprc", "atac_raw_celltype_auroc", "atac_raw_cell_auprc", "atac_raw_cell_auroc"
                            ]
    atac_res2.index=[f"{method}-{cluster}"]
    out_metrics = pd.concat([data_res, rna_res2,atac_res2],axis=1)
    if outf:
        # if 'imputation_rna' in rna_imp_path:
        #     out_metrics_path = rna_imp_path.split('imputation_rna')[0] + 'imputation-metrics.csv'
        # elif 'imputation-rna' in rna_imp_path:
        #     out_metrics_path = rna_imp_path.split('imputation-rna')[0] + 'imputation-metrics.csv'
        # else:
        #     out_metrics_path = rna_imp_path + 'imputation-metrics.csv'
        print(f"output to {outf}")
        out_metrics.to_csv(outf)

    else:
        print("output to stdout")
        return(out_metrics)
    


def mouse_brain_divide(func, adatas,
                        method,
                        cluster,
                        batch = 'batch',
                        label = 'cell_type',
                        outf = None
                        ):
    """
    10X mouse brain datasets calculated batch removal metrics for 
    WT(wild type) and AD brain separately, then took the average
     of metrics from two disease status group.
    """
    print("10x mouse data")
    out_dat1 = func(adatas[0], method, cluster,batch, label, None)
    out_dat2 = func(adatas[1], method, cluster,batch, label, None) 
    out_dat = out_dat1.copy()
    avg_col = out_dat1.columns.to_list()[3:]
    out_dat[avg_col] = (out_dat1[avg_col]+ out_dat2[avg_col])/2
    out_dat[["nCell"]] = out_dat1[["nCell"]] + out_dat2[["nCell"]]
    if outf:
        print(f"writing to {outf}")
        out_dat.to_csv(outf)
    else:    
        return(outf)

from scmmib.metrics import paired_graph_metrics,  mosaic_latent_metrics,\
    unpaired_latent_metrics
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import os

# we first provided a demo for evaluate unpair scRNA and scATAC diagonal integration output.
def test_unpair():
    # 1. load the latent files and metadata files, all demo files are deposited in the github folder
    metadata = "../manuscript_figure_script_and_data/stage2_res/metadata/SHARE_RNA+ATAC_raw_metadata.csv.gz"
    meta = pd.read_csv(metadata, index_col='barcode', header=0, dtype='category')
    myfiles = ["../test/SHARE-multiome-raw-scRNA+scATAC-bindSC-ATAC-latent.csv.gz",
               "../test/SHARE-multiome-raw-scRNA+scATAC-bindSC-RNA-latent.csv.gz"]
    latent_atac = pd.read_csv(myfiles[0], index_col=0, header=0)
    latent_rna = pd.read_csv(myfiles[1], index_col=0, header=0)
    latent_atac_reindex=latent_atac.reindex(meta.index)
    latent_rna_reindex= latent_rna.reindex(meta.index)
    # 2. generate the ann data format and calculate the metrics
    adata_unpaired = sc.AnnData(latent_rna_reindex, obs=meta, dtype='float32')
    adata_unpaired.obsm['RNA'] = latent_rna_reindex
    adata_unpaired.obsm['ATAC'] = latent_atac_reindex
    unpaired_latent_metrics(adata_unpaired, method = "bindSC", cluster = 'louvain', batch = None, label = 'cell_type', mods = ["RNA","ATAC"], outf=None, embed_acc=True) # outf=None,return stdout, or return the path or "outf" param
    # embed_acc determine wheter calculate the accuracy metrics for embed in "mods" params

# Then we provided a demo for method with graph output rather than embedding.
def test_graph():
    # 1. load the graph files and metadata files, all demo files are deposited in the github folde
    metadata = "../manuscript_figure_script_and_data/stage2_res/metadata/BMMC_RNA+ADT_p10_metadata.csv.gz"
    latent = pd.read_csv("../test/BMMC-CITE_seq-p10-CITE_seq-SeuratV4-multi-graph.csv.gz", index_col=0, header=0)
    meta = pd.read_csv(metadata, index_col='barcode', header=0, dtype='category')
    latent_reordered = latent.reindex(meta.index)
    
    # 2. generate the ann data format and calculate the metrics
    adata = sc.AnnData(latent_reordered, obs=meta, dtype='float32')
    method = "SeuratV4"
    adata.obsp[method] = latent_reordered
    paired_graph_metrics(adata, method = "bindSC", cluster = 'louvain', batch = 'batch', label = 'cell_type', mods = ["RNA","ATAC"], outf=None)

# Finally we showed a demo for Seurat v5 bridge mosaic scRNA and ADT integration.
def test_mosaic():
    # 1. load Seurat v5 generated latent embeddings of paired, unpaired RNA and unpaired ADT, as well as metadata
    metadata = "../manuscript_figure_script_and_data/stage2_res/metadata/BMMC_RNA+ADT_s2d1_s3d6_metadata.csv.csv.gz"
    myfiles = ["../test/BMMC-CITE_seq-s2d1_s3d6-scRNA+ADT-SeuratV5_adt_reduc_latent.csv","../test/BMMC-CITE_seq-s2d1_s3d6-scRNA+ADT-SeuratV5_multi_lap_latent.csv","../test/BMMC-CITE_seq-s2d1_s3d6-scRNA+ADT-SeuratV5_rna_reduc_latent.csv"]
    meta = pd.read_csv(metadata, index_col='barcode', header=0, dtype='category')
    paired="s2d1"
    unpaired="s3d6"
    batch="batch"
    # 2. match paired and unpaired cell information in metadata, and input the latents and metadatas for metrics evaluation
    # pair_cells = meta[meta[batch] == paired].index
    # unpair_cells = meta[meta[batch] == unpaired].index
    metadatas = [meta[meta[batch] == paired],meta[meta[batch] == unpaired],meta[meta[batch] == unpaired]]
    latent_pair=pd.read_csv(myfiles[1], index_col=0, header=0)
    latent_rna = pd.read_csv(myfiles[2], index_col=0, header=0)
    latent_mod2 = pd.read_csv(myfiles[0], index_col=0, header=0)
    latents = [latent_pair, latent_rna, latent_mod2]
    mosaic_latent_metrics(latents=latents,metadatas=metadatas,paired="s2d1", unpaired="s3d6", mod2="adt", batch="batch",label="cell_type",latent_path=myfiles[0], method='SeuratV5', writef=True)
    

test_unpair()
# test_graph()
# test_mosaic()

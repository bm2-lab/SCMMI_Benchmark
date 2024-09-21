from scmmib.metrics import paired_latent_metrics, paired_graph_metrics, scmmib_gc,\
    mosaic_cnk_latent_metrics, mosaic_latent_metrics,imputation_pair_rna_atac, \
    imputation_mosaic_rna_atac, imputation_mosaic_rna_adt, imputation_rna, \
    unpaired_latent_metrics, mouse_brain_divide, imputation_stabmap
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import os
from multiprocessing import Pool
import glob

def run_unpair():
    my_files = [i.rstrip() for i in open("/home/shaliu_fu/software/scmmib/test/unpair_rna_atac_append.txt","r")]
    pool = Pool(processes=8)
    result = []

    for i in range(0,len(my_files),2):
    # for i in range(0,4,2):  
        myf = my_files[i]
        outfile = myf.split('ATAC-latent')[0] + 'unpaired-metrics.csv'
        # if os.path.exists(outfile):
        #     continue
        batch="batch"
        label = "cell_type"
        method = myf.split('-')[-3]
        if "10x_mouse_brain" in myf:
            metadata="/home/wsg/BM/data/10x_mouse_brain/RNA+ATAC/2p5/metadata.csv"
            batch = "rep"
        elif "10x_PBMC" in myf:
            metadata="/home/wsg/BM/data/10x_PBMC/RNA+ATAC/RawData/metadata.csv"            
        elif "BMMC" in myf:
            
            label = 'cell_type'
            if "p10" in myf:
                metadata="/home/wsg/BM/data/BMMC/RNA+ATAC/p10/metadata.csv"
            elif "s3d10" in myf:
                metadata="/home/wsg/BM/data/BMMC/RNA+ATAC/s3d10/metadata.csv"
                batch=None
            else:
                print(myf)
                print("error")
        elif "HSPC" in myf:
            metadata="/home/wsg/BM/data/HSPC/RNA+ATAC/p10/metadata.csv"
        elif "SHARE" in myf:
            metadata="/home/wsg/BM/data/SHARE/RNA+ATAC/RawData/metadata.csv"
            batch=None
            label = 'cell_type'
        else:
            print(myf)
            print("error")

        meta = pd.read_csv(metadata, index_col='barcode', header=0, dtype='category')
        
        latent_atac = pd.read_csv(my_files[i], index_col=0, header=0)
        latent_rna = pd.read_csv(my_files[i+1], index_col=0, header=0)
        latent_atac_reindex=latent_atac.reindex(meta.index)
        latent_rna_reindex= latent_rna.reindex(meta.index)
        adata_unpaired = sc.AnnData(latent_rna_reindex, obs=meta, dtype='float32')
        adata_unpaired.obsm['RNA'] = latent_rna_reindex
        adata_unpaired.obsm['ATAC'] = latent_atac_reindex
        # unpaired_latent_metrics(adata_unpaired, method = method, cluster = 'louvain', batch = batch, label = label, mods = ["RNA","ATAC"], outf=outfile, embed_acc=True)
        result.append(pool.apply_async(unpaired_latent_metrics,(adata_unpaired, method, 'louvain', batch,label,["RNA","ATAC"], outfile, True, )))
    pool.close()
    pool.join()


def run_pair():
    # my_files = [i.rstrip() for i in open("/home/shaliu_fu/software/scmmib/test/pair_rna_atac_append.txt","r")]
    my_files = [i.rstrip() for i in open("/home/shaliu_fu/software/scmmib/test/pair_rna_adt_append.txt","r")]
    pool = Pool(processes=5)
    result = []
    for myf in my_files:
        outfile = myf.split('multi-latent')[0] + 'paired-metrics.csv'
        # if os.path.exists(outfile):
        #     continue
        batch="batch"
        label = None
        if "scATAC" in myf:
            if "10x_mouse_brain" in myf:
                metadata="/home/wsg/BM/data/10x_mouse_brain/RNA+ATAC/2p5/metadata.csv"
                batch = "rep"
            elif "10x_PBMC" in myf:
                metadata="/home/wsg/BM/data/10x_PBMC/RNA+ATAC/RawData/metadata.csv"            
            elif "BMMC" in myf:
                
                label = 'cell_type'
                if "p10" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ATAC/p10/metadata.csv"
                elif "s3d10" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ATAC/s3d10/metadata.csv"
                    batch=None
                else:
                    print(myf)
                    print("error")
            elif "HSPC" in myf:
                metadata="/home/wsg/BM/data/HSPC/RNA+ATAC/p10/metadata.csv"
            elif "SHARE" in myf:
                metadata="/home/wsg/BM/data/SHARE/RNA+ATAC/RawData/metadata.csv"
                batch=None
                label = 'cell_type'
            else:
                print(myf)
                print("error")
        elif "ADT" in myf:
            if "BMMC" in myf:
                label = 'cell_type'
                if "p10" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ADT/p10/metadata.csv"
                elif "s2d1" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ADT/s2d1/metadata.csv"
                    batch = None
                else:
                    print(myf)
                    print("error")
            elif "10x_kidney_cancer" in myf:
                metadata="/home/wsg/BM/data/10x_kidney_cancer/RNA+ADT/RawData/metadata.csv"

            elif "10x_NSCLC" in myf:    
                metadata="/home/wsg/BM/data/10x_NSCLC/RNA+ADT/RawData/metadata.csv"
                batch = "donor"
            elif "HSPC" in myf:
                metadata="/home/wsg/BM/data/HSPC/RNA+ADT/p10/metadata.csv"
            elif "lymph_node_A1" in myf:
                metadata="/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node_A1/metadata.csv"
                label="cell_type"
            elif "lymph_node" in myf:
                metadata="/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node/metadata.csv"
            elif "spleen" in myf:
                metadata="/home/wsg/BM/data/SPATIAL/RNA+ADT/spleen/metadata.csv"
            elif "thymus" in myf:
                metadata="/home/wsg/BM/data/SPATIAL/RNA+ADT/thymus/metadata.csv"
            else:
                print(myf)
                print("unknown dataset!!!")
        else:
            print(myf)
            print("unknown dataset!!!")
        if not os.path.exists(metadata):
            print("{0} not exist!".format(metadata))
        
        else:
            latent = pd.read_csv(myf, index_col=0, header=0)
            meta = pd.read_csv(metadata, index_col='barcode', header=0, dtype='category')
            latent_reordered = latent.reindex(meta.index)
            adata = sc.AnnData(latent_reordered, obs=meta, dtype='float32')
            method = myf.split('-')[-3]
            if metadata == "/home/wsg/BM/data/10x_mouse_brain/RNA+ATAC/2p5/metadata.csv":
                adata1 = adata[adata.obs['type']=="WT"].copy()
                adata2 = adata[adata.obs['type']=="AD"].copy()
                # result.append(pool.apply_async(mouse_brain_divide,(paired_latent_metrics,[adata1,adata2], method, 'louvain', batch,label, outfile,)))
                mouse_brain_divide(paired_latent_metrics,[adata1,adata2], method, 'louvain', batch,label, outfile)
            else:
                # paired_latent_metrics(adata_all, method = method, cluster = 'louvain', batch = batch, label = label,outf=paired_metrics_path)

                result.append(pool.apply_async(paired_latent_metrics,(adata, method, 'louvain', batch,label, outfile,)))
    pool.close()
    pool.join()

def run_graph():
    my_files = [i.rstrip() for i in open("/home/shaliu_fu/software/scmmib/test/graph_append.txt","r")]
    pool = Pool(processes=4)
    result = []
    for myf in my_files:
        outfile = myf.split('multi-graph')[0] + 'new-paired-metrics.csv'
        # if os.path.exists(outfile):
        #     continue
        batch="batch"
        label = None
        if "scATAC" in myf:
            if "10x_mouse_brain" in myf:
                metadata="/home/wsg/BM/data/10x_mouse_brain/RNA+ATAC/2p5/metadata.csv"
                batch = "rep"
            elif "10x_PBMC" in myf:
                metadata="/home/wsg/BM/data/10x_PBMC/RNA+ATAC/RawData/metadata.csv"            
            elif "BMMC" in myf:
                
                label = 'cell_type'
                if "p10" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ATAC/p10/metadata.csv"
                elif "s3d10" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ATAC/s3d10/metadata.csv"
                    batch=None
                else:
                    print(myf)
                    print("error")
            elif "HSPC" in myf:
                metadata="/home/wsg/BM/data/HSPC/RNA+ATAC/p10/metadata.csv"
            elif "SHARE" in myf:
                metadata="/home/wsg/BM/data/SHARE/RNA+ATAC/RawData/metadata.csv"
                batch=None
                label = 'cell_type'
            else:
                print(myf)
                print("error")
        elif "ADT" in myf:
            if "BMMC" in myf:
                label = 'cell_type'
                if "p10" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ADT/p10/metadata.csv"
                elif "s2d1" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ADT/s2d1/metadata.csv"
                    batch = None
                else:
                    print(myf)
                    print("error")
            elif "10x_kidney_cancer" in myf:
                metadata="/home/wsg/BM/data/10x_kidney_cancer/RNA+ADT/RawData/metadata.csv"

            elif "10x_NSCLC" in myf:    
                metadata="/home/wsg/BM/data/10x_NSCLC/RNA+ADT/RawData/metadata.csv"
                batch = "donor"
            elif "HSPC" in myf:
                metadata="/home/wsg/BM/data/HSPC/RNA+ADT/p10/metadata.csv"
            elif "lymph_node_A1" in myf:
                metadata="/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node_A1/metadata.csv"
                label="cell_type"
            elif "lymph_node" in myf:
                metadata="/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node/metadata.csv"
            elif "spleen" in myf:
                metadata="/home/wsg/BM/data/SPATIAL/RNA+ADT/spleen/metadata.csv"
            elif "thymus" in myf:
                metadata="/home/wsg/BM/data/SPATIAL/RNA+ADT/thymus/metadata.csv"
            else:
                print(myf)
                print("unknown dataset!!!")
        else:
            print(myf)
            print("unknown dataset!!!")
        if not os.path.exists(metadata):
            print("{0} not exist!".format(metadata))
        else:
            latent = pd.read_csv(myf, index_col=0, header=0)
            meta = pd.read_csv(metadata, index_col='barcode', header=0, dtype='category')
            latent_reordered = latent.reindex(meta.index)
            adata = sc.AnnData(latent_reordered, obs=meta, dtype='float32')
            method = myf.split('-')[-3]
            adata.obsp[method] = latent_reordered
            outfile = myf.split('multi-graph')[0] + 'new-paired-metrics.csv'
            print(outfile)  
            # paired_graph_metrics(adata, method = method, cluster = 'louvain', batch = batch, label = label, \
            #                               outf=outfile)
            if metadata == "/home/wsg/BM/data/10x_mouse_brain/RNA+ATAC/2p5/metadata.csv":
                adata1 = adata[adata.obs['type']=="WT"].copy()
                adata2 = adata[adata.obs['type']=="AD"].copy()
                result.append(pool.apply_async(mouse_brain_divide,(paired_graph_metrics,[adata1,adata2], method, 'louvain', batch,label, outfile,)))
                # mouse_brain_divide(paired_graph_metrics,[adata1,adata2], method, 'louvain', batch,label, outfile)
            else:
                result.append(pool.apply_async(paired_graph_metrics,(adata, method, 'louvain', batch,label, outfile,)))

    pool.close()
    pool.join()


def run_cnk():
    # my_files1: single file
    # my_files2: separate latent files.
    my_files1 = [i.rstrip() for i in open("/home/shaliu_fu/software/scmmib/test/cnk_append1.txt","r")]
    # my_files1 = [i.rstrip() for i in open("/home/shaliu_fu/software/scmmib/test/cnk_append1.txt","r")]
    my_files2 = []
    # my_files2 = [i.rstrip() for i in open("/home/shaliu_fu/software/scmmib/test/cnk_append2.txt","r")]
    # my_files2 = [i.rstrip() for i in open("/home/shaliu_fu/software/scmmib/test/mosaic_test.txt","r")]
    pool = Pool(processes=5)
    result = []
    for myf in my_files1:
    # for myf in ["/home/wsg/BM/results/task/mosaic_scRNA+ADT/robustness/BMMC/c5k_cnk/c5k_c5k/MIDAS/RUN_1/MIDAS_BMMC-CITE_seq-c5k_c5k-scRNA+ADT_latent.csv"]:
        # outfile = myf.split('latent')[0] + 'new-paired-metrics.csv'
        outfile = myf.split('latent')[0] + 'paired-metrics.csv'
        if os.path.exists(outfile):
            continue
        batch="data_size"
        label = 'cell_type'
        unpaired = "c5k"
        if "scATAC" in myf:
            mod2 = "atac"
            if "BMMC" in myf:
                if "c5k_c3k" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ATAC/c5k_c3k/metadata.csv"
                    paired = "c3k"
                elif "c5k_c5k" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ATAC/c5k_c5k/metadata.csv"
                    unpaired = "c5k_1"
                    paired = "c5k_2"
                elif "c5k_c10k" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ATAC/c5k_c10k/metadata.csv"
                    paired = "c10k"
                elif "c5k_c20k" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ATAC/c5k_c20k/metadata.csv"
                    paired = "c20k"
                elif "s1d1_s3d10" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ATAC/s1d1_s3d10/metadata.csv"
                    paired = "s1d1"
                    unpaired = "s3d10"
                    batch = "batch"
                else:
                    print(myf)
                    print("Unknown dataset!!")
            elif "SHARE" in myf:
                if "c5k_c3k" in myf:
                    metadata="/home/wsg/BM/data/SHARE/RNA+ATAC/c5k_c3k/metadata.csv"
                    paired = "c3k"
                elif "c5k_c5k" in myf:
                    metadata="/home/wsg/BM/data/SHARE/RNA+ATAC/c5k_c5k/metadata.csv"
                    unpaired = "c5k_1"
                    paired = "c5k_2"
                elif "c5k_c10k" in myf:
                    metadata="/home/wsg/BM/data/SHARE/RNA+ATAC/c5k_c10k/metadata.csv"
                    paired = "c10k"
                elif "c5k_c20k" in myf:
                    metadata="/home/wsg/BM/data/SHARE/RNA+ATAC/c5k_c20k/metadata.csv"
                    paired = "c20k"
                else:
                    print(myf)
                    print("Unknown data size!!")
            else:
                print(myf)
                print("Unknown datasets!!")

        elif "ADT" in myf:
            mod2 = "adt"
            if "BMMC" in myf:
                if "c5k_c3k" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ADT/c5k_c3k/metadata.csv"
                    paired = "c3k"
                elif "c5k_c5k" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ADT/c5k_c5k/metadata.csv"
                    unpaired = "c5k_1"
                    paired = "c5k_2"
                elif "c5k_c10k" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ADT/c5k_c10k/metadata.csv"
                    paired = "c10k"
                elif "c5k_c20k" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ADT/c5k_c20k/metadata.csv"
                    paired = "c20k"
                elif "s2d1_s3d6" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ADT/s2d1_s3d6/metadata.csv"
                    paired = "s3d6"
                    unpaired = "s2d1"
                    batch = "batch"
                else:
                    print(myf)
                    print("Unknown data size!!")
            else:
                print(myf)
                print("Unknown datasets!!")
        else:
            print(myf)
            print("unknown modality!!!")
        print(myf)
        print("start processing!")
        # method = myf.split('/')[-3]
        method = myf.split('/')[-1].split("_")[0]
        meta = pd.read_csv(metadata, index_col='barcode', header=0, dtype='category')
        pair_cells = meta[meta[batch] == paired].index
        unpair_cells = meta[meta[batch] == unpaired].index
        
        latent = pd.read_csv(myf, index_col=0, header=0)
        if method!="sciPENN" and method !="totalVI":
            metadatas = [meta[meta[batch] == paired],meta[meta[batch] == unpaired],meta[meta[batch] == unpaired]]
            latent_pair=latent.loc[pair_cells,:].copy()
            latent_rna = latent.loc[[x + '_rna' for x in unpair_cells],:].copy()
            latent_mod2 = latent.loc[[x + f'_{mod2}' for x in unpair_cells],:].copy()
            latents = [latent_pair, latent_rna, latent_mod2]
        else:
            metadatas = [meta[meta[batch] == paired],meta[meta[batch] == unpaired]]
            latent_pair=latent.loc[pair_cells,:].copy()
            latent_rna = latent.loc[[x + '_rna' for x in unpair_cells],:].copy()
            latents = [latent_pair, latent_rna]
        # mosaic_cnk_latent_metrics(latents,metadatas,paired, unpaired, 
                            # mod2, batch,label,myf,method,True)
        result.append(pool.apply_async(mosaic_cnk_latent_metrics,(latents,metadatas,paired, unpaired, 
                            mod2, batch,label,myf,method,True,)))
    
   
    for i in range(0,len(my_files2),3):
        myfiles2_reorder = [my_files2[i],my_files2[i+1],my_files2[i+2]]
        myf = myfiles2_reorder[2]
        outfile = myf.split('latent')[0] + 'paired-metrics.csv'
        # if os.path.exists(outfile):
        #     continue
        batch="data_size"
        label = "cell_type"
        unpaired = "c5k"
        if "scATAC" in myf:
            mod2 = "atac"
            if "BMMC" in myf:
                if "c5k_c3k" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ATAC/c5k_c3k/metadata.csv"
                    paired = "c3k"
                elif "c5k_c5k" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ATAC/c5k_c5k/metadata.csv"
                    unpaired = "c5k_1"
                    paired = "c5k_2"
                elif "c5k_c10k" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ATAC/c5k_c10k/metadata.csv"
                    paired = "c10k"
                elif "c5k_c20k" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ATAC/c5k_c20k/metadata.csv"
                    paired = "c20k"
                elif "s1d1_s3d10" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ATAC/s1d1_s3d10/metadata.csv"
                    paired = "s1d1"
                    unpaired = "s3d10"
                    batch = "batch"
                else:
                    print(myf)
                    print("Unknown dataset!!")
            elif "SHARE" in myf:
                if "c5k_c3k" in myf:
                    metadata="/home/wsg/BM/data/SHARE/RNA+ATAC/c5k_c3k/metadata.csv"
                    paired = "c3k"
                elif "c5k_c5k" in myf:
                    metadata="/home/wsg/BM/data/SHARE/RNA+ATAC/c5k_c5k/metadata.csv"
                    unpaired = "c5k_1"
                    paired = "c5k_2"
                elif "c5k_c10k" in myf:
                    metadata="/home/wsg/BM/data/SHARE/RNA+ATAC/c5k_c10k/metadata.csv"
                    paired = "c10k"
                elif "c5k_c20k" in myf:
                    metadata="/home/wsg/BM/data/SHARE/RNA+ATAC/c5k_c20k/metadata.csv"
                    paired = "c20k"
                else:
                    print(myf)
                    print("Unknown data size!!")
            else:
                print(myf)
                print("Unknown datasets!!")

        elif "ADT" in myf:
            mod2 = "adt"
            if "BMMC" in myf:
                if "c5k_c3k" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ADT/c5k_c3k/metadata.csv"
                    paired = "c3k"
                elif "c5k_c5k" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ADT/c5k_c5k/metadata.csv"
                    unpaired = "c5k_1"
                    paired = "c5k_2"
                elif "c5k_c10k" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ADT/c5k_c10k/metadata.csv"
                    paired = "c10k"
                elif "c5k_c20k" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ADT/c5k_c20k/metadata.csv"
                    paired = "c20k"
                elif "s2d1_s3d6" in myf:
                    metadata="/home/wsg/BM/data/BMMC/RNA+ADT/s2d1_s3d6/metadata.csv"
                    paired = "s3d6"
                    unpaired = "s2d1"
                    batch = "batch"
                else:
                    print(myf)
                    print("Unknown data size!!")
            else:
                print(myf)
                print("Unknown datasets!!")
        else:
            print(myf)
            print("unknown modality!!!")
        print(myf)
        print("start processing!")
        # method = myf.split('/')[-3]
        method = myf.split('-')[-1].split("_")[0]

        meta = pd.read_csv(metadata, index_col='barcode', header=0, dtype='category')
        pair_cells = meta[meta[batch] == paired].index
        unpair_cells = meta[meta[batch] == unpaired].index
        metadatas = [meta[meta[batch] == paired],meta[meta[batch] == unpaired],meta[meta[batch] == unpaired]]
        # latent = 

        latent_pair=pd.read_csv(myfiles2_reorder[1], index_col=0, header=0)
        latent_rna = pd.read_csv(myfiles2_reorder[2], index_col=0, header=0)
        latent_mod2 = pd.read_csv(myfiles2_reorder[0], index_col=0, header=0)
        latents = [latent_pair, latent_rna, latent_mod2]
        # mosaic_cnk_latent_metrics(latents,metadatas,paired, unpaired, 
                            # mod2, batch,label,myf,method,True)
        # result.append(pool.apply_async(mosaic_cnk_latent_metrics,(latents,metadatas,paired, unpaired, mod2, batch,label,myf,method,True,)))

    pool.close()
    pool.join()

def run_impute():
    # pair_file = [i.rstrip() for i in open("/home/shaliu_fu/software/scmmib/test/impute_pair_append.txt","r")]
    pair_file = []
    mosaic_rna_adt_file = [i.rstrip() for i in open("/home/shaliu_fu/software/scmmib/test/impute_mosaic_rna_adt.txt","r")]
    mosaic_rna_atac_file = [i.rstrip() for i in open("/home/shaliu_fu/software/scmmib/test/impute_mosaic_rna_atac.txt","r")]
    pool = Pool(processes=5)
    result = []
    for i in range(0,len(pair_file),2):
        rna_imp_path = pair_file[i+1]
        atac_imp_path = pair_file[i]
        if "BMMC" in rna_imp_path:
            data_path = "/home/wsg/BM/data/BMMC/RNA+ATAC/p10"
            metadata_path = glob.glob(data_path + '/metadata.csv')[0]
            rna_path = glob.glob(data_path + '/*-RNA-counts.h5ad')[0]
            atac_path = glob.glob(data_path + '/*-ATAC-peaks.h5ad')[0]
            method = rna_imp_path.split("-")[-3]
            outf = rna_imp_path.split('imputation-rna')[0] + 'imputation-metrics.csv'
            # result.append(pool.apply_async(imputation_pair_rna_atac,(metadata_path,rna_imp_path, rna_path,atac_imp_path, atac_path,method,outf, )))
            # imputation_pair_rna_atac(metadata_path,rna_imp_path, rna_path,atac_imp_path, atac_path,method, outf=outf )


        elif "SHARE" in rna_imp_path:
            data_path = "/home/wsg/BM/data/SHARE/RNA+ATAC/RawData"
            metadata_path = glob.glob(data_path + '/metadata.csv')[0]
            rna_path = glob.glob(data_path + '/*-RNA-counts.h5ad')[0]
            atac_path = glob.glob(data_path + '/*-ATAC-peaks.h5ad')[0]
            method = rna_imp_path.split("-")[-3]
            outf = rna_imp_path.split('imputation-rna')[0] + 'imputation-metrics.csv'
            # result.append(pool.apply_async(imputation_pair_rna_atac,(metadata_path,rna_imp_path, rna_path,atac_imp_path, atac_path,method,outf, )))
            # imputation_pair_rna_atac(metadata_path,rna_imp_path, rna_path,atac_imp_path, atac_path,method, outf=outf )
        else:
            print("Unknown dataset!")
    full_data = []
    noadt_data = []
    for line in mosaic_rna_adt_file:
        if "sciPENN" in line or "totalVI" in line:
            noadt_data.append(line)
        else:
            full_data.append(line)
    for myf in noadt_data:
        adt_imp_path = myf
        data_path = "/home/wsg/BM/data/BMMC/RNA+ADT/s2d1_s3d6"    
        metadata_path = glob.glob(data_path + '/metadata.csv')[0]
        adt_path = glob.glob(data_path + '/*-ADT-counts.h5ad')[0]
        method = adt_imp_path.split("/")[-1].split("_")[0]
        # result.append(pool.apply_async(imputation_mosaic_rna_adt,(metadata_path,None, None,adt_imp_path, adt_path,method,"s3d6","s2d1","batch",True, )))
        # out_res=imputation_mosaic_rna_adt(metadata_path,rna_imp_path=None, rna_path=None,adt_imp_path=adt_imp_path, adt_path=adt_path,method=method, paired="s3d6",unpaired="s2d1", batch="batch",outf=True )
    for i in range(0,len(full_data),2):
        rna_imp_path = full_data[i+1]
        adt_imp_path = full_data[i]
        data_path = "/home/wsg/BM/data/BMMC/RNA+ADT/s2d1_s3d6"    
        metadata_path = glob.glob(data_path + '/metadata.csv')[0]
        rna_path = glob.glob(data_path + '/*-RNA-counts.h5ad')[0]
        adt_path = glob.glob(data_path + '/*-ADT-counts.h5ad')[0]
        method = rna_imp_path.split("/")[-1].split("_")[0]
        outf = rna_imp_path.split('imputation_rna')[0] + 'imputation-metrics.csv'
        # imputation_mosaic_rna_adt(metadata_path,rna_imp_path, rna_path,adt_imp_path, adt_path,method,"s3d6","s2d1","batch",outf)
        # result.append(pool.apply_async(imputation_mosaic_rna_adt,(metadata_path,rna_imp_path, rna_path,adt_imp_path, adt_path,method,"s3d6","s2d1","batch",outf, )))

    share_res = []
    bmmc_res = []
    for line in mosaic_rna_atac_file:
        if "BMMC" in line:
            bmmc_res.append(line)
        elif "SHARE" in line:
            share_res.append(line)
        else:
            print("Unknown mosaic scRNA+scATAC datasets!")
    
    # print(share_res)
    for i in range(0,len(bmmc_res),2):
        rna_imp_path = bmmc_res[i+1]
        atac_imp_path = bmmc_res[i]
        data_path = "/home/wsg/BM/data/BMMC/RNA+ATAC/s1d1_s3d10"    
        metadata_path = glob.glob(data_path + '/metadata.csv')[0]
        rna_path = glob.glob(data_path + '/*-RNA-counts.h5ad')[0]
        atac_path = glob.glob(data_path + '/*-ATAC-peaks.h5ad')[0]
        method =  rna_imp_path.split("/")[-1].split("_")[0]
        outf = rna_imp_path.split('imputation_rna')[0] + 'imputation-metrics.csv'
        # result.append(pool.apply_async(imputation_mosaic_rna_atac,(metadata_path,rna_imp_path, rna_path,atac_imp_path, atac_path,method,"s1d1","s3d10","batch",outf, )))
        result.append(pool.apply_async(imputation_stabmap,(metadata_path,rna_imp_path, rna_path,atac_imp_path, atac_path,method,"s1d1","s3d10","batch",outf, )))
        # imputation_stabmap(metadata_path,rna_imp_path, rna_path,atac_imp_path, atac_path,method=method,paired="s1d1",unpaired="s3d10", batch="batch",outf=outf )
        # imputation_mosaic_rna_atac(metadata_path,rna_imp_path, rna_path,atac_imp_path, atac_path,method=rna_imp_path.split("/")[-3],paired="s1d1",unpaired="s3d10", batch="batch",outf=outf )
    for i in range(0,len(share_res),2):
        rna_imp_path = share_res[i+1]
        atac_imp_path = share_res[i]
        data_path = "/home/wsg/BM/data/SHARE/RNA+ATAC/c5k_c5k"    
        metadata_path = glob.glob(data_path + '/metadata.csv')[0]
        rna_path = glob.glob(data_path + '/*-RNA-counts.h5ad')[0]
        atac_path = glob.glob(data_path + '/*-ATAC-peaks.h5ad')[0]
        method =  rna_imp_path.split("/")[-1].split("_")[0]
        outf = rna_imp_path.split('imputation_rna')[0] + 'imputation-metrics.csv'
        # result.append(pool.apply_async(imputation_mosaic_rna_atac,(metadata_path,rna_imp_path, rna_path,atac_imp_path, atac_path, method,"c5k_2","c5k_1","data_size",outf, )))
        result.append(pool.apply_async(imputation_stabmap,(metadata_path,rna_imp_path, rna_path,atac_imp_path, atac_path, method,"c5k_2","c5k_1","data_size",outf, )))
        # imputation_mosaic_rna_atac(metadata_path,rna_imp_path, rna_path,atac_imp_path, atac_path,method=method,paired="c5k_2",unpaired="c5k_1", batch="data_size",outf=outf )
        # imputation_stabmap(metadata_path,rna_imp_path, rna_path,atac_imp_path, atac_path,method=rna_imp_path.split("/")[-3],paired="c5k_2",unpaired="c5k_1", batch="data_size",outf=outf )

    # scmvae_file = [i.rstrip() for i in open("/home/shaliu_fu/software/scmmib/test/impute_debug.txt","r")]
    # for i in range(0,len(scmvae_file),2):
    #     rna_imp_path = scmvae_file[i+1]
    #     data_path = "/home/wsg/BM/data/BMMC/RNA+ATAC/p10"    
    #     metadata_path = glob.glob(data_path + '/metadata.csv')[0]
    #     rna_path = glob.glob(data_path + '/*-RNA-counts.h5ad')[0]
    #     method =  rna_imp_path.split("-")[-3]
    #     outf = rna_imp_path.split('imputation-rna')[0] + 'imputation-metrics.csv'
    #     # imputation_rna(metadata_path,rna_imp_path, rna_path,method, outf=True )
    #     result.append(pool.apply_async(imputation_rna,(metadata_path,rna_imp_path, rna_path,method,outf, )))
    pool.close()
    pool.join()

def run_test():
    my_files2 = [i.rstrip() for i in open("/home/shaliu_fu/software/scmmib/test/mosaic_test.txt","r")]
    for i in range(0,len(my_files2),3):
        myfiles2_reorder = [my_files2[i],my_files2[i+1],my_files2[i+2]]
        myf = myfiles2_reorder[2]
        outfile = myf.split('latent')[0] + 'paired-metrics.csv'
        # if os.path.exists(outfile):
        #     continue
        label = "cell_type"
        metadata="/home/wsg/BM/data/BMMC/RNA+ADT/s2d1_s3d6/metadata.csv"
        paired = "s3d6"
        unpaired = "s2d1"
        batch = "batch"
        mod2 = 'adt'
        method = myf.split('/')[-3]

        meta = pd.read_csv(metadata, index_col='barcode', header=0, dtype='category')
        pair_cells = meta[meta[batch] == paired].index
        unpair_cells = meta[meta[batch] == unpaired].index
        metadatas = [meta[meta[batch] == paired],meta[meta[batch] == unpaired],meta[meta[batch] == unpaired]]
        # latent = 

        latent_pair=pd.read_csv(myfiles2_reorder[1], index_col=0, header=0)
        latent_rna = pd.read_csv(myfiles2_reorder[2], index_col=0, header=0)
        latent_mod2 = pd.read_csv(myfiles2_reorder[0], index_col=0, header=0)
        latents = [latent_pair, latent_rna, latent_mod2]
        mosaic_cnk_latent_metrics(latents,metadatas,paired, unpaired, 
                            mod2, batch,label,myf,method,True)



def run_mosaic():
    # my_files1: single file
    # my_files2: separate latent files.

    # my_files1 = [i.rstrip() for i in open("/home/shaliu_fu/software/scmmib/test/cnk_append1.txt","r")]
    my_files1 = [i.rstrip() for i in open("/home/shaliu_fu/software/scmmib/test/mosaic_append.txt","r")]
    # my_files2 = [i.rstrip() for i in open("/home/shaliu_fu/software/scmmib/test/mosaic_seurat_append.txt","r")]
    my_files2 = []
    # my_files2 = [i.rstrip() for i in open("/home/shaliu_fu/software/scmmib/test/mosaic_test.txt","r")]
    pool = Pool(processes=8)
    result = []
    for myf in my_files1:
        outfile = myf.split('latent')[0] + 'paired-metrics.csv'
        if os.path.exists(outfile):
            continue
        batch= "data_size"
        label = 'cell_type'
        
        if "scATAC" in myf:
            mod2 = "atac"
            if "BMMC" in myf:
                metadata="/home/wsg/BM/data/BMMC/RNA+ATAC/s1d1_s3d10/metadata.csv"
                paired = "s1d1"
                unpaired = "s3d10"
                batch = "batch"
            elif "SHARE" in myf:
                metadata="/home/wsg/BM/data/SHARE/RNA+ATAC/c5k_c5k/metadata.csv"
                unpaired = "c5k_1"
                paired = "c5k_2"
               
            else:
                print(myf)
                print("Unknown datasets!!")

        elif "ADT" in myf:
            mod2 = "adt"
            if "BMMC" in myf:        
                metadata="/home/wsg/BM/data/BMMC/RNA+ADT/s2d1_s3d6/metadata.csv"
                paired = "s3d6"
                unpaired = "s2d1"
                batch = "batch"
            else:
                print(myf)
                print("Unknown datasets!!")
        else:
            print(myf)
            print("unknown modality!!!")
        print(myf)
        print("start processing!")
        # method = myf.split('/')[-3] # position of algorithm name
        method = myf.split('/')[-1].split("_")[0]

        meta = pd.read_csv(metadata, index_col='barcode', header=0, dtype='category')
        pair_cells = meta[meta[batch] == paired].index
        unpair_cells = meta[meta[batch] == unpaired].index
        
        latent = pd.read_csv(myf, index_col=0, header=0)
        if method!="sciPENN" and method !="totalVI":
            metadatas = [meta[meta[batch] == paired],meta[meta[batch] == unpaired],meta[meta[batch] == unpaired]]
            latent_pair=latent.loc[pair_cells,:].copy()
            latent_rna = latent.loc[[x + '_rna' for x in unpair_cells],:].copy()
            latent_mod2 = latent.loc[[x + f'_{mod2}' for x in unpair_cells],:].copy()
            latents = [latent_pair, latent_rna, latent_mod2]
        else:
            metadatas = [meta[meta[batch] == paired],meta[meta[batch] == unpaired]]
            latent_pair=latent.loc[pair_cells,:].copy()
            latent_rna = latent.loc[[x + '_rna' for x in unpair_cells],:].copy()
            latents = [latent_pair, latent_rna]
        # mosaic_latent_metrics(latents,metadatas,paired, unpaired, 
        #                     mod2, batch,label,myf,method,True)
        result.append(pool.apply_async(mosaic_latent_metrics,(latents,metadatas,paired, unpaired, 
                            mod2, batch,label,myf,method,True,)))
    
   
    for i in range(0,len(my_files2),3):
        myfiles2_reorder = [my_files2[i],my_files2[i+1],my_files2[i+2]]
        myf = myfiles2_reorder[2]
        outfile = myf.split('latent')[0] + 'paired-metrics.csv'
        # if os.path.exists(outfile):
        #     continue
        batch="data_size"
        label = 'cell_type'
        
        if "scATAC" in myf:
            mod2 = "atac"
            if "BMMC" in myf:
                metadata="/home/wsg/BM/data/BMMC/RNA+ATAC/s1d1_s3d10/metadata.csv"
                paired = "s1d1"
                unpaired = "s3d10"
                batch = "batch"
            elif "SHARE" in myf:
                metadata="/home/wsg/BM/data/SHARE/RNA+ATAC/c5k_c5k/metadata.csv"
                unpaired = "c5k_1"
                paired = "c5k_2"
               
            else:
                print(myf)
                print("Unknown datasets!!")

        elif "ADT" in myf:
            mod2 = "adt"
            if "BMMC" in myf:        
                metadata="/home/wsg/BM/data/BMMC/RNA+ADT/s2d1_s3d6/metadata.csv"
                paired = "s3d6"
                unpaired = "s2d1"
                batch = "batch"
            else:
                print(myf)
                print("Unknown datasets!!")
        else:
            print(myf)
            print("unknown modality!!!")
        print(myf)
        print("start processing!")
        method = myf.split('-')[-1].split("_")[0]

        meta = pd.read_csv(metadata, index_col='barcode', header=0, dtype='category')
        pair_cells = meta[meta[batch] == paired].index
        unpair_cells = meta[meta[batch] == unpaired].index
        metadatas = [meta[meta[batch] == paired],meta[meta[batch] == unpaired],meta[meta[batch] == unpaired]]
        latent_pair=pd.read_csv(myfiles2_reorder[1], index_col=0, header=0)
        latent_rna = pd.read_csv(myfiles2_reorder[2], index_col=0, header=0)
        latent_mod2 = pd.read_csv(myfiles2_reorder[0], index_col=0, header=0)
        latents = [latent_pair, latent_rna, latent_mod2]
        # mosaic_latent_metrics(latents,metadatas,paired, unpaired, 
        #                     mod2, batch,label,myf,method,True)
        # result.append(pool.apply_async(mosaic_latent_metrics,(latents,metadatas,paired, unpaired, mod2, batch,label,myf,method,True,)))

    pool.close()
    pool.join()


if __name__ == "__main__":
    # run_graph()
    run_cnk()
    # run_mosaic()
    # run_impute()
    # run_unpair()
    # run_pair()
    # run_test()
    print("Complete!")
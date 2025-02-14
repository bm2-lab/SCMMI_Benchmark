## User manual of SCMMIB package

### paired_latent_metrics
```python
## for all paired integration metrics except for graph output
paired_latent_metrics(adata,
                          method = 'bindSC',
                          cluster = 'louvain',
                          batch = 'batch',
                          label = 'cell_type',
                          outf = './text.txt'
                         )
```
- adata: `anndata` object for joint embedding. The `obs` is celll metadata. The input latent embedding stores in `adata.X`.<br>
- method: name of the input algorithm. <br>
- cluster: clustering methods, deprecated. <br>
- batch: column name of batch used. Do not consider batch effect removal if `batch=None`.<br>
- label: column name of cell type used. Do not consider biological conservation if `label=None`<br>
- outf: output metric filename. Stdout if `outf = None`. <br>
<br>

### paired_graph_metrics
```python
## for all paired integration metrics of graph output

paired_graph_metrics(adata,
                          method = 'bindSC',
                          cluster = 'louvain',
                          batch = 'batch',
                          label = 'cell_type',
                          outf = './text.txt'
                         )
```
- adata: `anndata` object for joint graph. The `obs` is cell metadata. The input graph stores in `adata.obsp[method]`<br>
- method: name of the input algorithm. <br>
- cluster: clustering methods, deprecated. <br>
- batch: column name of batch used. Do not consider batch effect removal if `batch=None`.<br>
- label: column name of cell type used. Do not consider biological conservation if `label=None`<br>
- outf: output metric filename. Stdout if `outf = None`. <br>
<br>

### unpaired_latent_metrics
```python
## for all unpaired integration metrics
unpaired_latent_metrics(adata,
                            method,
                            cluster = 'louvain',
                            batch = 'batch',
                            label = 'cell_type',
                            mods = ['RNA', 'ATAC'],
                            outf = None,
                            embed_acc = True
                           ):
```
- adata: `anndata` object. the `obs` is celll metadata. Two latent embeddings are stored in obsm['RNA'] and obsm['ATAC'] or obsm['ADT']<br>
- method: name of the input algorithm. <br>
- cluster: clustering methods, deprecated. <br>
- batch: column name of batch used. Do not consider batch effect removal if `batch=None`.<br>
- label: column name of cell type used. Do not consider biological conservation if `label=None`<br>
- mods: the name of obsm data used for evaluation. <br>
- outf: output metric filename. Stdout if `outf = None`. <br>
<br>
<br>

### mosaic_latent_metrics
```python
## for mosaic integration metrics in accuracy and downsample robustness tasks.
mosaic_latent_metrics(latents,metadatas,paired="s1d1", unpaired="s3d10", 
                            mod2="atac", batch="batch",label="cell_type",
                            latent_path="", method='sciPENN', writef=True):
```
- latents: python list of np.ndarry or pandas Dataframe object for paired, RNA and the other modality. Detail see tutorials. <br>
- metadatas: python list of np.ndarry or pandas Dataframe object for cell metadata of paired, RNA and the other modality. <br>
- paired: batch name of paired cells in metadata batch column. <br>
- unpaired: batch name of unpaired cells in metadata batch column. <br>
- mod2: the second mosaic modality. "atac" or "adt" <br>
- batch: column name of batch used. Do not consider batch effect removal if `batch=None`.<br>
- label: column name of cell type used. Do not consider biological conservation if `label=None`<br>
- latent_path: name of input latent, must contain "latent" in filename. And two metrics files will generated in the samme path of different suffix. <br>
- method: name of the input algorithm. <br>
- writef: write two metrics file in `latent_path` if True. stdout if FALSE. <br>

<br>

### mosaic_cnk_latent_metrics
```python
# for mosaic integration metrics in paired sizes robustness tasks. only focus on unpair size
mosaic_cnk_latent_metrics(latents,metadatas,paired="s1d1", unpaired="s3d10", 
                            mod2="atac", batch="batch",label="cell_type",latent_path="",
                             method='sciPENN', writef=True):
```
- latents: python list of np.ndarry or pandas Dataframe object for paired, RNA and the other modality. <br>
- metadatas: python list of np.ndarry or pandas Dataframe object for cell metadata of paired, RNA and the other modality. <br>
- paired: batch name of paired cells in metadata batch column. <br>
- unpaired: batch name of unpaired cells in metadata batch column. <br>
- mod2: the second mosaic modality. "atac" or "adt" <br>
- batch: column name of batch used. Do not consider batch effect removal if `batch=None`.<br>
- label: column name of cell type used. Do not consider biological conservation if `label=None`<br>
- latent_path: name of input latent, must contain "latent" in filename. And two metrics files will generated in the samme path of different suffix. <br>
- method: name of the input algorithm. <br>
- writef: write two metrics file in `latent_path` if True. stdout if FALSE. <br>
<br>

### imputation_pair_rna_atac
```python
# For paired scRNA and scATAC imputation methods
imputation_pair_rna_atac(metadata_path,rna_imp_path, rna_path,atac_imp_path, atac_path, method, outf=None ):

```
- metadata_path: path of cell metadata.
- rna_imp_path: path of rna imputation csv file.
- rna_path: path of rna raw h5ad file. Used as rna gold standard.  
- atac_imp_path: path of atac imputation csv file.
- atac_path: path of atac raw h5ad file. Used as atac gold standard.  
- method: name of the input algorithm. <br>
- outf: output to given file name or pd.Dataframe() (None). <br>
<br>

### imputation_mosaic_rna_atac
```python
# For unpaired scRNA and scATAC mosaic imputation methods
imputation_mosaic_rna_atac(metadata_path,rna_imp_path, rna_path,atac_imp_path, atac_path,method,paired="s1d1",unpaired="s3d10",batch="batch", outf=None ):
```
- metadata_path: path of cell metadata.
- rna_imp_path: path of rna imputation csv file.
- rna_path: path of rna raw h5ad file. Used as rna gold standard.
- atac_imp_path: path of atac imputation csv file.
- atac_path: path of atac raw h5ad file. Used as atac gold standard.
- method: name of the input algorithm. <br>
- paired: batch name of paired cells in metadata batch column. <br>
- unpaired: batch name of unpaired cells in metadata batch column. <br>
- outf: output to given file name or pd.Dataframe() (None). <br>
<br>

### imputation_mosaic_rna_adt
```python
# For unpaired scRNA and ADT mosaic imputation methods
imputation_mosaic_rna_adt(metadata_path,rna_imp_path=None, rna_path=None,adt_imp_path=None, adt_path=None, method="sciPENN",paired="s3d6", unpaired="s2d1", batch="batch", outf=None ):
```
- metadata_path: path of cell metadata.
- rna_imp_path: path of rna imputation csv file.
- rna_path: path of rna raw h5ad file. Used as rna gold standard.
- adt_imp_path: path of atac imputation csv file.
- adt_path: path of atac raw h5ad file. Used as adt gold standard.
- method: name of the input algorithm. <br>
- paired: batch name of paired cells in metadata batch column. <br>
- unpaired: batch name of unpaired cells in metadata batch column. <br>
- outf: output to given file name or pd.Dataframe() (None). <br>
<br>

### imputation_rna
```python
imputation_rna(metadata_path,rna_imp_path, rna_path, method, outf=None ):
#  For scMVAE only. Only evaluate the scRNA modality.
```
- metadata_path: path of cell metadata.
- rna_imp_path: path of rna imputation csv file.
- rna_path: path of rna raw h5ad file. Used as rna gold standard.
- outf: output to given file name or pd.Dataframe() (None). <br>
<br>

### imputation_stabmap
```python
# For stabmap mosaic scRNA and scATAC only, which omit scATAC knn smoothing for too few scATAC peaks.

imputation_stabmap(metadata_path,rna_imp_path, rna_path,atac_imp_path, atac_path,method,paired="s1d1",unpaired="s3d10",batch="batch", outf=None ):

```
- metadata_path: path of cell metadata.
- rna_imp_path: path of rna imputation csv file.
- rna_path: path of rna raw h5ad file.
- atac_imp_path: path of atac imputation csv file.
- atac_path: path of atac raw h5ad file.
- method: name of the input algorithm. <br>
- paired: batch name of paired cells in metadata batch column. <br>
- unpaired: batch name of unpaired cells in metadata batch column. <br>
- outf: output to given file name or pd.Dataframe() (None). <br>
<br>

### mouse_brain_divide
```python
mouse_brain_divide(func, adatas,
                        method,
                        cluster = 'louvain',
                        batch = 'batch',
                        label = 'cell_type',
                        outf = None
                        ):
## A wrap function for mouse brain dataset and all functions above.

#    10X mouse brain datasets calculated batch removal metrics for 
#    WT(wild type) and AD brain separately, then took the average
#     of metrics from two disease status group.

# Example
mouse_brain_divide(paired_latent_metrics,[adata1,adata2], method, 'louvain', batch,label, outfile)
```

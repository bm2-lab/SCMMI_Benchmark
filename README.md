# SCMMI_Benchmark

---------------------------

### Description

SCMMI_benchmark is a benchmark pipeline for evaluating the usability, accuracy, robustness and scalability of single-cell multimodal integration algorithms, including 40 single-cell multi-modal integration methods involving modalities of DNA, RNA, protein and spatial multi-omics for paired integration, unpaired diagonal integration, and unpaired mosaic integration.

### Preprocessing

The scripts within the preprocess folder are designed to perform various manipulations on datasets, including single-batch sampling, downsampling and upsampling of cell numbers for scalability tasks, as well as downsampling of sequencing depth for RNA and ATAC for robustness tasks. Concurrently, the dataset is converted into several formats, comprising mtx, rds, h5ad, and h5Seurat, thereby accommodating the input requirements of various methodes. <br>
- data preparation and subsetting [code](scripts/preprocess/subset_data.R)
- input data format conversion [code](scripts/preprocess/format_convert.py)

### Assessment

The "assess" folder contains scripts to calculate single-embedding metrics and cross-embedding metrics for accuracy evaluation of multimodal integration.<br>

- single-embedding metrics [code](scripts/assess/single_embedding_metrics_benchmark.py)
- cross-embedding metrics [code](scripts/assess/cross_embedding_metrics_benchmark.py)

### Datasets
All datasets currently used in this benchmark study are listed below.
| dataset name       | multi-omics type | Batch                      | species | cell number | Tested tasks                                               |   |   |   |   |   |
|--------------------|------------------|----------------------------|---------|-------------|------------------------------------------------------------|---|---|---|---|---|
| BMMC Multiome      | RNA + ATAC       | 12 donors from 4 sites     | Human   | 69,249      | All accuracy metrics and robustness metrics                |   |   |   |   |   |
| BMMC CITE-seq      | RNA + ADT        | 12 donors from 4  sites    | Human   | 90,261      | All accuracy metrics and robustness metrics                |   |   |   |   |   |
| HSPC Multiome      | RNA + ATAC       | 4 donors  of 5 time points | Human   | 105,942     | Batch removal, cell alignment and imputation               |   |   |   |   |   |
| HSPC CITE-seq      | RNA + ADT        | 4 donors of 5 time points  | Human   | 70,988      | Batch removal, cell alignment and imputation               |   |   |   |   |   |
| SHARE-seq skin     | RNA + ATAC       | -                          | Mouse   | 34,774      | Biological conservation, cell alignment and imputation     |   |   |   |   |   |
| COVID19 CITE-seq   | RNA+ADT          | 143 donors                 | Human   | 781,123     | batch removal and paired and unpaired RNA+ADT Scalability, |   |   |   |   |   |
| Lymph node spatial | spatial+RNA+ADT  | 2 samples                  | Human   | 6843        | Biological conservation and batch removal                  |   |   |   |   |   |
| Thymus spatial     | spatial+RNA+ADT  | 4 samples                  | Mouse   | 17,824      | Batch removal                                              |   |   |   |   |   |
| Spleen SPOTS       | spatial+RNA+ADT  | 2samples                   | Mouse   | 5,336       | Batch removal                                              |   |   |   |   |   |
<br>


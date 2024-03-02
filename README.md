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
| dataset name       | multi-omics type | Batch                      | species | cell number | tissue                                  |
|--------------------|------------------|----------------------------|---------|-------------|-----------------------------------------|
| BMMC Multiome      | RNA + ATAC       | 12 donors from 4 sites     | Human   | 69,249      | bone marrow mononuclear cells           |
| BMMC CITE-seq      | RNA + ADT        | 12 donors from 4  sites    | Human   | 90,261      | bone marrow mononuclear cells           |
| HSPC Multiome      | RNA + ATAC       | 4 donors  of 5 time points | Human   | 105,942     | hematopoietic stem and progenitor cells |
| HSPC CITE-seq      | RNA + ADT        | 4 donors of 5 time points  | Human   | 70,988      | hematopoietic stem and progenitor cells |
| SHARE-seq skin     | RNA + ATAC       | -                          | Mouse   | 34,774      | skin                                    |
| COVID19 CITE-seq   | RNA+ADT          | 143 donors                 | Human   | 781,123     | peripheral blood immune cells           |
| Lymph node spatial | spatial+RNA+ADT  | 2 samples                  | Human   | 6843        | lymph node                              |
| Thymus spatial     | spatial+RNA+ADT  | 4 samples                  | Mouse   | 17,824      | thymus                                  |
| Spleen SPOTS       | spatial+RNA+ADT  | 2 samples                  | Mouse   | 5,336       | spleen                                  |
<br>


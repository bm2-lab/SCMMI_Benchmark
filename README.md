# SCMMI_Benchmark

---------------------------

### Description

SCMMI_benchmark is a benchmark pipeline for evaluating the usability, accuracy and robustness of single-cell multimodal integration algorithms, including 16 paired scRNA+scATAC integration methods, 13 unpaired scRNA+scATAC integration methods and 12 scRNA+ADT integration methods. This folder now contains accuracy assessment scripts for benchmark analysis.

### Preprocessing scripts for SCMMI_Benchmark manuscript
The "preprocess" folder includes scripts to prepare benchmark datasets for different algorithms from baseline datasets. <br>
- data preparation and subsetting [code](scripts/preprocess/subset_data.R)
- input data format conversion [code](scripts/preprocess/format_convert.py)


### Assessment scripts for SCMMI_Benchmark manuscript

The "assess" folder includes scripts to calculate single-embedding metrics and cross-embedding metrics for paired and unpaired multimodal integration.  <br>

- single-embedding metrics [code](scripts/assess/single_embedding_metrics_benchmark.py)
- cross-embedding metrics [code](scripts/assess/cross_embedding_metrics_benchmark.py)

### Datasets
The benchmark baseline datasets can be downloaded from [single-cell multimodal integration challenge competition](https://openproblems.bio/events/2021-09_neurips/).  <br>





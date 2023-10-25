# SCMMI_Benchmark

---------------------------

### Description

SCMMI_benchmark is a benchmark pipeline for evaluating the usability, accuracy and robustness of single-cell multimodal integration algorithms, including 16 paired scRNA+scATAC integration methods, 13 unpaired scRNA+scATAC integration methods and 12 scRNA+ADT integration. This folder now contains the accuracy assessment modules/scripts for benchmark analysis.


### Assessment scripts for SCMMI_Benchmark manuscript

The folder includes the scripts to calculate the single-embedding metrics and cross-embedding metrics for paired and unpaired multi-modal integration.  <br>

- [paired scRNA+scATAC metrics](scripts/assess/pair_metrics_benchmark.py)
- [paired scRNA+ADT metrics](scripts/assess/CITE_seq_metrics_benchmark.py)
- [unpaired scRNA+scATAC metrics](scripts/assess/unpair_metrics_benchmark.py)

### Datasets
The benchmark baseline datasets can be downloaded from [single-cell multimodal integration challenge competition](https://openproblems.bio/events/2021-09_neurips/).  <br>





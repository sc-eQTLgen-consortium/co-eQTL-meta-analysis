# Sc-gene networks in the brain

This repository contains the scripts needed for running the Snakemake (version 6.4.1) pipeline developed to create single cell gene networks. The pipeline takes a Seurat object containing pre-processed expression data and GWAS summary statistics and outputs an excel file containing a list of genes and their corresponding gene prioritization statistics. The produced co-expression matrix and prioritized genes can then be used to create gene networks.  

The scripts used in this pipeline can be found under **/scripts** and the pipeline can be run by using the `snakemake` command.

The **gene_pair_corrs.yaml** configuration file is used by the pipeline to specify file paths, parameters and settings used. 

### Dependencies
| Packages      | Version   |
|---------------|-----------|
| Downstreamer  | 2.2       |
| Matplotlib    | 3.3.2     |
| networkx      | 3.3       |
| NumPy         | 1.23.5    |
| optparse      | 1.7.3     |
| Pandas        | 2.0.3     |
| PascalX       | 0.0.3     |
| SciPy         | 1.10.1    |
| Seaborn       | 0.13.2    |
| Seurat        | 4.3.0.1   |
| statsmodels   | 0.14.2    |

*Author:* Anoek Kooijmans

*Date: 13-06-2024*

**License: GNU GPL v3**
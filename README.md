# Sc-gene networks in the brain

This documentation outlines the steps and components of the Snakemake pipeline developed to create single cell gene networks. The pipeline takes a Seurat object containing pre-processed expression data and GWAS summary statistics and outputs an excel file containing a list of genes and their corresponding gene prioritization statistics. The produced co-expression matrix and prioritized genes can then be used to create gene networks.  

The following text outlines the functions of the scripts implemented in this pipeline. The scripts can be found under `/scripts` and the pipeline can be run by 

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


## Step 1 - Calculate correlation

The first step of this analysis is to identify gene co-regulation patterns by calculating the correlation between each gene pair. Correlation analysis can be run parallel over different cell-types, cohorts, number of genes and correlation methods. 

| Argument    	| Description                                                             	|
|-------------	|-------------------------------------------------------------------------	|
| --celltype  	| Cell-type                                                               	|
| --cohort    	| Cohort-id                                                               	|
| --n         	| Number of genes to include. The top n expressed genes will be selected  	|
| --method    	| Correlation method. Either 'pearson' or 'spearman'                      	|
| --weight    	| Whether to perform weighted or unweighted correlation (True/False)      	|

### processRds.R

First, the Seurat file from the specified cell type is processed. The PFlog1pPF normalized expression counts are filtered on genes, using either a selected gene set, or the top n expressed genes. The counts are then filtered on a list of selected donors, saving a list of cell-barcode identifiers for each donor. Gene filtering is then performed over each donor, filtering out genes that have fewer than 10 cells showing detectable expression.  

| Argument     	| Description                                                           	|
|--------------	|-----------------------------------------------------------------------	|
| --genepath   	| File path to list of genes of interest. Can be set to 'nan            	|
| --donorpath  	| List of donors to include                                             	|
| --input      	| Seurat data for a specific cell type: '{ct}.Qced.Normalized.SCs.Rds'  	|
| --output     	| normalized-counts-{ct}-top-{n}.tsv.gz"                                	|
|              	| correlation-weight-{ct}-top-{n}.tsv.gz                                	|
|              	| donor_barcodes/barcodes-{donor}-{ct}-top-{n}.tsv.gz                   	|
|              	| gene_list/{ct}-top-{n}-gene-list.tsv.gz                               	|
|              	| donor_list/{ct}-top-{n}-donor-list.tsv.gz                             	|


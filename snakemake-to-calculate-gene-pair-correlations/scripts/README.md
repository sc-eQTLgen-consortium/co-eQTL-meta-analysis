## Step 1 - Calculate correlation

The first step of this analysis is to identify gene co-regulation patterns by calculating the correlation between each gene pair. Correlation analysis can be run parallel over different cell-types, cohorts, number of genes and correlation methods. 

| Argument    	| Description                                                             	|
|-------------	|-------------------------------------------------------------------------	|
| --celltype  	| Cell-type                                                               	|
| --cohort    	| Cohort-id                                                               	|
| --n         	| Number of genes to include. The top n expressed genes will be selected  	|
| --method    	| Correlation method. Either 'pearson' or 'spearman'                      	|
| --weight    	| Whether to perform weighted or unweighted correlation (True/False)      	|

#### processRds.R

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

#### calculate_correlations.py

Next, correlation is calculated on the normalized expression. The expression data is then loaded line by line to optimize computational time. This step is done in parallel over the donors, selecting the donor-specific correlation weights and counts from the expression matrix using the cell-barcodes saved per donor in the previous step. A correlation is calculated per gene, excluding genes with less than 10 non-zero counts. These are filled in as NaN. The correlation coefficient, and its corresponding correlation metrics (p-value and z-score) are saved.  

| Argument         	| Description                                               	|
|------------------	|-----------------------------------------------------------	|
| --count          	| Normalized counts file path                               	|
| --method         	| Correlation method                                        	|
| --weight         	| Whether to calculate weighted or unweighted correlation   	|
| --cell_barcodes  	| File path to list of barcodes per donor                   	|
| --output         	| {donor}-{ct}-top-{n}-{method}-{weight}.tsv.gz             	|

#### aggregate_metrics.py

The correlation metrics are aggregated over the donors resulting in a matrix that contains the correlation metric for each sample per gene-pair.  

| Argument  	| Description                                            	|
|-----------	|--------------------------------------------------------	|
| --input   	| {donor}-{ct}-top-{n}-{method}-{weight}.tsv.gz          	|
| --output  	| combined-corr-{ct}-top-{n}-{method}-{weight}.tsv.gz    	|
|           	| combined-std-{ct}-top-{n}-{method}-{weight}.tsv.gz     	|
|           	| combined-pval-{ct}-top-{n}-{method}-{weight}.tsv.gz    	|
|           	| combined-zscore-{ct}-top-{n}-{method}-{weight}.tsv.gz  	|

#### sample_quality_control.py

In this step donors with a cell count of less than 20 are removed from the gene-pair-by-donor correlation matrix created in the previous step. Next PCA is performed over the matrix and donors with a an absolute z-score higher than 3 on PC1 and PC2 were removed.  

| Argument              	| Description                                                             	|
|-----------------------	|-------------------------------------------------------------------------	|
| --input               	| combined-corr-{ct}-top-{n}-{method}-{weight}.tsv.gz                     	|
| --output              	| outliers-removed-combined-corr-{ct}-top-{n}-{method}-{weight}.tsv.gz    	|
| --donor_list          	| donor_list/{ct}-top-{n}-donor-list.tsv.gz                               	|
| --donor_list_updated  	| donor_list/{ct}-top-{n}-{method}-{weight}-donor-list-updated.tsv.gz     	|

#### create_matrix.py

The correlation coefficients are then meta-analyzed over donors into a single correlation coefficient per gene-pair. This is done by either taking the mean over the donors or calculating the Fisher's Z summary coefficient weighted by cell count per donor. Here, genes are only included if they were expressed in at least six donors.  This step outputs a gene-by-gene matrix of summarized correlation coefficients and p-values.  

| Arguement            	| Description                                                                             	|
|----------------------	|-----------------------------------------------------------------------------------------	|
| --meta_analysis      	| Meta-analysis method used to combine correlation over donors                            	|
| --donor_list         	| donor_list/{ct}-top-{n}-{method}-{weight}-donor-list-updated.tsv.gz                     	|
| --gene_list          	| gene_list/{ct}-top-{n}-gene-list.tsv.gz                                                 	|
| --gene_list_updated  	| gene_list/{meta_analysis}-corr-{ct}-top-{n}-{method}-{weight}-gene-list-updated.tsv.gz  	|
| --input              	| combined-corr-{ct}-top-{n}-{method}-{weight}.tsv.gz                                     	|
| --output             	| matrix-{meta-analysis}-corr-{ct}-top-{n}-{method}-{weight}.tsv.gz                       	|
|                      	| matrix-{meta-analysis}-pval-{ct}-top-{n}-{method}-{weight}.tsv.gz                       	|


## Step 2 - Gene prioritization 

Next, key genes in the gene-gene correlation matrix are prioritized using GWAS summary statistics with Downstreamer. Downstreamer identifies a set of key genes in the correlation matrix that are most associated with the provided GWAS genes. 

#### runRealGWAS.py

GWAS summary statistics are converted from p-values per variant to an aggregate p-value per gene, accounting for LD structures (1000 Genomes phase 3). Variant p-values are aggregated within a 25 kb window around the start and end of each gene.  

| Argument      	| Description                    	| Settings used                	|
|---------------	|--------------------------------	|------------------------------	|
| --refpanel    	| Reference panel path           	| 1kGP_high_coverage_Illumina  	|
| --gwas        	| GWAS summary statistics path   	| {sum_stats}.txt.gz           	|
| --annotation  	| Genome annotation file         	| genes_Ensembl94.txt          	|
| --threads     	| Number of threads for scoring  	| 20                           	|
| --outfile     	| Output file                    	| pascalX/{gwas}.txt           	|
| --rscol       	| Rs ids column index            	| 0                            	|
| --pscol       	| P-value column index           	| 1                            	|
| --window      	| Gene window                    	| 25000                        	|

#### eigendecomposeCorr.py

Eigen decomposition is performed on the meta-analyzed correlation matrix. Eigen decomposition enables the representation of the correlation matrix as a product of eigenvectors and eigenvalues. First, the gene-ids in the gene-gene correlation matrix are converted to Ensembl ids. Next eigen decomposition is performed and eigenvectors are selected that explain 80% of the cumulative explained variance.  

| Argument  	| Description                   	|                                                                    	|
|-----------	|-------------------------------	|--------------------------------------------------------------------	|
| --input   	| Gene-gene correlation matrix  	| matrix-{meta-analysis}-corr-{ct}-top-{n}-{method}-{weight}.tsv.gz  	|
|           	| Gene-annotation file          	| LimixAnnotationFile.txt                                            	|
| --output  	| Selected eigenvectors         	| EigVec-{meta_analysis}-{ct}-top-{n}-{method}-{weight}.txt          	|
|           	| Selected eigenvalues          	| EigVal-{meta_analysis}-{ct}-top-{n}-{method}-{weight}.txt          	|

#### reformat_eigenVec

The eigenvectors are then reformatted from `txt.gz` to `.cols.txt.gz`, `.rows.txt.gz` and `.datg`

| Argument  	| Description                                                       	|
|-----------	|-------------------------------------------------------------------	|
| --input   	| EigVec-{meta_analysis}-{ct}-top-{n}-{method}-{weight}.txt         	|
| --output  	| EigVec-{meta_analysis}-{ct}-top-{n}-{method}-{weight}.datg        	|
|           	| EigVec-{meta_analysis}-{ct}-top-{n}-{method}-{weight}.cols.txt.gz 	|
|           	| EigVec-{meta_analysis}-{ct}-top-{n}-{method}-{weight}.rows.txt.gz 	|

#### gene_prioritisation

Downstreamer prioritizes key genes using GWAS summary statistics. This is done by association analysis of the selected correlation eigenvectors with the GWAS gene p-values.  

| Argument                  	| Description                                                                                                           	| Settings used                                                	|
|---------------------------	|-----------------------------------------------------------------------------------------------------------------------	|--------------------------------------------------------------	|
| --mode                    	| Downstreamer mode                                                                                                     	| ENRICH                                                       	|
| --gwas                    	| GWAS gene p-values                                                                                                    	| pascalX/{gwas}.txt                                           	|
| --geneCorrelations        	| Gene expression matrix force normalized and split per chromosome                                                      	| permutationGeneCor/geneCorForceNormalchr_                    	|
| --output                  	| Output file path                                                                                                      	| {gwas}_keygenes_covCor                                       	|
| --genes                   	| File with gene information  col1: geneName (ensg) col2: chr col3: startPos col4: stopPos col5: geneType col6: chrArm  	| genes_Ensembl94_protein_coding.txt                           	|
| --expressionEigenVectors  	| Selected eigenvectors of gene-genecorrelation matrix                                                                  	| {ct}= EigVec-{meta_analysis}-{ct}-top-{n}-{method}-{weight}  	|
| --covariates              	| Corrects for potential inflation in correlation signal using median signal                                            	| medianGwasSignal.txt                                         	|
| --eh                      	| Exclude HLA locus during pathway enrichment (chr6 20mb - 40mb)                                                        	|                                                              	|
| -t                        	| Maximum number of calculation threads                                                                                 	| 8                                                            	|
| --forceNormalGenePvalues  	| Force normal gene p-values before pathway enrichment                                                                  	|                                                              	|
| --jblas                   	|                                                                                                                       	|                                                              	|
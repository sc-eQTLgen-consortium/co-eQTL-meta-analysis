# co-eQTL-meta-analysis

Caution: Pipeline is not yet ready for public use

## prerequisites

This pipeline makes use of the Snakemake. The pipeline was developed for Snakemake version 7.32.3. Python 3.12 is known to have issues with this Snakemake version, so python 3.11 is recommended.
The easiest way to get you up and running is to use conda (https://docs.anaconda.com/miniconda/):

```sh
conda create -n snakemake_env python==3.11
conda activate snakemake_env
conda install -n snakemake_env pip
pip install snakemake==7.32.3
```

## instructions

This pipeline consists of a number of major steps, listed below

### calculating gene-gene correlations per donor

we will need two major directories for files, and three supporting files:
- a place to store the software used in this pipeline
- a place to store intermediate and final results
- the singularity image with all the software
- the binary matrix jar file
- the limix annotation file that has gene info

#### setup the correlation pipeline steps

In the template yaml file, (https://github.com/sc-eQTLgen-consortium/co-eQTL-meta-analysis/blob/main/snakemake-to-calculate-gene-pair-correlations/gene_pair_corrs.yaml), I chose this directory: /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/coeqtl_redo_test/software

I next cd-ed into this directory and did a git clone of this repository:
```sh
mkdir -p /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/coeqtl_redo_test/software
cd /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/coeqtl_redo_test/software
git clone https://github.com/sc-eQTLgen-consortium/co-eQTL-meta-analysis.git
```

now let's go into the correlation calculation folder and get the full path
```sh
cd co-eQTL-meta-analysis/snakemake-to-calculate-gene-pair-correlations/
pwd
```

take note of this full path, as we will need it when setting up the yaml (the configuration file)

as mentioned previously, we will need a place to store intermediate and final files. In the example I chose this directory: /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/coeqtl_redo_test/output/

to make sure it exists, I did a mkdir:
```sh
mkdir -p /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/coeqtl_redo_test/output/
```

again, make note of this path, as we will need it for the yaml

we will also need to have the singularity image, which can be downloaded here !TODO. Make note of where this file is located, as we will need to supply its location in the yaml file.

next we need the binary matrix jar file, a java application for matrix operations. It can be downloaded here !TODO, once more make note of the location.

Finally we will need gene annotations to allows to do operations per chromosome. For this we will use the limix annotation file, available here !TODO. Store it somewhere and not the location.

#### correlation pipeline steps

In short, this pipeline does the following steps:

Step one is to generate files to allow snakemake to calculate all of the parallel jobs it will have to run.
This step uses the script: donor_ids.R

Step two is to separate the main seurat file into one file per individual (to allow for more parallelisation).
This step uses the script: create_donor_rds.R

Step three calculates the correlations for every gene pair within an individual.
This step uses the script: process_donor_rds.R

Step four is to aggregate all of these gene pair correlations files from one per individual into one file with all of the correlations.
This step uses the script: aggregate_metrics.py

Now we can start modifying the yaml. You can copy the contents of the yaml into a text editor, clone the github repo and make your own branch, or just edit in it nano. 

#### yaml parameters for the correlation pipeline

##### top_directory
directory that should be bound in singularity, all files should be a subdirectory of this one. If you need multiple binds, separate them with a comma.

##### output_directory
output directory of all the files, this is where the final output and the intermediate files will be saved. This is the output directory alluded to before

##### singularity_path
the location of the singularity image with all installed software, the one mentioned before

##### snakefile_path
the working directory for snakemake, location of this file and the snakefile. In the example this was /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/coeqtl_redo_test/software/co-eQTL-meta-analysis/snakemake-to-calculate-gene-pair-correlations/, and if you followed the instructions, this should be the output of the pwd command we did earlier

##### scripts_directory
the location of the scripts that are used, would usually be the scripts directory of the snakefile_path mentioned before

##### input_directory:
the directory that contains the RDS files of Seurat objects. If you followed the WG3 pipeline, this would be one of the output directories, and should house files named like NK.Qced.Normalized.SCs.Rds

##### qtl_celltype_directory
the directory containing the WG3-style eQTL input such as the qtlInput and PCs files, will be the super directory of the input_directory if you followed the WG3 pipeline to the letter, and should contain the L1 directory and smf file in that case.

##### qtl_celltype_subdirectory
the subdirectory, if applicable for the WG3-style eQTL input

##### limix_annotation_loc
the location of the limixannotation file

##### binary_matrix_loc
location of the binary matrix jar

##### relative_smf_path
the relative location of the .smf file from the qtl celltype and qtl celltype sub directories

##### seurat_object_prepend
prepend of name of RDS files. The part of the Seurat object name that goes before the cell type. In the case of following the wg3 pipeline, this is not present and can be kept empty

##### seurat_object_append
append of the RDS files. The part of the Seurat object name that goes after the cell type. In the case of following the wg3 pipeline, this would be .Qced.Normalized.SCs.Rds

##### seurat_assignment_column
name of the column in the RDS Seurat object, that denotes the sample assignment. In the case of the WG3 pipeline, this would be Assignment

##### weighting_method 
weighting method for the correlations, can be expression, zeroes or none

##### cohort
ID of cohort, this is up to you

##### cell_type
cell type(s) to calculate correlations for. The names need to match the file names you might have for the QTLinput and Seurat objects. If you are going to calculate specific gene-gene pairs, note the order of the cell types you supply

##### gene_list_path: 
gene lists to only calculate correlations for (instead of doing all genes), if not using a specific set, should be nan. If you want to use the same gene set for each celltype, supply one file. If you want a different geneset for each cell type, supply files in the same order as you supplied the respective celltypes the gene sets are for

#### running the correlation pipeline
In its current state, this pipeline runs the correlation calculation on the machine where the pipeline is started. As such, it is recommended to do this on a machine where you have compute available. The pipeline will also take a while, depending on the size of the dataset. As such it is also recommended to do this in a way where you can disconnect and reconnect without the pipeline stopping, such as sattach, screen on tmux. Here is an example using screen and SLURM:

create a screen session:

```sh
screen -S correlation_calculation
```

ask for resources. The more CPUs you ask, the faster things will go, but also the more memory you require. Larger datasets will also require more memory.

```sh
srun --cpus-per-task=4 --mem=64gb --nodes=1 --qos=priority --job-name=correlation_calculation --time=71:59:59 --tmp=1000gb --pty bash -i # this requests an interactive session
```

activate the conda environment with snakemake

```sh
conda activate snakemake_env
```

cd to our snakemake directory that we got with the pwd before, and start the pipeline (remember to set a correct number of cores)

```sh
cd /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/coeqtl_redo_test/software/co-eQTL-meta-analysis/snakemake-to-calculate-gene-pair-correlations/
snakemake -s gene_pair_corrs_Snakefile --cores 4 --rerun-incomplete
```

The first step should be very quick, as it will only create donor lists. Run the same command again to run the rest of the pipeline

```sh
snakemake -s gene_pair_corrs_Snakefile --cores 4 --rerun-incomplete
```

After this finishes, you can continue to the next step, where we do the co-eQTL mapping.

### performing co-eQTL mapping

Before we can start the mapping we need to generate some files. The underlying LIMIX pipeline works best when using bgen files. To reduce memory usage, it is also best to separate the genotype data per chromosome. To convert and split that data, it is recommended to use plink2, which is fast and efficient. To for example split a VCF of genotype data into per-chromosome bgen files, one can use a command like this:

```sh
CHROMS=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22')
for chrom in ${CHROMS[*]}
  do
    /groups/umcg-franke-scrna/tmp04/software/plink2_amd/plink2 \
    --vcf /groups/umcg-franke-scrna/tmp04/projects/venema-2022/ongoing/genotype/combined/imputed/all/sc-eqtlgen-pipeline-wg1/lpmcv2_imputed_hg38_info_filled_rsid.vcf.gz \
    --export bgen-1.2 \
    --out /groups/umcg-franke-scrna/tmp04/projects/venema-2022/ongoing/genotype/combined/imputed/all/sc-eqtlgen-pipeline-wg1/lpmcv2_imputed_hg38_info_filled_rsid_chr${chrom} \
    --chr ${chrom}
    done
```

An index needs to be generated for each of the bgen files. This index is created when the bgen file is opened for the first time. As the snakemake pipeline submits multiple parallel jobs, it might be that two processes try to open a file at the same time, while one of them might still be creating the index. As such, it is safest to just open each bgen file once, before running the pipeline, to make sure the index files are already fully generated. A small script to this, is included in the repository. If you have everything already set up from the previous step, you can go straight to running the script. If not, I will repeat the steps required here:

```sh
# make working directory
mkdir -p /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/coeqtl_redo_test/software
# go to working directory
cd /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/coeqtl_redo_test/software
# clone repo
git clone https://github.com/sc-eQTLgen-consortium/co-eQTL-meta-analysis.git
# activate snakemake environment]
conda activate snakemake_env
```

now we will go into the directory that has the mapping pipeline and get the path to that

```sh
cd snakemake-to-map-coeqtls/
pwd
```

make note of this directory. We will now install bgen reader into our snakemake environment to open the bgen files for indexing
```sh
pip install bgen_reader
```

with that sorted, we can now run the script for indexing the bgen files. The script is a subdirectory of the directory I told you to note earlier, and the script name is 'coeqtl_open_bgens.py'. It can use two parameters, the first is the directory of the bgen files, and the second is a regular expression to recognize the bgen files (*.bgen to just grab each bgen file for example). In my case the command looks like this:

```sh
python /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/coeqtl_redo_test/software/snakemake-to-map-coeqtls/scripts/coeqtl_open_bgens.py \
    --path /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/genotype_input/ \
    --regex EUR_imputed_hg38_varFiltered_chr*.bgen
```

next we are going to have to start creating some annotation files for the gene-pairs, and the variants that we want to test. Because the multiple testing burden is quite high when you consider gene pairs, we are only testing co-eQTLs for which one of the genes already shows an eQTL effect. Additionally, there is the option to include GWAS variants for those genes, if you have this type of information available. Finally we can only test genes that we included when we calculated the per-sample gene-gene correlations, and will only create annotations for those gene pairs, because annotation files for each pair of expressed genes would quickly balloon in size. This script is once again in the scripts folder, and is named 'coeqtl_make_limix_annotation_files_new.py'. The options and their names are listed below:

```sh
python /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/coeqtl_redo_test/software/snakemake-to-map-coeqtls/scripts/coeqtl_make_limix_annotation_files_new.py -h
# usage: coeqtl_make_limix_annotation_files_new.py [-h] [-g GWAS_LOC] [-l GENE_LOC] [-v VARIANT_LOC] [-q QTL_LOC] [-f FEATURES_OUT_LOC] [-o CO_LIMIX_ANNOTATION_PREPEND]
#
# options:
#   -h, --help            show this help message and exit
#   -g GWAS_LOC, --gwas_loc GWAS_LOC
#                         location of gwas file (string)
#   -l GENE_LOC, --gene_loc GENE_LOC
#                         location of gene list file (string)
#   -v VARIANT_LOC, --variant_loc VARIANT_LOC
#                         location of variant list file (string)
#   -q QTL_LOC, --qtl_loc QTL_LOC
#                         location of previous eQTL summary stats
#   -f FEATURES_OUT_LOC, --features_out_loc FEATURES_OUT_LOC
#                         location to write the to-test features
#   -o CO_LIMIX_ANNOTATION_PREPEND, --co_limix_annotation_prepend CO_LIMIX_ANNOTATION_PREPEND
#                         location of original feature annotations
```

We will unfortunately have to do this for each celltype separately for now, though you can loop the cell types in bash if you wanted to. In my case the command I used for my CD4+ T cells was the following:

```sh
/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/coeqtl_redo_test/software/snakemake-to-map-coeqtls/coeqtl_make_limix_annotation_files_new.py \
    --gwas_loc /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/sceqtlgen_coeqtl_map/GWAS_snp_gene_pairs_immune_related_disease.txt.gz \
    --gene_loc /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/calculate_correlation_metrics/gene_lists/F11_Decision_Tree_Geneswg3_wijst2018_Mono.Qced.Normalized.SCs.Rds.tsv.gz \
    --qtl_loc /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/sceqtlgen_coeqtl_map/eQTLs_finemapped_20240626/Mono.DS.wg3_Ye_wg3_wijst2018_wg3_sawcer_wg3_oneK1K_wg3_okada_wg3_Li_wg3_Franke_split_v3_wg3_Franke_split_v2_wg3_multiome_UT_wg3_idaghdour.csTop_qtl_results.txt \
    --features_out_loc /groups/umcg-franke-scrna/tmp04/projects/venema-2022/ongoing/qtl/coeqtl/input/replication_features_CD4_T_cells.tsv.gz \
    --co_limix_annotation_prepend /groups/umcg-franke-scrna/tmp04/projects/venema-2022/ongoing/qtl/coeqtl/input/replication_CD4_T_cells_co

```

Next the annotation file needs to split by chromosome, which is done by the script 'coeqtl_limix_anno_chr_files_new.R' This script expects the prepend of the --features_out_loc parameter of the previous script (without the cell type), the prepend of the --co_limix_annotation_prepend (without the cell type), and the cell type as used in the name of the files. The options and their names are listed below:

```sh
~/start_Rscript.sh /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/coeqtl_redo_test/software/snakemake-to-map-coeqtls/coeqtl_limix_anno_chr_files_new.R -h
# Usage: coeqtl_limix_anno_chr_files_new.R [options]
#
#
# Options:
# 	-c CHARACTER, --cell_type=CHARACTER
# 		cell type working on
#
# 	-a CHARACTER, --annotation_prepend=CHARACTER
# 		base name and location of annotation files
#
# 	-f CHARACTER, --features_test_prepend=CHARACTER
# 		base name and location of file with features to test
#
# 	-h, --help
# 		Show this help message and exit
```

lastly we need to make sure that the SMF file is annotated correctly. If you used the wg3 pipeline, you might have an SMF that contains the batch in the phenotype_id column. However if you subsequently followed the co-eQTL pipeline, it will have created gene-gene correlations per sample, not sample and batch. You can generate a new smf using the following R code

```r
smf <- read.table('/groups/umcg-franke-scrna/tmp04/projects/venema-2022/ongoing/qtl/eqtl/sc-eqtlgen/input/elmentaite_adult_martin_immune/cell_type_medhigh_inflammationsplit_mincor/smf.txt', header = T, sep = '\t')
smf$phenotype_id <- smf$genotype_id
smf <- smf[!duplicated(smf$genotype_id), ]
write.table(smf, '/groups/umcg-franke-scrna/tmp04/projects/venema-2022/ongoing/qtl/eqtl/sc-eqtlgen/input/elmentaite_adult_martin_immune/cell_type_medhigh_inflammationsplit_mincor/smf_doubleid.tsv', row.names = F, col.names = T, quote = F, sep = '\t')
```

Now we can start modifying the yaml in the snakemake-to-map-coeqtls folder. You can copy the contents of the yaml into a text editor, clone the github repo and make your own branch, or just edit in it nano. 

#### yaml parameters for the co-eQTL mapping pipeline


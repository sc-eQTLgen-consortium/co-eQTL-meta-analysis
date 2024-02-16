# co-eQTL-meta-analysis

Scripts usage in this pipeline:
Step one is to generate files to allow snakemake to calculate all of the parallel jobs it will have to run.
This step uses the script: donor_ids.R

Step two is to separate the main seurat file into one file per individual (to allow for more parallelisation).
This step uses the script: create_donor_rds.R

Step three calculates the correlations for every gene pair within an individual.
This step uses the script: process_donor_rds.R

Step four is to aggregate all of these gene pair correlations files from one per individual into one file with all of the correlations.
This step uses the script: aggregate_metrics.py

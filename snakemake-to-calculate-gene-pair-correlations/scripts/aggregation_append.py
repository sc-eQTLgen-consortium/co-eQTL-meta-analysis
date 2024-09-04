###
# Code to aggregate co-eQTL correlations into 1 file in a memory efficient manner
###


"""
IMPORTS
"""

import pandas as pd
import gzip
import sys
import re

"""
Input Variables
"""

list_of_donors = re.sub(r'[\[\]]','',sys.argv[1]).split(',')
celltype = sys.argv[2]
output = sys.argv[3]
top_path = sys.argv[4]
wg3_path = sys.argv[5]
chromosome = str(sys.argv[6])

print(f"First 10 donors of list: {list_of_donors[0:10]}")
print(f"Cell Type: {celltype}")
print(f"Output file: {output}-chr-{chromosome}.tsv.gz")
print(f"Input path: {top_path}")
print(f"Location of annotation file: {wg3_path}/input/LimixAnnotationFile.txt")

"""
Code Body
"""
anno_input = f"{wg3_path}/input/LimixAnnotationFile.txt"

print("Loading chromosome-gene annotation file.")
anno_file = pd.read_csv(anno_input,sep='\t')
first_run = True
chr_genes = set([i for i in anno_file.feature_id[anno_file.chromosome.astype(str)==chromosome]])
outdir = f"{output}-chr-{chromosome}.tsv.gz"
handleout = gzip.open(outdir, "wt")

del anno_file

gene_pair = []
corr = []

# Load in correlation file per line
for ind_ID in list_of_donors:
    print(f"Processing {ind_ID}")
    file = f"{top_path}/{ind_ID}-{celltype}-pearson-weighted.tsv.gz"
    gene_pair = []
    corr = []
    counter = 0

    with gzip.open(file, "rt") as f:
        header = None
        for line in f:
            counter += 1
            if header == None:
                header = line
            else:
                values = line.strip().split("\t")
                gene1, gene2 = values[0].split('_')
                if gene1 in chr_genes:
                    gene_pair.append(values[0])
                    corr.append(values[1])
            if counter % 1000000 == 0:
                print(f"{counter} gene pairs processed", end="\r")
    if first_run:
        first_header = gene_pair
        handleout.write("\t" + "\t".join(gene_pair))
        handleout.write(f"\n{ind_ID}")
        handleout.write("\t" + "\t".join(corr))
        first_run = False

    else:
        if not first_header == first_header:
            raise Exception("Error: gene pairs are not in the same order across files")

        handleout.write(f"\n{ind_ID}")
        handleout.write("\t" + "\t".join(corr))
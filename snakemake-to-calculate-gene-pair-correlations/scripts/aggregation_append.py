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
annotation_path = sys.argv[5]
chromosome = str(sys.argv[6])

print(f"First 10 donors of list: {list_of_donors[0:10]}")
print(f"Cell Type: {celltype}")
print(f"Output file: {output}-chr-{chromosome}.tsv.gz")
print(f"Input path: {top_path}")
print(f"Location of annotation file: {annotation_path}")

"""
Code Body
"""
anno_input = f"{annotation_path}"

print("Loading chromosome-gene annotation file.")
anno_file = pd.read_csv(anno_input,sep='\t')
first_run = True
chr_genes = set([i for i in anno_file.feature_id[anno_file.chromosome.astype(str)==chromosome]])
print(f"Writing file for chromosome {chromosome}.")
outdir = f"{output}-chr-{chromosome}.tsv.gz"
handleout = gzip.open(outdir, "wt")

for ind_ID in list_of_donors:
  sample = pd.read_csv(f"{top_path}/{ind_ID}-{celltype}-pearson-weighted.tsv.gz",sep='\t')
  gene_pairs = [i for i in sample.gene_pair.values]
  if first_run:
    indexes = set([i for i in range(len(gene_pairs)) if gene_pairs[i].split('_')[0] in chr_genes or gene_pairs[i].split('_')[1] in chr_genes])
    header = [sample.gene_pair.iloc[i] for i in range(len(sample.gene_pair)) if i in indexes]
    first_header = header
    first_run = False
    for gene_pair in header:
      handleout.write(f"\t{gene_pair}")

  corr_values = [sample.correlation.iloc[i] for i in range(len(sample.correlation)) if i in indexes]
  header = [sample.gene_pair.iloc[i] for i in range(len(sample.gene_pair)) if i in indexes]

  print(ind_ID)
  print(first_header)
  print(header)

  if not header == first_header:
    raise Exception("Error: gene pairs are not in the same order across files")

  handleout.write(f"\n{ind_ID}")
  for value in corr_values:
    handleout.write(f"\t{value}")



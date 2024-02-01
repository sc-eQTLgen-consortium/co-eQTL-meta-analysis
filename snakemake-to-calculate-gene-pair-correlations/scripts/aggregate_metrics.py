"""
Aggregates files from donors for a particular cell type
"""

import pandas as pd 
import sys
import json
import argparse
import gzip

parser = argparse.ArgumentParser(description="")
parser.add_argument("--celltype", required=True, type=str, help="Cell type")
parser.add_argument("--cohort", required=True, type=str, help="Cohort ID")
parser.add_argument("--n", required=True, type=str, help="Number of genes")
parser.add_argument("--method", required=True, type=str, help="Correlation method")
parser.add_argument("--weight", required=True, type=str, help="Whether correlation is weighted or unweighted")
parser.add_argument("--metric", required=True, type=str, help="Correlation metric")
parser.add_argument("--input", required=True, nargs="+", type=str, help="Input correlation file")
parser.add_argument("--output", required=True, type=str, help="Aggregated output file")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

# Load donor data
donor_df = pd.read_csv(args.input[0], sep = "\t", compression='gzip')
original_ids = donor_df['original_ids']
alt_ids = donor_df['alt_ids']
del(donor_df)

corr_dict = {}
for i in range(len(alt_ids)):
    file = f"{args.input[1]}{args.cohort}/{alt_ids[i]}-{args.metric}-{args.celltype}-top-{args.n}-{args.method}-{args.weight}.tsv.gz"
    gene_list = []
    corr_list = []
    with gzip.open(file, 'rt') as f:
        next(f)
        for line in f:
            values = line.strip("\n").split("\t")
            gene_pair,corr = values[0],values[1]
            gene_list.append(gene_pair)
            corr_list.append(corr)

    corr_dict[original_ids[i]] = corr_list

df = pd.DataFrame(data=corr_dict,index=gene_list,columns=corr_dict.keys())

x=[i.split('_') for i in df.index.values]
y=[i.sort() for i in x]
y=[f"{i[0]}_{i[1]}" for i in x]
df.index=y
df.sort_index(inplace=True)

df.to_csv(args.output, sep='\t', na_rep='NA', compression='gzip')
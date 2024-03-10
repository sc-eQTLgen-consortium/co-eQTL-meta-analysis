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
parser.add_argument("--input", required=True, nargs="+", type=str, help="Input correlation file")
parser.add_argument("--output", required=True, nargs="+", type=str, help="Aggregated output file")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def save_metric(metric_dict,gene_list,outFile):
    df = pd.DataFrame(data=metric_dict,index=gene_list,columns=metric_dict.keys())
    x=[i.split('_') for i in df.index.values]
    y=[i.sort() for i in x]
    y=[f"{i[0]}_{i[1]}" for i in x]
    df.index=y
    df.sort_index(inplace=True)

    df.to_csv(outFile, sep='\t', compression='gzip')
    
# Load donor data
donor_df = pd.read_csv(args.input[0], sep = "\t", compression='gzip')
original_ids = donor_df['original_ids']
alt_ids = donor_df['alt_ids']
del(donor_df)

corr_dict = {}
std_dict = {}
pval_dict = {}
zscore_dict = {}

for i in range(len(alt_ids)):
    file = f"{args.input[1]}{args.cohort}/{alt_ids[i]}-{args.celltype}-top-{args.n}-{args.method}-{args.weight}.tsv.gz"
    gene_list = []
    corr = []
    std = []
    pval = []
    zscore = []
    
    with gzip.open(file, 'rt') as f:
        for line in f:
            values = line.strip("\n").split("\t")
            gene_list.append(values[0])
            corr.append(values[1])
            std.append(values[2])
            pval.append(values[3])
            zscore.append(values[4])
             
    corr_dict[original_ids[i]] = corr
    std_dict[original_ids[i]] = std
    pval_dict[original_ids[i]] = pval
    zscore_dict[original_ids[i]] = zscore
    
save_metric(corr_dict,gene_list,args.output[0])
save_metric(std_dict,gene_list,args.output[1])
save_metric(pval_dict,gene_list,args.output[2])
save_metric(zscore_dict,gene_list,args.output[3])
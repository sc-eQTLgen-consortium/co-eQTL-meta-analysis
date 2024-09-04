"""
Splits list of gene pairs used for correlation into n chunks
"""
import numpy as np
import pandas as pd 
import gzip
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--gene_list", required=True, type=str, help="List of all expressed genes")
parser.add_argument("--n_chunks", required=True, type=int, help="Number of chunks")
parser.add_argument("--cohort", required=True, type=str, help="Cohort")
parser.add_argument("--cell_type", required=True, type=str, help="Cell type")
parser.add_argument("--out_dir", required=True, type=str, help="Output directory")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

n_chunks = args.n_chunks # set to number of chunks

print("Loading gene lists")
genes_to_test = pd.read_csv(args.gene_list, sep='\t', header=None).iloc[:,0].to_list() # all genes
genes_to_test.sort() #

total_pairs = len(genes_to_test) * (len(genes_to_test) - 1) // 2
pairs_per_chunk = total_pairs // n_chunks # number of gene pairs per chunk
remainder = total_pairs % n_chunks

current_chunk = 0
pair_count = 0

fileout = gzip.open(f"{args.out_dir}gene-pairs-{args.cell_type}-{current_chunk+1}.tsv.gz", "wt", 4)

print(f"Writing gene pairs to {n_chunks} files")

for i in range(0,len(genes_to_test)):
    genei = genes_to_test[i] # Get gene

    for j in range(i+1, len(genes_to_test)):
        genej = genes_to_test[j]
        fileout.write(f"{genei}_{genej}\n")

        pair_count += 1

        if current_chunk == 0 and pair_count >= pairs_per_chunk + remainder:
            fileout.close()
            current_chunk += 1
            remainder -= 1 
            fileout = gzip.open(f"{args.out_dir}gene-pairs-{args.cell_type}-{current_chunk+1}.tsv.gz", "wt", 4)
            pair_count = 0
            print(f"{current_chunk}/{n_chunks} written",end='\r')

        elif current_chunk > 0 and pair_count >= pairs_per_chunk and current_chunk < n_chunks - 1:
            fileout.close()
            current_chunk += 1
            remainder -= 1 
            fileout = gzip.open(f"{args.out_dir}gene-pairs-{args.cell_type}-{current_chunk+1}.tsv.gz", "wt", 4)
            pair_count = 0
            print(f"{current_chunk}/{n_chunks} written",end='\r')

print("Process complete")
fileout.close()
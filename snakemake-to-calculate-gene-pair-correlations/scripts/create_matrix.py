"""
Creates a gene-gene matrix from the aggregated correlation metric files
"""
import numpy as np
import pandas as pd 
import argparse
import gzip

parser = argparse.ArgumentParser(description="")
parser.add_argument("--n", required=True, type=str, help="Number of genes")
parser.add_argument("--aggrmethod", required=True, type=str, help="Aggregation method")
parser.add_argument("--input", required=True, nargs="+", type=str, help="Input correlation file")
parser.add_argument("--output", required=True, type=str, help="Aggregated output file")
#parser.add_argument("--metric", required=True, type=str, help="Correlation metric")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

m = np.empty((int(args.n), int(args.n)), dtype=np.float64)
np.fill_diagonal(m,1) 

positions = {}
index_counter = 0

# Read in correlation file line by line
with gzip.open(args.input[0], 'rt') as f:
    next(f) 
    for line in f:
        values = line.strip("\n").split("\t")
        gene1,gene2 = values[0].split("_")
        corr = [float(i) for i in values[1:]]

        if gene1 not in positions:
            positions[gene1] = index_counter
            index_counter += 1
        index1 = positions[gene1]

        if gene2 not in positions:
            positions[gene2] = index_counter
            index_counter += 1
        index2 = positions[gene2]

        if args.aggrmethod == 'mean':
            corr_combined = np.array(corr).mean()
        elif args.aggrmethod == 'fishersz':
            # Fishersz transformation
        else: 
            raise ValueError("Unsupported aggregation method: {}".format(args.aggrmethod))

        m[index1, index2] = corr_combined
        m[index2, index1] = corr_combined

f.close()

genes = positions.keys()
pd.DataFrame(m,index=genes,columns=genes).to_csv(args.output, sep="\t",header=True,index=True,compression='gzip')
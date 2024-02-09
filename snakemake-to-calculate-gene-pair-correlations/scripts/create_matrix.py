"""
Creates a gene-gene matrix from the aggregated correlation metric files
"""
import numpy as np
import pandas as pd 
import argparse
import gzip

parser = argparse.ArgumentParser(description="")
parser.add_argument("--n", required=True, type=str, help="Number of genes")
parser.add_argument("--meta_analysis", required=True, type=str, help="Meta-analysis method used to combine correlation over donors")
parser.add_argument("--donor_list", required=True, type=str, help="Donor list path")
parser.add_argument("--input", required=True, nargs="+", type=str, help="Input correlation file")
parser.add_argument("--output", required=True, type=str, help="Aggregated output file")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def fishers_z_transformation(corr, donor_path = args.donor_list):
    donor_list = pd.read_csv(donor_path, sep = "\t", compression='gzip',header=0)
    df = pd.DataFrame({"Sample": list(donor_list['alt_ids']),
                       "Correlation": list(corr),
                       "N": list(donor_list['count'])})
    df = df.set_index("Sample")
    df["Effect size (Y)"] = 0.5 * np.log( (1 + df["Correlation"]) / (1 - df["Correlation"]))
    df["Variance (V)"] = 1 / (df["N"] - 3)
    df["Weight (W)"] = 1 / df["Variance (V)"]
    df["WY"] = df["Weight (W)"] * df["Effect size (Y)"]
    M = df["WY"].sum() / df["Weight (W)"].sum()
    r = (np.exp(2 * M) - 1) / (np.exp(2 * M) + 1)
    return r

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
        corr = values[1:]

        if gene1 not in positions:
            positions[gene1] = index_counter
            index_counter += 1
        index1 = positions[gene1]

        if gene2 not in positions:
            positions[gene2] = index_counter
            index_counter += 1
        index2 = positions[gene2]
        
        corr_numeric = np.genfromtxt(corr, dtype=float, missing_values='NA', 
                                     filling_values=np.nan)
               
        if args.meta_analysis == 'fishersz':
            corr_combined = fishers_z_transformation(corr_numeric, donor_path = args.donor_list)
        elif args.meta_analysis == 'mean':
            corr_combined = np.nanmean(corr_numeric)
        else: 
            raise ValueError(f"Invalid meta-analysis method: {args.meta_analysis}. Supported methods: 'mean', 'fishersz'")

        m[index1, index2] = corr_combined
        m[index2, index1] = corr_combined

f.close()

genes = positions.keys()
pd.DataFrame(m,index=genes,columns=genes).to_csv(args.output, sep="\t",header=True,index=True,compression='gzip')
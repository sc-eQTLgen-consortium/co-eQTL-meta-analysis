"""
Aggregates files from donors for a particular cell type
"""

import pandas as pd 
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

# Load donor data
donor_df = pd.read_csv(args.input[0], sep = "\t", compression='gzip')
original_ids = donor_df['original_ids']
alt_ids = donor_df['alt_ids']
del(donor_df)

# Get file handles 
file_handles = []
for i in range(len(alt_ids)):
    file = f"{args.input[1]}{args.cohort}/{alt_ids[i]}-{args.celltype}-top-{args.n}-{args.method}-{args.weight}.tsv.gz"
    try:
        file_handle = gzip.open(file, "rt") 
        file_handles.append(file_handle)
    except FileNotFoundError:
        print(f"{file} not found.")

output_files = []
for output in args.output:
    # In order of corr, std, pval,zscore
    output_files.append(gzip.open(output, "wt"))

# Write header
for outfile in output_files:
    outfile.write("\t")
    outfile.write("\t".join([f"{alt_id}" for alt_id in alt_ids]) + "\n")

line_count = 0
safety = (int(args.n) * int(args.n)) - int(args.n)
all_eof = False
while not all_eof:
    lines = []
    eof = 0
    
    corr_list = []
    std_list = []
    pval_list = []
    zscore_list = []
    
    for file_handle in file_handles:
        line = file_handle.readline()
        values = line.strip("\n").split("\t")
        
        gene = values[0]
        corr_list.append(values[1])
        std_list.append(values[2])
        pval_list.append(values[3])
        zscore_list.append(values[4]) 

        if not line:
            eof += 1
            continue
        
    output_files[0].write(f"{gene}\t" + "\t".join([f"{corr}" for corr in corr_list]) + "\n")
    output_files[1].write(f"{gene}\t" + "\t".join([f"{std}" for std in std_list]) + "\n")
    output_files[2].write(f"{gene}\t" + "\t".join([f"{pval}" for pval in pval_list]) + "\n")
    output_files[3].write(f"{gene}\t" + "\t".join([f"{zscore}" for zscore in zscore_list]) + "\n")
    
    line_count += 1
    
    if eof == len(file_handles) or line_count > safety:
        all_eof = True
        break
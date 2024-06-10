"""
Aggregates files from donors for a particular cell type
"""

import pandas as pd 
import argparse
import gzip

parser = argparse.ArgumentParser(description="")
parser.add_argument("--celltype", required=True, type=str, help="Cell type")
parser.add_argument("--cohort", required=True, type=str, help="Cohort ID")
parser.add_argument("--n", required=True, type=int, help="Number of genes")
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

# Get file handles for each donor
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

# Set end of file indicator to False
all_eof = False

# Process every line for each input file simultaneously
max_n_lines = int((((args.n * args.n) / 2) - args.n) + 100)  # Doing 100 extra for safety
for _ in range(max_n_lines):
    eof_files = 0
    
    corr_list = []
    std_list = []
    pval_list = []
    zscore_list = []
    
    active_gene = None
    for file_handle in file_handles:
        line = file_handle.readline()
        if not line:
            eof_files += 1
                    
            corr_list.append("")
            std_list.append("")
            pval_list.append("")
            zscore_list.append("") 
            continue
        
        values = line.strip("\n").split("\t")
        
        gene = values[0]
        if active_gene is None:
            active_gene = gene
        else:
            if gene != active_gene:
                print("this is bad")  # edit this
                exit()
        
        corr_list.append(values[1])
        std_list.append(values[2])
        pval_list.append(values[3])
        zscore_list.append(values[4]) 

    if eof_files == len(file_handles):
        all_eof = True
        break
    # Write corr, std, pval and zscore to output files    
    output_files[0].write(f"{active_gene}\t" + "\t".join(corr_list) + "\n")
    output_files[1].write(f"{active_gene}\t" + "\t".join(std_list) + "\n")
    output_files[2].write(f"{active_gene}\t" + "\t".join(pval_list) + "\n")
    output_files[3].write(f"{active_gene}\t" + "\t".join(zscore_list) + "\n")

for outfile in output_files:
    outfile.close()
"""
Aggregates correlation per donors over n gene pairs
"""

import gzip
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--n_chunks", required=True, type=int, help="Number of chunks")
parser.add_argument("--donor", required=True, type=str, help="Donor id")
parser.add_argument("--cell_type", required=True, type=str, help="Cell type")
parser.add_argument("--output", required=True, type=str, help="Output file path")
parser.add_argument("--outdir", required=True, type=str, help="Output directory")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

n_chunks = args.n_chunks
fileHandles = []
for n in range(1,n_chunks+1):
    fileHandles.append(f"{args.outdir}/{args.donor}-{args.cell_type}-pearson-weighted-{n}.tsv.gz")

# Write header to file
fileOut = gzip.open(args.output, "wt", 4)
fileOut.write(f"gene_pair\tcorrelation\n")

for file_path in fileHandles:
    print(f"Processing {file_path}",end='\r')
    header = None

    with gzip.open(file_path,"rt") as f:
        for line in f:
            if header == None:
                gene_pair = line.strip().split("\t")
                header = True
            else:
                corr = line.strip().split("\t")

    # Write to output file
    for i in range(len(corr)):
        fileOut.write(f"{gene_pair[i]}\t{corr[i]}\n")

fileOut.close()
print("Done")
import gzip
import argparse
import pandas as pd 

parser = argparse.ArgumentParser(description="")
parser.add_argument("--outliers", required=True, type=str, help="Outlier table")
parser.add_argument("--outdir", required=True, type=str, help="Output directory")
parser.add_argument("--cohort", required=True, type=str, help="Cohort")
parser.add_argument("--celltype", required=True, type=str, help="Cell type")
parser.add_argument("--output", required=True, type=str, help="Output correlation file")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

# Load outiers df
outliers = gene_list = pd.read_csv(args.outliers, sep="\t", compression='gzip', header=None).iloc[:, 0].tolist()

# Load and write correlation matrix with outliers removed
fileout = gzip.open(args.output, "wt")
for i in range(1,23):
    file = f"{args.outdir}{args.cohort}/{args.cohort}-{args.celltype}-pearson-weighted-chr-{i}-final.tsv.gz"
    header = None
    try:
        with gzip.open(file,"rt") as f:
            for line in f:
                values = line.strip().split("\t")
                if header == None:
                    header = values[1:]
                    index = [idx for idx, id in enumerate(header) if id not in outliers]
                    header = [header[idx] for idx in index]
                    fileout.write("\t".join([""] + header) + "\n")
                else:
                    gene_pair = values[0]
                    corr = values[1:]
                    corr = [corr[idx] for idx in index]
                    fileout.write(gene_pair + "\t".join([""] + corr) + "\n")
    except FileNotFoundError:
        print(f"{file} not found")
fileout.close()
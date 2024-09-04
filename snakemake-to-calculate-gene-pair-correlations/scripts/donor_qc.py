"""
Plots the cell count, gene count and performs PCA over the donors to identify outliers
"""

import argparse
import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random
from pca import PCA

parser = argparse.ArgumentParser(description="")
parser.add_argument("--outdir", required=True, type=str, help="Output directory")
parser.add_argument("--donor_counts", required=True, type=str, help="Donor counts")
parser.add_argument("--color", required=True, type=str, help="Color")
parser.add_argument("--cohort", required=True, type=str, help="Cohort")
parser.add_argument("--celltype", required=True, type=str, help="Cell type")
parser.add_argument("--output", required=True, type=str, help="Output filepaths")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

## LOAD DATA
# Load cell counts
cell_count = pd.read_csv(args.donor_counts, sep="\t", header=0)
cell_count = cell_count.sort_values(by='Freq', ascending=False)
# Sort donor ids
sorted_ids = cell_count['Var1'].tolist()

subset = False
# Get gene counts
gene_dict = {}
for id in sorted_ids:
    file = f"{args.outdir}{args.cohort}/donor_gene_list/filtered-genes-{id}-{args.celltype}.tsv.gz"
    try:
        gene_list = pd.read_csv(file, sep="\t", compression='gzip', header=None).iloc[:, 0].tolist()
        gene_dict[id] = len(gene_list)
        if len(gene_list) > 4000:
            # Set to true if gene list is too large
            print("SUBSETTING DATA")
            subset = True
    except FileNotFoundError:
        print(f"{file} not found")


subset_size = 1000000
selected_lines = []
counter = 0

# Get correlations
corr_dict = {}
for i in range(1,23):
    file = f"{args.outdir}{args.cohort}/{args.cohort}-{args.celltype}-pearson-weighted-chr-{i}-final.tsv.gz"
    try:
        with gzip.open(file,"rt") as f:
            header = f.readline().strip().split("\t")[1:]
            for id in header:
                if id not in corr_dict:
                    corr_dict[id] = []

            for line in f:
                if subset:
                    if random.random() >= 0.5:
                        values = line.strip().split("\t")
                        corr = [float(val) if val != "NaN" else np.nan for val in values[1:]]
                        for idx, donor_id in enumerate(header):
                            corr_dict[donor_id].append(corr[idx])
                        counter+=1
                        if counter >= subset_size:
                            break

                else:
                    values = line.strip().split("\t")
                    corr = [float(val) if val != "NaN" else np.nan for val in values[1:]]
                    for idx, donor_id in enumerate(header):
                        corr_dict[donor_id].append(corr[idx])

    except FileNotFoundError:
        print(f"{file} not found")

corr_dict = {donor_id: corr_dict[donor_id] for donor_id in sorted_ids if donor_id in corr_dict}
corr_df = pd.DataFrame(corr_dict)
print(f"Correlation df shape: {corr_df.shape}")

## PLOTTING
# Create barplot cell counts per donor
plt.figure(figsize=(10, 6))
sns.set(style="white")
sns.barplot(x='Var1', y='Freq', data=cell_count, color=args.color)
plt.xlabel("Donors")
plt.ylabel("Cell counts")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f"{args.outdir}{args.cohort}/{args.cohort}-{args.celltype}-cell-counts.png", bbox_inches='tight')

# Plot barplot of gene counts per donor
plt.figure(figsize=(10, 6))
sns.set(style="white")
sns.barplot(x=list(gene_dict.keys()), y=list(gene_dict.values()), color=args.color)
plt.yscale('log')
plt.xlabel("Donors")
plt.ylabel("Gene counts")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f"{args.outdir}{args.cohort}/{args.cohort}-{args.celltype}-gene-counts.png", bbox_inches='tight')

# Plot boxplot of correlation distribution per donor
plt.figure(figsize=(10, 6))
sns.set(style="white")
sns.boxplot(data=corr_df, color=args.color, width=0.6, orient="v")
plt.xlabel("Donors")
plt.ylabel("Correlation coefficient")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f"{args.outdir}{args.cohort}/{args.cohort}-{args.celltype}-correlation-distribution.png", bbox_inches='tight')

##QUALITY CONTROL STEPS
# Removing samples with less than 20 cells
outlier_dict = {}
low_counts = cell_count[cell_count['Freq'] < 20]["Var1"].tolist()
print(f"low counts: {low_counts}")
print(corr_df.columns)
for x in low_counts:
    outlier_dict[x] = "< 20 cells"
    if x in corr_df.columns:
        corr_df.drop(columns=x, inplace=True)

# Performing pca to identify outliers
corr_df = corr_df.apply(pd.to_numeric, errors='coerce')
pca = PCA(corr_df)
pca.plot_pc1_vs_pc2(color="#808080",outfile=f"{args.outdir}{args.cohort}/{args.cohort}-{args.celltype}-pc1-pc2.png")
outliers = pca.return_outliers(threshold=3)

for x in outliers:
    outlier_dict[x] = "Outlier (z-score > 3)"
    corr_df.drop(columns=x, inplace=True)

outlier_df = pd.DataFrame.from_dict(outlier_dict, orient='index', columns=['Sample Status'])
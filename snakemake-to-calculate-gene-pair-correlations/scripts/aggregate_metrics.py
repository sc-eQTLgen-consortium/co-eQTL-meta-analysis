"""
Aggregates files from donors for a particular cell type
"""

import pandas as pd 
import sys
import json

list_of_donors = json.loads(sys.argv[1])
original_ids = json.loads(sys.argv[2])
cell_type = sys.argv[3]
metric = sys.argv[4]
cohort = sys.argv[5]
outfile = sys.argv[6]
infile = sys.argv[7]
min_samples = sys.argv[8]

print(f"First 10 donors of list: {list_of_donors[0:10]}")
print(f"Cell Type: {cell_type}")
print(f"Metric: {metric}")
print(f"Cohort: {cohort}")
print(f"Output directory: {outfile}")
print(f"Input directory: {infile}")
print(f"Input file: {infile}/{metric}-{cell_type}-ind-pearson-weighted.tsv.gz")
print(f"minimum number of correlations allowed per gene pair: {min_samples}")

df=pd.DataFrame()

for i in range(len(list_of_donors)):
  file=f"{infile}/{metric}-{cell_type}-{list_of_donors[i]}-pearson-weighted.tsv.gz"
  df[f"{original_ids[i]}"] = pd.read_csv(file, sep='\t')

print("Sorting all gene pairs to be in alphabetical order")
x=[i.split('_') for i in df.index.values]
y=[i.sort() for i in x]
y=[f"{i[0]}_{i[1]}" for i in x]
df.index=y
df.sort_index(inplace=True)

print(f"Dropping individuals with no values and gene pairs that do not have at least {min_samples} values")
df.dropna(how='all',axis=1)
df.dropna(thresh=min_samples,axis=0)

df.to_csv(outfile, sep='\t', na_rep='NA', compression='gzip')


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

print(f"First 10 donors of list: {list_of_donors[0:10]}")
print(f"Cell Type: {cell_type}")
print(f"Metric: {metric}")
print(f"Cohort: {cohort}")
print(f"Output directory: {outfile}")
print(f"Input directory: {infile}")
print(f"Input file: {infile}/{metric}-{cell_type}-ind-pearson-weighted.tsv.gz")

df=pd.DataFrame()

index_value=0
for i in range(len(list_of_donors)):
  file=f"{infile}/{metric}-{cell_type}-{list_of_donors[i]}-pearson-weighted.tsv.gz"
  ind=pd.read_csv(file, sep='\t')
  df[f"{original_ids[i]}"] = [val for val in ind['x']]
  try:
    if index_value == 0:
      index_value=ind.index.values
  except:
    continue


df.index=index_value
df.to_csv(outfile, sep='\t', na_rep='NA', compression='gzip')


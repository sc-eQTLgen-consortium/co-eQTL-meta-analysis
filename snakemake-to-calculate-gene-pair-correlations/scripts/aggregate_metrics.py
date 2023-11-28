"""
Aggregates files from donors for a particular cell type
"""

import pandas as pd 
import sys
import json

list_of_donors = json.loads(sys.argv[1])
cell_type = sys.argv[2]
metric = sys.argv[3]
outfile = sys.argv[4]
infile = sys.argv[5]

print(f"First 10 donors of list: {list_of_donors[0:10]}")
print(f"Cell Type: {cell_type}")
print(f"Metric: {metric}")
print(f"Output directory: {outfile}")
print(f"Input directory: {infile}")
print(f"Input file: {infile}/{metric}_{cell_type}_ind_pearson_weighted.tsv.gz")

df=pd.DataFrame()

index_value=0
for i in list_of_donors:
  file=f"{infile}/{metric}_{cell_type}_{i}_pearson_weighted.tsv.gz"
  ind=pd.read_csv(file, sep='\t')
  df[f"{i}"] = [val for val in ind['x']]
  try:
    if index_value == 0:
      index_value=ind.index.values
  except:
    continue


df.index=index_value
df.to_csv(outfile, sep='\t', na_rep='NA', compression='gzip')


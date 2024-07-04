"""
This script creates empty files when file does not exist already
"""

import sys

list_of_donors = sys.argv[1].split(',')
output_path = sys.argv[2]
cell_type = sys.argv[3]

print(list_of_donors[0:10])
print(f"out: {output_path}")
print(f"cell type: {cell_type}")

for donor in list_of_donors:
  quiet = open(f"{output_path}/counts/normalized-counts-{donor}-{cell_type}.mtx", 'w+')
  quiet = open(f"{output_path}/donor_weight/correlation-weight-{donor}-{cell_type}.tsv.gz", 'w+')
  quiet = open(f"{output_path}/donor_gene_list/filtered-genes-{donor}-{cell_type}.tsv.gz", 'w+')


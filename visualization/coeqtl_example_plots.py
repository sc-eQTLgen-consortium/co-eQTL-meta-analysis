### Plot co-eQTL examples

import pandas as pd
import os
import sys
import glob
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib
import numpy as np
from matplotlib.ticker import FuncFormatter
import pyranges as pr
from tqdm import tqdm


#####
# Getting gene info from GTF file
#####


gtf = pr.read_gtf("Downloads/genes.gtf.gz")
gtf_df = gtf.df

genes = ['AP1S3','RNASET2','RNASET2','SH3YL1']
snps = ['2:223824522:A:G','6:166974681:A:G','6:166980633:C:T','2:228088:G:A']

genes = ['AP1S3']
snps = ['2:223824522:A:G']

var_name = snps[0]
gene_name = genes[0]
gene_df = gtf_df[gtf_df.gene_name==gene_name].copy()
chromosome = int("".join("".join(gene_df.Chromosome.unique()).split("chr")))
if chromosome not in range(1,23):
  raise Exception(f"Chromosome {chromosome} not in normal range, check value, MT,X,Y is fine")

if len(gene_df[gene_df.Feature=="gene"].Start.unique()) != 1:
  print("Multiple gene start locations detected!")

if len(gene_df[gene_df.Feature=="gene"].End.unique()) != 1:
  print("Multiple gene end locations detected!")

if len(gene_df[gene_df.Feature=="gene"].Strand.unique()) != 1:
  print("Multiple gene end locations detected!")

start = int(gene_df[gene_df.Feature=="gene"].Start.unique()[0])
end = int(gene_df[gene_df.Feature=="gene"].End.unique()[0])
strand = gene_df[gene_df.Feature=="gene"].Strand.unique()[0]

exons = {}
for exon_id in tqdm(gene_df.exon_id.unique()):
  if exon_id is np.nan:
    continue
  exon_start = gene_df[np.array(gene_df.exon_id==exon_id) & np.array(gene_df.Feature=="exon")].Start.unique()
  exon_end = gene_df[np.array(gene_df.exon_id==exon_id) & np.array(gene_df.Feature=="exon")].End.unique()
  if len(exon_start) == 1:
    exons[exon_id] = [int(exon_start[0]), int(exon_end[0])]
  else:
    raise Exception(f"Error: There are multiple exon start sites for exon: {exon_id}")


def format_xticks(value, _):
    return f'{value / 1000:.0f}'


### Plot a single gene
fig = plt.figure(figsize=(12, 8)) 
ax = fig.subplots()
gene_size = end-start
y=4
font_weight = 'bold'
ax.hlines(y=y, xmin=start, xmax=end, colors='black', linestyles='solid', linewidth=3)
ax.vlines(x=start, ymin=y-1, ymax=y+1, colors='black', linestyles='solid', linewidth=1)
ax.vlines(x=end, ymin=y-1, ymax=y+1, colors='black', linestyles='solid', linewidth=1)
ax.text((start + end) / 2, y + 1, gene_name, ha='center', va='bottom', fontsize=7, fontweight=font_weight, fontstyle='italic')

### Plot exons if necessary
for exon in exons:
  ax.hlines(y=y, xmin=exons[exon][0], xmax=exons[exon][1], colors='blue', linestyles='solid', linewidth=20)


### Plotting gene start

real_start = start
if strand == "-":
  real_start = end

ax.vlines(x=real_start, ymin=y+2, ymax=y+4, colors='green', linestyles='solid', linewidth=1)
ax.hlines(y=y+4, xmin=eval(f"real_start{strand}1000"), xmax=real_start, colors='green', linestyles='solid', linewidth=1)
path = [[eval(f"real_start{strand}1000"), y+3.5],[eval(f"real_start{strand}1000"), y+4.5],[eval(f"real_start{strand}3000"), y+4]]
ax.add_patch(patches.Polygon(path,facecolor="green"))


##########
### Plot variant(s)
##########

var_name = snps[0]
gene_name = genes[0]
var_position = int(var_name.split(":")[1])
dist_from_gene = min([abs(start-var_position),abs(end-var_position)])

if var_position >= start and var_position <= end:
  dist_from_gene = 0


ax.vlines(x=var_position, ymin=y-0.5, ymax=y+0.5, colors='red', linestyles='solid', linewidth=1.5)
ax.text(var_position, y - 2, var_name, ha='center', va='bottom', fontsize=7, fontweight=font_weight, fontstyle='italic')


## add extra variants to plot if wanted

extra_vars = ['6:166982169:CTCTTTCCCCTTCTTCCT:C','6:166963191:TG:T','6:166956409:C:A']
offset = 0.5
for var_name in tqdm(extra_vars):
  new_dist = min([abs(start-var_position),abs(end-var_position)])
  if new_dist > dist_from_gene:
    dist_from_gene = new_dist
  offset += 0.5
  var_position = int(var_name.split(":")[1])
  ax.vlines(x=var_position, ymin=y-0.5, ymax=y+0.5, colors='green', linestyles='solid', linewidth=1.5)
  ax.text(var_position, y - (2+offset), var_name, ha='center', va='bottom', fontsize=7, fontweight=font_weight, fontstyle='italic')



### Tweak plot size
start_of_plot = start-(dist_from_gene+25000)
end_of_plot = end+(dist_from_gene+25000)
ax.set_xlim((start_of_plot), (end_of_plot))
#ax.hlines(y=y, xmin=start-(dist_from_gene+90000), xmax=end+(dist_from_gene+90000), colors='black', linestyles='solid', linewidth=1)
ax.set_xlabel(f'Position on chromosome {chromosome} (kb)', fontsize=10)
#ax.set_xticks(ax2.get_xticks())
#ax.set_xticks([ax2.get_xticks()[0], ax2.get_xticks()[0] + abs(((ax2.get_xticks()[0] - ax2.get_xticks()[-1]) / 2)), ax2.get_xticks()[-1]])
ax.set_ylim(0, 50)
ax.set_yticks([])

# Create a formatter
formatter = FuncFormatter(format_xticks)

# Apply formatter to all right-side plots
ax.xaxis.set_major_formatter(formatter)

# Add line for DNA backbone
ax.hlines(y=y, xmin=start-(dist_from_gene+25000), xmax=end+(dist_from_gene+25000), colors='black', linestyle='solid', linewidth=1)


####################
# Save plot as PDF #
####################

plt.savefig(f"Downloads/{gene_name}_{var_name.replace(':','_')}_gene_var_location.pdf")
plt.cla()
plt.close()

print(f"Plot start at position: {start_of_plot} and ends at position: {end_of_plot}")







####################
# Determine SNP to gene distance + if it is in an exon
####################


gtf = pr.read_gtf("Downloads/genes.gtf.gz")
gtf_df = gtf.df
sel = pd.read_excel("Downloads/coeqtl_examples_full_table_annotations_selection.xlsx")
sel.drop_duplicates(subset=["snp_id","eGene"],inplace=True,ignore_index=True)
man_sel = pd.DataFrame()
man_sel["eGene"] = ["CMTM8","JAZF1","CD164"]
man_sel["snp_id"] = ["3:32304282:C:G","7:28117268:C:T","6:109380015:T:C"]

results_table = pd.read_csv("Downloads/coeqtl_examples_full_table_annotations.tsv.gz",sep='\t')
man_sel = results_table.copy()
man_sel = man_sel[["eGene","snp_id","snp_egene"]]
man_sel.drop_duplicates(subset="snp_egene",inplace=True,ignore_index=True)

gtf_df_sub = gtf_df[gtf_df.gene_name.isin(man_sel.eGene)]

dist_from_gene_col = []
in_exon_col = []
for row in tqdm(range(len(man_sel.eGene))):
  gene_name = man_sel.eGene.iloc[row]
  if "_" in gene_name:
    dist_from_gene_col.append(np.nan)
    in_exon_col.append(np.nan)
    continue
  gene_df = gtf_df_sub[gtf_df_sub.gene_name==gene_name].copy()
  chromosome = int("".join("".join(gene_df.Chromosome.unique()).split("chr")))
  start = int(gene_df[gene_df.Feature=="gene"].Start.unique()[0])
  end = int(gene_df[gene_df.Feature=="gene"].End.unique()[0])
  strand = gene_df[gene_df.Feature=="gene"].Strand.unique()[0]
  var_name = man_sel.snp_id.iloc[row]
  var_position = int(var_name.split(":")[1])
  dist_from_gene = [abs(start-var_position),abs(end-var_position)]
  dist_from_gene = dist_from_gene.index(min(dist_from_gene))
  dist_from_gene = [var_position-start,var_position-end][dist_from_gene]
  if var_position >= start and var_position <= end:
    dist_from_gene = 0
  dist_from_gene = eval(f"{dist_from_gene} * {strand}1")
  dist_from_gene_col.append(dist_from_gene)
  in_exon = False
  if dist_from_gene == 0:
    exons = {}
    for exon_id in gene_df.exon_id.unique():
      if exon_id is np.nan:
        continue
      exon_start = gene_df[np.array(gene_df.exon_id==exon_id) & np.array(gene_df.Feature=="exon")].Start.unique()
      exon_end = gene_df[np.array(gene_df.exon_id==exon_id) & np.array(gene_df.Feature=="exon")].End.unique()
      if len(exon_start) == 1:
        exons[exon_id] = [int(exon_start[0]), int(exon_end[0])]
      else:
        raise Exception(f"Error: There are multiple exon start sites for exon: {exon_id}")
    for exon in exons:
      if var_position >= exons[exon][0] and var_position <= exons[exon][1]:
        in_exon = True
      if in_exon:
        break
  #tqdm.write(f"{gene_name} {strand} {var_name} {dist_from_gene} in exon: {in_exon}")
  in_exon_col.append(int(in_exon))
  

man_sel["dist_from_gene"] = dist_from_gene_col
man_sel["in_exon"] = in_exon_col

my_dict = {}
for i in tqdm(range(len(man_sel.snp_egene))):
  snp_egene = man_sel.snp_egene.iloc[i]
  dist = man_sel.dist_from_gene.iloc[i]
  in_exon = man_sel.in_exon.iloc[i]
  my_dict[snp_egene] = {"dist_from_gene":float(dist), "in_exon":float(in_exon)}


results_table["dist_from_gene"] = [my_dict[i]["dist_from_gene"] for i in tqdm(results_table.snp_egene)]
results_table["in_exon"] = [my_dict[i]["in_exon"] for i in tqdm(results_table.snp_egene)]



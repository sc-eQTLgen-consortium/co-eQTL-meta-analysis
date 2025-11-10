
##########
# GWAS co-eQTL coloc.
##########

import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import pyranges as pr
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


##########
# co-eQTL plotting
##########

df = pd.read_csv("/scratch/hb-sceqtlgen/GRN-Project/finemapped_coeqtl_results/joint_co_eQTL_fine_mapping_CD4_T_SF_coloc_20250411/meta_joint_susie_finemap_CD4_T_formatted.tsv.gz",sep='\t')
#df = pd.read_csv("/scratch/hb-sceqtlgen/GRN-Project/finemapped_coeqtl_results/joint_co_eQTL_fine_mapping_Mono_SF_coloc_20250411/meta_joint_susie_finemap_Mono_formatted.tsv.gz",sep='\t')
#df = pd.read_csv("/scratch/hb-sceqtlgen/GRN-Project/finemapped_coeqtl_results/joint_co_eQTL_fine_mapping_CD8_T_SF_coloc_20250411/meta_joint_susie_finemap_CD8_T_formatted.tsv.gz",sep='\t')
#df = pd.read_csv("/scratch/hb-sceqtlgen/GRN-Project/finemapped_coeqtl_results/joint_co_eQTL_fine_mapping_B_SF_coloc_20250411/meta_joint_susie_finemap_B_formatted.tsv.gz",sep='\t')

gene_pair = "ID2_IL7R"
snp_id = "5:35877403:G:A"
cell_type = "CD4_T"

region_df = df[df.feature_id==gene_pair].copy()
chromosome = int(region_df.snp_chromosome.unique())
region_start = min(list(region_df.snp_position.values))
region_end = max(list(region_df.snp_position.values))

snp_df = region_df[region_df.snp_id==snp_id]

gtf = pr.read_gtf("/scratch/hb-sceqtlgen/GRN-Project/finemapped_coeqtl_results/genes.gtf.gz")
gtf_df = gtf.df

gtf_df = gtf_df[(gtf_df.Feature=="gene") & (gtf_df.Chromosome==f"chr{chromosome}") & (gtf_df.gene_type=="protein_coding")]
gtf_df = gtf_df[((gtf_df.Start>=region_start)&(gtf_df.Start<=region_end)) | ((gtf_df.End>=region_start)&(gtf_df.End<=region_end))]
gtf_df.sort_values(by="Start",inplace=True,ignore_index=True)



# Plotting code

fig, ax = plt.subplots()

x_vals = list(region_df.snp_position.values)
y_vals = list(-np.log10(region_df.p_value.values))
snp_y = list(-np.log10(snp_df.p_value.values)) # allows highlighting of specific SNP

ax.scatter(x_vals, y_vals, s=1, color='gray')
ax.scatter(list(snp_df.snp_position.values), snp_y, s=1, color='blue')

step = int((region_end - region_start) / 3)
ax.set_xticks([i for i in range(region_start, region_end+1, step)])
ax.set_xticklabels([round(i/1000000,1) for i in range(region_start, region_end+1, step)])

step = 5
ax.set_yticks([i for i in range(0, int(round(max(y_vals),-1))+1, step)])
ax.set_yticklabels([int(i)for i in range(0, int(round(max(y_vals),-1))+1, step)])

ax.set_xlabel(f"Position chromosome {chromosome} (Mb)")
ax.set_ylabel("-log10 P-value")
ax.set_title(f"co-eQTL {cell_type}: {gene_pair}")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)



# Plot surrounding genes


prev_gene_end = 0
y = max(y_vals)*0.05
y_offset = y*0.25
for row in range(len(gtf_df.gene_name)):
  edge_of_plot_start = False #variable to show gene is cut off
  edge_of_plot_end = False #variable to show gene is cut off
  gene_start = gtf_df.Start.iloc[row]
  if gene_start < region_start:
    gene_start = region_start
    edge_of_plot_start = True
  gene_end = gtf_df.End.iloc[row]
  if gene_end > region_end:
    gene_end = region_end
    edge_of_plot_end = True
  gene_name = gtf_df.gene_name.iloc[row]
  if gene_start <= prev_gene_end:
    y += -1
  else:
    y = -2
  prev_gene_end = gene_end
  ax.hlines(y=y+y_offset, xmin=gene_start, xmax=gene_end, colors='dimgrey', linestyles='solid', linewidth=1)
  if not edge_of_plot_start:
    ax.vlines(x=gene_start, ymin=y, ymax=y+(2*y_offset), colors='dimgrey', linestyles='solid', linewidth=0.2)
  if not edge_of_plot_end:
    ax.vlines(x=gene_end, ymin=y, ymax=y+(2*y_offset), colors='dimgrey', linestyles='solid', linewidth=0.2)
  ax.text((gene_start + gene_end) / 2, y-(2*y_offset), gene_name, ha='center', va='bottom', fontsize=2, fontweight='bold', fontstyle='italic')


#fig.set_size_inches(12, 6)
fig.set_size_inches(6, 2)
plt.savefig(f"pval_coloc_{gene_pair}_{cell_type}_resize.pdf")
plt.cla()
plt.close()

print(f"file saved to: {os.getcwd()}/pval_coloc_{gene_pair}_{cell_type}_resize.pdf")


##########
# GWAS plotting
##########

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

## Things to be aware of:
# Load correct coeqtl file, check which cell type it needs to be
# Edit gene pair / celltype variables


#df_coeqtl = pd.read_csv("Downloads/meta_joint_susie_finemap_B_formatted.tsv.gz",sep='\t')
df_coeqtl = pd.read_csv("Downloads/meta_joint_susie_finemap_CD4_T_formatted.tsv.gz",sep='\t')

#df = pd.read_csv("Downloads/GWAS_coloc_files/crohns_RNASET2-EFO_0000384.tsv.gz",sep='\t')
#df = pd.read_csv("Downloads/GWAS_coloc_files/GCST90267986_Cystatin_AP1S3.h.tsv.gz",sep='\t')
df = pd.read_csv("Downloads/finngen_R12_L12_DERMATITISNAS.tsv",sep='\t')

gene_pair = "ID2_IL7R"
snp_id = "5:35877403:G:A"
region_df = df_coeqtl[df_coeqtl.feature_id==gene_pair].copy()
chromosome = int(region_df.snp_chromosome.unique())
region_start = min(list(region_df.snp_position.values))
region_end = max(list(region_df.snp_position.values))

##For most GWAS catalog sum stats
#region_df = df[(df.hm_chrom==chromosome) & (df.hm_pos>=region_start) & (df.hm_pos<=region_end)].copy()
#region_df['hm_variant_id'] = [i.replace('_',':') for i in region_df.hm_variant_id]
#region_df=region_df[region_df.hm_variant_id.isin(df_coeqtl.snp_id)]
#snp_df = region_df[region_df.hm_variant_id==snp_id]

##For Finngen file
region_df = df[(df["#chrom"]==chromosome) & (df.pos>=region_start) & (df.pos<=region_end)].copy()
region_df['hm_variant_id'] = df[["#chrom","pos","ref","alt"]].astype(str).agg(':'.join, axis=1)
region_df=region_df[region_df.hm_variant_id.isin(df_coeqtl.snp_id)]
snp_df = region_df[region_df.hm_variant_id==snp_id]



# Plotting code

fig, ax = plt.subplots()

#x_vals = list(region_df.hm_pos.values)
x_vals = list(region_df.pos.values)

#y_vals = list(-np.log10(region_df.p_value.values))
y_vals = list(-np.log10(region_df.pval.values))

#snp_y = list(-np.log10(snp_df.p_value.values))
snp_y = list(-np.log10(snp_df.pval.values))

# to re-scale y values
#scales = {'original':{'lower':min(y_vals), 'upper':max(y_vals)}, 'desired':{'lower':0, 'upper':15}}

#y_scaled = [scales['desired']['lower'] + (x - scales['original']['lower']) * (scales['desired']['upper'] - scales['desired']['lower']) / (scales['original']['upper'] - scales['original']['lower']) for x in y_vals]
#snp_scaled = [scales['desired']['lower'] + (x - scales['original']['lower']) * (scales['desired']['upper'] - scales['desired']['lower']) / (scales['original']['upper'] - scales['original']['lower']) for x in snp_y]

ax.scatter(x_vals, y_vals, s=1, color='gray')
if len(snp_df) > 0:
  ax.scatter(list(snp_df.pos.values), snp_y, s=1, color='blue')

step = int((region_end - region_start) / 3)
ax.set_xticks([i for i in range(region_start, region_end+1, step)])
ax.set_xticklabels([round(i/1000000,1) for i in range(region_start, region_end+1, step)])

step = 5
ax.set_yticks([i for i in range(0, int(round(max(y_vals),-1))+1, step)])
ax.set_yticklabels([int(i)for i in range(0, int(round(max(y_vals),-1))+1, step)])

ax.set_xlabel(f"Position chromosome {chromosome} (Mb)")
ax.set_ylabel("-log10 P-value")
ax.set_title(f"GWAS: Dermatitis finngen {gene_pair}")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

#fig.set_size_inches(12, 6)
fig.set_size_inches(6, 2)
plt.savefig(f"Downloads/pval_gwas_coloc_{gene_pair}_resize_2.pdf")
plt.cla()
plt.close()











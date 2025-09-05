### co-eQTL GWAS enrichment

## LD done on Habrok
# see files at: /scratch/hb-sceqtlgen/GRN-Project/sceqtlgen_gwas/

import pandas as pd
import random
import os
from tqdm import tqdmx	
import time
import scipy.stats as stats

#######################
### Steps in script ###
#######################

## Load variants and filter
# filter out none eQTL SNPs from co-eQTLs
# filter out co-eQTL SNPs from eQTLs
# filter out HLA region
# filter based on LD
## calculate enrichment for co-eQTLs over eQTLs
# all traits
# all traits - cell count
# disease traits
# immune traits
## same as above but per cell type

######################
### Required files ###
######################

#set of eQTL variants:       sceqtlgen_sig_vars_no_dc.txt.gz
#set of coeQTL variants:     5DS_Meta_Analysis_Sign_coeQTLs133p_val_mt_eGene_.csv
#LD between all variants:    sceqtlgen_all_vars_tested_ld_out_01LD.tsv.gz
#GWAS sum stats:             gwas_filtered_updated_08032024.tsv.gz
#LD between QTL & GWAS vars: sceqtlgen_actually_all_vs_gwas_ld_out.tsv.gz


################################
### Load variants and filter ###
################################

sig_vars_eqtl = pd.read_csv("/scratch/hb-sceqtlgen/GRN-Project/sceqtlgen_gwas/sceqtlgen_sig_vars_no_dc.txt.gz",sep='\t')
sig_vars_eqtl.sort_values(by="p_value",inplace=True)
sig_vars_eqtl.drop_duplicates(subset="snp_id",inplace=True)

ld_all = pd.read_csv("/scratch/hb-sceqtlgen/GRN-Project/sceqtlgen_gwas/LD_output/sceqtlgen_all_vars_tested_ld_out_01LD.tsv.gz",sep='\t')
ld_all = ld_all[ld_all.snp.isin(sig_vars_eqtl.snp_id.values)] ## ld_all contains LD for all SNPs tested, here we are only interested in sig eQTL sNPs
ld_all["chr_pos"] = [f"{i.split(':')[0]}:{i.split(':')[1]}" for i in ld_all.snp]

coeqtl_df = pd.read_csv("/scratch/hb-sceqtlgen/GRN-Project/coeqtl_results/updated_results_feb_2025/5DS_Meta_Analysis_Sign_coeQTLs133p_val_mt_eGene_.csv")

sig_vars_coeqtl = coeqtl_df[["snp_id","p_value"]]
sig_vars_coeqtl["chr_pos"] = [f"{sig_vars_coeqtl.snp_id.iloc[i].split(':')[0]}:{sig_vars_coeqtl.snp_id.iloc[i].split(':')[1]}" for i in range(len(sig_vars_coeqtl.snp_id))]
sig_vars_coeqtl.sort_values(by="p_value",inplace=True)
sig_vars_coeqtl.drop_duplicates(subset="snp_id",inplace=True)
sig_vars_coeqtl=sig_vars_coeqtl[sig_vars_coeqtl.snp_id.isin(sig_vars_eqtl.snp_id)] # filter out none eQTL SNPs

#filter out HLA region
hla_snps = [i for i in ld_all.chr_pos if i.split(":")[0] == "6" and int(i.split(":")[1]) > 25000000 and int(i.split(":")[1]) < 35000000]
sig_vars_eqtl = sig_vars_eqtl[~sig_vars_eqtl.chr_pos.isin(hla_snps)]
sig_vars_coeqtl = sig_vars_coeqtl[~sig_vars_coeqtl.chr_pos.isin(hla_snps)]

# filter out co-eQTL SNPs to have a set of eqtl only snps
sig_vars_eqtl = sig_vars_eqtl[~sig_vars_eqtl.snp_id.isin(sig_vars_coeqtl.snp_id)]


## Filter with p-value ordering (files are already ordered by p-value)

coeqtl_snps = [i for i in sig_vars_coeqtl.chr_pos]
eqtl_snps = [i for i in sig_vars_eqtl.chr_pos]

# first find independent significant eQTL SNPs
ld_all_dict = {}
for i in tqdm(range(len(ld_all.chr_pos))):
  ld_all_dict[ld_all.chr_pos.iloc[i]] = ld_all.ld_vars.iloc[i]


seen_snps = set() ## used to track snps in LD
final_snps = [] ## list of snps that are not in LD with each other
for snp in tqdm(coeqtl_snps):
  snp_chr_pos = ":".join(snp.split(":")[0:2])
  if snp_chr_pos not in seen_snps:
    final_snps.append(snp)
    chrom = snp.split(":")[0]
    temp=[seen_snps.add(f"{chrom}:{pos}") for pos in ld_all_dict[snp].split(";")] ## adds LD snps to set, assigned to a temp variable to avoid statements in console


sig_vars_coeqtl = sig_vars_coeqtl[sig_vars_coeqtl.chr_pos.isin(final_snps)]

seen_snps = set()
final_snps = []
for snp in tqdm(eqtl_snps):
  snp_chr_pos = ":".join(snp.split(":")[0:2])
  if snp_chr_pos not in seen_snps:
    final_snps.append(snp)
    chrom = snp.split(":")[0]
    temp=[seen_snps.add(f"{chrom}:{pos}") for pos in ld_all_dict[snp].split(";")]


sig_vars_eqtl = sig_vars_eqtl[sig_vars_eqtl.chr_pos.isin(final_snps)]

##################################
### GWAS enrichment all traits ###
##################################

results_dict = {}
results_dict["celltypes_combined"] = {}
results_dict["celltypes_combined"]["N_coeqtl_SNPs"] = sig_vars_coeqtl.shape[0]
results_dict["celltypes_combined"]["N_eqtl_SNPs"] = sig_vars_eqtl.shape[0]


gwas = pd.read_csv("/scratch/hb-sceqtlgen/GRN-Project/danLD/gwas_filtered_updated_08032024.tsv.gz",sep='\t')

ld_all_gwas = pd.read_csv("/scratch/hb-sceqtlgen/GRN-Project/sceqtlgen_gwas/LD_output/sceqtlgen_actually_all_vs_gwas_ld_out.tsv.gz",sep='\t')
ld_all_gwas['chr_pos'] = ld_all_gwas.snp

ld_all_gwas["chromosome"]=[int(i.split(":")[0]) for i in ld_all_gwas.snp]
gwas_snps=set([i.replace("_",":") for i in gwas.chr_pos])


sig_vars_eqtl["is_gwas"] = ["yes" if i in gwas_snps else "no" for i in sig_vars_eqtl.chr_pos.values]
ld_vars = ld_all_gwas[ld_all_gwas.chr_pos.isin(sig_vars_eqtl.chr_pos)]
ld_vars.sort_values(by="chr_pos",ignore_index=True, inplace=True)
ld_vars = [i for i in ld_vars.ld_vars]
sig_vars_eqtl.sort_values(by="chr_pos", ignore_index=True, inplace=True)
sig_vars_eqtl["is_gwas_ld"] = ["yes" if ";" in ld_vars[i] else "no" for i in range(len(sig_vars_eqtl.chr_pos))]
sig_vars_eqtl["gwas"] = ["yes" if (sig_vars_eqtl.is_gwas_ld.iloc[i]=="yes" or sig_vars_eqtl.is_gwas.iloc[i]=="yes") else "no" for i in range(len(sig_vars_eqtl.is_gwas_ld))]

sig_vars_coeqtl["is_gwas"] = ["yes" if i in gwas_snps else "no" for i in sig_vars_coeqtl.chr_pos.values]
ld_vars = ld_all_gwas[ld_all_gwas.chr_pos.isin(sig_vars_coeqtl.chr_pos)]
ld_vars.sort_values(by="chr_pos",ignore_index=True, inplace=True)
ld_vars = [i for i in ld_vars.ld_vars]
sig_vars_coeqtl.sort_values(by="chr_pos", ignore_index=True, inplace=True)
sig_vars_coeqtl["is_gwas_ld"] = ["yes" if ";" in ld_vars[i] else "no" for i in range(len(sig_vars_coeqtl.chr_pos))]
sig_vars_coeqtl["gwas"] = ["yes" if (sig_vars_coeqtl.is_gwas_ld.iloc[i]=="yes" or sig_vars_coeqtl.is_gwas.iloc[i]=="yes") else "no" for i in range(len(sig_vars_coeqtl.is_gwas_ld))]


data = [[sum([i for i in sig_vars_coeqtl.gwas=="yes"]),len(sig_vars_coeqtl.snp_id)-sum([i for i in sig_vars_coeqtl.gwas=="yes"])],[sum([i for i in sig_vars_eqtl.gwas=="yes"]), len(sig_vars_eqtl.snp_id)-sum([i for i in sig_vars_eqtl.gwas=="yes"])]]
odd_ratio, p_value = stats.fisher_exact(data) 
print('odd ratio is : ' + str(odd_ratio)) 
print('p_value is : ' + str(p_value)) 
print(f'confidence interval: {stats.contingency.odds_ratio(data).confidence_interval()}')


results_dict["celltypes_combined"]["Fishers_Table_all"] = data
results_dict["celltypes_combined"]["Odds_ratio_all"] = odd_ratio
results_dict["celltypes_combined"]["P_value_all"] = p_value
results_dict["celltypes_combined"]["CI_all"] = stats.contingency.odds_ratio(data).confidence_interval()



## GWAS enrichment without cell counts

## removing cell count GWAS traits
# Hematological measurement: EFO:0004503
# basophil count: EFO:0005090
# basophil percentage of leukocytes: EFO:0007992
# basophil percentage of granulocytes: EFO:0007995

og_gwas_snps = set([i.replace("_",":") for i in gwas.chr_pos]) # og means before removing traits

efo_terms_to_drop = ["EFO_0004503","EFO_0005090","EFO_0007992","EFO_0007995"]
gwas_no_cell_count = gwas[~gwas.EFO_term.isin(efo_terms_to_drop)]
gwas_no_cell_count_snps = set([i.replace("_",":") for i in gwas_no_cell_count.chr_pos])
## SNPs need to be removed from LD files as we are not interested in these GWAS SNPs anymore
gwas_snps_to_remove = pd.DataFrame()
gwas_snps_to_remove["snp"]=[i for i in og_gwas_snps if i not in gwas_no_cell_count_snps]
gwas_snps_to_remove["chromosome"]=[i.split(":")[0] for i in gwas_snps_to_remove.snp]
gwas_snps_to_remove["pos"]=[i.split(":")[1] for i in gwas_snps_to_remove.snp]
## Removing cell count GWAS SNPs from LD file
new_ld_gwas = pd.DataFrame()
new_ld_snp = []
new_ld_ld_vars = []
for chrom in range(1,23):
  print(chrom)
  subset_gwas_snps_to_remove = gwas_snps_to_remove[gwas_snps_to_remove.chromosome==str(chrom)]
  subset_ld_all_gwas = ld_all_gwas[ld_all_gwas.chromosome==int(chrom)]
  pos_set = set(subset_gwas_snps_to_remove.pos.values)
  for i in tqdm(range(len(subset_ld_all_gwas.snp))):
    new_ld_snp.append(subset_ld_all_gwas.snp.iloc[i])
    new_ld_ld_vars.append(";".join([i for i in subset_ld_all_gwas.ld_vars.iloc[i].split(";") if i not in pos_set]))


new_ld_gwas["snp"]=new_ld_snp
new_ld_gwas["ld_vars"]=new_ld_ld_vars
ld_all_gwas_no_cell_count = new_ld_gwas
ld_all_gwas_no_cell_count['chr_pos'] = ld_all_gwas_no_cell_count.snp

gwas_no_cell_count_snps = set([i.replace("_",":") for i in gwas_no_cell_count.chr_pos])


sig_vars_eqtl["is_gwas_no_cell_count"] = ["yes" if i in gwas_no_cell_count_snps else "no" for i in sig_vars_eqtl.chr_pos.values]
ld_vars = ld_all_gwas_no_cell_count[ld_all_gwas_no_cell_count.chr_pos.isin(sig_vars_eqtl.chr_pos)]
ld_vars.sort_values(by="chr_pos",ignore_index=True, inplace=True)
ld_vars = [i for i in ld_vars.ld_vars]
sig_vars_eqtl.sort_values(by="chr_pos", ignore_index=True, inplace=True)
sig_vars_eqtl["is_gwas_ld_no_cell_count"] = ["yes" if ";" in ld_vars[i] else "no" for i in range(len(sig_vars_eqtl.chr_pos))]
sig_vars_eqtl["gwas_no_cell_count"] = ["yes" if (sig_vars_eqtl.is_gwas_ld_no_cell_count.iloc[i]=="yes" or sig_vars_eqtl.is_gwas_no_cell_count.iloc[i]=="yes") else "no" for i in range(len(sig_vars_eqtl.is_gwas_ld))]

sig_vars_coeqtl["is_gwas_no_cell_count"] = ["yes" if i in gwas_no_cell_count_snps else "no" for i in sig_vars_coeqtl.chr_pos.values]
ld_vars = ld_all_gwas_no_cell_count[ld_all_gwas_no_cell_count.chr_pos.isin(sig_vars_coeqtl.chr_pos)]
ld_vars.sort_values(by="chr_pos",ignore_index=True, inplace=True)
ld_vars = [i for i in ld_vars.ld_vars]
sig_vars_coeqtl.sort_values(by="chr_pos", ignore_index=True, inplace=True)
sig_vars_coeqtl["is_gwas_ld_no_cell_count"] = ["yes" if ";" in ld_vars[i] else "no" for i in range(len(sig_vars_coeqtl.chr_pos))]
sig_vars_coeqtl["gwas_no_cell_count"] = ["yes" if (sig_vars_coeqtl.is_gwas_ld_no_cell_count.iloc[i]=="yes" or sig_vars_coeqtl.is_gwas_no_cell_count.iloc[i]=="yes") else "no" for i in range(len(sig_vars_coeqtl.is_gwas_ld))]


data = [[sum([i for i in sig_vars_coeqtl.gwas_no_cell_count=="yes"]),len(sig_vars_coeqtl.snp_id)-sum([i for i in sig_vars_coeqtl.gwas_no_cell_count=="yes"])],[sum([i for i in sig_vars_eqtl.gwas_no_cell_count=="yes"]), len(sig_vars_eqtl.snp_id)-sum([i for i in sig_vars_eqtl.gwas_no_cell_count=="yes"])]]
odd_ratio, p_value = stats.fisher_exact(data) 
print('odd ratio is : ' + str(odd_ratio)) 
print('p_value is : ' + str(p_value)) 
print(f'confidence interval: {stats.contingency.odds_ratio(data).confidence_interval()}')


results_dict["celltypes_combined"]["Fishers_Table_all_no_cell_count"] = data
results_dict["celltypes_combined"]["Odds_ratio_all_no_cell_count"] = odd_ratio
results_dict["celltypes_combined"]["P_value_all_no_cell_count"] = p_value
results_dict["celltypes_combined"]["CI_all_no_cell_count"] = stats.contingency.odds_ratio(data).confidence_interval()




## disease specific

dis = gwas[gwas.disease==1]
dis_snps = set([i.replace("_",":") for i in dis.chr_pos])

ld_all_gwas["chromosome"] = [int(i.split(":")[0]) for i in ld_all_gwas.snp]
ld_coeqtl_snps = ld_all_gwas[ld_all_gwas.snp.isin(sig_vars_coeqtl.chr_pos)]
ld_eqtl_snps = ld_all_gwas[ld_all_gwas.snp.isin(sig_vars_eqtl.chr_pos)]

sig_dis_snps = []
for i_vars in range(len(ld_coeqtl_snps.ld_vars)):
  temp_var = ld_coeqtl_snps.snp.iloc[i_vars]
  chrom = ld_coeqtl_snps.chromosome.iloc[i_vars]
  if any([f"{chrom}:{i}" in dis_snps for i in ld_coeqtl_snps.ld_vars.iloc[i_vars].split(";")]):
    sig_dis_snps.append(f"yes")
  else:
    sig_dis_snps.append(f"no")

ld_coeqtl_snps["disease"] = sig_dis_snps


other_dis_snps = []
for i_vars in range(len(ld_eqtl_snps.ld_vars)):
  temp_var = ld_eqtl_snps.snp.iloc[i_vars]
  chrom = ld_eqtl_snps.chromosome.iloc[i_vars]
  if any([f"{chrom}:{i}" in dis_snps for i in ld_eqtl_snps.ld_vars.iloc[i_vars].split(";")]):
    other_dis_snps.append(f"yes")
  else:
    other_dis_snps.append(f"no")

ld_eqtl_snps["disease"] = other_dis_snps


sig_vars_coeqtl["disease"] = ["".join(ld_coeqtl_snps[ld_coeqtl_snps.chr_pos==i].disease) for i in sig_vars_coeqtl.chr_pos]
sig_vars_eqtl["disease"] = ["".join(ld_eqtl_snps[ld_eqtl_snps.chr_pos==i].disease) for i in sig_vars_eqtl.chr_pos]

data = [[sum([i=="yes" for i in sig_vars_coeqtl.disease]), len([i=="yes" for i in sig_vars_coeqtl.disease])-sum([i=="yes" for i in sig_vars_coeqtl.disease])],[sum([i=="yes" for i in sig_vars_eqtl.disease]), len([i=="yes" for i in sig_vars_eqtl.disease])-sum([i=="yes" for i in sig_vars_eqtl.disease])]]
odd_ratio, p_value = stats.fisher_exact(data) 
print('odd ratio is : ' + str(odd_ratio)) 
print('p_value is : ' + str(p_value)) 
print(f'confidence interval: {stats.contingency.odds_ratio(data).confidence_interval()}')

results_dict["celltypes_combined"]["Fishers_Table_disease"] = data
results_dict["celltypes_combined"]["Odds_ratio_disease"] = odd_ratio
results_dict["celltypes_combined"]["P_value_disease"] = p_value
results_dict["celltypes_combined"]["CI_disease"] = stats.contingency.odds_ratio(data).confidence_interval()


## immune specific

dis = gwas[gwas.immune_related==1]
dis_snps = set([i.replace("_",":") for i in dis.chr_pos])

ld_all_gwas["chromosome"] = [int(i.split(":")[0]) for i in ld_all_gwas.snp]
ld_coeqtl_snps = ld_all_gwas[ld_all_gwas.snp.isin(sig_vars_coeqtl.chr_pos)]
ld_eqtl_snps = ld_all_gwas[ld_all_gwas.snp.isin(sig_vars_eqtl.chr_pos)]

sig_dis_snps = []
for i_vars in range(len(ld_coeqtl_snps.ld_vars)):
  temp_var = ld_coeqtl_snps.snp.iloc[i_vars]
  chrom = ld_coeqtl_snps.chromosome.iloc[i_vars]
  if any([f"{chrom}:{i}" in dis_snps for i in ld_coeqtl_snps.ld_vars.iloc[i_vars].split(";")]):
    sig_dis_snps.append(f"yes")
  else:
    sig_dis_snps.append(f"no")

ld_coeqtl_snps["disease"] = sig_dis_snps


other_dis_snps = []
for i_vars in range(len(ld_eqtl_snps.ld_vars)):
  temp_var = ld_eqtl_snps.snp.iloc[i_vars]
  chrom = ld_eqtl_snps.chromosome.iloc[i_vars]
  if any([f"{chrom}:{i}" in dis_snps for i in ld_eqtl_snps.ld_vars.iloc[i_vars].split(";")]):
    other_dis_snps.append(f"yes")
  else:
    other_dis_snps.append(f"no")

ld_eqtl_snps["disease"] = other_dis_snps


sig_vars_coeqtl["disease"] = ["".join(ld_coeqtl_snps[ld_coeqtl_snps.chr_pos==i].disease) for i in sig_vars_coeqtl.chr_pos]
sig_vars_eqtl["disease"] = ["".join(ld_eqtl_snps[ld_eqtl_snps.chr_pos==i].disease) for i in sig_vars_eqtl.chr_pos]

data = [[sum([i=="yes" for i in sig_vars_coeqtl.disease]), len([i=="yes" for i in sig_vars_coeqtl.disease])-sum([i=="yes" for i in sig_vars_coeqtl.disease])],[sum([i=="yes" for i in sig_vars_eqtl.disease]), len([i=="yes" for i in sig_vars_eqtl.disease])-sum([i=="yes" for i in sig_vars_eqtl.disease])]]
odd_ratio, p_value = stats.fisher_exact(data) 
print('odd ratio is : ' + str(odd_ratio)) 
print('p_value is : ' + str(p_value)) 
print(f'confidence interval: {stats.contingency.odds_ratio(data).confidence_interval()}')


results_dict["celltypes_combined"]["Fishers_Table_immune"] = data
results_dict["celltypes_combined"]["Odds_ratio_immune"] = odd_ratio
results_dict["celltypes_combined"]["P_value_immune"] = p_value
results_dict["celltypes_combined"]["CI_immune"] = stats.contingency.odds_ratio(data).confidence_interval()

from pprint import pp
pp(results_dict)


new_df = pd.DataFrame()
new_df["ALL"] = [results_dict['celltypes_combined']['P_value_all'], results_dict['celltypes_combined']['Odds_ratio_all'], results_dict['celltypes_combined']['CI_all'][1], results_dict['celltypes_combined']['CI_all'][0]]
new_df.index = ["p_value","estimate","conf_upper","conf_lower"]

new_df["ALL_no_cell_count"]=[results_dict['celltypes_combined']['P_value_all_no_cell_count'], results_dict['celltypes_combined']['Odds_ratio_all_no_cell_count'], results_dict['celltypes_combined']['CI_all_no_cell_count'][1], results_dict['celltypes_combined']['CI_all_no_cell_count'][0]]

new_df["ALL_disease"] = [results_dict['celltypes_combined']['P_value_disease'], results_dict['celltypes_combined']['Odds_ratio_disease'], results_dict['celltypes_combined']['CI_disease'][1], results_dict['celltypes_combined']['CI_disease'][0]]
new_df["ALL_immune"] = [results_dict['celltypes_combined']['P_value_immune'], results_dict['celltypes_combined']['Odds_ratio_immune'], results_dict['celltypes_combined']['CI_immune'][1], results_dict['celltypes_combined']['CI_immune'][0]]







########################################
#### Separate cell type enrichments ####
########################################

celltypes = ["B","CD4_T","CD8_T","Mono","NK"]

path = "/scratch/hb-sceqtlgen/GRN-Project/sceqtlgen_gwas/eQTLs_finemapped_20250212/"
files = os.listdir(path)
coeqtl_df = pd.read_csv("/scratch/hb-sceqtlgen/GRN-Project/coeqtl_results/updated_results_feb_2025/5DS_Meta_Analysis_Sign_coeQTLs133p_val_mt_eGene_.csv", index_col=0)


results_dict = {}

for ct in celltypes:
  results_dict[ct] = {}
  file = "".join([i for i in files if f"{ct}.Ds" in i])
  sig_vars_eqtl =pd.read_csv(f"{path}{file}",sep='\t')
  sig_vars_eqtl = sig_vars_eqtl[["snp_id","p_value"]]
  sig_vars_eqtl.sort_values(by="p_value",inplace=True)
  sig_vars_eqtl.drop_duplicates(subset="snp_id",inplace=True,ignore_index=True)
  sig_vars_eqtl["chr_pos"] = [f"{i.split(':')[0]}:{i.split(':')[1]}" for i in sig_vars_eqtl.snp_id]
  ld_all = pd.read_csv("/scratch/hb-sceqtlgen/GRN-Project/sceqtlgen_gwas/LD_output/sceqtlgen_all_vars_tested_ld_out_01LD.tsv.gz",sep='\t')
  ld_all = ld_all[ld_all.snp.isin(sig_vars_eqtl.snp_id.values)]
  ld_all["chr_pos"] = [f"{i.split(':')[0]}:{i.split(':')[1]}" for i in ld_all.snp]
  coeqtl_df_ct = coeqtl_df[coeqtl_df.cell_type==ct]
  sig_vars_coeqtl = coeqtl_df[["snp_id","p_value"]]
  sig_vars_coeqtl["chr_pos"] = [f"{sig_vars_coeqtl.snp_id.iloc[i].split(':')[0]}:{sig_vars_coeqtl.snp_id.iloc[i].split(':')[1]}" for i in range(len(sig_vars_coeqtl.snp_id))]
  sig_vars_coeqtl.sort_values(by="p_value",inplace=True)
  sig_vars_coeqtl.drop_duplicates(subset="snp_id",inplace=True,ignore_index=True)
  sig_vars_coeqtl=sig_vars_coeqtl[sig_vars_coeqtl.snp_id.isin(sig_vars_eqtl.snp_id)] # filter out none eQTL SNPs
  #filter out HLA region
  hla_snps = [i for i in ld_all.chr_pos if i.split(":")[0] == "6" and int(i.split(":")[1]) > 25000000 and int(i.split(":")[1]) < 35000000]
  sig_vars_eqtl = sig_vars_eqtl[~sig_vars_eqtl.chr_pos.isin(hla_snps)]
  sig_vars_coeqtl = sig_vars_coeqtl[~sig_vars_coeqtl.chr_pos.isin(hla_snps)]
  # filter out co-eQTL SNPs
  sig_vars_eqtl = sig_vars_eqtl[~sig_vars_eqtl.snp_id.isin(sig_vars_coeqtl.snp_id)]
  # filter based on LD
  coeqtl_snps = [i for i in sig_vars_coeqtl.chr_pos]
  eqtl_snps = [i for i in sig_vars_eqtl.chr_pos]
  # first find independent significant eQTL SNPs
  ld_all_dict = {}
  for i in tqdm(range(len(ld_all.chr_pos))):
    ld_all_dict[ld_all.chr_pos.iloc[i]] = ld_all.ld_vars.iloc[i]
  seen_snps = set()
  final_snps = []
  for snp in tqdm(coeqtl_snps):
    snp_chr_pos = ":".join(snp.split(":")[0:2])
    if snp_chr_pos not in seen_snps:
      final_snps.append(snp)
      chrom = snp.split(":")[0]
      temp = [seen_snps.add(f"{chrom}:{pos}") for pos in ld_all_dict[snp].split(";")]
  sig_vars_coeqtl = sig_vars_coeqtl[sig_vars_coeqtl.chr_pos.isin(final_snps)]
  seen_snps = set()
  final_snps = []
  for snp in tqdm(eqtl_snps):
    snp_chr_pos = ":".join(snp.split(":")[0:2])
    if snp_chr_pos not in seen_snps:
      final_snps.append(snp)
      chrom = snp.split(":")[0]
      temp=[seen_snps.add(f"{chrom}:{pos}") for pos in ld_all_dict[snp].split(";")]
  sig_vars_eqtl = sig_vars_eqtl[sig_vars_eqtl.chr_pos.isin(final_snps)]
  #
  #
  ####### GWAS enrichment
  gwas_snps = set([i.replace("_",":") for i in gwas.chr_pos])
  sig_vars_eqtl["is_gwas"] = ["yes" if i in gwas_snps else "no" for i in sig_vars_eqtl.chr_pos.values]
  ld_vars = ld_all_gwas[ld_all_gwas.chr_pos.isin(sig_vars_eqtl.chr_pos)]
  ld_vars.sort_values(by="chr_pos",ignore_index=True, inplace=True)
  ld_vars = [i for i in ld_vars.ld_vars]
  sig_vars_eqtl.sort_values(by="chr_pos", ignore_index=True, inplace=True)
  sig_vars_eqtl["is_gwas_ld"] = ["yes" if ";" in ld_vars[i] else "no" for i in range(len(sig_vars_eqtl.chr_pos))]
  sig_vars_eqtl["gwas"] = ["yes" if (sig_vars_eqtl.is_gwas_ld.iloc[i]=="yes" or sig_vars_eqtl.is_gwas.iloc[i]=="yes") else "no" for i in range(len(sig_vars_eqtl.is_gwas_ld))]
  sig_vars_coeqtl["is_gwas"] = ["yes" if i in gwas_snps else "no" for i in sig_vars_coeqtl.chr_pos.values]
  ld_vars = ld_all_gwas[ld_all_gwas.chr_pos.isin(sig_vars_coeqtl.chr_pos)]
  ld_vars.sort_values(by="chr_pos",ignore_index=True, inplace=True)
  ld_vars = [i for i in ld_vars.ld_vars]
  sig_vars_coeqtl.sort_values(by="chr_pos", ignore_index=True, inplace=True)
  sig_vars_coeqtl["is_gwas_ld"] = ["yes" if ";" in ld_vars[i] else "no" for i in range(len(sig_vars_coeqtl.chr_pos))]
  sig_vars_coeqtl["gwas"] = ["yes" if (sig_vars_coeqtl.is_gwas_ld.iloc[i]=="yes" or sig_vars_coeqtl.is_gwas.iloc[i]=="yes") else "no" for i in range(len(sig_vars_coeqtl.is_gwas_ld))]
  data = [[sum([i for i in sig_vars_coeqtl.gwas=="yes"]),len(sig_vars_coeqtl.snp_id)-sum([i for i in sig_vars_coeqtl.gwas=="yes"])],[sum([i for i in sig_vars_eqtl.gwas=="yes"]), len(sig_vars_eqtl.snp_id)-sum([i for i in sig_vars_eqtl.gwas=="yes"])]]
  results_dict[ct]["Fishers_Table_all"] = data
  odd_ratio, p_value = stats.fisher_exact(data) 
  results_dict[ct]["Odds_Ratio_all"] = odd_ratio 
  results_dict[ct]["P_value_all"] = p_value 
  results_dict[ct]["confidence_interval_all"] = stats.contingency.odds_ratio(data).confidence_interval()
  print(f"\n{ct} traits done.\n")
  #
  #
  ####### GWAS enrichment no cell count traits
  sig_vars_eqtl["is_gwas"] = ["yes" if i in gwas_no_cell_count_snps else "no" for i in sig_vars_eqtl.chr_pos.values]
  ld_vars = ld_all_gwas_no_cell_count[ld_all_gwas_no_cell_count.chr_pos.isin(sig_vars_eqtl.chr_pos)]
  ld_vars.sort_values(by="chr_pos",ignore_index=True, inplace=True)
  ld_vars = [i for i in ld_vars.ld_vars]
  sig_vars_eqtl.sort_values(by="chr_pos", ignore_index=True, inplace=True)
  sig_vars_eqtl["is_gwas_ld"] = ["yes" if ";" in ld_vars[i] else "no" for i in range(len(sig_vars_eqtl.chr_pos))]
  sig_vars_eqtl["gwas"] = ["yes" if (sig_vars_eqtl.is_gwas_ld.iloc[i]=="yes" or sig_vars_eqtl.is_gwas.iloc[i]=="yes") else "no" for i in range(len(sig_vars_eqtl.is_gwas_ld))]
  sig_vars_coeqtl["is_gwas"] = ["yes" if i in gwas_snps else "no" for i in sig_vars_coeqtl.chr_pos.values]
  ld_vars = ld_all_gwas[ld_all_gwas.chr_pos.isin(sig_vars_coeqtl.chr_pos)]
  ld_vars.sort_values(by="chr_pos",ignore_index=True, inplace=True)
  ld_vars = [i for i in ld_vars.ld_vars]
  sig_vars_coeqtl.sort_values(by="chr_pos", ignore_index=True, inplace=True)
  sig_vars_coeqtl["is_gwas_ld"] = ["yes" if ";" in ld_vars[i] else "no" for i in range(len(sig_vars_coeqtl.chr_pos))]
  sig_vars_coeqtl["gwas"] = ["yes" if (sig_vars_coeqtl.is_gwas_ld.iloc[i]=="yes" or sig_vars_coeqtl.is_gwas.iloc[i]=="yes") else "no" for i in range(len(sig_vars_coeqtl.is_gwas_ld))]
  data = [[sum([i for i in sig_vars_coeqtl.gwas=="yes"]),len(sig_vars_coeqtl.snp_id)-sum([i for i in sig_vars_coeqtl.gwas=="yes"])],[sum([i for i in sig_vars_eqtl.gwas=="yes"]), len(sig_vars_eqtl.snp_id)-sum([i for i in sig_vars_eqtl.gwas=="yes"])]]
  results_dict[ct]["Fishers_Table_all_no_cell_count"] = data
  odd_ratio, p_value = stats.fisher_exact(data) 
  results_dict[ct]["Odds_Ratio_all_no_cell_count"] = odd_ratio 
  results_dict[ct]["P_value_all_no_cell_count"] = p_value 
  results_dict[ct]["confidence_interval_all_no_cell_count"] = stats.contingency.odds_ratio(data).confidence_interval()
  print(f"\n{ct} traits done.\n")  
  #
  #
  ## disease specific
  dis = gwas[gwas.disease==1]
  dis_snps = set([i.replace("_",":") for i in dis.chr_pos])
  ld_all_gwas["chromosome"] = [int(i.split(":")[0]) for i in ld_all_gwas.snp]
  ld_coeqtl_snps = ld_all_gwas[ld_all_gwas.snp.isin(sig_vars_coeqtl.chr_pos)]
  ld_eqtl_snps = ld_all_gwas[ld_all_gwas.snp.isin(sig_vars_eqtl.chr_pos)]
  sig_dis_snps = []
  for i_vars in range(len(ld_coeqtl_snps.ld_vars)):
    temp_var = ld_coeqtl_snps.snp.iloc[i_vars]
    chrom = ld_coeqtl_snps.chromosome.iloc[i_vars]
    if any([f"{chrom}:{i}" in dis_snps for i in ld_coeqtl_snps.ld_vars.iloc[i_vars].split(";")]):
      sig_dis_snps.append(f"yes")
    else:
      sig_dis_snps.append(f"no")
  ld_coeqtl_snps["disease"] = sig_dis_snps
  other_dis_snps = []
  for i_vars in range(len(ld_eqtl_snps.ld_vars)):
    temp_var = ld_eqtl_snps.snp.iloc[i_vars]
    chrom = ld_eqtl_snps.chromosome.iloc[i_vars]
    if any([f"{chrom}:{i}" in dis_snps for i in ld_eqtl_snps.ld_vars.iloc[i_vars].split(";")]):
      other_dis_snps.append(f"yes")
    else:
      other_dis_snps.append(f"no")
  ld_eqtl_snps["disease"] = other_dis_snps
  sig_vars_coeqtl["disease"] = ["".join(ld_coeqtl_snps[ld_coeqtl_snps.chr_pos==i].disease) for i in sig_vars_coeqtl.chr_pos]
  sig_vars_eqtl["disease"] = ["".join(ld_eqtl_snps[ld_eqtl_snps.chr_pos==i].disease) for i in sig_vars_eqtl.chr_pos]
  data = [[sum([i=="yes" for i in sig_vars_coeqtl.disease]), len([i=="yes" for i in sig_vars_coeqtl.disease])-sum([i=="yes" for i in sig_vars_coeqtl.disease])],[sum([i=="yes" for i in sig_vars_eqtl.disease]), len([i=="yes" for i in sig_vars_eqtl.disease])-sum([i=="yes" for i in sig_vars_eqtl.disease])]]
  results_dict[ct]["Fishers_Table_disease"] = data
  odd_ratio, p_value = stats.fisher_exact(data) 
  results_dict[ct]["Odds_Ratio_disease"] = odd_ratio 
  results_dict[ct]["P_value_disease"] = p_value 
  results_dict[ct]["confidence_interval_disease"] = stats.contingency.odds_ratio(data).confidence_interval()
  print(f"\n{ct} disease done.\n")
  #
  #
  ## immune specific
  dis = gwas[gwas.immune_related==1]
  dis_snps = set([i.replace("_",":") for i in dis.chr_pos])
  ld_all_gwas["chromosome"] = [int(i.split(":")[0]) for i in ld_all_gwas.snp]
  ld_coeqtl_snps = ld_all_gwas[ld_all_gwas.snp.isin(sig_vars_coeqtl.chr_pos)]
  ld_eqtl_snps = ld_all_gwas[ld_all_gwas.snp.isin(sig_vars_eqtl.chr_pos)]
  sig_dis_snps = []
  for i_vars in range(len(ld_coeqtl_snps.ld_vars)):
    temp_var = ld_coeqtl_snps.snp.iloc[i_vars]
    chrom = ld_coeqtl_snps.chromosome.iloc[i_vars]
    if any([f"{chrom}:{i}" in dis_snps for i in ld_coeqtl_snps.ld_vars.iloc[i_vars].split(";")]):
      sig_dis_snps.append(f"yes")
    else:
      sig_dis_snps.append(f"no")
  ld_coeqtl_snps["disease"] = sig_dis_snps
  other_dis_snps = []
  for i_vars in range(len(ld_eqtl_snps.ld_vars)):
    temp_var = ld_eqtl_snps.snp.iloc[i_vars]
    chrom = ld_eqtl_snps.chromosome.iloc[i_vars]
    if any([f"{chrom}:{i}" in dis_snps for i in ld_eqtl_snps.ld_vars.iloc[i_vars].split(";")]):
      other_dis_snps.append(f"yes")
    else:
      other_dis_snps.append(f"no")
  ld_eqtl_snps["disease"] = other_dis_snps
  sig_vars_coeqtl["disease"] = ["".join(ld_coeqtl_snps[ld_coeqtl_snps.chr_pos==i].disease) for i in sig_vars_coeqtl.chr_pos]
  sig_vars_eqtl["disease"] = ["".join(ld_eqtl_snps[ld_eqtl_snps.chr_pos==i].disease) for i in sig_vars_eqtl.chr_pos]
  data = [[sum([i=="yes" for i in sig_vars_coeqtl.disease]), len([i=="yes" for i in sig_vars_coeqtl.disease])-sum([i=="yes" for i in sig_vars_coeqtl.disease])],[sum([i=="yes" for i in sig_vars_eqtl.disease]), len([i=="yes" for i in sig_vars_eqtl.disease])-sum([i=="yes" for i in sig_vars_eqtl.disease])]]
  results_dict[ct]["Fishers_Table_immune"] = data
  odd_ratio, p_value = stats.fisher_exact(data) 
  results_dict[ct]["Odds_Ratio_immune"] = odd_ratio 
  results_dict[ct]["P_value_immune"] = p_value 
  results_dict[ct]["confidence_interval_immune"] = stats.contingency.odds_ratio(data).confidence_interval()
  print(f"\n{ct} done.\n")


from pprint import pp
pp(results_dict)

for ct in celltypes:
  new_df[f"{ct}_all_traits"]=[results_dict[ct]['P_value_all'], results_dict[ct]['Odds_Ratio_all'], results_dict[ct]['confidence_interval_all'][1], results_dict[ct]['confidence_interval_all'][0]]
  new_df[f"{ct}_all_traits_no_cell_count"]=[results_dict[ct]['P_value_all_no_cell_count'], results_dict[ct]['Odds_Ratio_all_no_cell_count'], results_dict[ct]['confidence_interval_all_no_cell_count'][1], results_dict[ct]['confidence_interval_all_no_cell_count'][0]]
  new_df[f"{ct}_disease"]=[results_dict[ct]['P_value_disease'], results_dict[ct]['Odds_Ratio_disease'], results_dict[ct]['confidence_interval_disease'][1], results_dict[ct]['confidence_interval_disease'][0]]
  new_df[f"{ct}_immune"]=[results_dict[ct]['P_value_immune'], results_dict[ct]['Odds_Ratio_immune'], results_dict[ct]['confidence_interval_immune'][1], results_dict[ct]['confidence_interval_immune'][0]]




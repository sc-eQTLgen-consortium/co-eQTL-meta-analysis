

###
# MAKE FILE OF ALL FEATURES TO TEST
##
import pandas as pd
import numpy as np
import os
import time
from pprint import pp

now = time.time()
ct = 'Mono'
tmp = 'tmp01'
dataset = 'wg3_multiome'
ongoing = 'Depracated'
sc_dir = 'sceqtlgen_coeqtl_map'

#gwas = pd.read_csv(f'/groups/umcg-franke-scrna/{tmp}/projects/sc-eqtlgen-consortium-pipeline/{ongoing}/wg3/{dataset}/{sc_dir}/GWAS_snp_gene_pairs_immune_related_disease.txt.gz',sep='\t')

genes = pd.read_csv(f'/groups/umcg-franke-scrna/{tmp}/projects/sc-eqtlgen-consortium-pipeline/{ongoing}/wg3/{dataset}/calculate_correlation_metrics/gene_lists/F11_Decision_Tree_Genes{dataset}_{ct}.Qced.Normalized.SCs.Rds.tsv.gz',sep='\t',header=None)

df = pd.read_csv(f'/groups/umcg-franke-scrna/{tmp}/projects/sc-eqtlgen-consortium-pipeline/{ongoing}/wg3/{dataset}/{sc_dir}/eQTLs_finemapped_20240626/{ct}.DS.wg3_Ye_wg3_wijst2018_wg3_sawcer_wg3_oneK1K_wg3_okada_wg3_Li_wg3_Franke_split_v3_wg3_Franke_split_v2_wg3_multiome_UT_wg3_idaghdour.csTop_qtl_results.txt',sep='\t')


gene_list = [i for i in genes[0]]
egenes = np.unique([i for i in df.feature_id])
filtered_egenes= [i for i in egenes if i in gene_list]
df = df[df.feature_id.isin(filtered_egenes)]
df.index=[i for i in range(len(df.feature_id))]
#gwas=gwas[gwas.feature_id.isin(filtered_egenes)]

features1=[i for i in df.feature_id]
#features2=[i for i in gwas.feature_id]
#all_feats=features1+features2
all_feats=features1  #Not using GWAS variants
features1=0
features2=0
snps1=[i for i in df.snp_id]
snps2=[i for i in gwas.snp_id]
all_snps=snps1+snps2
snps1=0
snps2=0

all_egene_snps=pd.DataFrame()
all_egene_snps['feature_id']=all_feats
all_egene_snps['snp_id']=all_snps
all_feats=0
all_snps=0

unq_ids = [f"{all_egene_snps.feature_id.iloc[i]}_{all_egene_snps.snp_id.iloc[i]}" for i in range(len(all_egene_snps.snp_id))]

now=time.time()
final_features = []
for i in range(len(unq_ids)):
  if time.time() - now > 10:
    now=time.time()
    print(round(i/len(unq_ids),3),end='\r')
  temp_gene = unq_ids[i].split('_')[0]
  temp_snp = unq_ids[i].split('_')[1]
  for j in gene_list:
    if j == temp_gene:
      continue
    temp_feature = f"{temp_gene}_{j};{temp_snp}"
    final_features.append(temp_feature)


new_df=pd.DataFrame()
new_df['snp_id']=[i.split(';')[1] for i in final_features]
final_features = [i.split(';')[0] for i in final_features]
final_features = ['_'.join(sorted(i.split('_'))) for i in final_features]
new_df['feature_id']=[i.split(';')[0] for i in final_features]
output=f'all_features_to_test_{ct}.tsv'

new_df.to_csv(output,index=None,sep='\t')



###
# MAKE ANNOTATION FILES
###

df = df.drop_duplicates(subset='feature_id')
unq_ids = [f"{df.feature_id.iloc[i]}_{df.snp_id.iloc[i]}" for i in range(len(df.snp_id))]
chromosome = [i for i in df.feature_chromosome]
start = [i for i in df.feature_start]
end = [i for i in df.feature_end]
ENSG = [i for i in df.ENSG]
biotype = [i for i in df.biotype]

final_features = []
new_chr = []
new_start = []
new_end = []
new_ensg = []
new_biotype = []
for i in range(len(unq_ids)):
  if time.time() - now > 10:
    now=time.time()
    print(round(i/len(unq_ids),3),end='\r')
  temp_gene = unq_ids[i].split('_')[0]
  temp_snp = unq_ids[i].split('_')[1]
  for j in gene_list:
    if j == temp_gene:
      continue
    temp_feature = f"{temp_gene}_{j};{temp_snp}"
    final_features.append(temp_feature)
    new_chr.append(chromosome[i])
    new_start.append(start[i])
    new_end.append(end[i])
    new_ensg.append(ENSG[i])
    new_biotype.append(biotype[i])


final_features = [i.split(';')[0] for i in final_features]
final_features = ['_'.join(sorted(i.split('_'))) for i in final_features]


new_df=pd.DataFrame()
new_df['feature_id']=final_features
new_df['chromosome']=new_chr
new_df['start']=new_start
new_df['end']=new_end
new_df['ENSG']=new_ensg
new_df['biotype']=new_biotype

new_df.sort_values(by='feature_id',inplace=True,ignore_index=True)

features=[i for i in new_df.feature_id]
pGene = 0
feats_to_split = []
other_feats = []
for i in range(len(features)):
  if pGene == features[i]:
    feats_to_split.append(i)
  else:
    other_feats.append(i)
  pGene = features[i]


# feats_to_split and other_feats are lists of indexes to quickly subset the big df
annotation2 = new_df.iloc[feats_to_split]
annotation1 = new_df.iloc[other_feats]

output=f'LimixAnnotationFile_co_{ct}.txt'
annotation1.to_csv(output,index=None,sep='\t')
output=f'LimixAnnotationFile_co2_{ct}.txt'
annotation2.to_csv(output,index=None,sep='\t')

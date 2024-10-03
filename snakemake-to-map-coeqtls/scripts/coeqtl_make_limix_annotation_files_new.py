"""
This script is creating gene annotation files

authors: Dan Kaptijn, Marc-Jan Bonder, Roy Oelen
"""

#############
# libraries #
#############
import pandas as pd
import numpy as np
import os
import time
from pprint import pp
import argparse

#############
# functions #
#############



#############
# main code #
#############

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gwas_loc', type = str, help = 'location of gwas file (string)')
parser.add_argument('-l', '--gene_loc', type = str, help = 'location of gene list file (string)')
parser.add_argument('-v', '--variant_loc', type = str, help = 'location of variant list file (string)')
parser.add_argument('-q', '--qtl_loc', type = str, help = 'location of previous eQTL summary stats')
parser.add_argument('-f', '--features_out_loc', type = str, help = 'location to write the to-test features')
parser.add_argument('-o', '--co_limix_annotation_prepend', type = str, help = 'location of original feature annotations')
args = parser.parse_args()


####################
# to-test features #
####################

# read the limix annotation file (we'll need it later)
limix_orig = pd.read_csv(args.limix_annotation_loc, sep = '\t', header = 0)

# get genes
genes = pd.read_csv(args.gene_loc, sep = '\t', header = None)
# get the genes
gene_list = [i for i in genes[0]]
# set those as the genes to test
features1 = gene_list
features2 = []
# and the variants to test
snps1 = []
snps2 = []

# get the eQTL summary stats
eqtls = pd.read_csv(args.qtl_loc, sep = '\t')
# get genes with an eQTL effect
egenes = np.unique([i for i in eqtls.feature_id])
# filter genes to egenes
filtered_egenes= [i for i in egenes if i in gene_list]
# filter eQTL sum stats to genes in the data
eqtls = eqtls[eqtls.feature_id.isin(filtered_egenes)]
# set the index
eqtls.index=[i for i in range(len(eqtls.feature_id))]
# if we had eQTL sumstats, we would instead use those as the gene1
features1=[i for i in eqtls.feature_id]
# and the variants
snps1=[i for i in eqtls.snp_id]

# get variants from GWAS
gwas = None
if args.gwas_loc is not None:
    gwas = pd.read_csv(args.gwas_loc, sep = '\t')
    # subset to genes we have
    gwas=gwas[gwas.feature_id.isin(gene_list)]
    # if we had GWAS sumstats with gene IDs, we will test those as gene2
    features2=[i for i in gwas.feature_id]
    # and the variants
    snps2=[i for i in gwas.snp_id]

# or from a specific variant file
variants = None
if args.variant_loc is not None:
    snps2 = pd.read_csv(args.variant_loc, sep = '\t', header = None)
    # get the genes
    snps2 = [i for i in snps2[0]]

# get both those lists
all_feats = features1 + features2

# get all variants
all_snps = snps1 + snps2

# create the snp-gene pairs table
all_egene_snps=pd.DataFrame()
all_egene_snps['feature_id']=all_feats
all_egene_snps['snp_id']=all_snps

# get the unique combinations
unq_ids = [f"{all_egene_snps.feature_id.iloc[i]}_{all_egene_snps.snp_id.iloc[i]}" for i in range(len(all_egene_snps.snp_id))]

# I don't know exactly what is going on here, I should ask Dan
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

# create a new dataframe
new_df=pd.DataFrame()
# have the variants in there
new_df['snp_id']=[i.split(';')[1] for i in final_features]
# and the genes
final_features = [i.split(';')[0] for i in final_features]
# then join them on all the other genes
final_features = ['_'.join(sorted(i.split('_'))) for i in final_features]
# add that to the new dataframe
new_df['feature_id']=[i.split(';')[0] for i in final_features]
# write this to a file
new_df.to_csv(args.features_out_loc, index = None, sep = '\t')


####################
# annotation files #
####################

# grab the variables
limix_orig = eqtls.drop_duplicates(subset='feature_id')
# and subset to genes that we have
limix_orig=limix_orig[limix_orig.feature_id.isin(gene_list)]
# extract all info
unq_ids = [f"{limix_orig.feature_id.iloc[i]}_{limix_orig.snp_id.iloc[i]}" for i in range(len(limix_orig.snp_id))]
chromosome = [i for i in limix_orig.feature_chromosome]
start = [i for i in limix_orig.feature_start]
end = [i for i in limix_orig.feature_end]
ENSG = [i for i in limix_orig.ENSG]
biotype = [i for i in limix_orig.biotype]

# create a new annotation
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


new_limix=pd.DataFrame()
new_limix['feature_id']=final_features
new_limix['chromosome']=new_chr
new_limix['start']=new_start
new_limix['end']=new_end
new_limix['ENSG']=new_ensg
new_limix['biotype']=new_biotype

new_limix.sort_values(by='feature_id',inplace=True,ignore_index=True)

features=[i for i in new_limix.feature_id]
pGene = 0
feats_to_split = []
other_feats = []
for i in range(len(features)):
  if pGene == features[i]:
    feats_to_split.append(i)
  else:
    other_feats.append(i)
  pGene = features[i]


# feats_to_split and other_feats are lists of indexes to quickly subset the big limix_orig
annotation2 = new_limix.iloc[feats_to_split]
annotation1 = new_limix.iloc[other_feats]

output=args.co_limix_annotation_prepend + '.txt'
annotation1.to_csv(output,index=None,sep='\t')
output=args.co_limix_annotation_prepend + '2.txt'
annotation2.to_csv(output,index=None,sep='\t')

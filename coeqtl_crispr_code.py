
###########################
# co-eQTL CRISPR analysis #
###########################

## location = Downloads folder on laptop or Vaxtron


import pandas as pd
import numpy as np
import re
from tqdm import tqdm
import random

import statistics as stats
import scipy.stats as sp


## Load in data

coeqtl_df = pd.read_csv("/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/GRN-reconstruction/coeqtl_results/5DS_Meta_Analysis_Sign_coeQTLs133p_val_mt_eGene_.csv.gz",index_col=0)


#get headers of the columns of the file

df_crispr = pd.read_csv('/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/GRN-reconstruction/crispr_validation/cripsr_matrix_replogle.csv.gz')

df_crispr = df_crispr[[i for i in df_crispr.columns if i != "Unnamed: 1"]]

# drop first two rows
df_crispr = df_crispr.iloc[2:]
df_crispr = df_crispr.rename(columns={'gene_transcript': 'gene_name'})

col_names = []
for col in df_crispr.columns[1:]:
    col = col.split('_')
    col_names.append(col[1])


col_names.insert(0, 'gene_name')
df_copy = df_crispr.copy()
df_copy = df_crispr.set_axis(col_names, axis=1)


### Creating stacked df from looking up crispr specific values
df_crispr = df_copy.copy()
df_crispr.index = [i for i in df_crispr.gene_name]
df_crispr = df_crispr[df_crispr.columns[1:]]
df_crispr = df_crispr.stack().reset_index()
df_crispr.columns = ["eGene_measured","coeGene_target","crispr_score"]
df_crispr["abs_crispr_score"] = [abs(i) for i in df_crispr.crispr_score]
df_crispr["feature_id_e_co"] = [f"{df_crispr.eGene_measured.iloc[i]}_{df_crispr.coeGene_target.iloc[i]}" for i in tqdm(range(len(df_crispr.eGene_measured)))]

### Filtering crispr data to co-eQTL genes
coeqtl_df["feature_id_e_co"] = [f"{coeqtl_df.eGene.iloc[i]}_{coeqtl_df.coeGene.iloc[i]}" for i in range(len(coeqtl_df.eGene))]
df_crispr_coeqtl = df_crispr[df_crispr.feature_id_e_co.isin(coeqtl_df.feature_id_e_co)].copy()
df_crispr_coeqtl.drop_duplicates(subset="feature_id_e_co",inplace=True,ignore_index=True)
coegene_list = df_crispr_coeqtl.coeGene_target.unique()
df_crispr = df_crispr[df_crispr.coeGene_target.isin(coegene_list)]
coeqtl_set = set(coeqtl_df.feature_id_e_co)

is_coeqtl = []
for i in df_crispr.feature_id_e_co:
  if i in coeqtl_set:
    is_coeqtl.append(1)
  else:
    is_coeqtl.append(0)

df_crispr["coeQTL"]=is_coeqtl


coeqtl_egenes = set(coeqtl_df.eGene.unique())

df_crispr = df_crispr[df_crispr.eGene_measured.isin(coeqtl_egenes)]



#Determine z-scores of coeqtls (sig) and non-coeqtls (nsig) test for significant difference
sig = list(df_crispr[df_crispr['coeQTL'] == 1 ]['crispr_score'])
abs_sig = [abs(i) for i in sig]
nSig =list(df_crispr[df_crispr['coeQTL'] == 0]['crispr_score'])
abs_nSig = [abs(i) for i in nSig]

print(np.mean((sig)))
print(np.mean((nSig)))
print(f"size difference: {len(sig)} vs {len(nSig)}")
sp.ttest_ind((sig),(nSig))

print(np.mean((abs_sig)))
print(np.mean((abs_nSig)))
print(f"size difference: {len(abs_sig)} vs {len(abs_nSig)}")
sp.ttest_ind((abs_sig),(abs_nSig))


df_crispr.columns = ['eGene', 'coeGene', 'crispr_score', 'abs_crispr_score', 'feature_id_e_co', 'coeQTL']
def crispr_coegene_sig_nsig_unique(df):
    #coegenes = list(df['coeGene'].unique())
    """
    Return unique crispr_scores based on significance in lists
    Return mean crispr_scores per coegene in lists
    df = dataframe of coeqtls with crispr scores with annotated significance
    Columns: coeQTL, coeGene, eGene, crispr_score
    """
    sig_coegenes = list(df[df['coeQTL'] == 1]['coeGene'].unique())
    mean_sig_list = []
    mean_nsig_list = []
    unique_sig_list = []
    unique_nsig_list = []
    coegene_list = []
    for coegene in tqdm(sig_coegenes):
        coegene_list.append(coegene)
        df_one_coegene = df.query('coeGene == @coegene')
        #Significant coegenes    
        sig = [abs(i) for i in list(df_one_coegene.query('coeQTL == 1')['crispr_score'].unique())]
        mean_sig_list.append(np.mean(sig))
        for x in sig:
            unique_sig_list.append(abs(x))
        #Non signficiant coegenes
        #nsig = [abs(i) for i in list(df_one_coegene.query('coeQTL == 0')['crispr_score'].unique())]# here n is not the same for sig an not_sig
        nsig = random.sample([abs(i) for i in list(df_one_coegene.query('coeQTL == 0')['crispr_score'].unique())],len(sig)) # here n is the same for sig an not_sig
        mean_nsig_list.append(np.mean(nsig))
        for y in nsig:
            unique_nsig_list.append(abs(y))
    return unique_sig_list, unique_nsig_list, mean_sig_list, mean_nsig_list, coegene_list



unique_sig_list, unique_nsig_list, mean_sig_list, mean_nsig_list, coegene_list = crispr_coegene_sig_nsig_unique(df_crispr)


list_of_pvalues = []
for cycle in tqdm(range(100)):
  unique_sig_list, unique_nsig_list, mean_sig_list, mean_nsig_list, coegene_list = crispr_coegene_sig_nsig_unique(df_crispr)
  list_of_pvalues.append(sp.wilcoxon((mean_sig_list), (mean_nsig_list))[1])
  

print(np.mean(mean_sig_list))
print(np.mean(mean_nsig_list))
#Assuming the scores are not independent due to the experiment

##T test
sp.ttest_rel((mean_sig_list), (mean_nsig_list))

#Wilcoxon signed rank test
sp.wilcoxon((mean_sig_list), (mean_nsig_list))

## boxplot of CRISPR results
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns

df = pd.DataFrame({"Not_coeqtl": mean_nsig_list, "Sig_coeqtl":mean_sig_list})

my_palette = {"Sig_coeqtl": "darkorange", "Not_coeqtl": "dodgerblue"}
ax = sns.boxplot(data=df, palette=my_palette)
ax.set(ylabel="absolute CRISPR z-score")
ax.set_title("N=625, P=8.54x10-7")
fig = ax.get_figure()
fig.savefig("replogle_crispr_boxplot_egene_coeqtl.pdf")
plt.cla()



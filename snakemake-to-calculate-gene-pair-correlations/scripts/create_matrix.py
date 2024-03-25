"""
Creates a gene-gene matrix from the aggregated correlation metric files
"""
import numpy as np
import pandas as pd 
import argparse
import gzip
import scipy.stats 

parser = argparse.ArgumentParser(description="")
parser.add_argument("--meta_analysis", required=True, type=str, help="Meta-analysis method used to combine correlation over donors")
parser.add_argument("--donor_list", required=True, type=str, help="Donor list path")
parser.add_argument("--gene_list", required=True, type=str, help="Gene list path")
parser.add_argument("--input", required=True, nargs="+", type=str, help="Input correlation file")
parser.add_argument("--output", required=True, nargs="+", type=str, help="Aggregated output file")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def fishers_z_transformation(corr,N):
    
    corr = np.array(corr, dtype=float)
    N = np.array(N, dtype=float)
    corr,N = remove_nan(corr,N)
    
    return get_fishers_z_metrics(corr,N)
    
def get_fishers_z_metrics(corr,N):
    
    Y = calc_effect_size(corr)
    V = calc_variance_within(N)
    W = calc_weight(V)

    sum_W = np.sum(W)
    sum_WY = np.sum(W * Y)
    sum_WY2 = np.sum(W * Y * Y)
    sum_W2 = np.sum(W * W)

    M = sum_WY / sum_W
    var = 1 / sum_W
    se = np.sqrt(var)
    z = M / se
    pval = calc_p_value(z)
    
    return M, pval

def remove_nan(corr,N):
    mask = ~np.isnan(corr)
    return corr[mask], N[mask]

def calc_effect_size(corr):
    effect_size = 0.5 * np.log((1 + corr) / (1 - corr))
    return effect_size

def calc_variance_within(N):
    return 1 / (N - 3)

def calc_weight(V):
    return 1 / V

def calc_p_value(z):
    pval = scipy.stats.norm.sf(abs(z))*2
    pval = np.where(pval == 0.0, 2.2250738585072014e-308, pval)
    return pval

donor_list = pd.read_csv(args.donor_list, sep = "\t", compression='gzip',header=0)
gene_list = pd.read_csv(args.gene_list, sep = "\t", compression='gzip',header=None).iloc[:, 0].tolist()
n = len(gene_list)

m_corr = np.empty((n, n), dtype=np.float64)
np.fill_diagonal(m_corr,1.0)
m_pval = np.empty((n, n), dtype=np.float64)
np.fill_diagonal(m_pval,1.0)

positions = {}
index_counter = 0
header = None

# Read in correlation file line by line
with gzip.open(args.input[0], 'rt') as f:
    for line in f:
        
        values = line.strip("\n").split("\t")
        if header == None:
            header = values
            donor_list = donor_list.set_index('alt_ids').reindex(index=header[1:]).reset_index()
            N = donor_list['count'].tolist()
            continue
        gene1,gene2 = values[0].split("_")
        corr = values[1:]

        if gene1 not in positions:
            positions[gene1] = index_counter
            index_counter += 1
        index1 = positions[gene1]

        if gene2 not in positions:
            positions[gene2] = index_counter
            index_counter += 1
        index2 = positions[gene2]
        
        if args.meta_analysis == 'fishersz':
            corr_combined,pval = fishers_z_transformation(corr,N)
        elif args.meta_analysis == 'mean':
            #corr_combined = np.nanmean(corr_numeric)
            pass
        else: 
            raise ValueError(f"Invalid meta-analysis method: {args.meta_analysis}. Supported methods: 'mean', 'fishersz'")

        m_corr[index1, index2] = corr_combined
        m_corr[index2, index1] = corr_combined
        m_pval[index1, index2] = pval
        m_pval[index2, index1] = pval       

f.close()

genes = positions.keys()
pd.DataFrame(m_corr,index=genes,columns=genes).to_csv(args.output[0], sep="\t",header=True,index=True,compression='gzip')
pd.DataFrame(m_pval,index=genes,columns=genes).to_csv(args.output[1], sep="\t",header=True,index=True,compression='gzip')
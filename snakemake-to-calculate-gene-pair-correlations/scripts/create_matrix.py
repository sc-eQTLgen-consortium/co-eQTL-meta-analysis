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
parser.add_argument("--cohort", required=True, type=str, help="Cohort id")
parser.add_argument("--celltype", required=True, type=str, help="Cell type")
parser.add_argument("--n", required=True, type=str, help="Number of genes")
parser.add_argument("--input", required=True, nargs="+", type=str, help="Input correlation file")
parser.add_argument("--output", required=True, nargs="+", type=str, help="Aggregated output file")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")
    
def fishers_z_transformation(corr,N):
    
    N = np.array(N, dtype=float)
    corr,N = remove_nan(corr,N)
    
    if np.isnan(corr).any() or np.isnan(N).any():
        print("Warning: NaN values encountered in correlation or sample count.")    
    
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

positions = {}
index_counter = 0
header = None
# Read in correlation file line by line
with gzip.open(args.input[0], 'rt') as f:
    for line in f:
        
        values = line.strip("\n").split("\t")
        if header == None:
            header = values
            donor_list = donor_list[donor_list['alt_ids'].isin(header)]
            donor_list.set_index('alt_ids', inplace=True)
            donor_list = donor_list.reindex(index=header[1:])
            N = donor_list['count'].tolist()

            # Filter out genes that aren't expressed in more than 5 samples
            gene_count_dict = {}
            for x in donor_list.index.tolist():
                print(x)
                #TODO: change file location when adding to snakemake
                file = f"donor_gene_list/filtered-genes-{x}-{args.celltype}-top-{args.n}.tsv.gz"
                genes = pd.read_csv(file, sep="\t", compression='gzip', header=None).iloc[:, 0].tolist()

                for gene in genes:
                    gene_count_dict[gene] = gene_count_dict.get(gene, 0) + 1
            
            gene_count_df = pd.DataFrame(gene_count_dict.items(), columns=['Gene', 'Count'])
            gene_count_df['Gene'] = gene_count_df['Gene'].astype(str)
            gene_count_df = gene_count_df[gene_count_df['Count'] > 5]
            gene_list = gene_count_df["Gene"].tolist()
            n = len(gene_list)
            gene_df = pd.DataFrame(gene_list, columns=['GeneName'])
            gene_df.to_csv(args.gene_list, sep='\t', index=False, header=False, quoting=None)
            
            # initialise matrices
            m_corr = np.empty((n, n), dtype=np.float64)
            np.fill_diagonal(m_corr,1.0)
            m_pval = np.empty((n, n), dtype=np.float64)
            np.fill_diagonal(m_pval,1.0)
            continue

        gene1,gene2 = values[0].split("_")

        if gene1 not in gene_list or gene2 not in gene_list:
            continue

        corr = values[1:]
        corr = np.array(corr, dtype=float)

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
            corr_combined = np.nanmean(corr_numeric)
        else: 
            raise ValueError(f"Invalid meta-analysis method: {args.meta_analysis}. Supported methods: 'mean', 'fishersz'")

        m_corr[index1, index2] = corr_combined
        m_corr[index2, index1] = corr_combined
        m_pval[index1, index2] = pval
        m_pval[index2, index1] = pval       

f.close()

genes = list(positions.keys())

while np.isnan(m_corr).any():

    nan_dict = {}

    for i in range(len(m_corr)):
        corr = m_corr[i, :]
        if np.isnan(corr).any():
            nan_dict[i] = np.count_nonzero(np.isnan(corr))
    
    print(sum(nan_dict.values()))
    top_nan = max(nan_dict, key=nan_dict.get)
    print(top_nan)

    m_corr = np.delete(m_corr, top_nan, 0)
    m_corr = np.delete(m_corr, top_nan, 1)

    m_pval = np.delete(m_pval, top_nan, 0)
    m_pval = np.delete(m_pval, top_nan, 1)
 
    genes.pop(top_nan)
    
pd.DataFrame(m_corr, index=genes, columns=genes).replace("", np.nan).to_csv(args.output[0], sep="\t", header=True, index=True, compression='gzip', na_rep="nan")
pd.DataFrame(m_pval, index=genes, columns=genes).replace("", np.nan).to_csv(args.output[1], sep="\t", header=True, index=True, compression='gzip', na_rep="nan")     
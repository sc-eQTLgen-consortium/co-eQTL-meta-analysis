"""
Calculates correlation per donor 
"""
import numpy as np
import pandas as pd 
import argparse
import gzip
import os
from scipy import stats
from scipy.special import betainc
from scipy.io import mmread

parser = argparse.ArgumentParser(description="")
parser.add_argument("--counts", required=True, type=str, help="Normalized counts file path")
parser.add_argument("--method", required=True, type=str, help="Correlation method")
parser.add_argument("--weight", required=True, nargs = "+", type=str, help="Whether to calculate weighted or unweighted correlation (True/False), if true provide filepath to weights")
parser.add_argument("--gene_list_donor", required=True, type=str, help="List of genes that pass filters in this donor")
parser.add_argument("--gene_list", required=True, type=str, help="Genes to include")
parser.add_argument("--output", required=True, type=str, help="output file name")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def calculate_correlation(x,y,weight):
    # Convert data types
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)
    weight = np.array(weight, dtype=float)
    # Calculate correlation with wtd_cor
    return wtd_cor(x=x,y=y,weight=weight)

def wtd_cor(x,y,weight):
    # Normalise weight
    weight = weight / np.mean(weight)
    return onecor_wtd(x=x,y=y,weight=weight)

def onecor_wtd(x,y,weight):
    # Remove NaN
    mask = ~np.isnan(x) & ~np.isnan(y) & ~np.isnan(weight)
    x = x[mask]
    y = y[mask]
    weight = weight[mask]
    
    return lm_wfit(x=stdz(x, weight), y=stdz(y, weight), w=weight)

def stdz(x, weight):
    # Standardize variable
    x = x - wtd_mean(x, weight)
    x = x / np.sqrt(wtd_var(x, weight))

    return x

def wtd_mean(x, weight):
    # Weighted mean
    return np.sum(x * weight) / np.sum(weight)

def wtd_var(x, weight):
    # Weighted variance
    sw = np.sum(weight)
    xbar = np.sum(weight * x) / sw
    wtd_var = np.sum(weight * ((x - xbar) ** 2)) / (sw - 1)
    return wtd_var

def lm_wfit(x, y, w):
    # Weight variables; if unweighted, variables will be multiplied by (sqrt(1))
    wts = np.sqrt(w)
    x = x * wts
    y = y * wts
    # Add intercept column which is a vector of 1's multiplied by the weight vector
    # A matrix is created with the intercept (column 1) and x (column 2)
    X = np.vstack((wts, x)).T
    return lm(X=X, y=y)

def lm(X, y):
    # Store the sample size (n) and degrees of freedom (df)
    n = X.shape[0]
    df = X.shape[1]
    
    inv_m = inverse(X)

    betas = fit(X=X, y=y, inv_m=inv_m)
    
    #y_hat = predict(X=X, betas=betas)
    
    #rss = calc_rss(y=y, y_hat=y_hat)
    
    #std = calc_std(rss=rss, n=n, df=df, inv_m=inv_m)
    
    # Calculate t-value
    #t_value = betas / std

    #p_value = calc_p_value(t_value=t_value, n=n)   
    
    #zscore = calc_z_score(pval=p_value[1],corr=betas[1])
    
    # Only return x-vector statistics
    return betas[1]

def inverse(X): 
    # The first part of solving for beta is taking the inverse of X
    # If columns are colinear, produces a pseudo-inverse matrix
    X_square = X.T.dot(X)
    try:
        return np.linalg.inv(X_square)
    except np.linalg.LinAlgError:
        print("Warning: using pseudo-inverse")
        return np.linalg.pinv(X_square)

def fit(X, y, inv_m=None):
    # Calculate model beta's
    if inv_m is None:
        inv_m = inverse(X)
    return inv_m.dot(X.T).dot(y)
    
def predict(X, betas):
    # Get the predicted value of y
    return np.dot(X, betas)

def calc_rss(y, y_hat):
    # Calculate the residual sum of squares (the difference between y and y_hat ^ 2 summed together)
    res = calc_residuals(y=y, y_hat=y_hat)
    res_squared = res * res
    return np.sum(res_squared)

def calc_residuals(y, y_hat):
    # Calculate difference between y and y_hat
    return y - y_hat

def calc_std(rss, n, df, inv_m):
    # Calculate the standard error
    return np.sqrt(rss / (n - df) * np.diag(inv_m))

def calc_p_value(t_value, n):
    # Calculate p-value
    v = n - 2
    p_value = betainc(v / 2, 1 / 2, v / ((t_value * t_value) + v))
    p_value[p_value < 5e-324] = 5e-324
    return p_value

def calc_z_score(pval,corr):
    # Calculate z-score 
    zscore = np.abs(stats.norm.ppf(pval/2))
    if corr < 0:
        zscore *= -1
    return zscore

def check_gz_nonempty(file_loc):
    # check if a file is empty
    with gzip.open(file_loc, 'rb') as f:
        # try to open the file
        try:
            file_contents = f.read(1)
            return len(file_contents) > 0
        except:
            return False

print("Loading genes")
if check_gz_nonempty(args.gene_list_donor) is False:
    print(' '.join([args.gene_list_donor, 'contains no genes, skipped']))
    sys.exit()

gene_list = pd.read_csv(args.gene_list_donor, sep = "\t", header=None).iloc[:, 0].tolist()

if check_gz_nonempty(args.gene_list) is False:
    print(' '.join([args.gene_list, 'contains no genes, skipped']))
    sys.exit()

genes = pd.read_csv(args.gene_list, sep='\t', header=None).iloc[:,0].to_list()
# make sure these are all strings
genes = [str(element) for element in genes]
# sort the genes, because we sorted the counts by genes in the process_donor.R step
genes.sort()
print("Loading weights")
weight = pd.read_csv(args.weight[1], sep = "\t",header=0)
weight = [i for i in weight.weight]

print("Loading normalised count matrix")
norm_counts = mmread(args.counts)
values = norm_counts.toarray()

all_genes = len(genes)
n_genes = len(gene_list)
print(f"Number of genes filtered: {all_genes - n_genes} out of {all_genes} genes")

# If spearman, rank values
if args.method == 'spearman':
    ranked_values = []
    for i in range(n_genes):
        ranked_values.append(list(stats.rankdata(values[i])))
    values = ranked_values
    
fileout = gzip.open(args.output, "wt", 4)
fileout.write("gene_pair\tcorrelation\n")
print("Calculating correlation")
for i in range(0,all_genes):
    genei = genes[i]
    if genei in gene_list:
        gene_index = gene_list.index(genei)
        x = values[gene_index]
    
    for j in range(i+1,all_genes):
        genej = genes[j]
        if genej in gene_list:
            gene_index = gene_list.index(genej)
            y = values[gene_index]

        if genei not in gene_list or genej not in gene_list:
            fileout.write(f"{genei}_{genej}\t{np.nan}\n")
            continue

        corr = calculate_correlation(x,y,weight)
        fileout.write(f"{genei}_{genej}\t{corr}\n")



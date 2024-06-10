"""
Calculates correlation per donor 
"""
import numpy as np
import pandas as pd 
import argparse
import gzip
from scipy import stats
from scipy.special import betainc

parser = argparse.ArgumentParser(description="")
parser.add_argument("--counts", required=True, type=str, help="Normalised counts file path")
parser.add_argument("--method", required=True, type=str, help="Correlation method")
parser.add_argument("--weight", required=True, nargs = "+", type=str, help="Whether to calculate weighted or unweighted correlation (True/False), if true provide filepath to weights")
parser.add_argument("--cell_barcodes", required=True, type=str, help="Donor names and corresponding indices")
parser.add_argument("--gene_list", required=True, type=str, help="Donor genes to include")
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

    # The first part of solving for beta is taking the inverse of X
    # If columns are colinear, produces a pseud-inverse matrix
    inv_m = inverse(X)
    
    # Calculate model beta's
    betas = fit(X=X, y=y, inv_m=inv_m)
    
    # Get the predicted value of y
    y_hat = predict(X=X, betas=betas)
    
    # Calculate the residual sum of squares (the difference between y and y_hat ^ 2 summed together)
    rss = calc_rss(y=y, y_hat=y_hat)
    
    # Calculate the standard error
    std = calc_std(rss=rss, n=n, df=df, inv_m=inv_m)
    
    # Get the t-value
    t_value = betas / std

    p_value = calc_p_value(t_value=t_value, n=n)   
    
    # Calculate z-score 
    zscore = calc_z_score(pval=p_value[1],corr=betas[1])
    
    # Only return x-vector statistics
    return betas[1], std[1], p_value[1], zscore
    
def inverse(X): 
    X_square = X.T.dot(X)
    try:
        return np.linalg.inv(X_square)
    except np.linalg.LinAlgError:
        print("Warning: using pseudo-inverse")
        return np.linalg.pinv(X_square)

def fit(X, y, inv_m=None):
    if inv_m is None:
        inv_m = inverse(X)
    return inv_m.dot(X.T).dot(y)
    
def predict(X, betas):
    return np.dot(X, betas)

def calc_rss(y, y_hat):
    res = calc_residuals(y=y, y_hat=y_hat)
    res_squared = res * res
    return np.sum(res_squared)

def calc_residuals(y, y_hat):
    return y - y_hat

def calc_std(rss, n, df, inv_m):
    return np.sqrt(rss / (n - df) * np.diag(inv_m))

def calc_p_value(t_value, n):
    """
    last row == stats.t.sf(np.abs(t_value), df=n - 2) * 2 but this
    is faster somehow.

    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.t.html
    https://en.wikipedia.org/wiki/Student%27s_t-distribution (Cumulative distribution function)
    """
    v = n - 2
    p_value = betainc(v / 2, 1 / 2, v / ((t_value * t_value) + v))
    p_value[p_value < 5e-324] = 5e-324
    return p_value

def calc_z_score(pval,corr):
    zscore = np.abs(stats.norm.ppf(pval/2))
    if corr < 0:
        zscore *= -1
    return zscore

# Load barcodes
barcode = set(pd.read_csv(args.cell_barcodes, sep = "\t", compression='gzip',header=None).iloc[:, 0].tolist())
# Load genes
gene_list = pd.read_csv(args.gene_list, sep = "\t", compression='gzip',header=None).iloc[:, 0].tolist()

# Load weight
weight = pd.read_csv(args.weight[1], sep = "\t", compression='gzip',header=0)
weight_dict = dict(zip(weight.index, weight["weight"]))

# Calculate correlation per donor
genes = []
values = []
processed_pairs = set()
with gzip.open(args.counts, 'rt') as f:
    header = None
    mask = None
    
    for line in f:
        values_list = line.split("\t")
        if header == None:
            mask = [i for i, value in enumerate(values_list) if value in barcode]
            header = [values_list[i] for i in mask]
            continue
        genes.append(values_list[0])
        all_values = values_list[1:]
        values.append([all_values[i] for i in mask])
    f.close()

all_genes = len(genes)
n_genes = len(gene_list)
print(f"Number of genes filtered: {all_genes - n_genes}")
weight = [1] * n_genes
if args.weight[0]:
    weight = [weight_dict[barcode] for barcode in header]

if args.method == 'spearman':
    ranked_values = []
    for i in range(n_genes):
        ranked_values.append(list(stats.rankdata(values[i])))
    values = ranked_values
    
fileout = gzip.open(args.output, "wt")

# Calculate correlation
for i in range(0,n_genes):
    x = values[i]
    genei = genes[i]
    
    for j in range(i+1,n_genes):
        y = values[j]
        genej = genes[j]
        
        if genei not in gene_list or genej not in gene_list:
            fileout.write(f"{genei}_{genej}\t{np.nan}\t{np.nan}\t{np.nan}\t{np.nan}\n")
            continue

        corr,std,pval,zscore = calculate_correlation(x,y,weight)
        fileout.write(f"{genei}_{genej}\t{corr}\t{std}\t{pval}\t{zscore}\n")
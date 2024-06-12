"""
Perform eigendecomposition on a n x n matrix
"""

import pandas as pd
import numpy as np
import argparse
import gzip
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input", required=True, nargs="+", type=str, help="Correlation matrices")
parser.add_argument("--output", required=True, nargs="+", type=str, help="Output file name")

args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def eigen_decomposition(matrix):
    if not np.all(np.isfinite(matrix)) or np.any(np.isnan(matrix)):
        nan_rows, nan_cols = np.where(np.isnan(matrix))
        inf_rows, inf_cols = np.where(np.isinf(matrix))
        
        raise ValueError(f"Correlation matrix contains NaNs or infs at rows {nan_rows} and columns {nan_cols}, and infs at rows {inf_rows} and columns {inf_cols}")
    
    eigenvalues, eigenvectors =np.linalg.eig(matrix)
    order = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    return eigenvalues, eigenvectors

def plot_exp_var(exp_var, cum_var, eigenvalues, save_path=args.output[0]):
    plt.figure(figsize=(10, 4))

    # Plot cumulative explained variance
    plt.subplot(1, 2, 1)
    sns.set(style="white")
    line = plt.plot(range(1, len(cum_var) + 1), cum_var, linestyle='-', color="#808080", label='Cumulative Explained Variance')
    plt.axhline(y=0.8, color='r', linestyle='--', alpha=0.7)
    plt.xlabel('Components')
    plt.ylabel('Explained Variance')
    plt.title('')
    plt.legend(loc='upper left')
    #plt.ylim(0, 1.1)
    comp_to_include = [i for i in cumulative_variances if i < 0.8]
    plt.text(len(exp_var) - 50, 0.82, f'n-components = {len(comp_to_include)}', ha='right', va='bottom', fontsize=10,color="r")

    # Plot scatter plot of eigenvalues
    plt.subplot(1, 2, 2)
    plt.scatter(range(1, len(eigenvalues) + 1), eigenvalues, color='#808080', label='Eigenvalues',s=8)
    plt.axhline(y=0, color='r', linestyle='--', alpha=0.7, label='Zero Line')
    plt.xlabel('Components')
    plt.ylabel('Eigenvalues')
    plt.title('Scatter Plot of Eigenvalues')
    plt.legend()
    
    plt.savefig(save_path, bbox_inches="tight")
    
corr_matrix =pd.read_csv(args.input[0], sep = "\t", compression='gzip',header=0,index_col=0)
print(f"len: {len(corr_matrix)}")
file = "/groups/umcg-biogen/tmp03/output/2024-01-11-scMetaBrainConsortium/2023-11-22-Workgroup4-GeneNetworks/Anoek-CalcCorrMatrix/input/LimixAnnotationFile.txt"
annotation_df = pd.read_csv(file, sep = "\t",header=0)
id_mapping = dict(zip(annotation_df['feature_id'], annotation_df['ENSG']))

corr_matrix.columns = corr_matrix.columns.map(id_mapping)
corr_matrix.index = corr_matrix.index.map(id_mapping)
print(list(corr_matrix.index))
unmapped_genes = [gene for gene in corr_matrix.columns if gene not in list(annotation_df['ENSG'])]

print("Unmapped genes: ")
print(unmapped_genes)
print(len(corr_matrix))
corr_matrix = corr_matrix.drop(index=unmapped_genes, columns=unmapped_genes)
print(len(corr_matrix))

if len(list(corr_matrix.index)) != len(set(list(corr_matrix.index))):
    print("Duplicates found!")

if len(list(corr_matrix.columns)) != len(set(list(corr_matrix.columns))):
    print("Duplicates found!")

def find_duplicates(lst):
    duplicates = []
    seen = set()
    for item in lst:
        if item in seen:
            if item not in duplicates:
                duplicates.append(item)
        else:
            seen.add(item)
    return duplicates

# Example usage:
my_list = list(corr_matrix.index)
duplicates = find_duplicates(my_list)
print(duplicates)


# Example usage:
my_list = list(corr_matrix.columns)
duplicates = find_duplicates(my_list)
print(duplicates)




eigenvalues, eigenvectors = eigen_decomposition(corr_matrix)
explained_variance = (eigenvalues / eigenvalues.sum())
cumulative_variances = np.cumsum(explained_variance)
#plot_exp_var(explained_variance,cumulative_variances,eigenvalues)

comp_to_include = [i for i in cumulative_variances if i < 0.8]
n_comp = len(comp_to_include)
n_comp = 50

column_names = [f"PC{i}" for i in range(1, n_comp+1)]
eigenvectors_df = pd.DataFrame(eigenvectors[:,0:n_comp], columns=column_names)
eigenvectors_df.index = corr_matrix.index
eigenvectors_df.to_csv(args.output[1], sep="\t",header=True)

column_names = [f"PC{i}" for i in range(1, len(eigenvalues)+1)]
print(len(column_names))
print(len(eigenvalues))
eigenvalues_df = pd.DataFrame({"Comp": column_names, "eigenValue": eigenvalues})
print(eigenvalues_df)
eigenvalues_df.to_csv(args.output[2], sep="\t",header=True,index=False)

"""
Detects possible sample outliers
"""
import pandas as pd
import numpy as np
#import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import gzip
from itertools import combinations
from scipy.stats import pearsonr

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input", required=True, type=str, help="Sample by gene pair correlation matrix")
parser.add_argument("--output", required=True, nargs = "+", type=str, help="Output sample by gene pair correlation matrix")
parser.add_argument("--donor_list", required=True, type=str, help="Donor list path")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")

def create_corr_matrix(df):

    corr_matrix = pd.DataFrame(index=df.columns, columns=df.columns)
    for sample1_sample2 in combinations(df.columns,2):
        x,y = sample1_sample2
        gene_pairs = df[[x,y]].dropna()
    
        if len(gene_pairs) > 3:
            corr,pval = pearsonr(gene_pairs[x],gene_pairs[y])
            corr_matrix.at[x,y] = corr
            corr_matrix.at[y,x] = corr
        else:
            corr_matrix.at[x,y] = 0.0
            corr_matrix.at[y,x] = 0.0
            
    # Set diagonal to 1
    for sample in df.columns:
        corr_matrix.at[sample, sample] = 1.0
    
    return corr_matrix

class PCA:
    def __init__(self, data):
        self.data = data
        self.X = self.initialize_data()
        self.eigenvalues, self.eigenvectors = self.eigen_decomposition()
        self.explained_var = self.calculate_explained_variance()
        self.transformed_data = self.transform_data()
        self.z_scores = self.calculate_z_scores()

    def initialize_data(self):
        self.data.replace('NA', np.nan, inplace=True)
        X = self.data.values
        return X

    def eigen_decomposition(self):
        self.eigenvalues, self.eigenvectors =np.linalg.eig(self.X)
        order = self.eigenvalues.argsort()[::-1]
        self.eigenvalues = self.eigenvalues[order]
        self.eigenvectors = self.eigenvectors[:, order]
        return self.eigenvalues, self.eigenvectors
    
    def calculate_explained_variance(self):
        self.explained_var = (self.eigenvalues / self. eigenvalues.sum())*100
        return self.explained_var
    
    def transform_data(self):
        self.transformed_data = np.dot(self.X, self.eigenvectors)
        return self.transformed_data 
    
    def calculate_z_scores(self):
        
        pc1 = self.transformed_data[:, 0]
        pc2 = self.transformed_data[:, 1]
    
        mean_pc1, std_dev_pc1 = np.mean(pc1), np.std(pc1)
        mean_pc2, std_dev_pc2 = np.mean(pc2), np.std(pc2)
    
        z_scores_pc1 = (pc1 - mean_pc1) / std_dev_pc1
        z_scores_pc2 = (pc2 - mean_pc2) / std_dev_pc2
    
        return np.column_stack((z_scores_pc1, z_scores_pc2))
    
    def return_outliers(self,threshold=3):
        outliers = (np.abs(self.z_scores) > threshold).any(axis=1)
        return self.data.columns[outliers]

    def plot_pc1_vs_pc2(self,color,outfile):

        data = pd.DataFrame({'PC1': self.transformed_data[:, 0], 
            'PC2': self.transformed_data[:, 1] })

        outliers = (np.abs(self.z_scores) > 3).any(axis=1)        

        sns.set(style="white")
        plt.figure(figsize=(8, 8))
        
        plt.xlim(data['PC1'].min()-1, data['PC1'].max()+1)
        plt.ylim(data['PC2'].min()-1, data['PC2'].max()+1)

        sns.scatterplot(data=data[~outliers], x='PC1', y='PC2', color=color, alpha=0.5, label='Sample')
        sns.scatterplot(data=data[outliers], x='PC1', y='PC2', color='red', alpha=0.5, label='Outlier')

        for i in self.data.columns[outliers]:
            pc1_value = data.loc[data.index[outliers], 'PC1']
            pc2_value = data.loc[data.index[outliers], 'PC2']
            plt.text(pc1_value, pc2_value, i, fontsize=8, color='black', ha='right', va='bottom')

        plt.xlabel(f'PC1 ({self.explained_var[0]:.2f}%)')
        plt.ylabel(f'PC2 ({self.explained_var[1]:.2f}%)')
        num_points = data.shape[0]
        cumulative_variance = self.explained_var[0] + self.explained_var[1]
        plt.text(data['PC1'].max()+1, data['PC2'].min()-0.8, f'Number of samples: {num_points}', ha='right', va='bottom', fontsize=10)
        plt.text(data['PC1'].max()+1, data['PC2'].min()-0.9, f'Cumulative Variance: {cumulative_variance:.2f}%', ha='right', va='bottom', fontsize=10)
        sns.despine()
        plt.savefig(outfile, bbox_inches="tight")

print("Loading sample by gene pair matrix")
df = pd.read_csv(args.input, sep="\t", compression='gzip', header=0, index_col=0)
print(df.head())
outlier_dict = {}

print(f"\nRemoving samples with less than 20 cells")
donor_list = pd.read_csv(args.donor_list, sep = "\t", compression='gzip',header=0)
donor_list["count"] = donor_list["count"].astype(int)
filtered_counts = donor_list[donor_list['count'] < 20]
low_counts = filtered_counts["alt_ids"].tolist()

for x in low_counts:
    outlier_dict[x] = "Cell count < 20 cells"
    df.drop(columns=x, inplace=True)

print("\nRemoving samples with no expressed genes")
no_expressed_genes = df.columns[df.isnull().all()]

for x in no_expressed_genes:
    outlier_dict[x] = "No expressed genes"
    df.drop(columns=x, inplace=True)    

print("\nRemoving samples with no shared expressed genes")
corr_matrix = create_corr_matrix(df)
column_sums = corr_matrix.sum(axis=0)
zero_samples = (column_sums == 1).tolist()

for x in corr_matrix.columns[zero_samples]:
    outlier_dict[x] = "No shared expressed genes across samples"
    df.drop(columns=x, inplace=True)    
print("\nPerforming PCA")
corr_matrix = corr_matrix.apply(pd.to_numeric)

#PCA(corr_matrix).plot_pc1_vs_pc2(color="#808080",outfile=args.output[0])

print("\nRemoving outliers")
outliers = PCA(corr_matrix).return_outliers(threshold=3)

for x in outliers:
    outlier_dict[x] = "Outlier (z-score > 3)"
    
for sample in outlier_dict.keys():
    if sample in df.columns:
        df.drop(columns=sample, inplace=True)
        
# Repeat process to check if more outliers would be detected
corr_matrix = create_corr_matrix(df)
    
print("\nRepeating PCA to check if outliers were removed")
corr_matrix = corr_matrix.apply(pd.to_numeric)
#PCA(corr_matrix).plot_pc1_vs_pc2(color="#808080",outfile=args.output[1])
outliers = PCA(corr_matrix).return_outliers(threshold=3)

outlier_df = pd.DataFrame.from_dict(outlier_dict, orient='index', columns=['Sample Status'])
outlier_df.to_csv(args.output[2],sep="\t",compression="gzip",index=True)
print("\nThe following samples were removed from the matrix:")
print(outlier_df)

print("\nSaving updated matrix")
df.to_csv(args.output[3],sep="\t",compression="gzip",index=True,na_rep="nan")
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations

class PCA:
    def __init__(self,data):
        self.data = data
        self.corr_matrix = self.create_corr_matrix()
        self.X = self.initialize_data()
        self.eigenvalues, self.eigenvectors = self.eigen_decomposition()
        self.explained_var = self.calculate_explained_variance()
        self.transformed_data = self.transform_data()
        self.z_scores = self.calculate_z_scores()

    def create_corr_matrix(self):
        df = self.data
        corr_matrix = pd.DataFrame(index=df.columns, columns=df.columns)
        for sample1_sample2 in combinations(df.columns,2):
            x,y = sample1_sample2
            gene_pairs = df[[x,y]].dropna()

            # Get correlation of each sample pair
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

    def initialize_data(self):
        # Remove NA
        self.corr_matrix.replace('NA', np.nan, inplace=True)
        X = self.corr_matrix.values.astype(float)
        return X

    def eigen_decomposition(self):
        # Convert input matrix to its eigenvectors and eigenvalues
        self.eigenvalues, self.eigenvectors =np.linalg.eig(self.X)
        order = self.eigenvalues.argsort()[::-1]
        self.eigenvalues = self.eigenvalues[order]
        self.eigenvectors = self.eigenvectors[:, order]
        return self.eigenvalues, self.eigenvectors

    def calculate_explained_variance(self):
        # Get explained variance
        self.explained_var = (self.eigenvalues / self. eigenvalues.sum())*100
        return self.explained_var

    def transform_data(self):
        # Matrix dot multiplication to transform data
        self.transformed_data = np.dot(self.X, self.eigenvectors)
        return self.transformed_data

    def calculate_z_scores(self):
        # Get z-scores thresholds for pc1 and pc2
        pc1 = self.transformed_data[:, 0]
        pc2 = self.transformed_data[:, 1]

        mean_pc1, std_dev_pc1 = np.mean(pc1), np.std(pc1)
        mean_pc2, std_dev_pc2 = np.mean(pc2), np.std(pc2)

        z_scores_pc1 = (pc1 - mean_pc1) / std_dev_pc1
        z_scores_pc2 = (pc2 - mean_pc2) / std_dev_pc2

        return np.column_stack((z_scores_pc1, z_scores_pc2))

    def return_outliers(self,threshold=3):
        # Get sample outliers
        outliers = (np.abs(self.z_scores) > threshold).any(axis=1)
        return self.data.columns[outliers]

    def plot_pc1_vs_pc2(self, color, outfile):
        data = pd.DataFrame({'PC1': self.transformed_data[:, 0], 'PC2': self.transformed_data[:, 1]})

        outliers = (np.abs(self.z_scores) > 3).any(axis=1)

        sns.set_style("whitegrid", {'grid.linestyle': '--', 'grid.color': "#E8E8E8"})
        plt.figure(figsize=(6, 6))

        plt.xlim(data['PC1'].min() - 1, data['PC1'].max() + 1)
        plt.ylim(data['PC2'].min() - 1, data['PC2'].max() + 1)

        sns.scatterplot(data=data[~outliers], x='PC1', y='PC2', color=color, alpha=0.5, label='Donor')
        sns.scatterplot(data=data[outliers], x='PC1', y='PC2', color='red', alpha=0.5, label='Outlier')

        for i in self.data.columns[outliers]:
            pc1_value = data.loc[data.index[outliers], 'PC1']
            pc2_value = data.loc[data.index[outliers], 'PC2']
            plt.text(pc1_value, pc2_value, i, fontsize=8, color='black', ha='right', va='bottom')

        # Calculate z-score thresholds
        mean_pc1, std_dev_pc1 = np.mean(self.transformed_data[:, 0]), np.std(self.transformed_data[:, 0])
        mean_pc2, std_dev_pc2 = np.mean(self.transformed_data[:, 1]), np.std(self.transformed_data[:, 1])
        z_score_threshold_pc1 = 3 * std_dev_pc1
        z_score_threshold_pc2 = 3 * std_dev_pc2

        plt.axvline(x=mean_pc1 + z_score_threshold_pc1, color='red', linestyle='--', linewidth=1,alpha=0.5)
        plt.axvline(x=mean_pc1 - z_score_threshold_pc1, color='red', linestyle='--', linewidth=1,alpha=0.5)
        plt.axhline(y=mean_pc2 + z_score_threshold_pc2, color='red', linestyle='--', linewidth=1,alpha=0.5)
        plt.axhline(y=mean_pc2 - z_score_threshold_pc2, color='red', linestyle='--', linewidth=1,alpha=0.5)

        plt.xlabel(f'PC1 ({self.explained_var[0]:.2f}%)')
        plt.ylabel(f'PC2 ({self.explained_var[1]:.2f}%)')
        num_points = data.shape[0]
        cumulative_variance = self.explained_var[0] + self.explained_var[1]
        plt.text(data['PC1'].max() + 0.8, data['PC2'].min() - 0.8, f'Number of donors: {num_points}', ha='right', va='bottom', fontsize=10)
        plt.text(data['PC1'].max() + 0.8, data['PC2'].min() - 0.9, f'Cumulative Variance: {cumulative_variance:.2f}%', ha='right', va='bottom', fontsize=10)
        plt.savefig(outfile, bbox_inches="tight")
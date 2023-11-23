#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#calculating correlations for single donor, to allow for parallelisation in snakemake pipeline

library(Seurat)
library(matrixStats)
library(weights)
library(stringr)

# functions
weighted_corr <- function(normalized_counts_mtr_d, raw_counts_mtr_d,method){
	  
  if(method =='spearman'){
    for(i in 1:nrow(normalized_counts_mtr_d)){
      normalized_counts_mtr_d[i,] <- rank(normalized_counts_mtr_d[i,])
    }
  } else {
      normalized_counts_mtr_d = normalized_counts_mtr_d
  }
  
  corMat = wtd.cor(t(normalized_counts_mtr_d),weight = (colSums(raw_counts_mtr_d!=0)))
  corVec <- corMat$correlation[upper.tri(corMat$correlation)]
  pvVec <- corMat$p.value[upper.tri(corMat$p.value)]
  pvVec[pvVec<1.0e-150] = 1.0e-150
  ZscMat <- abs(qnorm(pvVec/2))
  ZscMat[corVec<0] <- ZscMat[corVec<0]*-1
  mat_pval <- as.data.frame(cbind(corVec, pvVec,ZscMat))
  expanded <- expand.grid(rownames(corMat$correlation), rownames(corMat$correlation))
  mat <- matrix(paste(expanded$Var1, expanded$Var2, sep="_"), nrow=length(rownames(corMat$correlation)), ncol=length(rownames(corMat$correlation)))
  gene_p <- mat[upper.tri(mat)]
	    
  rownames(mat_pval) <- gene_p
  return(mat_pval)
}


get_corr_values_donor_from_rds_file <- function(donor,dir_with_seurat,cohort_id,output_dir,gene_pairs,cell_type){

  cell_type_data <- readRDS(paste0(dir_with_seurat,cell_type,'_',donor,'.rds'))
  cell_type_data_S <- cell_type_data[rownames(genes),]
        
  # Getting a number of counts per cell
  raw_counts_mtr <- as.data.frame(as.matrix(cell_type_data_S@assays$RNA@counts ))
        
  number_of_counts <- as.data.frame(colSums(raw_counts_mtr))
  ncells_per_donor <- as.data.frame(table(cell_type_data$Assignment))
	  
  raw_counts_mtr_d <- as.data.frame(as.matrix(cell_type_data_S[,cell_type_data_S$Assignment ==donor]@assays$RNA@counts ))
  normalized_counts_mtr_d <- as.data.frame(as.matrix(cell_type_data_S[,cell_type_data_S$Assignment ==donor]@assays$data@data ))
	    
  rm(cell_type_data_S)
  rm(cell_type_data)
	      
  raw_counts_mtr_d_non_zero <- as.data.frame(rowSums(raw_counts_mtr_d != 0))
  colnames(raw_counts_mtr_d_non_zero) <- 'counts'
	        
  raw_counts_mtr_d_non_zero$mean <- rowMeans(raw_counts_mtr_d)  # save
  raw_counts_mtr_d_non_zero$UMI <- rowSums(raw_counts_mtr_d)  # save
  genes_with_more_than_10_non_zero <- rownames(raw_counts_mtr_d_non_zero[raw_counts_mtr_d_non_zero$counts >10,]) # save
  raw_counts_mtr_d_non_zero <- raw_counts_mtr_d_non_zero[rownames(raw_counts_mtr_d_non_zero) %in%genes_with_more_than_10_non_zero, ]
  raw_counts_mtr_d_UMI <- as.data.frame(colSums(raw_counts_mtr_d != 0))
  colnames(raw_counts_mtr_d_UMI) <- 'UMI'

  normalized_counts_mtr_d <- normalized_counts_mtr_d[rownames(normalized_counts_mtr_d) %in%genes_with_more_than_10_non_zero, ]

  #getting PCC and SCC  
  
  print("output directory:")
  method = 'pearson'
  print(gsub(' ','',paste0(output_dir,cohort_id,'/','zscore_',cell_type,"_",donor,"_",method,'_weighted.tsv.gz')))
  print("calculating Pearson weighted correlation")
  mat_pval_PCC <- weighted_corr(normalized_counts_mtr_d, raw_counts_mtr_d, method)
  mat_pval_PCC$x <- rownames(mat_pval_PCC)
  mat_pval_PCC <- merge(gene_pairs, mat_pval_PCC, by='x' , all.x=T)
  rownames(mat_pval_PCC) <- mat_pval_PCC$x

  corr_w_gz = gzfile(gsub(' ','',paste0(output_dir,cohort_id,'/','corr_',cell_type,"_",donor,"_",method,'_weighted.tsv.gz')),'w')
  write.table(paste(rownames(mat_pval_PCC),mat_pval_PCC$corVec,sep='\t'), corr_w_gz,sep='\t', row.names = F, quote = F)
  close(corr_w_gz)

  pval_w_gz = gzfile(gsub(' ','',paste0(output_dir,cohort_id,'/','pval_',cell_type,"_",donor,"_",method,'_weighted.tsv.gz')),'w')
  write.table(paste(rownames(mat_pval_PCC),mat_pval_PCC$pvVec,sep='\t'), pval_w_gz,sep='\t', row.names = F, quote = F)
  close(pval_w_gz)

  zsco_w_gz = gzfile(gsub(' ','',paste0(output_dir,cohort_id,'/','zscore_',cell_type,"_",donor,"_",method,'_weighted.tsv.gz')),'w')
  write.table(paste(rownames(mat_pval_PCC),mat_pval_PCC$ZscMat,sep='\t'), zsco_w_gz,sep='\t', row.names = F, quote = F)
  close(zsco_w_gz)

  print("done.")  
}


cell_type = args[1]
cohort_id = args[2]
donor = args[3]
dir_with_seurat = args[4]
path_for_donor_rds = args[5]
output_dir = args[6]

# selection of genes 

sc_data <- readRDS(paste0(dir_with_seurat,paste(cell_type,'.Qced.Normalized.SCs.Rds',sep='')))

mono_expressing_genes <- as.data.frame(rowSums(sc_data@assays$data@data))
colnames(mono_expressing_genes) <- 'sum_of_exp'
mono_expressing_genes$zeros <- rowSums(sc_data@assays$data@data==0)
mono_expressing_genes2 <- as.data.frame(mono_expressing_genes[order(mono_expressing_genes$sum_of_exp, decreasing = T),])

genes.use <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA|^MT|^ENSG|^RPL",rownames(mono_expressing_genes2),value=TRUE, invert=TRUE)

mono_expressing_genes3 <- mono_expressing_genes2[rownames(mono_expressing_genes2) %in% genes.use, ]
genes <- mono_expressing_genes3[1:1000,]

print("Calculating all possible gene pairs unidirectionally")
n=1
gene_pairs = c()
gene_names = rownames(genes)
for (i in seq(length(gene_names)-1)){ 
  gene_pairs = c(gene_pairs, paste(rep(gene_names[n],(length(gene_names)-n)),gene_names[(n+1):length(gene_names)],sep='_')) 
  n = n+1
}

print("Starting correlation functions.")
get_corr_values_donor_from_rds_file(donor,path_for_donor_rds,cohort_id,output_dir,gene_pairs,cell_type)



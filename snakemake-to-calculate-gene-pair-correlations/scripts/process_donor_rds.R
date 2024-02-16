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


get_corr_values_donor_from_rds_file <- function(donor,dir_with_seurat,cohort_id,output_dir,gene_pairs,cell_type,gene_names){

  cell_type_data <- readRDS(paste0(dir_with_seurat,cell_type,'-',donor,'.rds'))
  cell_type_data_S <- cell_type_data[rownames(gene_names),]
        
  # Getting a number of counts per cell  
  raw_counts_mtr_d <- as.data.frame(as.matrix(cell_type_data_S@assays$RNA@counts ))
  normalized_counts_mtr_d <- as.data.frame(as.matrix(cell_type_data_S@assays$data@data ))
	    
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

  corr_w_gz = gzfile(gsub(' ','',paste0(output_dir,cohort_id,'/','corr-',cell_type,"-",donor,"-",method,'-weighted.tsv.gz')),'w')
  write.table(paste(rownames(mat_pval_PCC),mat_pval_PCC$corVec,sep='\t'), corr_w_gz,sep='\t', row.names = F, quote = F)
  close(corr_w_gz)

  pval_w_gz = gzfile(gsub(' ','',paste0(output_dir,cohort_id,'/','pval-',cell_type,"-",donor,"-",method,'-weighted.tsv.gz')),'w')
  write.table(paste(rownames(mat_pval_PCC),mat_pval_PCC$pvVec,sep='\t'), pval_w_gz,sep='\t', row.names = F, quote = F)
  close(pval_w_gz)

  zsco_w_gz = gzfile(gsub(' ','',paste0(output_dir,cohort_id,'/','zscore-',cell_type,"-",donor,"-",method,'-weighted.tsv.gz')),'w')
  write.table(paste(rownames(mat_pval_PCC),mat_pval_PCC$ZscMat,sep='\t'), zsco_w_gz,sep='\t', row.names = F, quote = F)
  close(zsco_w_gz)

  print("done.")  
}


cell_type = args[1]
cohort_id = args[2]
donor = args[3]
gene_list_file = args[4]
path_for_donor_rds = args[5]
output_dir = args[6]
alt_gene_list = args[7]

# selection of genes 

print(paste("Cell type:",cell_type))
print(paste("Cohort:",cohort_id))
print(paste("Donor:",donor))
print(paste("gene list file:",gene_list_file))
print(paste("donor rds path:",path_for_donor_rds))
print(paste("output directory:",output_dir))
print(paste("alternative gene list path:",alt_gene_list))

print("Calculating all possible gene pairs unidirectionally")
if(alt_gene_list == 'nan'){
  gene_names = as.data.frame(read.table(gene_list_file,header=T)[,1])
  colnames(gene_names) = c('gene')
  rownames(gene_names) = gene_names$gene
} else {
  gene_names = as.data.frame(read.table(alt_gene_list, header=T)[,1])
  colnames(gene_names) = c('gene')
  rownames(gene_names) = gene_names$gene
}
print(paste("Number of genes in list:",length(gene_names$gene)))

gene_pairs = c()
n = 1
for (i in seq(length(gene_names$gene)-1)){ 
  gene_pairs = c(gene_pairs, paste(rep(gene_names$gene[n],(length(gene_names$gene)-n)),gene_names$gene[(n+1):length(gene_names$gene)],sep='_')) 
  n = n+1
}

print("Starting correlation functions.")
get_corr_values_donor_from_rds_file(donor,path_for_donor_rds,cohort_id,output_dir,gene_pairs,cell_type,gene_names)



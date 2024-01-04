# output donor tables
# to be done before running snakemake pipeline in order to have the donor information

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

cell_type = args[1]
cohort_id = args[2]
dir_with_seurat = args[3]
donor_rds_dir = args[4]
gene_list_out = args[5]
alt_gene_list = args[6]

library(Seurat)

output_dir <- donor_rds_dir

sc_data <- readRDS(paste0(dir_with_seurat,paste(cell_type,'.Qced.Normalized.SCs.Rds',sep='')))

print("seurat object loaded, continuing with analysis")

donortab  <- as.data.frame(table(sc_data$Assignment))
donor_list <- list(donortab[donortab$Freq >10,]$Var1)
names(donor_list) <- c('original_labels')
donor_list$filt_labels <- gsub(pattern='_', replacement='', x=donor_list$original_labels)
donor <- unique(sc_data$Assignment)[2]

write.table(donor_list, paste0(output_dir,'/',cohort_id,'_',cell_type,'_donor_list.tsv'),sep='\t', row.names = F, quote = F)
write.table(as.data.frame(table(sc_data$Assignment)), paste0(output_dir,'/',cohort_id,'_',cell_type,'_donor_counts.tsv'),sep='\t', row.names = F, quote = F)

print("donors filtered.")

if (alt_gene_list == 'nan'){
  print("Creating list of most expressed genes.")
  
  expressing_genes <- as.data.frame(rowSums(sc_data@assays$data@data))
  expressing_genes <- cbind(rownames(sc_data@assays$data@data),expressing_genes)
  print("data frame loaded")
  
  colnames(expressing_genes) <- c('gene_names','sum_of_exp')
  expressing_genes$sum_of_exp <- as.numeric(expressing_genes$sum_of_exp)
  expressing_genes <- expressing_genes[order(expressing_genes$sum_of_exp, decreasing = T),]
  
  ###warning following code removes genes with ensemble ID
  #genes.use <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA|^MT|^ENSG|^RPL",rownames(expressing_genes),value=TRUE, invert=TRUE)
  #print("Ribosomal and mitochondrial genes removed.")
  #expressing_genes <- expressing_genes[rownames(expressing_genes) %in% genes.use,]
  
  ### outputting donor rds
  print("Selecting top expressed genes")
  genes <- expressing_genes[1:3000,]
  
  write.table(rownames(genes), paste(gene_list_out,cohort_id,'_',cell_type,'_genes.tsv',sep=''), sep='\t',row.names=F,quote=F)

} else{
  print("A gene list was provided, no new gene list will be generated.")
}

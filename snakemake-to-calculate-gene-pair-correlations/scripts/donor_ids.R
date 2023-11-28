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

donortab  <- as.data.frame(table(sc_data$Assignment))
donor_list <- donortab[donortab$Freq >10,]$Var1
donor <- unique(sc_data$Assignment)[2]

expressing_genes <- as.data.frame(rowSums(sc_data@assays$data@data))
colnames(expressing_genes) <- 'sum_of_exp'
expressing_genes$zeros <- rowSums(sc_data@assays$data@data==0)
expressing_genes2 <- as.data.frame(expressing_genes[order(expressing_genes$sum_of_exp, decreasing = T),])

genes.use <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA|^MT|^ENSG|^RPL",rownames(expressing_genes2),value=TRUE, invert=TRUE)

expressing_genes3 <- expressing_genes2[rownames(expressing_genes2) %in% genes.use, ]

### outputting donor rds
print("Selecting top 1000 genes")
if(alt_gene_list=='nan'){
  genes <- expressing_genes3[1:1000,]
} else {
  genes_to_use = read.table(alt_gene_list, header=T)[,1]
  genes = expressing_genes3[rownames(expressing_genes3) %in% genes_to_use,]
}

write.table(rownames(genes), paste(gene_list_out,cohort_id,'_',cell_type,'_genes.tsv',sep=''), sep='\t',row.names=F,quote=F)
write.table(sort(donor_list), paste0(output_dir,'/',cohort_id,'_',cell_type,'_donor_list.tsv'),sep='\t', row.names = F, quote = F)
write.table(as.data.frame(table(sc_data$Assignment)), paste0(output_dir,'/',cohort_id,'_',cell_type,'_donor_counts.tsv'),sep='\t', row.names = F, quote = F)


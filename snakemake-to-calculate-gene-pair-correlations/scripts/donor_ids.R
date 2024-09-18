# output donor tables
# to be done before running snakemake pipeline in order to have the donor information

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

cell_type = args[1]
cohort_id = args[2]
seurat_object_path = args[3]
donor_rds_dir = args[4]
gene_list_out = args[5]
alt_gene_list = args[6] # THIS DOES NOT SEEM TO BE USED!
smf = args[7]
qtl_input_path = args[8]
seurat_assignment_column = args[9]

print(paste0("donor_ids.R run with", "\n",
"param1: cell_type: ",cell_type, "\n",
"param2: cohort: ",cohort_id, "\n",
"param3: seurat object: ",seurat_object_path, "\n",
"param4: output directory: ",donor_rds_dir, "\n",
"param5: standard_gene_list: ",gene_list_out, "\n",
"param6: alternative_gene_list: ",alt_gene_list, "\n"  # THIS DOES NOT SEEM TO BE USED!,
"param7: sample mapping file: ",smf, "\n",
"param8: wg3 QTL inputs: ",qtl_input_path, "\n",
"param9: sample assignment Seurat column: ",seurat_assignment_column))

library(Seurat)
library(stringr)

output_dir <- donor_rds_dir

sc_data <- readRDS(seurat_object_path)

print("seurat object loaded, continuing with analysis")

print("Filtering RDS file with smf and Pcs files")
smf = read.csv(smf,sep='\t')
Pcs = read.csv(paste(qtl_input_path,cell_type,".qtlInput.Pcs.txt",sep=''),sep='\t')

# Filter out donors that are not in the smf file
x = sc_data@meta.data[[seurat_assignment_column]] %in% smf$genotype_id
if (sum(x) != length(sc_data@meta.data[[seurat_assignment_column]])){
  cells.use = colnames(sc_data[,x])
  subset_file = subset(sc_data, cells=cells.use)
  sc_data = subset_file
}

# Filter out donors that are not in Pcs file
Pcs_donors = str_split(Pcs$X,';',simplify=T)[,1]

Pcs_donors = Pcs_donors[Pcs_donors %in% sc_data@meta.data[[seurat_assignment_column]]]
x=sc_data@meta.data[[seurat_assignment_column]] %in% Pcs_donors
if (sum(x) != length(sc_data@meta.data[[seurat_assignment_column]])){
  cells.use = colnames(sc_data[,x])
  subset_file = subset(sc_data, cells=cells.use)
  sc_data = subset_file
}


donortab  <- as.data.frame(table(sc_data[[seurat_assignment_column]]))
donor_list <- list(donortab[donortab$Freq >10,]$Var1)
names(donor_list)=c('original_labels')
donor_list$filt_labels <- gsub(pattern='_', replacement='', x=donor_list$original_labels)
donor <- unique(sc_data[[seurat_assignment_column]])[2]

print("donors filtered.")

expressing_genes <- as.data.frame(rowSums(sc_data@assays$data@data))
expressing_genes = cbind(rownames(sc_data@assays$data@data),expressing_genes)
print("data frame loaded")

colnames(expressing_genes) <- c('gene_names','sum_of_exp')
expressing_genes$sum_of_exp <- as.numeric(expressing_genes$sum_of_exp)
expressing_genes <- expressing_genes[order(expressing_genes$sum_of_exp, decreasing = T),]

### outputting donor rds
print("Selecting top expressed genes")
genes <- expressing_genes[1:3000,]

write.table(rownames(genes), paste(gene_list_out,cohort_id,'-',cell_type,'-genes.tsv',sep=''), sep='\t',row.names=F,quote=F)
write.table(donor_list, paste0(output_dir,'/',cohort_id,'-',cell_type,'-donor-list.tsv'),sep='\t', row.names = F, quote = F)
write.table(as.data.frame(table(sc_data[[seurat_assignment_column]])), paste0(output_dir,'/',cohort_id,'-',cell_type,'-donor-counts.tsv'),sep='\t', row.names = F, quote = F)


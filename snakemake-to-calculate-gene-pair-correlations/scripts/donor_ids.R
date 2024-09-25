#!/usr/bin/env Rscript

# output donor tables
# to be done before running snakemake pipeline in order to have the donor information


####################
# Functions        #
####################


normalize_mj <- function(seurat_object) {
  # get the count matrix where we have the correct cell type
  count_matrix <- GetAssayData(seurat_object, slot = "counts")
  # ignore genes that are never expressed
  count_matrix <-  count_matrix[which(rowSums(count_matrix) != 0), ]
  # create new object to store the counts in
  norm_count_matrix <- count_matrix
  # do mean sample-sum normalization
  sample_sum_info = colSums(norm_count_matrix)
  mean_sample_sum = mean(sample_sum_info)
  sample_scale = sample_sum_info / mean_sample_sum
  # divide each column by sample_scale
  norm_count_matrix@x <- norm_count_matrix@x / rep.int(sample_scale, diff(norm_count_matrix@p))
  if ('layers' %in% slotNames(seurat_object[['RNA']])) {
    print('using Seurat v5 style \'layer\'')
    seurat_object[['data']] <- CreateAssayObject(data = norm_count_matrix)

  } else {
    print('using Seurat v3/4 style \'slot\'')
    seurat_object[['data']] <- CreateAssayObject(data = norm_count_matrix)
  }
  return(seurat_object)
}



####################
# Main Code        #
####################

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

cat(paste0("donor_ids.R run with", "\n",
"param1: cell_type: ",cell_type, "\n",
"param2: cohort: ",cohort_id, "\n",
"param3: seurat object: ",seurat_object_path, "\n",
"param4: output directory: ",donor_rds_dir, "\n",
"param5: standard_gene_list: ",gene_list_out, "\n",
"param6: alternative_gene_list: ",alt_gene_list, "\n",  # THIS DOES NOT SEEM TO BE USED!,
"param7: sample mapping file: ",smf, "\n",
"param8: wg3 QTL inputs: ",qtl_input_path, "\n",
"param9: sample assignment Seurat column: ",seurat_assignment_column, '\n'))

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


donortab  <- as.data.frame(table(sc_data@meta.data[[seurat_assignment_column]]))
donor_list <- list(donortab[donortab$Freq >10,]$Var1)
names(donor_list)=c('original_labels')
donor_list$filt_labels <- gsub(pattern='_', replacement='', x=donor_list$original_labels)
#donor <- unique(sc_data[[seurat_assignment_column]])[2]

print("donors filtered.")

# add MJ normalization if not already present
if (!('data' %in% names(sc_data))) {
  cat("\nadding MJ style normalization\n")
  sc_data <- normalize_mj(sc_data)
}

# get the expressing genes
expressing_genes <- NULL
if ('data' %in% names(sc_data)) {
  # use either the data slot that was created in WG3
  print('using data slot/layer')
  if ('layers' %in% slotNames(sc_data[['data']])) {
    print('using Seurat v5 style \'layer\'')
    expressing_genes <- as.data.frame(rowSums(sc_data@layers$data@layers$data))
    # get the dataframe of feature names
    feature_names_df <- data.frame(sc_data@assays$data@features)
    # get the feature names for this layer
    feature_names <- rownames(feature_names_df[feature_names_df$data == T, ])
    # add as a column to the expressed genes
    expressing_genes <- cbind(feature_names, expressing_genes)
  } else {
    print('using Seurat v3/4 style \'slot\'')
    expressing_genes <- as.data.frame(rowSums(sc_data@assays$data@data))
    expressing_genes <- cbind(rownames(sc_data@assays$data@data), expressing_genes)
  }
} else if ('RNA' %in% names(sc_data)) {
  # or the regular RNA slot used in Seurat
  print('using data slot/layer')
  if ('layers' %in% slotNames(sc_data[['RNA']])) {
    print('using Seurat v5 style \'layer\'')
    expressing_genes <- as.data.frame(rowSums(sc_data@assays$RNA@layers$data))
    # get the dataframe of feature names
    feature_names_df <- data.frame(sc_data@assays$RNA@features)
    # get the feature names for this layer
    feature_names <- rownames(feature_names_df[feature_names_df$data == T, ])
    # add as a column to the expressed genes
    expressing_genes <- cbind(feature_names, expressing_genes)
  } else {
    print('using Seurat v3/4 style \'slot\'')
    expressing_genes <- as.data.frame(rowSums(sc_data@assays$RNA@data))
    expressing_genes <- cbind(rownames(sc_data@assays$RNA@data), expressing_genes)
  }
} else {
  stop(paste('seurat object contains neither \'data\' nor \'RNA\' attributes'))
}
print("data frame loaded")

colnames(expressing_genes) <- c('gene_names','sum_of_exp')
expressing_genes$sum_of_exp <- as.numeric(expressing_genes$sum_of_exp)
expressing_genes <- expressing_genes[order(expressing_genes$sum_of_exp, decreasing = T),]

### outputting donor rds
print("Selecting top expressed genes")
genes <- expressing_genes[1:3000,]

write.table(genes[['gene_names']], paste(gene_list_out,cohort_id,'-',cell_type,'-genes.tsv',sep=''), sep='\t',row.names=F,quote=F, col.names = F)
write.table(donor_list, paste0(output_dir,'/',cohort_id,'-',cell_type,'-donor-list.tsv'),sep='\t', row.names = F, quote = F)
write.table(as.data.frame(table(sc_data[[seurat_assignment_column]])), paste0(output_dir,'/',cohort_id,'-',cell_type,'-donor-counts.tsv'),sep='\t', row.names = F, quote = F)


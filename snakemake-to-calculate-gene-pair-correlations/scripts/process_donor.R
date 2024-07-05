#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

cell_type = args[1]
cohort_id = args[2]
dir_with_seurat = args[3]
donor_rds_dir = args[4]
genes_to_use = args[5]
smf = args[6]

library(Seurat)
library(matrixStats)
library(weights)
library(stringr)
library(Matrix)

print(paste("cell type:",cell_type))
print(paste("cohort id:",cohort_id))
print(paste("seurat file: ",dir_with_seurat,cell_type,".Qced.Normalized.SCs.Rds",sep=''))
print(paste("output directory: ",donor_rds_dir,cohort_id,'/','donor_rds','/',sep=''))
print(paste("Sample mapping file: ",smf,sep=''))
print(paste("Pcs file: ",dir_with_seurat,cell_type,".qtlInput.Pcs.txt",sep=''))
print(paste("Using gene list: ", genes_to_use, sep=''))

# selection of genes 

print("Loading big rds file")
sc_data <- readRDS(paste0(dir_with_seurat,paste(cell_type,'.Qced.Normalized.SCs.Rds',sep='')))

print("Filtering RDS file with smf and Pcs files")
smf = read.csv(smf,sep='\t')
Pcs = read.csv(paste(dir_with_seurat,cell_type,".qtlInput.Pcs.txt",sep=''),sep='\t')

# Filter out donors that are not in the smf file
x = sc_data@meta.data$Assignment %in% smf$genotype_id
if (sum(x) != length(sc_data@meta.data$Assignment)){
  cells.use = colnames(sc_data[,x])
  subset_file = subset(sc_data, cells=cells.use)
  sc_data = subset_file
}

# Filter out donors that are not in the Pcs file
if (length(Pcs$X) != 0){
  Pcs_donors = str_split(Pcs$X,';',simplify=T)[,1]
} else {
  Pcs_donors = str_split(rownames(Pcs),';',simplify=T)[,1]
}
Pcs_donors = Pcs_donors[Pcs_donors %in% sc_data@meta.data$Assignment]
x=sc_data@meta.data$Assignment %in% Pcs_donors
if (sum(x) != length(sc_data@meta.data$Assignment)){
  cells.use = colnames(sc_data[,x])
  subset_file = subset(sc_data, cells=cells.use)
  sc_data = subset_file
}

cat("\nUsing provided gene list")
genes_to_use = read.table(genes_to_use, header=F)[,1]

sc_data_filtered <- sc_data[genes_to_use,]
rm(sc_data)

donortab  <- as.data.frame(table(sc_data_filtered$Assignment))
donor_list <- donortab[donortab$Freq >10,]$Var1

# Filter low expressed genes per donor
for(donor in donor_list){
  cat("\nProcessing donor:", donor)
  donor_rds <- sc_data_filtered[,sc_data_filtered$Assignment == donor ]
  #cell_barcodes <- colnames(donor_rds)
  #barcodesOut = paste0(donor_rds_dir,cohort_id,"/donor_barcodes/barcodes-",donor,"-",cell_type,".tsv.gz")
  #write.table(cell_barcodes, file=gzfile(barcodesOut),sep="\t",row.names = FALSE, col.names=FALSE, quote = FALSE) 
  
  # Gene filtering
  raw_counts <- as.data.frame(as.matrix(donor_rds@assays$RNA@counts))

  raw_counts_df <- as.data.frame(rowSums(raw_counts != 0))
  raw_counts_df$gene_name <- rownames(raw_counts_df)
  colnames(raw_counts_df) <- c('counts','gene_names')
  raw_counts <- raw_counts[str_order(rownames(raw_counts)),]
  filtered_genes <- raw_counts_df$gene_names[raw_counts_df$counts > 10] 
  cat(" Number of genes after filtering:",length(filtered_genes))

  donor_rds <- donor_rds[filtered_genes,]

  genesOut  <- paste0(donor_rds_dir,cohort_id,"/donor_gene_list/filtered-genes-",donor,"-",cell_type,".tsv.gz")
  write.table(filtered_genes, gzfile(genesOut),sep="\t",row.names = FALSE, col.names=FALSE, quote = FALSE)

  # Extract weights per donor
  raw_counts <- raw_counts[rownames(raw_counts) %in% filtered_genes, ]
  weights <- data.frame(weight = colSums(raw_counts != 0))
  weightOutput <- paste0(donor_rds_dir,cohort_id,"/donor_weight/correlation-weight-",donor,"-",cell_type,".tsv.gz")
  write.table(weights, gzfile(weightOutput),sep="\t",row.names = TRUE, quote = FALSE)

  # Save norm counts
  norm_sparse <- donor_rds@assays$data@data
  norm_sparse=norm_sparse[str_order(rownames(norm_sparse)),]
  countOutput <- paste0(donor_rds_dir,cohort_id,"/counts/normalized-counts-",donor,"-",cell_type,".mtx")
  writeMM(norm_sparse, file=countOutput)
}

print("Process Complete.")

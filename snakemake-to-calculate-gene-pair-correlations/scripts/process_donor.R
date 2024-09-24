#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

cell_type = args[1]
cohort_id = args[2]
seurat_object_path = args[3]
donor_rds_dir = args[4]
genes_to_use_loc = args[5]
smf = args[6]
qtl_input_path = args[7]
seurat_assignment_column = args[8]
weight_method = args[9]

library(Seurat)
library(matrixStats)
library(weights)
library(stringr)
library(Matrix)

cat(paste("param1: cell type:",cell_type, "\n", 
"param2: cohort id: ",cohort_id, "\n", 
"param3: seurat file: ",seurat_object_path, "\n", 
"param4: output directory: ",donor_rds_dir,cohort_id,'/','donor_rds','/', "\n", 
"param5: Using gene list: ", genes_to_use_loc, "\n", 
"param6: Sample mapping file: ", smf, "\n", 
"param7: wg3 QTL inputs: ", qtl_input_path, "\n", 
"param8: sample assignment Seurat column: ", seurat_assignment_column, "\n", 
"param9: weighting method: ", weight_method, "\n", 
"inferred PCs file: ", qtl_input_path,cell_type,".qtlInput.Pcs.txt", "\n", 
sep=''))

# selection of genes
print("Loading big rds file")
sc_data <- readRDS(seurat_object_path)

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

# Filter out donors that are not in the Pcs file
Pcs_donors = str_split(Pcs$X,';',simplify=T)[,1]
Pcs_donors = Pcs_donors[Pcs_donors %in% sc_data@meta.data[[seurat_assignment_column]]]
x=sc_data@meta.data[[seurat_assignment_column]] %in% Pcs_donors
if (sum(x) != length(sc_data@meta.data[[seurat_assignment_column]])){
  cells.use = colnames(sc_data[,x])
  subset_file = subset(sc_data, cells=cells.use)
  sc_data = subset_file
}

# use gene list if supplied
genes_to_use <- NULL
if (!is.na(genes_to_use_loc) & 
 !is.null(genes_to_use_loc) & 
 genes_to_use_loc != '' & 
 genes_to_use_loc != 'nan' & 
 genes_to_use_loc != 'None') {
 cat("\nUsing provided gene list\n")
  genes_to_use = read.table(genes_to_use_loc, header=F)[,1]
  # check if there are any genes overlapping
  if (length(intersect(genes_to_use, rownames(sc_data))) == 0) {
    stop(paste('no genes overlap between the gene list', genes_to_use_loc, 'and the object', seurat_object_path))
  }
  sc_data_filtered <- sc_data[genes_to_use,]
} else {
  cat("\nno gene list provided (or nan/None), using all genes\n")
  sc_data_filtered <- sc_data
}
rm(sc_data)

donortab  <- as.data.frame(table(sc_data_filtered@meta.data[[seurat_assignment_column]]))
donor_list <- donortab[donortab$Freq >10,]$Var1

# Filter low expressed genes per donor
for(donor in donor_list){
  cat("\nProcessing donor:", donor)
  donor_rds <- sc_data_filtered[,sc_data_filtered[[seurat_assignment_column]] == donor ]
  #cell_barcodes <- colnames(donor_rds)
  #barcodesOut = paste0(donor_rds_dir,cohort_id,"/donor_barcodes/barcodes-",donor,"-",cell_type,".tsv.gz")
  #write.table(cell_barcodes, file=gzfile(barcodesOut),sep="\t",row.names = FALSE, col.names=FALSE, quote = FALSE) 
  
  # Gene filtering
  raw_counts <- NULL
  if ('layers' %in% slotNames(donor_rds[['RNA']])) {
    print('using Seurat v5 style \'layer\'')
    raw_counts <- donor_rds@assays$RNA@layers$counts
    # get the dataframe of feature names
    feature_names_df <- data.frame(donor_rds@assays$RNA@features)
    # get the feature names for this layer
    feature_names <- rownames(feature_names_df[feature_names_df$counts == T, ])
    # set as the rownames
    rownames(raw_counts) <- feature_names
  } else {
    print('using Seurat v3/4 style \'slot\'')
    raw_counts <- as.data.frame(as.matrix(donor_rds@assays$RNA@counts))
  }
  

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
  raw_counts <- raw_counts[rownames(raw_counts) %in% filtered_genes, , drop = F]
  weights <- NULL
  # do weighting depending on the method
  if (weight_method == 'expression') {
    # this is based on the total expression
    weights <- data.frame(weight = colSums(raw_counts))
  } else if(weight_method == 'zeroes') {
    # this is based on the number of non-zero counts
    weights <- data.frame(weight = colSums(raw_counts != 0))
  } else if(weight_method == 'none') {
    # this just dummies everything to 1, so it will not weight at all
    weights <- data.frame(weight = rep(1, times = ncol(raw_counts)))
  } else {
    # if we don't have a valid weighting method, we'll crash
    stop(paste('invalid weighting method, valid options are \'expression\', \'zeroes\', or \'none\', you supplied', weight_method))
  }
  weightOutput <- paste0(donor_rds_dir,cohort_id,"/donor_weight/correlation-weight-",donor,"-",cell_type,".tsv.gz")
  write.table(weights, gzfile(weightOutput),sep="\t",row.names = TRUE, quote = FALSE)

  # Save norm counts
  norm_sparse <- NULL
  if ('data' %in% names(donor_rds)) {
    # use either the data slot that was created in WG3
    print('using data slot/layer')
    if ('layers' %in% slotNames(donor_rds[['data']])) {
      print('using Seurat v5 style \'layer\'')
      norm_sparse <- donor_rds@assays$data@layers$data
      # get the dataframe of feature names
      feature_names_df <- data.frame(donor_rds@assays$data@features)
      # get the feature names for this layer
      feature_names <- rownames(feature_names_df[feature_names_df$data == T, , drop = F])
      # set as the rownames
      rownames(norm_sparse) <- feature_names
    } else {
      print('using Seurat v3/4 style \'slot\'')
      norm_sparse <- donor_rds@assays$data@data
    }
  } else if ('RNA' %in% names(donor_rds)) {
   # use either the data slot that was created in WG3
    print('using RNA slot/layer')
    if ('layers' %in% slotNames(donor_rds[['RNA']])) {
      print('using Seurat v5 style \'layer\'')
      norm_sparse <- donor_rds@assays$RNA@layers$data
      # get the dataframe of feature names
      feature_names_df <- data.frame(donor_rds@assays$RNA@features)
      # get the feature names for this layer
      feature_names <- rownames(feature_names_df[feature_names_df$data == T, , drop = F])
      # set as the rownames
      rownames(norm_sparse) <- feature_names
    } else {
      print('using Seurat v3/4 style \'slot\'')
      norm_sparse <- donor_rds@assays$RNA@data
    }
  } else {
  stop(paste('seurat object contains neither \'data\' nor \'RNA\' attributes'))
  }
  norm_sparse=norm_sparse[str_order(rownames(norm_sparse)),]
  countOutput <- paste0(donor_rds_dir,cohort_id,"/counts/normalized-counts-",donor,"-",cell_type,".mtx")
  writeMM(norm_sparse, file=countOutput)
}

print("Process Complete.")

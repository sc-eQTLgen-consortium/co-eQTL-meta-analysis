#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Marc-Jan Bonder, Dan Kaptijn, Roy Oelen
# Name: process_donor.R
# Function: get donor-specific mtx files, filtered feature lists and weights for each donor
############################################################################################################################

####################
# libraries        #
####################

library(Seurat)
library(matrixStats)
library(weights)
library(stringr)
library(Matrix)
library(optparse)

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

# make command line options
option_list <- list(
  make_option(c("-c", "--cell_type"), type="character", default=NULL,
              help="cell type working on", metavar="character"),
  make_option(c("-r", "--cohort_id"), type="character", default=NULL,
              help="name of cohort", metavar="character"),
  make_option(c("-s", "--seurat_object_path"), type="character", default=NULL,
              help="location to Seurat object (rds)", metavar="character"),
  make_option(c("-d", "--donor_rds_dir"), type="character", default=NULL,
              help="where to store filtered gene lists, counts and weights", metavar="character"),
  make_option(c("-g", "--genes_to_use_loc"), type="character", default=NULL,
              help="file with list of genes to include", metavar="character"),
  make_option(c("-m", "--smf"), type="character", default=NULL,
              help="location of sample mapping file", metavar="character"),
  make_option(c("-q", "--qtl_input_path"), type="character", default=NULL,
              help="folder that contains QTL output from WG3, so qtlInput.txt.gz etc.", metavar="character"),
  make_option(c("-a", "--seurat_assignment_column"), type="character", default=NULL,
              help="name of column in Seurat metadata that has the assignment of the cell to a donor", metavar="character"),
  make_option(c("-w", "--weight_method"), type="character", default=NULL,
              help="which method to use for weighting: expression/zeroes/none", metavar="character"), 
  make_option(c("-n", "--skip_filter"), type="character", default=0,
              help="skip filtering step of Seurat object with smf and QTL inputs: 0/1, default is 0", metavar="numeric")            
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# read the parameters
cell_type = opt[['cell_type']]
cohort_id = opt[['cohort_id']]
seurat_object_path = opt[['seurat_object_path']]
donor_rds_dir = opt[['donor_rds_dir']]
genes_to_use_loc = opt[['genes_to_use_loc']]
smf = opt[['smf']]
qtl_input_path = opt[['qtl_input_path']]
seurat_assignment_column = opt[['seurat_assignment_column']]
weight_method = opt[['weight_method']]
skip_filter_int = opt[['skip_filter']]

# make the filter skip a boolean
skip_filter <- F
if (skip_filter_int == 1) {
  skip_filter <- T
} else if(skip_filter < 0 | skip_filter > 1) {
  stop("skip_filter needs to be '0' or '1'\n")
}

# args = commandArgs(trailingOnly=TRUE)
# cell_type = args[1]
# cohort_id = args[2]
# seurat_object_path = args[3]
# donor_rds_dir = args[4]
# genes_to_use_loc = args[5]
# smf = args[6]
# qtl_input_path = args[7]
# seurat_assignment_column = args[8]
# weight_method = args[9]

cat(paste("cell type:",cell_type, "\n", 
"cohort id: ",cohort_id, "\n", 
"seurat file: ",seurat_object_path, "\n", 
"output directory: ",donor_rds_dir,cohort_id,'/','donor_rds','/', "\n", 
"Using gene list: ", genes_to_use_loc, "\n", 
"Sample mapping file: ", smf, "\n", 
"wg3 QTL inputs: ", qtl_input_path, "\n", 
"sample assignment Seurat column: ", seurat_assignment_column, "\n", 
"weighting method: ", weight_method, "\n", 
"inferred PCs file: ", qtl_input_path,cell_type,".qtlInput.Pcs.txt", "\n", 
"skippping filter: ", skip_filter, "\n", 
sep=''))

# selection of genes
print(paste("Loading big rds file", seurat_object_path, '...'))
sc_data <- readRDS(seurat_object_path)

# init smf and PC variables
smf <- NULL
Pcs <- NULL

# read sample mapping and PC files
if (!skip_filter) {
  print("Filtering RDS file with smf and Pcs files")
  smf = read.table(smf, sep='\t')
  Pcs = read.table(paste(qtl_input_path, '/', cell_type, ".qtlInput.Pcs.txt", sep=''), sep='\t', row.names = 1)
}

# Filter out donors that are not in the smf file
if (!skip_filter) {
  x = sc_data@meta.data[[seurat_assignment_column]] %in% smf$genotype_id
  if (sum(x) != length(sc_data@meta.data[[seurat_assignment_column]])){
    # get which we removed for debugging purposes
    donors_removed <- setdiff(unique(sc_data@meta.data[[seurat_assignment_column]]), smf$genotype_id)
    warning(paste('removed samples not present in smf:', paste(donors_removed, collapse = ',')))
    # get the cellbarcodes to include
    cells.use = colnames(sc_data[,x])
    # subset the Seurat object
    subset_file = subset(sc_data, cells=cells.use)
    sc_data = subset_file
  }
}

# Filter out donors that are not in the Pcs file
if (!skip_filter) {
  #Pcs_donors = str_split(Pcs$X,';',simplify=T)[,1]
  Pcs_donors = str_split(rownames(Pcs),';',simplify=T)[,1]
  Pcs_donors = Pcs_donors[Pcs_donors %in% sc_data@meta.data[[seurat_assignment_column]]]
  x=sc_data@meta.data[[seurat_assignment_column]] %in% Pcs_donors
  if (sum(x) != length(sc_data@meta.data[[seurat_assignment_column]])){
    # get which we removed for debugging purposes
    donors_removed <- setdiff(unique(sc_data@meta.data[[seurat_assignment_column]]), Pcs_donors)
    warning(paste('removed samples not present in PC file:', paste(donors_removed, collapse = ',')))
    # get the cellbarcodes to include
    cells.use = colnames(sc_data[,x])
    # subset the Seurat object
    subset_file = subset(sc_data, cells=cells.use)
    sc_data = subset_file
  }
}

# add MJ normalization if not already present
if (!('data' %in% names(sc_data))) {
  cat("\nadding MJ style normalization\n")
  sc_data <- normalize_mj(sc_data)
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
  cat(paste("\nProcessing donor:", donor, "\n"))
  donor_rds <- sc_data_filtered[,sc_data_filtered[[seurat_assignment_column]] == donor ]
  #cell_barcodes <- colnames(donor_rds)
  #barcodesOut = paste0(donor_rds_dir,cohort_id,"/donor_barcodes/barcodes-",donor,"-",cell_type,".tsv.gz")
  #write.table(cell_barcodes, file=gzfile(barcodesOut),sep="\t",row.names = FALSE, col.names=FALSE, quote = FALSE) 
  
  # Gene filtering
  raw_counts <- NULL
  if ('layers' %in% slotNames(donor_rds[['RNA']])) {
    print('using Seurat v5 style \'layer\'')
    raw_counts <- donor_rds@assays$RNA@layers$counts
    # check if there are no dimensions
    if (is.null(dim(raw_counts))) {
      raw_counts <- Matrix(raw_counts, byrow = TRUE, nrow = 1, sparse = TRUE, sparse = T)
    }
    # get the dataframe of feature names
    feature_names_df <- data.frame(donor_rds@assays$RNA@features)
    # get the feature names for this layer
    feature_names <- rownames(feature_names_df[feature_names_df$counts == T, , drop = F])
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
  cat(paste(" Number of genes after filtering:",length(filtered_genes), "\n"))

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
  # set rownames correctly
  rownames(weights) <- colnames(raw_counts)
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
      # check if there are no dimensions
      if (is.null(dim(norm_sparse))) {
         norm_sparse <- Matrix(norm_sparse, byrow = T, nrow = 1, sparse = T)
      }
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
      # check if there are no dimensions
      if (is.null(dim(norm_sparse))) {
         norm_sparse <- Matrix(norm_sparse, byrow = T, nrow = 1, sparse = T)
      }
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
  norm_sparse=norm_sparse[str_order(rownames(norm_sparse)), , drop = F]
  countOutput <- paste0(donor_rds_dir,cohort_id,"/counts/normalized-counts-",donor,"-",cell_type,".mtx")
  writeMM(norm_sparse, file=countOutput)
}

print("Process Complete.")

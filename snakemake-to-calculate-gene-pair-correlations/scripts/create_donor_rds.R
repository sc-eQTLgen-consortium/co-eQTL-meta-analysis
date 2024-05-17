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

print(paste("cell type:",cell_type))
print(paste("cohort id:",cohort_id))
print(paste("seurat file: ",dir_with_seurat,cell_type,".Qced.Normalized.SCs.Rds",sep=''))
print(paste("output directory: ",donor_rds_dir,cohort_id,'/','donor_rds','/',sep=''))
print(paste("Sample mapping file: ",smf,sep=''))
print(paste("Pcs file: ",dir_with_seurat,cell_type,".qtlInput.Pcs.txt",sep=''))

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
Pcs_donors = str_split(Pcs$X,';',simplify=T)[,1]
Pcs_donors = Pcs_donors[Pcs_donors %in% sc_data@meta.data$Assignment]
x=sc_data@meta.data$Assignment %in% Pcs_donors
if (sum(x) != length(sc_data@meta.data$Assignment)){
  cells.use = colnames(sc_data[,x])
  subset_file = subset(sc_data, cells=cells.use)
  sc_data = subset_file
}

expressing_genes <- as.data.frame(rowSums(sc_data@assays$data@data))
expressing_genes = cbind(rownames(sc_data@assays$data@data),expressing_genes)
print("data frame loaded")

colnames(expressing_genes) <- c('gene_names','sum_of_exp')
expressing_genes <- expressing_genes[order(expressing_genes$sum_of_exp, decreasing = T),]

### outputting donor rds
if(genes_to_use=="nan"){
  print("As no gene list provided, top 1000 genes will be calculated after removal of ribosomal and mitochondrial genes.")
  #warning following code removes genes with ensemble ID
  genes.use <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA|^MT|^ENSG|^RPL",rownames(expressing_genes),value=TRUE, invert=TRUE)
  expressing_genes <- expressing_genes[rownames(expressing_genes) %in% genes.use,]
  genes <- expressing_genes[1:1000,]
} else{
  print("Using gene list provided.")
  genes_to_use = read.table(genes_to_use, header=T)[,1]
  genes = expressing_genes[rownames(expressing_genes) %in% genes_to_use,]
}

sc_data_S <- sc_data[rownames(genes),]
  
# Getting a number of counts per cell
raw_counts_mtr <-  as.data.frame(as.matrix(sc_data_S@assays$RNA@counts ))
  
number_of_counts <- as.data.frame(colSums(raw_counts_mtr))
ncells_per_donor <- as.data.frame(table(sc_data$Assignment))

donortab  <- as.data.frame(table(sc_data$Assignment))
donor_list <- donortab[donortab$Freq >10,]$Var1

print(paste("Saving rds file per donor to ",donor_rds_dir,cohort_id,'/','donor_rds',sep=''))

for(donor in donor_list){
  donor_filename = gsub(pattern='_',replacement='',x=donor)
  donor_rds <- sc_data_S[,sc_data_S$Assignment ==donor ]
  saveRDS(donor_rds, paste(donor_rds_dir,cohort_id,'/',"donor_rds",'/',cell_type,'-',donor_filename,'.rds',sep=""))
}

print("done.")


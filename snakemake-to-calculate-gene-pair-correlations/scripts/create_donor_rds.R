#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

cell_type = args[1]
cohort_id = args[2]
dir_with_seurat = args[3]
donor_rds_dir = args[4]

library(Seurat)
library(matrixStats)
library(weights)
library(stringr)

print(paste("cell type:",cell_type))
print(paste("cohort id:",cohort_id))
print(paste("seurate file: ",dir_with_seurat,cell_type,".Qced.Normalized.SCs.Rds",sep=''))
print(paste("output directory: ",donor_rds_dir,cohort_id,'/','donor_rds','/',sep=''))

# selection of genes 

print("Loading big rds file")
sc_data <- readRDS(paste0(dir_with_seurat,paste(cell_type,'.Qced.Normalized.SCs.Rds',sep='')))

mono_expressing_genes <- as.data.frame(rowSums(sc_data@assays$data@data))
colnames(mono_expressing_genes) <- 'sum_of_exp'
mono_expressing_genes$zeros <- rowSums(sc_data@assays$data@data==0)
mono_expressing_genes2 <- as.data.frame(mono_expressing_genes[order(mono_expressing_genes$sum_of_exp, decreasing = T),])

genes.use <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA|^MT|^ENSG|^RPL",rownames(mono_expressing_genes2),value=TRUE, invert=TRUE)

mono_expressing_genes3 <- mono_expressing_genes2[rownames(mono_expressing_genes2) %in% genes.use, ]

### outputting donor rds
print("Selecting top 1000 genes")
genes <- mono_expressing_genes3[1:1000,]
  
sc_data_S <- sc_data[rownames(genes),]
  
# Getting a number of counts per cell
raw_counts_mtr <-  as.data.frame(as.matrix(sc_data_S@assays$RNA@counts ))
  
number_of_counts <- as.data.frame(colSums(raw_counts_mtr))
ncells_per_donor <- as.data.frame(table(sc_data$Assignment))

donortab  <- as.data.frame(table(sc_data$Assignment))
donor_list <- donortab[donortab$Freq >10,]$Var1

print(paste("Saving rds file per donor to ",donor_rds_dir,cohort_id,'/','donor_rds',sep=''))

for(donor in donor_list){
  donor_rds <- sc_data_S[,sc_data_S$Assignment ==donor ]
  saveRDS(donor_rds, paste(donor_rds_dir,cohort_id,'/',"donor_rds",'/',cell_type,'_',donor,'.rds',sep=""))
}

print("done.")


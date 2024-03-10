.libPaths("/usr/local/lib/R/site-library")

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("yaml"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("weights"))

option_list <- list(
  make_option("--celltype", type = "character", help = "Cell type"),
  make_option("--cohort", type = "character", help = "Cohort ID"),
  make_option("--genelist", type = "character", help = "List of genes. If not including a specific set of genes set to: nan"),
  make_option("--n", type = "integer", help = "Number of genes"),
  make_option("--donors", type = "character", help = "List of donors"),
  make_option("--genepath", type = "character", help = "Output gene list filepath"),
  make_option("--donorpath", type = "character", help = "Output donor list filepath"),
  make_option("--input", type = "character", help = "Input rds filepath"),
  make_option("--output", type = "character", help = "Output directory")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

print("Options in effect:")
for (name in names(args)) {
  print(paste0("  --", name, " ", args[[name]]))
}

cat("\n\nLoading Rds file:", args$input)
sc_data <- readRDS(args$input)

cat("\nExtracting expressed genes")
expr_genes <- data.frame(sum_expr = rowSums(sc_data@assays$data@data),gene_name = rownames(sc_data@assays$data@data))
expr_genes <- expr_genes[order(expr_genes$sum_expr, decreasing = T),]

if(args$genelist == "nan"){
  cat("\nNo gene list provided.") 
  cat("\nUsing",args$n,"most expressed genes (exluding ribosomal and mitochondrial genes)")
  #warning following code removes genes with ensemble ID
  pattern <- "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA|^MT|^ENSG|^RPL"
  genes.use <- grep(pattern, expr_genes$gene_name ,value=TRUE, invert=TRUE)
  expr_genes <- expr_genes[expr_genes$gene_name %in% genes.use,]
  genes <- expr_genes[1:args$n,]
} else {
  cat("\nUsing provided gene list")
  gene_list = read.table(args$genelist, header=F)[,1]
  genes <- expr_genes[expr_genes$gene_name %in% gene_list,]
}

sc_data_filtered <- sc_data[rownames(genes),]

rm(sc_data,expr_genes)

cat("\n\nUsing provided donor list")
donor_list <- read.table(args$donors, header=F)$V1
donor_count <- as.data.frame(table(sc_data_filtered$Assignment))
donor_count <- donor_count[donor_count$Var1 %in% donor_list, ]
# Geting alternative id (only if ids contain "_")
colnames(donor_count) <- c("original_ids", "count")
donor_count$alt_ids <- gsub(pattern='_', replacement='', x=donor_count$original_ids)

write.table(donor_count, file = gzfile(args$donorpath), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nDonor data saved to:",args$donorpath)

rm(donor_count)

for(donor in donor_list){
  cat("\nProcessing donor:", donor)
  donor_rds <- sc_data_filtered[,sc_data_filtered$Assignment == donor ]
  cell_barcodes <- colnames(donor_rds)  
  barcodesOut = paste0(args$output,args$cohort,"/donor_barcodes/barcodes-",donor,"-",args$celltype,"-top-",args$n,".tsv.gz")
  write.table(cell_barcodes, gzfile(barcodesOut),sep="\t",row.names = FALSE, col.names=FALSE, quote = FALSE)

  # Gene filtering
  raw_counts <- as.data.frame(as.matrix(donor_rds@assays$RNA@counts))
  n_counts <- as.data.frame(as.matrix(donor_rds@assays$data@data))

  raw_counts_df <- as.data.frame(rowSums(raw_counts != 0))
  raw_counts_df$gene_name <- rownames(raw_counts_df)
  colnames(raw_counts_df) <- c('counts','gene_names')
  filtered_genes <- raw_counts_df$gene_names[raw_counts_df$counts > 10] 
  # save list of genes per donor
  genesOut  <- paste0(args$output,args$cohort,"/donor_gene_list/filtered-genes-",donor,"-",args$celltype,"-top-",args$n,".tsv.gz")
  write.table(filtered_genes, gzfile(genesOut),sep="\t",row.names = FALSE, col.names=FALSE, quote = FALSE)

  # Extract weights per donor
  raw_counts <- raw_counts[rownames(raw_counts) %in% filtered_genes, ]
  weights <- data.frame(weight = colSums(raw_counts != 0))
  weightOutput <- paste0(args$output,args$cohort,"/donor_weight/correlation-weight-",donor,"-",args$celltype,"-top-",args$n,".tsv.gz")
  write.table(weights, gzfile(weightOutput),sep="\t",row.names = TRUE, quote = FALSE)
}
cat("\nCleanup")
rm(donor_rds,cell_barcodes,n_counts,raw_counts,raw_counts_df,filtered_genes,weights)
cat("\nGetting raw counts")
raw_counts <- as.data.frame(as.matrix(sc_data_filtered@assays$RNA@counts))
raw_counts <- raw_counts[rowSums(raw_counts) != 0, ]
cat("\nFilter n_counts")
n_counts <- as.data.frame(as.matrix(sc_data_filtered@assays$data@data))
n_counts <- n_counts[rownames(n_counts) %in% rownames(n_counts), ]

cat("Cleanup")
rm(sc_data_filtered)
cat("\nSave genes")
# Save gene list
write.table(rownames(n_counts), file = gzfile(args$genepath), sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
cat(paste("\nGene list saved to:", args$genepath))
cat("\nSaving counts")
countOutput <- paste0(args$output,args$cohort,"/normalized-counts-",args$celltype,"-top-",args$n,".tsv.gz")
cat("\n\nSaving normalized counts:",countOutput)
write.table(n_counts, gzfile(countOutput),sep="\t",row.names = TRUE, quote = FALSE)

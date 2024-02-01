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

write.table(rownames(genes), file = gzfile(args$genepath), sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
cat(paste("\nGene list saved to:", args$genepath))

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

rm(gene_list,genes,donor_count)

cat("\n\nSaving rds file per donor: n = ", length(donor_list))
for(donor in donor_list){
  donor_filename <- gsub(pattern='_',replacement='',x=donor)
  donor_rds <- sc_data_filtered[,sc_data_filtered$Assignment == donor ]
  donor_output <- paste0(args$output,args$cohort,"/donor_rds/",donor,"-",args$celltype,"-top-",args$n,".rds")
  saveRDS(donor_rds, donor_output)
  rm(donor_rds)
}
.libPaths("/usr/local/lib/R/site-library")

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("yaml"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("weights"))

option_list <- list(
  make_option("--celltype", type = "character", help = "Cell type"),
  make_option("--cohort", type = "character", help = "Cohort ID"),
  make_option("--n", type = "integer", help = "Number of genes"),
  make_option("--donor", type = "character", help = "Donor"),
  make_option("--method", type = "character", help = "Correlation method"),
  make_option("--weight", type = "character", help = "Weight for correlation method. Either 'weighted' or 'unweighted'"),
  make_option("--genelist", type = "character", help = "List of genes"),
  make_option("--input", type = "character", help = "Input donor rds file"),
  make_option("--output", type = "character", help = "Output directory"),
  make_option("--metadata", type = "logical", help = "Include metadata", default = FALSE)
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

print("Options in effect:")
for (name in names(args)) {
  print(paste0("  --", name, " ", args[[name]]))
}

# Define functions
calculate_correlation <- function(n_counts, raw_counts, method, weighted) {
  
  # Check if input is correct
  if (!(method %in% c('pearson', 'spearman'))) {
    stop("Invalid method. Use 'pearson' or 'spearman'.")}

  if (!is.logical(weighted)) {
    stop("Invalid 'weighted' argument. Use TRUE or FALSE.")}

  # Calculate correlation
  if (method == 'spearman'){
    for (i in 1:nrow(n_counts)) { n_counts[i,] <- rank(n_counts[i,]) }}
  
  if (weighted) { weights <- as.data.frame(colSums(raw_counts != 0))[, 1]
  } else { weights <- rep(1, ncol(raw_counts)) }
  
  corMat <- wtd.cor(t(n_counts), weight=weights)
    
  return(corMat)
}

get_cor_stats <- function(corMat) {
  
  # Get correlation, p-value and z-score
  corr <- corMat$correlation[upper.tri(corMat$correlation)]
  
  pval <- corMat$p.value[upper.tri(corMat$p.value)]
  pval[pval < 1.0e-150] <- 1.0e-150
  
  zscore <- abs(qnorm(pval/2))
  zscore[corr < 0] <- zscore[corr < 0] * -1
  
  cor_df <- as.data.frame(cbind(corr, pval, zscore))
  
  # Get gene-pairs
  expanded <- expand.grid(rownames(corMat$correlation), rownames(corMat$correlation))
  
  mat <- matrix(paste(expanded$Var1, expanded$Var2, sep = "_"), 
                nrow = length(rownames(corMat$correlation)), 
                ncol = length(rownames(corMat$correlation)))
  
  gene_p <- mat[upper.tri(mat)]
  
  rownames(cor_df) <- gene_p
  
  return(cor_df)
}

get_paired_zeros <- function(counts, gene_pairs) {
  
  paired_zeros <- c()
  total_zeros <- c()
  
  for (x in gene_pairs) {
    g <- strsplit(x, "_")[[1]]
    gene1 <- counts[g[1],]
    gene2 <- counts[g[2],]
    
    count1 <- which(gene1 %in% 0)
    count2 <- which(gene2 %in% 0)
    
    paired_zeros <- append(paired_zeros, length(intersect(count1,count2)))
    total_zeros <- append(total_zeros, length(count1) + length(count2))
  }
  df <- cbind.data.frame(paired_zeros,total_zeros)
  rownames(df) <- gene_pairs
  return(df)
}

cat("\n\nLoading list of genes")
gene_list <- read.table(gzfile(args$genelist), header = FALSE, sep = "\t")$V1
cat("\nNumber of genes loaded:", length(gene_list))

cat("\n\nLoading RDS file:",args$input)
sc_data <- readRDS(args$input)

cat("\nExtracting counts")  
n_counts <- as.data.frame(as.matrix(sc_data@assays$data@data))
raw_counts <- as.data.frame(as.matrix(sc_data@assays$RNA@counts))
rm(sc_data)

# Fitering on counts > 10 non-zero values
raw_counts_df <- as.data.frame(rowSums(raw_counts != 0))
raw_counts_df$gene_name <- rownames(raw_counts_df)
colnames(raw_counts_df) <- c('counts','gene_names')
filtered_genes <- raw_counts_df$gene_names[raw_counts_df$counts > 10] ## TODO: add variable for number to be filtered?

raw_counts_filtered <- raw_counts[rownames(raw_counts) %in% filtered_genes, ]
n_counts_filtered <- n_counts[rownames(n_counts) %in% filtered_genes, ]

gene_pairs <- combn(gene_list, 2, FUN = function(x) paste(x, collapse = "_"), simplify = TRUE)
cat("\nNumber of unique gene pairs:", length(gene_pairs))

if (args$metadata) {
  cat("Saving metadata")
  paired_zero_counts <- get_paired_zeros(raw_counts,gene_pairs)
  gene_count_df <- as.data.frame(list(counts_sum = rowSums(raw_counts),
                                      counts_mean = rowMeans(raw_counts),
                                      counts_sd = apply(raw_counts, 1, sd),
                                      zero_counts = rowSums(raw_counts == 0)))
  write.table(paired_zero_counts, file = gzfile(paste0(args$output,args$cohort,"/",args$donor,"-gene-pair-metadata.tsv.gz")), sep = "\t", quote = FALSE, row.names = TRUE) 
  write.table(gene_count_df, file = gzfile(paste0(args$output,args$cohort,"/",args$donor,"-gene-count-metadata.tsv.gz")), sep = "\t", quote = FALSE, row.names = TRUE)
}

cat("\n\nCalculating correlation:",args$weight,args$method)
weighted <- ifelse(args$weight == 'weighted', TRUE, FALSE)
corMat <- calculate_correlation(n_counts_filtered, raw_counts_filtered, args$method, weighted)

rm(n_counts,n_counts_filtered,raw_counts,raw_counts_filtered)

cat("\nExtracting correlation for unique gene pairs\n")
cor_df <- get_cor_stats(corMat)
cor_df$x <- rownames(cor_df)
cor_df <- merge(gene_pairs, cor_df, by='x' , all.x=T)
rownames(cor_df) <- cor_df$x
rm(corMat)

metrics = c("corr","pval","zscore")
cat("\n\nSaving correlation for each metric\n\n")
for (x in metrics){
    out_path <- paste0(args$output,args$cohort,"/",args$donor,"-",x,"-",args$celltype,"-top-",args$n,"-",args$method,"-",args$weight,".tsv.gz")
    write.table(cor_df[x], gzfile(out_path),sep="\t",row.names = TRUE, quote = FALSE)
    cat(x, "\nSaved to", out_path,"\n")
  }
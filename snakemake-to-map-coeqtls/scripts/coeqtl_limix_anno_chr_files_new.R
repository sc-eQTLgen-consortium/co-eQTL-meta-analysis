#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Dan Kaptijn, Marc-Jan Bonder, Roy Oelen
# Name: coeqtl_limix_anno_chr_files_new.R
# Function: create feature annotation files for co-eQTL mapping
# Example usage: 
# ~/start_Rscript.sh coeqtl_limix_anno_chr_files_new.R \
#   --cell_type CD4_T_cells \
#   --annotation_prepend /groups/umcg-franke-scrna/tmp04/projects/venema-2022/ongoing/qtl/coeqtl/input/replication_ \
#   --features_test_prepend /groups/umcg-franke-scrna/tmp04/projects/venema-2022/ongoing/qtl/coeqtl/input/replication_features_
#
############################################################################################################################

####################
# libraries        #
####################

library(optparse)

####################
# Functions        #
####################

####################
# Main Code        #
####################

# make command line options
option_list <- list(
  make_option(c("-c", "--cell_type"), type="character", default=NULL,
              help="cell type working on", metavar="character"),
  make_option(c("-a", "--annotation_prepend"), type="character", default=NULL,
              help="base name and location of annotation files", metavar="character"),
  make_option(c("-f", "--features_test_prepend"), type="character", default=NULL,
              help="base name and location of file with features to test", metavar="character")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# read the parameters
ct <- opt[['cell_type']]
annotation_prepend <- opt[['annotation_prepend']]
features_test_prepend <- opt[['features_test_prepend']]

# paste together the full paths, and read the files
annot <- read.delim(paste0(annotation_prepend, ct,"_co.tsv.gz"))
e_annot <- read.delim(paste0(annotation_prepend, ct, "_co2.tsv.gz"))
fvf <-  read.delim(paste0(features_test_prepend, ct, ".tsv.gz"))

# do each chromosome
for (chr in 22:1){
  print(paste("chromsome",chr))
  annotRel <- annot[which(annot$chromosome==chr),]
  e_annotRel <- e_annot[which(e_annot$chromosome==chr),]
  
  annotRel <- rbind(annotRel,e_annotRel)
  
  fvfRel = fvf[which(fvf$feature_id %in% annotRel$feature_id),]
  
  e_annotRel = annotRel[which(duplicated(annotRel$feature_id)),]
  annotRel = annotRel[which(!duplicated(annotRel$feature_id)),]
  
  write.table(annotRel,file = gzfile(paste0(annotation_prepend, "co_",ct,"_chr",chr,".txt.gz")),sep="\t",quote=F,row.names=F)
  write.table(e_annotRel,file = gzfile(paste0(annotation_prepend, "co2_",ct,"_chr",chr,".txt.gz")),sep="\t",quote=F,row.names=F)
  write.table(fvfRel,file = gzfile(paste0(features_test_prepend, "",ct,"_chr",chr,".txt.gz")),sep="\t",quote=F,row.names=F)
}

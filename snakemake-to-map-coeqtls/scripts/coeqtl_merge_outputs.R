#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: coeqtl_merge_outputs.R
# Function: merge nominal coeqtl output files
# Example usage: 
# ~/start_Rscript.sh coeqtl_merge_outputs.R \
#   --input_dir /groups/umcg-franke-scrna/tmp02/projects/venema-2022/ongoing/qtl/coeqtl/output/CD4_T_cells/qtl/ \
#   --output_file /groups/umcg-franke-scrna/tmp02/projects/venema-2022/ongoing/qtl/coeqtl/output/CD4_T_cells/qtl/CD4_T_cells.tsv.gz
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
  make_option(c("-i", "--input_dir"), type="character", default=NULL,
              help="input directory of per-chunk co-eQTL mapping outputs", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default=NULL,
              help="output file", metavar="character")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# get the values
if (is.null(opt[['input_dir']])) {
  error('no input directory supplied')
}
if (is.null(opt[['output_file']])) {
  error('no output file supplied')
}
input_dir <- opt[['input_dir']]
output_file <- opt[['output_file']]

# we'll save the input of all chunks
table_per_chunk <- list()
# list the directories
dirs_input <- list.dirs(input_dir, recursive = F, full.names = F)
# subset to the ones with the a pattern we expect
dirs_chunks <- dirs_input[grepl('\\d+_\\d+_\\d+', dirs_input)]
# check each chunk
for (chunk in dirs_chunks) {
  # get the full path to the file
  chunk_path <- paste(input_dir, '/', chunk, '/', 'qtl_results_all.txt.gz', sep = '')
  # check if the file is not empty
  if (file.size(chunk_path) != 0L) {
    # then read the file if it's not empty
    chunk_contents <- read.table(chunk_path, sep = '\t', header = T)
    # add chunk
    chunk_contents[['chunk']] <- chunk
    # and add to the list
    table_per_chunk[[chunk]] <- chunk_contents
  }
}
# now combine all the chunks
all_qtl <- do.call('rbind', table_per_chunk)
# add MTC
all_qtl[['BH']] <- p.adjust(all_qtl[['p_value']], method = 'BH')
all_qtl[['bonferroni']] <- p.adjust(all_qtl[['p_value']], method = 'bonferroni')
# order by P
all_qtl <- all_qtl[order(all_qtl[['p_value']]), ]

# gz file ends with .gz
if (grepl('.gz$', output_file)) {
  output_file <- gzfile(output_file)
}
# write output
write.table(all_qtl, output_file, row.names = F, col.names = T, quote = F, sep = '\t')

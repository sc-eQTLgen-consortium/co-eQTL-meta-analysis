"""
This script is creating chunking files for the co-eQTL mapping pipeline

authors: Roy Oelen

example usage:

python coeqtl_make_limix_chunking_file.py \
    --qtl_loc /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/sceqtlgen_coeqtl_map/eQTLs_finemapped_20240626/CD4_T.Ds.wg3_Ye_wg3_wijst2018_wg3_Trynka_wg3_sawcer_wg3_oneK1K_wg3_okada_wg3_Nawijn_wg3_Li_wg3_Franke_split_v3_wg3_Franke_split_v2_wg3_multiome_UT_wg3_idaghdour.5.csTop_qtl_results.txt \
    --out_loc /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/sceqtlgen_coeqtl_map/input/cd4_chunking_all.txt.gz

"""

#############
# libraries #
#############
import pandas as pd
import numpy as np
import os
import time
from pprint import pp
import argparse

#############
# functions #
#############



#############
# main code #
#############

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-q', '--qtl_loc', type = str, help = 'location of previous eQTL summary stats')
parser.add_argument('-o', '--out_loc', type = str, help = 'location of original feature annotations')
parser.add_argument('-c', '--chromosome_column', type = str, help = 'column in the QTL output that has the chromosome of the features', default = 'feature_chromosome')
parser.add_argument('-s', '--start_column', type = str, help = 'column in the QTL output that has the start position of the features', default = 'feature_start')
parser.add_argument('-e', '--end_column', type = str, help = 'column in the QTL output that has the end position of the features', default = 'feature_end')
parser.add_argument('-p', '--significance_column', type = str, help = 'column in the QTL output that has the significance to filter the features on (leave empty to not filter)', default = None)
parser.add_argument('-v', '--significance_cutoff', type = float, help = 'significance cutoff to filter QTL output on (leave empty to not filter)', default = None)
args = parser.parse_args()


####################
# to-test features #
####################

# we need to make sure some of these columns are strings
string_columns = [args.chromosome_column, args.start_column, args.end_column]
# make that into the dtypes dict
dtypes = {x : 'str'  for x in string_columns}
# get the eQTL summary stats
eqtls = pd.read_csv(args.qtl_loc, sep = '\t', dtype = dtypes)
# subset if those parameters were supplied
if args.significance_column is not None and args.significance_cutoff is not None:
    eqtls = eqtls[eqtls[args.significance_column] < args.significance_cutoff]

# subset to the columns we care about
eqtl_features = eqtls[args.chromosome_column] + ':' + eqtls[args.start_column] + '-' + eqtls[args.end_column]
# only keep unique entries
eqtl_features.drop_duplicates(inplace = True)
# check if the output file ends in gz
if args.out_loc.endswith('.gz'):
    eqtl_features.to_csv(args.out_loc, index = None, header = False, sep = '\t', compression = 'gzip')
else:
    eqtl_features.to_csv(args.out_loc, index = None, header = False, sep = '\t')

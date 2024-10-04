"""
This script is for opening each bgen file, so an index is created

authors: Roy Oelen

example usage:

python coeqtl_open_bgens.py \
    --path /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/genotype_input/ \
    --regex EUR_imputed_hg38_varFiltered_chr*.bgen

"""

#############
# libraries #
#############

import glob
import argparse
from bgen_reader import read_bgen


#############
# functions #
#############

def open_and_close_bgen(full_bgen_path):
    # say we'll start to open the file
    print(''.join(['starting to open ', full_bgen_path]))
    # read the file
    bgen_handle = read_bgen(full_bgen_path, verbose = True)
    # and close the file
    bgen_handle.bgen_file_close()
    # and say we closed it
    print(''.join(['closed ', full_bgen_path]))


#############
# main code #
#############

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-l', '--path', type = str, help = 'location of bgen files (string)')
parser.add_argument('-r', '--regex', type = str, help = 'regex of bgen files to open (string)')
args = parser.parse_args()

# get the arguments
path_bgens = args.path
regex_bgens = args.regex

# if there was no path supplied, we'll assume it is the current directory
if args.path is None:
    path_bgens = './'

# if there was no regex, we'll assume it was everything ending in .bgen
if args.regex is None:
    regex_bgens = '*.bgen'

# paste together what to match
path_to_match = '/'.join([path_bgens, regex_bgens])

# now actually use glob to get the files
bgen_files = glob.glob(path_to_match)

# say what we'll open
print('will try to open these files:')
print(','.join(bgen_files))

# open and close each file
for f in bgen_files:
    open_and_close_bgen(f)

# say we opened them all
print('opened and closed each file')
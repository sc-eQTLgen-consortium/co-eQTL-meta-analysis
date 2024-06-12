#!/usr/bin/env python3 

import argparse
import csv
from PascalX import genescorer

parser = argparse.ArgumentParser(description='Arguments to Pascal scoring.')
parser.add_argument("-r", "--refpanel", help = "Reference panel path")
parser.add_argument("-g", "--gwas", help = "GWAS summary statistics path")
parser.add_argument("-a", "--annotation", help = "Genome annotation file")
parser.add_argument("-t", "--threads", help = "Number of threads for scoring")
parser.add_argument("-o", "--outfile", help = "Output file", nargs = '?', const=".")
parser.add_argument("-s", "--rscol", help = "rs id column")
parser.add_argument("-p", "--pcol", help = "p-value column")
parser.add_argument("-w", "--window", help = "gene window")

args = parser.parse_args()

Scorer = genescorer.chi2sum(window=int(args.window))
Scorer.load_refpanel(args.refpanel, keepfile=None, parallel=int(args.threads), qualityT = None)
#Scorer.load_refpanel("/groups/umcg-fg/tmp01/projects/downstreamer/downstreamer_bundle/reference_datasets/1000G_EUR_noFIN/PascalX_get1KGGRCh37/EUR.1KG.GRCh37", keepfile=None, parallel=4)
print("Reference panel loaded.")

Scorer.load_genome(args.annotation,ccol=1,cid=0,csymb=0,cstx=2,cetx=3,chrStart=0,NAgeneid='n/a',header=True) # csymb = cid to get ensmbl
#Scorer.load_genome("/groups/umcg-fg/tmp01/projects/downstreamer/depict2_bundle/reference_datasets/human_b37/ensgR75_protein_coding.txt",ccol=1,cid=0,csymb=0,cstx=2,cetx=3,chrStart=0,NAgeneid='n/a',header=True) # csymb = cid to get ensmbl
print("Genome annotation loaded.")
print(args.annotation)

print(args.gwas)
print(int(args.rscol))
print(int(args.pcol))

Scorer.load_GWAS(args.gwas,rscol=int(args.rscol),pcol=int(args.pcol), header=True)#, delimiter = "\t"
#Scorer.load_GWAS("test_noempty_nonas.txt",rscol=0,pcol=1, header=True)
print("GWAS dataset loaded.")

print("Starting scoring...")
R = Scorer.score_all(parallel=int(args.threads), method="saddle",autorescore=False, nobar=True)
#R = Scorer.rescore(R,method="ruben", parallel=int(args.threads))
#R = Scorer.score_chr(chrs=[21], method="auto",autorescore=True)

#Scorer.save_scores(args.output)

print("outfile")
print(args.outfile)

with open(args.outfile,'w') as out:
    csv_out=csv.writer(out, delimiter='\t')
    csv_out.writerow(['gene','pvalue','nsnps', 'min_pvalue'])
    for row in R[0]:
        csv_out.writerow(row)

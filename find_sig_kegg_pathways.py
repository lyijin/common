#!/usr/bin/env python3

"""
> find_sig_kegg_pathways.py <

Based a background of genes <--> pathways, when given a list of genes of
interest, perform a Fisher's exact test + correction to identify pathways
that are enriched/depleted in genes from the list.
"""
import argparse
import csv
import os

import calc_fishers_exact
import correct_p_values

parser = argparse.ArgumentParser(description="""
Based a background of genes <--> pathways, when given a list of genes of
interest, perform a Fisher's exact test + correction to identify pathways
that are enriched/depleted in genes from the list.""")
parser.add_argument('genes_of_interest', metavar='txt_file',
                    type=argparse.FileType('r'), 
                    help='txt file with genes of interest, one per line.')
parser.add_argument('--universe', metavar='txt_file',
                    type=argparse.FileType('r'), 
                    help='txt file of (smaller) universe, one per line.')
group = parser.add_mutually_exclusive_group()
group.add_argument('--spis', action='store_true',
                    help='Background list: S. pistillata.')
group.add_argument('--smic', action='store_true',
                    help='Background list: S. microadriaticum.')
group.add_argument('--aiptasia', action='store_true',
                    help='Background list: Aiptasia.')

args = parser.parse_args()

# parse universe (e.g. subset of genes that are methylated)
gene_universe = []
if args.universe:
    tsv_reader = csv.reader(args.universe, delimiter='\t')
    for row in tsv_reader:
        if not row: continue
        
        # assume gene name is in the first column
        gene_universe.append(row[0])
        
    gene_universe = list(set(gene_universe))

# parse background list
if args.spis:
    background_list = os.path.dirname(os.path.abspath(__file__)) \
                      + '/kegg_backgrounds/spis.kegg.tsv'
elif args.smic:
    background_list = os.path.dirname(os.path.abspath(__file__)) \
                      + '/kegg_backgrounds/smic.kegg.tsv'
elif args.aiptasia:
    background_list = os.path.dirname(os.path.abspath(__file__)) \
                      + '/kegg_backgrounds/aiptasia.kegg.tsv'

# ko_background: a dict to store genes that belong to a specific pathway
#   ko_background['koXXXXX'] = [geneA, geneD, geneG, ...]
ko_background = {}
all_background_genes = []
tsv_reader = csv.reader(open(background_list), delimiter='\t')
for row in tsv_reader:
    if not row: continue
    
    gene = row[0]
    
    # ignore genes not in the universe, if it exists
    if args.universe:
        if gene not in gene_universe:
            continue
    
    gene_kos = row[1].split(';')
    all_background_genes.append(gene)
    
    for g in gene_kos:
        if g not in ko_background:
            ko_background[g] = []
            
        ko_background[g].append(gene)

all_background_genes = list(set(all_background_genes))

# parse genes of interest
genes_of_interest = []
tsv_reader = csv.reader(args.genes_of_interest, delimiter='\t')
for row in tsv_reader:
    if not row: continue
    
    # assume gene name is in the first column
    genes_of_interest.append(row[0])

if args.universe:
    genes_of_interest = list(set(genes_of_interest) & set(gene_universe))
else:
    genes_of_interest = list(set(genes_of_interest))

# perform p value calculations
p_values = {}
for ko in ko_background:
    k = len(set(genes_of_interest) & set(ko_background[ko]))
    m = len(ko_background[ko])
    n = len(set(genes_of_interest) & set(all_background_genes))
    N = len(all_background_genes)
    direction = 'Enriched' if k/n > m/N else 'Depleted'
    
    annot = '\t'.join([str(x) for x in [ko, k, m, n, N, direction]])
    p_values[annot] = calc_fishers_exact.fishers_exact_cdf(k, m, n, N)
    
p_values = correct_p_values.correct_p_values(p_values, 'BH')

for p in p_values:
    print (p, p_values[p], sep='\t')
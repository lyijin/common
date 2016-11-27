#!/usr/bin/env python3

"""
> map_genes_to_scaffolds.py <

Script parses gff file, noting which genes map back to which scaffold. When
given a list of genes, return a table containing the genes and the scaffold
names.

Note that --include operates on a FILE, and its use is compulsory.
"""
import argparse

parser = argparse.ArgumentParser(description="""
Script parses gff file, noting which genes map back to which scaffold. When
given a list of genes, return a table containing the genes and the scaffold
names.""")
parser.add_argument('gff_file', metavar='gff_file',
                    type=argparse.FileType('r'), 
                    help='GFF3 file to define gene/scaffold relationship')
parser.add_argument('--include', metavar='gene_file', required=True,
                    type=argparse.FileType('r'), nargs='+',
                    help='get the scaffolds for genes in the given files')
parser.add_argument('--space', metavar='n', type=int, 
                    help='''ugly hack: split gene name by space,
                    then pick n-th piece as gene name''')
                    
args = parser.parse_args()

# start of script - get gff data
import parse_gff
gff_data = parse_gff.parse_gff(args.gff_file, 'gene')

# create dictionary that maps genes back to scaffold
scaffold_info = {}
for g in gff_data:
    for h in gff_data[g]['gene']:
        scaffold_info[h] = g

# read in gene lists
for i in args.include:
    for gene in i:
        gene = gene.strip()
        
        if args.space:
            gene = gene.split(' ')[args.space]
        
        print (gene, scaffold_info[gene], sep='\t')

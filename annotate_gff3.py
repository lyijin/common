#!/usr/bin/env python3

"""
> annotate_gff3.py <

Script expands the 'Name=' field in the gff3 field to include putative function
for the gene.

This is mostly helpful for genome browsers, as users can skip the step
of noting down interesting gene numbers --> check what their functions are.
"""
import argparse
import csv
import re

parser = argparse.ArgumentParser(description="""
Script expands the 'Name=' field in the gff3 field to include putative function
for the gene.""")

parser.add_argument('genome_gff3', metavar='gff3_file',
                    type=argparse.FileType('r'),
                    help='genome GFF3 annotation.')
parser.add_argument('annot_file', metavar="cov_filename",
                    type=argparse.FileType('r'),
                    help='file in the form *_tabulated_expanded_gene_name.tsv.')
parser.add_argument('--regex', '-r', type=str,
                    help='''hack to resolve nomenclature differences between 
                            the two provided files. The same regex is applied
                            to both identifiers to check for equality. For
                            instance, to match Pdae1 and PdaeGene1, use the 
                            regex "\d+". Invalid regexes would result in 
                            horrible, uncaught exceptions.''')

args = parser.parse_args()

# create a dictionary to store the expanded gene name
expanded_gene_name = {}

tsv_reader = csv.reader(args.annot_file, delimiter='\t')
next(tsv_reader)        # skip header
for row in tsv_reader:
    if not row: continue
    if len(row) < 3: continue
    
    gene = row[0]
    annot = row[3]
    
    if args.regex:
        g = re.search('({})'.format(args.regex), gene).group(1)
        expanded_gene_name[g] = ' // ' + annot
    else:
        expanded_gene_name[gene] = ' // ' + annot

# read the gff3 file
tsv_reader = csv.reader(args.genome_gff3, delimiter='\t')
for row in tsv_reader:
    # target gff3 lines for 'gene', specifically; otherwise, print original line
    if len(row) < 3:
        print (*row, sep='\t')
        continue
    
    if row[2] != 'gene':
        print (*row, sep='\t')
        continue
    
    gff_annots = row[8]
    gene_name = re.search('Name=(.*?)(;|$)', gff_annots).group(1)
    
    if args.regex:
        g = re.search('({})'.format(args.regex), gene_name).group(1)
        if g not in expanded_gene_name:
            print (*row, sep='\t')
            continue
        expansion = expanded_gene_name[g]
    else:
        if gene_name not in expanded_gene_name:
            print (*row, sep='\t')
            continue
        expansion = expanded_gene_name[gene_name]
    
    gff_annots = gff_annots.replace('Name=' + gene_name, 
                                    'Name=' + gene_name + expansion)
    
    print (*row[:8], gff_annots, sep='\t')

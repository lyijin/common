#!/usr/bin/env python3

"""
> map_genes_to_full_annots.py <

Script parses file (typically called *_tabulated_annots.tsv) that assigns 
protein annotation to genes. When given a list of genes, return a table 
containing the genes and the full annotations.

Script can _sort_of_ handle partial gene annotations as well. If a gene
annotation contains separators ("_", "|" or " "), a synonyms dictionary is
created to handle potential partial matches. This makes it useful in handling
either word in "Smic1234 SmicGene1234", but not for "adi_v1.00002" (as the 
left word "adi_" is not unique).

--full option prints EVERY THING, including blast results and GO annotations.
By default, only annotations are printed out (1st to 4th column of the file).

--simple only prints the gene name, with a header line on top. Mutually 
exclusive with the --full flag.

--noheader suppresses the printing of the header.
"""
import argparse
import csv
import re

parser = argparse.ArgumentParser(description="""
Script parses file (typically called *_tabulated_annots.tsv) that assigns 
protein annotation to genes. When given a list of genes, return a table 
containing the genes and the full annotations.""")
parser.add_argument('annot_file', metavar='tsv_file',
                    type=argparse.FileType('r'), 
                    help='tsv file containing gene annotations.')
parser.add_argument('include', metavar='gene_file',
                    type=argparse.FileType('r'), nargs='+',
                    help='genes that need annotating.')
parser.add_argument('--noheader', action='store_true',
                    help='Skips printing of header.')
group = parser.add_mutually_exclusive_group()
group.add_argument('--full', action='store_true',
                    help='print all details associated with annotation.')
group.add_argument('--simple', action='store_true',
                    help='only print gene names.')

args = parser.parse_args()

def trim_line(line, flag):
    # depending on which flag is set, print stuff out of different lengths
    if flag == 'full':
        return '\t'.join(line)
    elif flag == 'simple':
        return '\t'.join(line[:1])
    else:
        return '\t'.join(line[:4])

# start of script - read tabulated annot into memory
annot_data = {}
synonyms = {}       # values in synonyms are keys in annot_data

tsv_reader = csv.reader(args.annot_file, delimiter='\t')
header = next(tsv_reader)
for row in tsv_reader:
    if not row: continue
    
    annot_data[row[0]] = row
    if '_' or ' ' in row[0]:
        partial_annots = re.split(' |_|\|', row[0])
        for p in partial_annots:
            synonyms[p] = row[0]

# read in gene lists and print full annots
if args.full:
    flag = 'full'
elif args.simple:
    flag = 'simple'
else:
    flag = 'default'

if not args.noheader:
    print (trim_line(header, flag))

for i in args.include:
    for gene in i:
        gene = gene.strip()
        
        if gene in annot_data:
            print (trim_line(annot_data[gene], flag))
        elif synonyms:
            if gene in synonyms:
                print (trim_line(annot_data[synonyms[gene]], flag))
            else:
                print (gene)
        else:
            print (gene)

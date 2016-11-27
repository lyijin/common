#!/usr/bin/env python3

"""
> map_go_to_genes.py <

Script parses file (typically called *_tabulated_annots.tsv) that assigns 
GO terms to genes. When given a file containing GO terms somewhere along the 
line, return the same file with gene names after the GO terms.

The returned GO terms have no concept of significance!
"""
import argparse
import csv
import re

parser = argparse.ArgumentParser(description="""
Script parses file (typically called *_tabulated_annots.tsv) that assigns 
GO terms to genes. When given a file containing GO terms somewhere along the 
line, return the same file with gene names after the GO terms.""")
parser.add_argument('annot_file', metavar='annot_file',
                    type=argparse.FileType('r'), 
                    help='tsv file containing gene annotations.')
parser.add_argument('go_file', metavar='text_file',
                    type=argparse.FileType('r'),
                    help='file with GO terms that need annotating.')
parser.add_argument('--shorten', '-s', metavar='n', type=int, default=0,
                    help='shorten gene names to n-th word of name.')

args = parser.parse_args()

# store GO:gene relationship in a dict
#   go_data['GO:XXXXXXX'] = ['Gene1', 'Gene2', ...]
go_data = {}

tsv_reader = csv.reader(args.annot_file, delimiter='\t')
for row in tsv_reader:
    if not row: continue
    if len(row) < 16: continue
    if row[15][:3] != 'GO:' : continue
    
    go_terms = row[15].split(',')
    for g in go_terms:
        if g not in go_data:
            go_data[g] = []
        
        if args.shorten:
            gene_name = row[0].split(' ')[args.shorten]
        else:
            gene_name = row[0]
        
        go_data[g].append(gene_name)

with args.go_file as f:
    for line in f:
        go_term = re.search('(GO:\d{7})', line)
        line = line.rstrip()
        if go_term:
            # if multiple GO terms are present, take first one
            go_term = go_term.group(1)
            if go_term in go_data:
                print (line, ','.join(go_data[go_term]), sep='\t')
            else:
                print (line, 'no_data' , sep='\t')
        else:
            print (line)

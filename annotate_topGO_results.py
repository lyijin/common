#!/usr/bin/env python3

"""
> annotate_topGO_results.py <

For each line in the topGO results, this script adds in the genes that contain
the GO term as an additional column in the results file.

Script _can_ accept multiple topGO results files, and the mandatory --annot 
switch takes in the annotation file that assigns GO terms to genes.
"""
import argparse
import csv
import os.path
import re

def natural_sort(input_list):
    tryint = lambda x: int(x) if x.isdigit() else x
    chunked_text = lambda x: [tryint(y) for y in re.split('([0-9]+)', x)]
    sorted_list = sorted(input_list, key=chunked_text)

    return sorted_list

parser = argparse.ArgumentParser(description="""
For each line in the topGO results, this script adds in the genes that contain
the GO term as an additional column in the results file.""")
parser.add_argument('topgo_files', metavar='tsv_files',
                    type=argparse.FileType('r'), nargs='+',
                    help='TAB-separated topGO output files.')
parser.add_argument('--annot', metavar='annotation_tsv_file',
                    type=argparse.FileType('r'), required=True,
                    help='file containing GO terms assigned to genes.')

args = parser.parse_args()

# read the annot file, and create dict of GO ID to genes
go_to_genes = {}
tsv_reader = csv.reader(args.annot, delimiter='\t')
for row in tsv_reader:
    if not row: continue
    if len(row) < 2: continue
    if row[1][:3] != 'GO:': continue
    
    for x in row[1].split(','):
        if x not in go_to_genes:
            go_to_genes[x] = []
        
        go_to_genes[x].append(row[0])

# for each topGO results file, add in the genes containing the GO term
for tg in args.topgo_files:
    tsv_reader = csv.reader(tg, delimiter='\t')
    output_filename = '{}.annot{}'.format(os.path.splitext(tg.name)[0],
                                              os.path.splitext(tg.name)[1])
    
    with open(output_filename, 'w') as of:
        for row in tsv_reader:
            no_go = False
            if not row: no_go = True
            if len(row) < 6: no_go = True
            if row[1][:3] != 'GO:': no_go = True
            
            if no_go:
                print ('\t'.join(row), file=of)
                continue
            
            genelist = ','.join(natural_sort(go_to_genes[row[1]]))
            print ('\t'.join(row + [genelist]), file=of)

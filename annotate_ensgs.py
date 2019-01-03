#!/usr/bin/env python3

"""
> annotate_ensgs.py <

Given a text input, checks whether there are human genes in the form of
ENSGxxxxxxxxxxx (ENSEMBL ID)--if there are, add the symbol and full name of
the gene, then print to standard output. If no input filename is provided,
stdin is expected.

Why not use BiomaRt? It's a massive hassle to pipe output to/from R, and that
library requires internet access every time it's run. On the other hand, this
script uses an offline data source for gene names that can be wget/curl from
  ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt
"""
import argparse
import csv
import gzip
import os.path
import re
import sys

# read the annot file, and create dict that contains details of all ENSGs
ensgs_dict = {}

hgnc_file = os.path.dirname(__file__) + '/hgnc_complete_set.txt.gz'
tsv_reader = csv.reader(gzip.open(hgnc_file, 'rt'), delimiter='\t')
for row in tsv_reader:
    if not row: continue
    if len(row) < 20: continue
    if 'ENSG' not in row[19]: continue
    
    symbol = row[1]
    name = row[2]
    ensg = row[19]
    
    ensgs_dict[ensg] = {'symbol': symbol, 'name': name}

def get_annot(ensg, warnings=False):
    if ensg in ensgs_dict:
        return ensgs_dict[ensg]['symbol'], ensgs_dict[ensg]['name']
    else:
        if warnings:
            print (f'WARNING: {ensg} does not have a gene symbol/name assigned.',
                   file=sys.stderr)
        return '', ''

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    Given a text input, checks whether there are human genes in the form of
    ENSGxxxxxxxxxxx (ENSEMBL ID)--if there are, add the symbol and full name of
    the gene, then print to standard output. If no input filename is provided,
    stdin is expected.""")
    parser.add_argument('text_file', metavar='text_file', 
                        type=argparse.FileType('r'), nargs='?', default=sys.stdin,
                        help='text file containing ENSGs.')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='prints diagnostic stuff to stderr.')
    args = parser.parse_args()

    # for each topGO results file, add in the genes containing the GO term
    ae_warning = False
    with args.text_file as f:
        for line in f:
            try:
                ensg = re.search('ENSG\d+', line).group(0)
                print (line.strip(), *get_annot(ensg, warnings=args.verbose),
                       sep='\t')
            
            except AttributeError:
                # line has no 'ENSGxxxx' in it
                if not ae_warning and args.verbose:
                    print ('WARNING: there is at least a line with no ENSG in it.',
                           file=sys.stderr)
                    ae_warning = True
                print (line.strip())

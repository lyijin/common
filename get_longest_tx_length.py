#!/usr/bin/env python3

"""
> get_longest_tx_length.py <

Script takes in cDNA/ncRNA FASTA file (multiple; can be compressed) from
Ensembl, then gets the longest transcript length for each gene. Outputs a 
two-column table with gene and max length.

Transcript lengths are useful for conversion of CPM (count-per-million) values
to TPM (transcript-per-million).

cDNA/ncRNA files can be downloaded from
https://www.ensembl.org/info/data/ftp/index.html

While cDNAs are slightly shorter than full transcripts, it's a good-enough
approximation.
"""
import argparse
from pathlib import Path
import re

import parse_fasta

parser = argparse.ArgumentParser(description="""
Script takes in cDNA/ncRNA FASTA file (multiple; can be compressed) from
Ensembl, then gets the longest transcript length for each gene. Outputs a 
two-column table with gene and max length.
""")
parser.add_argument('fasta_file', metavar='fasta_file', nargs='+',
                    type=Path, help='cDNA file downloaded from Ensembl.')

args = parser.parse_args()

# define dict to store longest tx lengths
longest_tx = {}

# iterate through the fasta files provided, store longest tx length
for fa in args.fasta_file:
    sequences = parse_fasta.get_all_sequences(fa, 'fasta')
    for seq in sequences:
        gene = re.search(r'gene:(\w+\d+)\.\d+ ', seq)[1]
        
        if gene not in longest_tx:
            longest_tx[gene] = 0
        
        if longest_tx[gene] < len(sequences[seq]):
            longest_tx[gene] = len(sequences[seq])

# print results
print ('gene', 'max_tx_length', sep='\t')
for g in sorted(longest_tx):
    print (g, longest_tx[g], sep='\t')

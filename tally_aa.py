#!/usr/bin/env python3

docstring = """
> tally_aa.py <

Script to tally the proportion of amino acids in a FASTA/FASTQ file.

X = any degenerate aa
* = stop codon
"""
import argparse
import collections
from pathlib import Path

import parse_fasta

parser = argparse.ArgumentParser(
    description=docstring, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('fasta_files', metavar='fasta_files',
                    type=Path, nargs='+',
                    help='FASTA files that will be tallied.')
parser.add_argument('--fastq', action='store_true', default=False,
                    help='Tallied files are FASTQ, not FASTA.')

args = parser.parse_args()

valid_aa = 'ACDEFGHIKLMNPQRSTVWY*X'

# header row
print ('File', *valid_aa, sep='\t')

for f in args.fasta_files:
    base_composition = collections.Counter()
    seqs = parse_fasta.get_all_sequences(f, 'fastq' if args.fastq else 'fasta',
                                         sequences_only=True)
    
    for s in seqs:
        base_composition += collections.Counter(s.upper())
    
    # print results out
    print (f.name, *[base_composition[x] for x in valid_aa], sep='\t')

#!/usr/bin/env python3

"""
> tally_dinuc.py <

Script to tally the proportion of dinucleotides in a FASTA/FASTQ file.
Supports opening gzipped files with the appropriate suffix (.gz).

Valid dinucleotides are [ACGT]/[ACGT]. Script will also tally GC% in valid
dinucleotides (i.e. excluding bases that are before or after degenerate bases).
"""
import argparse
import collections
from pathlib import Path

import parse_fasta

parser = argparse.ArgumentParser(description="""
Script to tally the proportion of dinucleotides in a FASTA/FASTQ file.
Supports opening gzipped files with the appropriate suffix (.gz).""")

parser.add_argument(
    'fasta_files', metavar='fasta_files', type=Path, nargs='+',
    help='FASTA files that will be tallied.')
parser.add_argument(
    '--fastq', action='store_true', default=False,
    help='Tallied files are FASTQ, not FASTA.')

args = parser.parse_args()

for f in args.fasta_files:
    # print headers
    print (f.name)
    print ('dinuc', 'observed', 'expected', sep='\t')
    
    base_composition = collections.Counter()
    seqs = parse_fasta.get_all_sequences(f, 'fastq' if args.fastq else 'fasta',
                                         sequences_only=True)
    
    for s in seqs:
        s = s.upper()
        dinucs = [''.join(x) for x in zip(s[::2], s[1::2])] + \
                 [''.join(y) for y in zip(s[1::2], s[2::2])]
        base_composition += collections.Counter(dinucs)
    
    # tally stuff
    all_dinucs = 0
    mononucs = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for i in ['A', 'C', 'G', 'T']:
        for j in ['A', 'C', 'G', 'T']:
            # get sum of all dinucleotides
            all_dinucs += base_composition[i + j]
            
            # tally mono-nucleotides to calculate GC%
            mononucs[i] += base_composition[i + j]
            mononucs[j] += base_composition[i + j]
    
    # note: all_mononucs is approx 2 x all_dinucs, because each base is 
    # doubly-counted (see code 3 and 4 lines above)
    all_mononucs = sum(mononucs.values())
    
    # print results
    for i in ['A', 'C', 'G', 'T']:
        for j in ['A', 'C', 'G', 'T']:
            exp = mononucs[i] * mononucs[j] / (all_mononucs ** 2) * all_dinucs
            print (i + j, base_composition[i + j], int(exp), sep='\t')
    
    print ('sum', all_dinucs, sep='\t')
    gc_pct = (mononucs['C'] + mononucs['G']) / sum(mononucs.values()) * 100
    print ('GC%', round(gc_pct, 2), sep='\t')
    print ()

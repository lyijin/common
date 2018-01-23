#!/usr/bin/env python3

"""
> tally_dinuc.py <

Script to tally the proportion of dinucleotides in a FASTA/FASTQ file.
Supports opening gzipped files with the --gzip / -g flag.

Valid dinucleotides are [ACGT]/[ACGT]. Script will also tally GC% in valid
dinucleotides (i.e. excluding bases that are before or after degenerate bases).
"""
import argparse
import collections

import parse_fasta

parser = argparse.ArgumentParser(description="""
Script to tally the proportion of dinucleotides in a FASTA/FASTQ file.
Supports opening gzipped files with the --gzip / -g flag.""")

parser.add_argument('fasta_files', metavar='fasta_files',
                    type=argparse.FileType('r'), nargs='+',
                    help='FASTA files that will be tallied.')
parser.add_argument('--fastq', action='store_true', default=False,
                    help='Tallied files are FASTQ, not FASTA.')
parser.add_argument('--gzip', '-g', action='store_true', default=False,
                    help='Tallied files are gzip-compressed.')

args = parser.parse_args()

for f in args.fasta_files:
    # print headers
    print (f.name)
    print ('dinuc', 'observed', 'expected', sep='\t')
    
    base_composition = collections.Counter()
    seqs = parse_fasta.get_all_sequences(f, 'fastq' if args.fastq else 'fasta',
                                         gzip_compressed=args.gzip,
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
    cg_pct = (mononucs['C'] + mononucs['G']) / sum(mononucs.values()) * 100
    print ('CG%', round(cg_pct, 2), sep='\t')
    print ()

#!/usr/bin/env python3

"""
> tally_ACGTN.py <

Script to tally the proportion of A/C/G/T/N bases in a FASTA/FASTQ file.
Supports opening gzipped files with the --gzip / -g flag.

N = any degenerate base (R/Y/...), not just "N".
"""
import argparse
import collections

import parse_fasta

parser = argparse.ArgumentParser(description="""
Script to tally the proportion of A/C/G/T/N bases in a FASTA/FASTQ file.
Supports opening gzipped files with the --gzip / -g flag.""")

parser.add_argument('fasta_files', metavar='fasta_files',
                    type=argparse.FileType('r'), nargs='+',
                    help='FASTA files that will be tallied.')
parser.add_argument('--fastq', action='store_true', default=False,
                    help='Tallied files are FASTQ, not FASTA.')
parser.add_argument('--gzip', '-g', action='store_true', default=False,
                    help='Tallied files are gzip-compressed.')

args = parser.parse_args()

# header row
print ('File', 'A', 'C', 'G', 'T', 'N', 'ACGT', 'ACGTN', 'GC%', sep='\t')

for f in args.fasta_files:
    base_composition = collections.Counter()
    seqs = parse_fasta.get_all_sequences(f, 'fastq' if args.fastq else 'fasta',
                                         gzip_compressed=args.gzip,
                                         sequences_only=True)
    
    for s in seqs:
        base_composition += collections.Counter(s.upper())
    
    a = base_composition['A']
    c = base_composition['C']
    g = base_composition['G']
    t = base_composition['T']
    acgt = a + c + g + t
    acgtn = sum(base_composition.values())
    non_acgt = acgtn - acgt

    gc_pct = round((c + g) / acgt * 100, 3)
    
    # print results out
    print (f.name, a, c, g, t, non_acgt, acgt, acgtn, gc_pct, sep='\t')
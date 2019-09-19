#!/usr/bin/env python3

"""
> tally_per-base_stats.py <

Script to tally per-base proportions of A/C/G/T/N bases in a file containing
a sequence every line. Produces a table specific for every detected sequence
length.

N = any degenerate base (R/Y/...), not just "N".

Assumes standard input by default, but can be overridden by a filename.
"""
import argparse
import collections
import sys

parser = argparse.ArgumentParser(description="""
Script to tally per-base proportions of A/C/G/T/N bases in a file containing
a sequence every line. Produces a table specific for every detected sequence
length.""")

parser.add_argument('input_file', metavar='input_file',
                    type=argparse.FileType('r'), nargs='?',
                    default=sys.stdin,
                    help='File containing sequences, one per line.')
parser.add_argument('-n', action='store_true',
                    help='print whole numbers instead of proportions.')

args = parser.parse_args()

# read file into dict that is stratified by sequence lengths
sequences = {}
for line in args.input_file:
    seq = line.strip()
    if len(seq) not in sequences:
        sequences[len(seq)] = []
    
    sequences[len(seq)].append(seq)

# print per-base values
for seqlen in sorted(sequences):
    # compile a dict that contains base compositions, which makes it easier to
    # print stuff out in a tabular fashion later
    bc = {}
    for n in range(seqlen):
        temp = collections.Counter([x[n] for x in sequences[seqlen]])
        bc[n] = {'A': temp['A'], 'C': temp['C'],
                 'G': temp['G'], 'T': temp['T']}
        # define 'N' as "non-ACGT"
        bc[n]['N'] = sum(temp.values()) - sum(bc[n].values())
    
    # print table out
    print (f'## Sequence length = {seqlen} ##')
    print ('', *[x + 1 for x in range(seqlen)], sep='\t')
    for base in 'ACGTN':
        if args.n:
            # print raw frequencies
            print (base, *[bc[x][base] for x in range(seqlen)],
                   sep='\t')
        else:
            # print proportions
            print (base,
                   *[round(bc[x][base] / sum(bc[x].values()), 4) for x in range(seqlen)],
                   sep='\t')
    print ()

#!/usr/bin/env python3

docstring = """
> tally_per-aa_stats.py <

Script to tally per-amino acid proportions in a file containing a peptide 
sequence every line. Produces a table specific for every detected sequence
length.

X = any degenerate aa
* = stop codon

Assumes standard input by default, but can be overridden by a filename.
"""
import argparse
import collections
import sys

parser = argparse.ArgumentParser(
    description=docstring, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('input_file', metavar='input_file',
                    type=argparse.FileType('r'), nargs='?',
                    default=sys.stdin,
                    help='File containing sequences, one per line.')
parser.add_argument('-n', action='store_true',
                    help='print whole numbers instead of proportions.')

args = parser.parse_args()

valid_aa = 'ACDEFGHIKLMNPQRSTVWY*X'

# read file into dict that is stratified by sequence lengths
sequences = {}
for line in args.input_file:
    seq = line.strip()
    if len(seq) not in sequences:
        sequences[len(seq)] = []
    
    sequences[len(seq)].append(seq)

# print per-aa values
for seqlen in sorted(sequences):
    # compile a dict that contains aa compositions, which makes it easier to
    # print stuff out in a tabular fashion later
    aa_comp = {}
    for n in range(seqlen):
        temp = collections.Counter([x[n] for x in sequences[seqlen]])
        aa_comp[n] = {x: temp[x] for x in valid_aa}
    
    # print table out
    print (f'## Sequence length = {seqlen} ##')
    print ('', *[x + 1 for x in range(seqlen)], sep='\t')
    for aa in valid_aa:
        if args.n:
            # print raw frequencies
            print (aa, *[aa_comp[x][aa] for x in range(seqlen)],
                   sep='\t')
        else:
            # print proportions
            print (aa,
                   *[round(aa_comp[x][aa] / sum(aa_comp[x].values()), 4) for x in range(seqlen)],
                   sep='\t')
    print ()

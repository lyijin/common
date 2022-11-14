#!/usr/bin/env python3

docstring = """
> diff_fasta.py <

When given two FASTA files, run a hacky diff on sequences shared in common,
and calculate the proportion of sequences that are identical.

Supports reading compressed FASTAs too.
""".strip()

import argparse
import math
from pathlib import Path
import statistics

import parse_fasta

def trim_annot(orig_annot):
    # based by assuming Illumina annots having e.g., EAS139:136:FC706VJ:2:2104:15343:197393
    # trim annots to first 7 colon-delimited fields
    orig_annot = orig_annot.split(' ')[0]
    orig_annot = ':'.join(orig_annot.split(':')[:7])
    return orig_annot

parser = argparse.ArgumentParser(
    description=docstring, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('fasta1', metavar='fasta_file', type=Path,
    help='first FASTA file')
parser.add_argument('fasta2', metavar='fasta_file', type=Path,
    help='second FASTA file')
parser.add_argument('-f', '--fuzzy_annots', action='store_true',
    help='brute-force matching assuming Illumina naming conventions.')

args = parser.parse_args()

# very inefficient code--reads all sequences from both files into memory
fasta1_seqs = parse_fasta.get_all_sequences(args.fasta1, 'fasta')
fasta2_seqs = parse_fasta.get_all_sequences(args.fasta2, 'fasta')

# print length details about both input files
fasta1_mean_len = statistics.mean(len(x) for x in fasta1_seqs.values())
fasta2_mean_len = statistics.mean(len(x) for x in fasta2_seqs.values())
print ('First file:', f'{args.fasta1.name} '
       f'(n = {len(fasta1_seqs)}, mean {round(fasta1_mean_len, 2)} bp)', sep='\t')
print ('Second file:', f'{args.fasta2.name} '
       f'(n = {len(fasta2_seqs)}, mean {round(fasta2_mean_len, 2)} bp)', sep='\t')
print ()

if args.fuzzy_annots:
    trimmed_labels = [trim_annot(x) for x in fasta1_seqs]
    assert len(set(trimmed_labels)) == len(trimmed_labels), \
        (f'For the first fasta file, there are {len(set(trimmed_labels))} '
         f'unique annots amongst {len(trimmed_labels)} overall annots. All '
         'annotations must be unique.')
    
    trimmed_labels = [trim_annot(x) for x in fasta2_seqs]
    assert len(set(trimmed_labels)) == len(trimmed_labels), \
        (f'For the second fasta file, there are {len(set(trimmed_labels))} '
         f'unique annots amongst {len(trimmed_labels)} overall annots. All '
         'annotations must be unique.')
    
    fasta1_seqs = {trim_annot(x):fasta1_seqs[x] for x in fasta1_seqs}
    fasta2_seqs = {trim_annot(x):fasta2_seqs[x] for x in fasta2_seqs}

# get sequence annots that are in common
common_annots = set(fasta1_seqs) & set(fasta2_seqs)

# check that there are common annots in both files
assert len(common_annots), \
    'There were no annotations in common across both files. Try -f to force fuzzy matching.'

# start tallying stats across files
num_identical_seq = 0
num_same_length = 0
list_nonidentical_len = []
print ('Total common annotations:  ', len(common_annots), sep='\t')
for ca in common_annots:
    if fasta1_seqs[ca] == fasta2_seqs[ca]:
        num_identical_seq += 1

    if len(fasta1_seqs[ca]) == len(fasta2_seqs[ca]):
        num_same_length += 1
    else:
        # append absolute length difference into a list
        list_nonidentical_len.append(abs(len(fasta1_seqs[ca]) - len(fasta2_seqs[ca])))

# dump stats out to stdout
print ('Identical sequences:        ', f'{num_identical_seq} '
       f'({round(num_identical_seq / len(common_annots) * 100, 2)}% of total)', sep='\t')
print ('Sequences with same length: ', f'{num_same_length} '
       f'({round(num_same_length / len(common_annots) * 100, 2)}% of total)', sep='\t')

if list_nonidentical_len:
    print ('\nFor sequences of dissimilar length,')
    print ('Mean length difference:     ', 
           f'{round(statistics.mean(list_nonidentical_len), 2)} bp', sep='\t')
    print ('SE length difference:       ', 
           f'{round(statistics.stdev(list_nonidentical_len) / math.sqrt(len(list_nonidentical_len)), 2)} bp',
           sep='\t')

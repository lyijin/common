#!/usr/bin/env python3

doc = f"""
> tabulate_peptide_freqs.py <

Tabulates frequencies from the any/all randomised segments for
`*.peptide_freq.tsv` files. Prints tabulated counts to stdout.
""".strip()

import argparse
import functools
import gzip
from pathlib import Path
import sys

import pandas as pd

parser = argparse.ArgumentParser(
    description=doc, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    'nucleotide_freq_tsv', type=Path, nargs='+',
    help='nucleotide frequency table files (can be compressed).')
parser.add_argument(
    '-c', '--cdr', type=int, default=0,
    help='focus on a particular CDR (default = all) that are delimited by "|".')
parser.add_argument(
    '-t', '--threshold', type=int, default=0,
    help='only output sequences with FPM > threshold (default = 0).')

args = parser.parse_args()

# read data
giant_dict = {}
for tsv_file in args.nucleotide_freq_tsv:
    df = pd.read_table(
        tsv_file if tsv_file.name[-3:] != '.gz' else gzip.open(tsv_file.name, 'rt'),
        header=None,
        names=['seq', 'count'],
        usecols=[0, 1])
    
    # manipulate sequence if user wants a particular CDR
    if args.cdr > 0:
        # the '- 1' is because python is 0-based and CDRs are numbered in a
        # 1-based fashion
        df['seq'] = df['seq'].str.split('|', expand=True)[args.cdr - 1]
        
        # then sum across rows with identical sequences
        df = df.groupby(['seq'], as_index=False)['count'].sum()
    
    # convert counts into FPM
    # (as to why i don't parse FPMs directly from the table, it's because those
    # numbers are rounded--and summing across rows will cause rounding errors)
    df['count'] = df['count'] / df['count'].sum() * 1000000
    
    # remove rows below threshold (when set)
    if args.threshold: df = df[df['count'] > args.threshold]
    
    # rename column header to filename
    short_fname = tsv_file.name.split('.peptide_freq')[0]
    df.columns = ['seq', short_fname]
    
    giant_dict[tsv_file.name] = df

# merge data
dfs = functools.reduce(
    lambda x, y: pd.merge(x, y, on='seq', how='outer'),
    [giant_dict[z] for z in giant_dict])
dfs = dfs.fillna(0)
dfs = dfs.sort_values(by=['seq'])

# print data
dfs.to_csv(sys.stdout, sep='\t', float_format='%.3f', index=False)

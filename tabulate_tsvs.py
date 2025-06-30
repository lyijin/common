#!/usr/bin/env python3

"""
> tabulate_tsvs.py <

Script to merge similar tsvs together, by checking common keys in the first
column (by default, can be toggled with --key).

Headers are assumed to be NOT PRESENT by default!

Design considerations: it's possible for keys and retained columns to be
interspersed, i.e. key1, col1, key2, col2, col3, ... but for output
purposes, keys always appear first, then columns next, i.e. key1, key2, 
col1_file1, col2_file1, col3_file1, col1_file2, col2_file2, col3_file2, etc.

Autodetects whether input files are gzip-compressed from `.gz` suffix.
"""
import argparse
from pathlib import Path
import sys
import time

import pandas as pd

parser = argparse.ArgumentParser(description="""
Script to merge similar tsvs together, by checking common keys in the first
column (by default, can be toggled with --key).""")

parser.add_argument('tsv_files', metavar='tsv_filenames', type=Path, nargs='+',
                    help='tab-separated filenames.')
parser.add_argument('--header', action='store_true',
                    help='tsv files have a header line (default: no headers!).')
parser.add_argument('--key', '-k', metavar='key_columns',
                    type=int, nargs='+', default=[0],
                    help='(0-based) columns as keys (default: 0).')
parser.add_argument('--col', '-c', metavar='retained_columns', type=int, nargs='+',
                    help='(0-based) columns as values (default: all except -k).')
parser.add_argument('--comment', metavar='comment_chars', type=str,
                    help='skip lines starting with comment chars.')
parser.add_argument('--how', default='left',
                    help='specify how to join tsvs: left/right/outer/inner (default: left).')
parser.add_argument('--fillna',
                    help='replace NaN values with specified string (default: blank).')
parser.add_argument('--mem_efficient', '-m', action='store_true',
                    help='memory-efficient mode, uses PyArrow strings. Depends on PyArrow.')
parser.add_argument('-v', action='store_true',
                    help='verbose mode, prints extra details to stderr.')
args = parser.parse_args()

# sanity checks
assert args.how in ['left', 'right', 'outer', 'inner'], \
    "args.how has to be one of 'left', 'right', 'outer', or 'inner'."
assert len(args.key) == len(set(args.key)), \
    'possible duplication of values in --key parameter.'
assert len(args.col) == len(set(args.col)), \
    'possible duplication of values in --col parameter.'

def diag_print(*args):
    """
    Helper function to print diagnostic-level verbose output to stderr.
    """
    print (f'[{time.asctime()}]', *args, file=sys.stderr)

if args.v:
    diag_print(f'Files used: {", ".join(x.name for x in args.tsv_files)}')
    
    if args.header:
        diag_print('Headers are PRESENT.')
    else:
        diag_print('Headers are NOT PRESENT.')

# if -c are used, define maxwidth so that the read_csv function below
# only reads maxwidth columns
maxwidth = 0
if args.col:
    maxwidth = max(args.key + args.col)

# create the header line for printing later
num_keys = len(args.key)
header = '\t' * num_keys

# go through the tsv files provided as args
for n, tsv_file in enumerate(args.tsv_files):
    # read data
    if args.v:
        diag_print(f'Reading file #{n + 1}: {tsv_file.name}...')
    
    df = pd.read_table(
        tsv_file,
        header=0 if args.header else None, 
        usecols=None if not maxwidth else range(0, maxwidth + 1),
        dtype='string[pyarrow]' if args.mem_efficient else 'string',
        comment=args.comment)
    
    # manipulate dataframe if args.col is not None, or when args.key isn't [0]
    if args.col:
        # in df, keys must appear before columns
        df = df.iloc[:, args.key + args.col]
    elif args.key != [0] and not args.col:
        num_cols = len(df.columns)
        non_keys = sorted(list(set(range(num_cols)) - set(args.key)))
        df = df.iloc[:, args.key + non_keys]
    
    # now that the keys are on the LHS, append filenames to non-key columns
    # to avoid errors during merging
    colnames = df.columns.values
    colnames = list(colnames[:len(args.key)]) + \
        [f'{x}_{tsv_file.name}' for x in colnames[len(args.key):]]
    df.columns = colnames
    
    # extend header line
    header += tsv_file.name
    header += '\t' * (len(df.columns) - num_keys)
    
    # merge data
    if args.v:
        diag_print(f'Merging file #{n + 1}: {tsv_file.name}...')
    
    if n == 0:
        # first past the post, no need to merge
        combined_df = df
    else:
        # remember, keys are always on the leftmost columns
        key_cols = range(len(args.key))
        combined_df = pd.merge(combined_df, df,
                               how=args.how,
                               left_on=list(combined_df.columns[key_cols]),
                               right_on=list(df.columns[key_cols]))

# fill NAs with an optional value if specified
if args.fillna:
    combined_df = combined_df.fillna(args.fillna)    

if args.v:
    diag_print(f'Union of all files produces {len(combined_df)} rows.')

# remove the final tab character before printing, then the whole thing
print (header[:-1])
print (combined_df.to_csv(sep='\t', header=args.header, index=False), end='')

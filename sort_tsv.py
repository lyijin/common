#!/usr/bin/env python3

"""
> sort_tsv.py <

Script to sort tsv by mean values of specified columns. Tries its hardest
to coerce strings into values, but might lead to unpredictable behaviour.
User should try to make sure these columns do not have strings.

By default, keys are assumed to be the first column, and values are all other
columns except the key. If the table is large, it might be easier to rule out
text-containing columns using -k <text_column1> <text_column2> than specifying
all value-containing columns.

Headers are assumed to be PRESENT, and a single line!

Autodetects whether input files are gzip-compressed from `.gz` suffix.
"""
import argparse
import gzip
from pathlib import Path
import sys
import time

import pandas as pd

def benchmark_print(message):
    """
    Prints a timestamped message to stderr, to aid visual assessment of
    chokepoints in code. 
    """
    print (f'[{time.asctime()}] {message}', file=sys.stderr)

parser = argparse.ArgumentParser(description="""
Script to sort tsv by mean values of specified columns. Breaks horribly if
the columns do not contain values, onus is on user to make sure these columns
do not have strings.""")

parser.add_argument('tsv_file', metavar='tsv_filename',
                    type=Path, nargs='?', default=sys.stdin,
                    help='tab-separated file.')
parser.add_argument('--noheader', action='store_true',
                    help='tsv file has no header line (default: has headers!).')
parser.add_argument('-k', '--key', metavar='key_columns',
                    type=int, nargs='+', default=[0],
                    help='(0-based) columns as keys (default: 0).')
parser.add_argument('-c', '--col', metavar='retained_columns',
                    type=int, nargs='+',
                    help='(0-based) columns as values (default: all except -k).')
parser.add_argument('-a', '--ascending', action='store_true',
                    help='sort in ascending manner (default: descending).')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='verbose mode, prints extra details to stderr.')
args = parser.parse_args()

if args.verbose:
    benchmark_print(f'File used: {args.tsv_file.name}')
        
    if args.noheader:
        benchmark_print('Headers: NOT present')
    else:
        benchmark_print('Headers: present')

# read data. dtype is forced to str to prevent floats getting increased
# precision when printed out
df = pd.read_table(args.tsv_file,
                   header=None if args.noheader else 0,
                   dtype=str)

# define which columns are value-containing. if user provided values for -c,
# respect that and don't auto-infer
if not args.col:
    num_cols = len(df.columns)
    args.col = sorted(list(set(range(num_cols)) - set(args.key)))

if args.verbose:
    benchmark_print(f'Value-containing columns: {args.col}')

# as all columns is read as str, they have to be coerced into numeric
temp = df.iloc[:, args.col].copy()
temp_nonnumeric_cols = temp.select_dtypes('object').columns.to_list()
if temp_nonnumeric_cols:
    for tnc in temp_nonnumeric_cols:
        temp[tnc] = pd.to_numeric(temp[tnc], errors='coerce')

# calculate means, store it in a new column. column has unrealistically long
# name to avoid collisions with columns in sorted table
df['meanmeanmeanmeanmean'] = temp.mean(axis=1)

# sort table in descending order (or ascending, if the flag is set)
df = df.sort_values(by='meanmeanmeanmeanmean', ascending=args.ascending)

# drop the mean column to restore original table
df = df.drop(columns='meanmeanmeanmeanmean')

# rename columns with 'Unnamed' back into '' (stops the printing of 'Unnamed')
df_cols = df.columns.to_list()
df_cols = ['' if 'Unnamed: ' in x else x for x in df_cols]
df.columns = df_cols

# print content
print (df.to_csv(sep='\t', header=not args.noheader, index=False), end='')

if args.verbose:
    benchmark_print(f'Finished sorting')

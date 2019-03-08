#!/usr/bin/env python3

"""
> tabulate_tsvs.py <

Script to merge similar tsvs together, by checking common keys in the first
column (by default, can be toggled with --key).

Headers are assumed to be NOT PRESENT by default!

Design considerations: it's possible for keys and retained columns to be
interspersed, i.e. key1, col1, key2, col2, col3, ... but for output
purposes, keys always appear first, then columns next, i.e. key1, key2, 
col1_file1, col2_file1, col3_file1, col1_file2, col2_file2, col3_file2, ...
"""
import argparse
import sys

import pandas as pd

parser = argparse.ArgumentParser(description="""
Script to merge similar tsvs together, by checking common keys in the first
column.""")

parser.add_argument('tsv_files', metavar='tsv_filenames',
                    type=argparse.FileType('r'), nargs='+',
                    help='tab-separated filenames.')
parser.add_argument('--header', action='store_true',
                    help='tsv files have a header line (default: no headers!).')
parser.add_argument('--key', '-k', metavar='key_columns',
                    type=int, nargs='+', default=[0],
                    help='(0-based) columns as keys (default: 0).')
parser.add_argument('--col', '-c', metavar='retained_columns',
                    type=int, nargs='+',
                    help='(0-based) columns as values (default: all except -k).')
parser.add_argument('--how', default='left',
                    help='specify how to join tsvs: left (default)/right/outer/inner.')
parser.add_argument('-v', action='store_true',
                    help='verbose mode, prints extra details to stderr.')
args = parser.parse_args()

# sanity checks
assert args.how in ['left', 'right', 'outer', 'inner'], \
    "args.how has to be one of 'left', 'right', 'outer', or 'inner'."

if args.v:
    print ('Files used: {}'.format(', '.join([x.name for x in args.tsv_files])),
        file=sys.stderr)
        
    if args.header:
        print ('Headers are PRESENT.', file=sys.stderr)
    else:
        print ('Headers are NOT PRESENT.', file=sys.stderr)

# if -c are used, define maxwidth so that the read_csv function below
# only reads maxwidth columns
maxwidth = 0
if args.col:
    maxwidth = max(args.key + args.col)

# read data
giant_dict = {}
for tsv_file in args.tsv_files:
    giant_dict[tsv_file.name] = \
        pd.read_csv(tsv_file, sep='\t',
                    header=0 if args.header else None, 
                    usecols=None if not maxwidth else range(0, maxwidth + 1),
                    dtype=str)
    
    # manipulate dataframe if args.col is not None, or when args.key isn't [0]
    if args.col:
        # in df, keys must appear before columns
        giant_dict[tsv_file.name] = \
            giant_dict[tsv_file.name].iloc[:, args.key + args.col]
    elif args.key != [0] and not args.col:
        num_cols = len(giant_dict[tsv_file.name].columns)
        non_keys = sorted(list(set(range(num_cols)) - set(args.key)))
        giant_dict[tsv_file.name] = \
            giant_dict[tsv_file.name].iloc[:, args.key + non_keys]
    
    if args.v:
        print (f'\rReading file #{args.tsv_files.index(tsv_file) + 1}',
               end='', file=sys.stderr)

# merge data
combined_data = giant_dict[args.tsv_files[0].name]

# remember, keys are always on the leftmost columns
key_cols = range(len(args.key))
for n in range(1, len(args.tsv_files)):
    left = combined_data
    right = giant_dict[args.tsv_files[n].name]
    
    combined_data = pd.merge(left, right,
                             how=args.how,
                             left_on=list(left.columns[key_cols]),
                             right_on=list(right.columns[key_cols]))

if args.v:
    print (f'\nUnion of all files produces {len(combined_data)} rows.',
           file=sys.stderr)

# print combined data out
num_keys = len(args.key)

# header line contain filenames that were merged
header = '\t' * num_keys
for tsv_file in args.tsv_files:
    header += tsv_file.name
    header += '\t' * (len(giant_dict[tsv_file.name].columns) - num_keys)

# remove the final tab character before printing
print (header[:-1])

# print content
print (combined_data.to_csv(sep='\t', header=args.header, index=False), end='')

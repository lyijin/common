#!/usr/bin/env python3

"""
> fishers_exact_two_cols.py <

Given a file with multiple columns (and at least two rows), carry out a Fisher's
exact test and B-H correct the calculated p-value between 1st and 2nd column
by default (i.e. --col 0 1).

To test multiple columns, use the --col/--c argument, supplying pairs in the
form of --col a b --col c d. Note that column numberings are zero-based
(i.e. column_0<TAB>column_1<TAB>...).
"""
import argparse
import collections
import csv
import sys

import scipy.stats

import correct_p_values

parser = argparse.ArgumentParser(description="""
Given a file with two columns (and at least two rows), carry out a Fisher's
exact test and B-H correct the calculated p-value.""")

parser.add_argument('tsv_file', metavar="tsv_file",
                    type=argparse.FileType('r'), nargs='?',
                    default=sys.stdin, help="TAB-separated file with values.")
parser.add_argument('--col', '-c', metavar="column_pairings",
                    type=int, nargs=2, action='append',
                    help="pairs of zero-based column numbers.")

args = parser.parse_args()
if not args.col: args.col = [[0, 1]]

tsv_reader = csv.reader(args.tsv_file, delimiter='\t')

# get the sum of each columns across all rows
all_cols = sorted(list(set([y for x in args.col for y in x])))
col_sum = {x: 0 for x in all_cols}

for row in tsv_reader:
    if not row: continue
    
    for a in all_cols:
        if row[a]:
            col_sum[a] += int(row[a])

# Fisher's exact test: let the four values in the 2x2 table be:
#   w  |  x
# -----------
#   y  |  z
for c in args.col:
    p_values = collections.OrderedDict()
    args.tsv_file.seek(0)
    for n, row in enumerate(tsv_reader):
        if not row: continue
        
        # if no data exists for w/y, assume they're both 0
        try:
            w = int(row[c[0]])
        except:
            w = 0
        x = col_sum[c[0]] - w
        
        try:
            y = int(row[c[1]])
        except:
            y = 0
        z = col_sum[c[1]] - y
                
        row_id = '{}_{}_{}'.format(n, w, y)
        p_values[row_id] = scipy.stats.fisher_exact([[w, x], [y, z]])[-1]
        
    corrected_p = correct_p_values.correct_p_values(p_values)
    
    # printing
    for p in p_values:
        orig_row = '\t'.join(p.split('_')[-2:])
        print (orig_row, corrected_p[p], sep='\t')
    print ()

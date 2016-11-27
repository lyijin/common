#!/usr/bin/env python3

"""
> fishers_exact_within_row.py <

Given a file with four or more columns, carry out a Fisher's exact test and 
(optionally) B-H correct the calculated p-value between four columns within
the row.

By default, columns 0, 1, 2 and 3 provide values for the variables w, x, y, z
in the following Fisher's exact table (--col 0 1 2 3):
    w  |  x
  -----------
    y  |  z

To test multiple columns, use the --col/--c argument, supplying quads in the
form of --col a b c d --col e f g h. Note that column numberings are zero-based
(i.e. column_0<TAB>column_1<TAB>...).
"""
import argparse
import collections
import csv
import sys

import scipy.stats

import correct_p_values

def check_int(string):
    if string:
        return int(string)
    else:
        return 0

parser = argparse.ArgumentParser(description="""
Given a file with four or more columns, carry out a Fisher's exact test and 
(optionally) B-H correct the calculated p-value between four columns within
the row.""")

parser.add_argument('tsv_file', metavar="tsv_file",
                    type=argparse.FileType('r'), nargs='?',
                    default=sys.stdin, help="TAB-separated file with values.")
parser.add_argument('--col', '-c', metavar="column_pairings",
                    type=int, nargs=4, action='append',
                    help="pairs of zero-based column numbers.")
parser.add_argument('--corr', '-r', metavar="correction type", type=str, 
                    default='BH', help="P value correction type, default 'BH'.")
parser.add_argument('--directionality', '-d', action='store_true',
                    help="""return P values that indicate enrichment/depletion.
                            +ve P values are enriched/not determined,
                            -ve P values are depleted.""")

args = parser.parse_args()
if not args.col: args.col = [[0, 1, 2, 3]]

tsv_reader = csv.reader(args.tsv_file, delimiter='\t')

# Fisher's exact test: let the four values in the 2x2 table be:
#   w  |  x
# -----------
#   y  |  z
directionality = {}
for c in args.col:
    p_values = collections.OrderedDict()
    args.tsv_file.seek(0)
    for n, row in enumerate(tsv_reader):
        if not row: continue
        
        # if no data exists for w/x/y/z, assume they're 0
        w = check_int(row[c[0]])
        x = check_int(row[c[1]])
        y = check_int(row[c[2]])
        z = check_int(row[c[3]])
        
        row_id = '_'.join([str(xx) for xx in [n, w, x, y, z]])
        
        if args.directionality:
            directionality[row_id] = 1
            try:
                if w/x < y/z:
                    directionality[row_id] = -1
            except:
                pass
        p_values[row_id] = scipy.stats.fisher_exact([[w, x], [y, z]])[-1]
        
    corrected_p = correct_p_values.correct_p_values(p_values, args.corr)
    
    # printing
    for p in p_values:
        orig_row = '\t'.join(p.split('_')[1:])
        if args.directionality:
            print (orig_row, corrected_p[p] * directionality[p], sep='\t')
        else:
            print (orig_row, corrected_p[p], sep='\t')
    print ()

#!/usr/bin/env python3

docstring = """
> fishers_exact_within_row.py <

Given a file with four or more columns, carry out a Fisher's exact test and 
(optionally) B-H correct the calculated p-value between four columns within
the row.

By default, columns 0, 1, 2 and 3 provide values for the variables w, x, y, z
in the following Fisher's exact table (--col 0 1 2 3):
    w  |  x  |  m
  -----------
    y  |  z  |
  -----------
    n  |  o  |  N

There's another mode --kmnN where values w, x, y, z are inferred from
  k = w             (# of successes)
  m = w + x         (row sum)
  n = w + y         (col sum)
  N = w + x + y + z (overall total)

... and another mode --wxno where top row and column sums are provided.
Inferring of w, x, y, z follows similar approaches where
  w = w
  x = x
  n = w + y
  o = x + z

To test multiple columns, use the --col/--c argument, supplying quads in the
form of --col a b c d --col e f g h. Note that column numberings are zero-based
(i.e. column_0<TAB>column_1<TAB>...).

Requires Python version > 3.7: dictionaries became ordered by default.
""".strip()

import argparse
import csv
import sys

import scipy.stats

def check_int(string):
    if string:
        try:
            return int(string)
        except ValueError:
            return 0
    else:
        return 0

parser = argparse.ArgumentParser(
    description=docstring, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('tsv_file', metavar='tsv_file',
                    type=argparse.FileType('r'), nargs='?',
                    default=sys.stdin, help='TAB-separated file with values.')
parser.add_argument('--col', '-c', metavar='column_pairings',
                    type=int, nargs=4, action='append',
                    help='pairs of zero-based column numbers.')
group = parser.add_mutually_exclusive_group()
group.add_argument('--kmnN', action='store_true',
                   help='provided values are NOT wxyz, they are kmnN.')
group.add_argument('--wxno', action='store_true',
                   help='provided values are NOT wxyz, they are wxno.')
parser.add_argument('--corr', '-r', metavar='correction type', type=str, 
                    default='BH', help="P value correction type, default 'BH'."
                                       "Use 'none' to suppress corretion.")
parser.add_argument('--directionality', '-d', action='store_true',
                    help='return P values that indicate enrichment/depletion.'
                         '+ve P values are enriched/not determined,'
                         '-ve P values are depleted.')

args = parser.parse_args()
if not args.col: args.col = [[0, 1, 2, 3]]

tsv_reader = csv.reader(args.tsv_file, delimiter='\t')

# Fisher's exact test: let the four values in the 2x2 table be:
#     w  |  x  |  m
#   -----------
#     y  |  z  |
#   -----------
#     n  |  o  |  N

directionality = {}
for c in args.col:
    p_values = {}
    for row_no, row in enumerate(tsv_reader):
        if not row: continue
        
        if len(row) < 4:
            # urgh row is not complete. assume wxyz are all 0
            w = x = y = z = 0
        else:
            # if no data exists for w/x/y/z, assume they're 0
            w = check_int(row[c[0]])
            x = check_int(row[c[1]])
            y = check_int(row[c[2]])
            z = check_int(row[c[3]])
        
        # deal with non-wxyz inputs
        if args.kmnN:
            # infer actual w, x, y, z values from the provided k, m, n, N
            # w == k, ignore
            x = x - w
            y = y - w
            z = z - (w + x + y)
        elif args.wxno:
            # infer actual w, x, y, z values from the provided w, x, n, o
            # w == w, ignore
            # x == x, ignore
            y = y - w
            z = z - x
        
        row_id = '_'.join([str(xx) for xx in [row_no, w, x, y, z]])
        
        if args.directionality:
            directionality[row_id] = 1
            try:
                if w/x < y/z:
                    directionality[row_id] = -1
            except:
                pass
        
        p_values[row_id] = scipy.stats.fisher_exact([[w, x], [y, z]])[-1]
    
    if args.corr == 'none':
        # no correction
        corrected_p = p_values
    else:
        # this is another script found in github:lyijin/common. download and
        # place this in the same folder as this script
        import correct_p_values
        
        corrected_p = correct_p_values.correct_p_values(p_values, args.corr)
    
    # printing
    for p in p_values:
        orig_row = '\t'.join(p.split('_')[1:])
        if args.directionality:
            print (orig_row, corrected_p[p] * directionality[p], sep='\t')
        else:
            print (orig_row, corrected_p[p], sep='\t')
    print ()

#!/usr/bin/env python3

"""
> plot_venn.py <

Given two/three input lists and optionally two/three labels, plot a Venn
diagram and provide a text list (to standard output) of the intersections
between the lists.

Lists are assumed to be one-per-line, and if there are multiple columns,
only the first one is considered.

Requires matplotlib, matplotlib_venn and numpy.
"""
import argparse
import collections
import csv

from matplotlib import pyplot as plt
import numpy as np

def get_venn2_subsets(a, b):
    """
    Given two sets, returns a dictionary containing:
      {'10': stuff present in a but not b,
       '11': stuff present in both a, b ...}
    """
    dict_venn2 = {'10': a - b,
                  '01': b - a,
                  '11': a & b}
    return dict_venn2

def get_venn3_subsets(a, b, c):
    """
    Given three sets, returns a dictionary containing:
      {'100': stuff present in a but not b/c,
       '011': stuff present in b/c but not a, ...}
    """
    dict_venn3 = {'100': a - (b | c),
                  '010': b - (a | c),
                  '001': c - (a | b),
                  '110': (a & b) - c,
                  '101': (a & c) - b,
                  '011': (b & c) - a,
                  '111': a & b & c}
    return dict_venn3

parser = argparse.ArgumentParser(description="""
Given two/three input lists and optionally two/three labels, plot a Venn
diagram and provide a text list (to standard output) of the intersections
between the lists.""")

parser.add_argument('tsv_files', metavar='tsv_files',
                    type=argparse.FileType('r'), nargs='+',
                    help='two or three input lists.')
parser.add_argument('-u', '--unweighted', action='store_true',
                    help='plot unweighted circles.')
parser.add_argument('-l', '--labels', metavar='labels',
                    type=str, nargs='+',
                    help='labels for the lists (default: filenames).')
parser.add_argument('--summary', action='store_true',
                    help='print (text) summary of overlaps.')
parser.add_argument('-s', '--savefig', metavar='pdf_filename', 
                    nargs='?', const='venn_diagram.pdf',
                    help="saves pdf output (default: 'venn_diagram.pdf').")
args = parser.parse_args()

# sanity checking
assert 2 <= len(args.tsv_files) <= 3, \
    'There can only be two or three input files!'

if args.labels:
    assert (len(args.tsv_files) == len(args.labels)), \
        'The number of labels must correspond to the number of files!'
else:
    args.labels = [x.name for x in args.tsv_files]

# read info from files
set_info = collections.OrderedDict()
for n, f in enumerate(args.tsv_files):
    tsv_reader = csv.reader(f, delimiter='\t')
    
    set_info[args.labels[n]] = []
    for row in tsv_reader:
        if not row: continue
        
        set_info[args.labels[n]].append(row[0])
    
    set_info[args.labels[n]] = set(set_info[args.labels[n]])

# plotting
plt.figure(figsize=(8,8))
if len(set_info) == 3:
    from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
    
    venn_data = get_venn3_subsets(*[set_info[x] for x in set_info])
    venn_sizes = {x:len(venn_data[x]) for x in venn_data}
    
    if args.unweighted:
        venn3_unweighted(venn_sizes, set_labels=args.labels)
        venn3_circles({'100': 1, '010': 1, '001': 1, '110': 1, 
                       '101': 1, '011': 1, '111': 1})
    else:
        venn3(venn_sizes, set_labels=args.labels)
        venn3_circles(venn_sizes)
else:
    from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
    
    venn_data = get_venn2_subsets(*[set_info[x] for x in set_info])
    venn_sizes = {x:len(venn_data[x]) for x in venn_data}
    
    if args.unweighted:
        venn2_unweighted(venn_sizes, set_labels=args.labels)
        venn2_circles({'10': 1, '01': 1, '11': 1})
    else:
        venn2(venn_sizes, set_labels=args.labels)
        venn2_circles(venn_sizes)

if args.savefig:
    # without bbox_inches, the saved figure has truncated axes.
    plt.savefig(args.savefig, bbox_inches='tight')

# print parsed output. fancy lambda basically enforces the order of keys
# being printed out, i.e. '111', '110', '101', '011', '100', '010', '001'.
for v in sorted(venn_data, key=lambda x: (x.count('1'), int(x)), reverse=True):
    print (', '.join([args.labels[n] for n, x in enumerate(v) if x == '1']),
           end=':\n')
    if args.summary:
        print (len(venn_data[v]), 'items.')
    else:
        print ('\n'.join(venn_data[v]))
    
    print ()

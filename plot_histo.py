#!/usr/bin/env python3

"""
> plot_histo.py <

Uses seaborn to plot a generic histogram of the input data.

Assumes that data is one-per-line, and if there are multiple columns, only the
first one is considered (i.e. same assumptions as plot_venn.py).
"""
import argparse
import csv

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

parser = argparse.ArgumentParser(description="""
Uses seaborn to plot a generic histogram of the input data.""")

parser.add_argument('tsv_file', metavar='tsv_file',
                    type=argparse.FileType('r'), nargs='+',
                    help='one/more input lists.')
parser.add_argument('--labels', '-l', metavar=('xlabel', 'ylabel'),
                    type=str, nargs=2,
                    help='labels for the x- and y-axes (default: none).')
parser.add_argument('--xlim', '-x', metavar=('min', 'max'),
                    type=int, nargs=2,
                    help='minimum and maximum limits of data (default: none).')
parser.add_argument('--bins', '-b', metavar='bin_number',
                    type=int, help='number of bins ' + \
                        '(default: matplotlib-controlled).')
parser.add_argument('--kde', '-k', action='store_false', default=True,
                    help='do not plot gaussian kde line.')
parser.add_argument('--bw', type=float,
                    help='change kde bandwidth to an int/float.')
parser.add_argument('--norm_hist', '-n', action='store_false', default=True,
                    help='do not normalise (i.e. plot absolute numbers).')
parser.add_argument('--savefig', '-s', metavar='pdf_filename', 
                    nargs='?', const='histogram.pdf',
                    help="save pdf output (default: 'histogram.pdf').")
args = parser.parse_args()

def int_or_float(s):
    try:
        return int(s)
    except ValueError:
        return float(s)
        
# read data
data = {}
for t in args.tsv_file:
    data[t.name] = []
    tsv_reader = csv.reader(t, delimiter='\t')
    for row in tsv_reader:
        if not row: continue
        
        data[t.name].append(int_or_float(row[0]))

# filter based on xlim
if args.xlim:
    for d in data:
        temp = []
        for x in data[d]:
            if args.xlim[0] < x < args.xlim[1]:
                temp.append(x)
            elif x <= args.xlim[0]:
                temp.append(args.xlim[0])
            elif x >= args.xlim[1]:
                temp.append(args.xlim[1])
                
        data[d] = temp

# calculate bins
if args.xlim and args.bins:
    desired_bins = np.linspace(args.xlim[0], args.xlim[1], args.bins + 1)
elif args.bins:
    desired_bins = args.bins
else:
    desired_bins = None

# modify kde bandwidth if args.bw is called
kde_kws = {}
if args.kde and args.bw:
    kde_kws = {'bw': args.bw}

# seaborn
sns.set_style('white')
sns.set_style('ticks')
#sns.set_palette('husl', len(args.tsv_file) + 1)
fig, ax = plt.subplots(figsize=(8, 5))       # smaller text
for d in sorted(data):
    sns.distplot(np.array(data[d]), bins=desired_bins, kde=args.kde, 
                 hist_kws={'alpha': 0.3}, kde_kws=kde_kws, 
                 norm_hist=args.norm_hist, label=d)

#if len(args.tsv_file) > 1:
plt.legend(title='File', loc=9, ncol=1)

if args.xlim:
    xlim_range = max(args.xlim) - min(args.xlim)
    ax.set_xlim(args.xlim[0] - .02 * xlim_range, args.xlim[1] + .02 * xlim_range)

sns.despine(offset=10, trim=True)

if args.labels:
    sns.axlabel(args.labels[0], args.labels[1])

if args.savefig:
    fig = plt.gcf()
    
    # without bbox_inches, the saved figure has truncated axes.
    fig.savefig(args.savefig, bbox_inches='tight')

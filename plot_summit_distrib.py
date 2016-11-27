#!/usr/bin/env python3

"""
> plot_summit_distrib.py <

Given histone modification summit positions across a stretch of sequence, 
calculate whether these positions are distributed differently from start/ends
of genes, than if locations were chosen at random.

This should be run after generate_dist_to_genes.py, as it produces an output
file for the latter script to run on.

Expected distributions within a list is calculated via Monte Carlo - the 
"all" list with n positions will have len(obs) positions picked randomly, 
and this is repeated 10,000 times to generate the expected distribution.
"""
import argparse
import csv
import random
import re
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import seaborn as sns

import natural_sort
import parse_fasta

parser = argparse.ArgumentParser(description="""
Given histone modification summit positions across a stretch of sequence, 
calculate whether these positions are distributed differently from start/ends
of genes, than if locations were chosen at random.""")

parser.add_argument('dist_tsv', metavar='tsv_file',
                    type=argparse.FileType('r'),
                    help='''file containing distances from 
                            generate_dist_to_genes.py.''')
parser.add_argument('bismark_cov', metavar='cov_file',
                    type=argparse.FileType('r'),
                    help='file with summit coordinates.')
parser.add_argument('--restrict', metavar='genome_feature',
                    type=str, default='genic_and_intergenic',
                    help='''restricts analysis to certain genomic features.
                            Valid values are ['genic_and_intergenic', 'genic',
                            'exonic', 'intronic', 'intergenic'].''')
parser.add_argument('--window', metavar='n', type=int, default=4000,
                    help='''only consider summit positions that are within 
                            n bp upstream/downstream of a gene's start/end.''')
parser.add_argument('--dist', metavar='type_of_distance',
                    type=str, default='upstream_to_start',
                    help='''calculates distance upstream/downstream to
                            start/ends of genes. Valid values are 
                            ['upstream_to_start', 'downstream_from_start',
                             'upstream_to_end', 'downstream_from_end'].''')
parser.add_argument('-s', '--savefig', metavar='pdf_file', 
                    nargs='?', const='meth_distrib.pdf',
                    help="save pdf output (default: 'meth_distrib.pdf').")
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints diagnostic stuff to stderr.')
args = parser.parse_args()

# sanity check
assert args.restrict in ['genic_and_intergenic', 'genic', 'exonic',
    'intronic', 'intergenic'], 'invalid --restrict keyword used!'
assert args.dist in ['upstream_to_start', 'downstream_from_start',
    'upstream_to_end', 'downstream_from_end'], 'invalid --dist keyword used!'

WINDOW = args.window

if args.verbose:
    print ('[{}] Start.'.format(time.asctime()), file=sys.stderr)

# read the file that stores how far each position is from gene start/
# gene ends, and save that to a dict
distances = {}
csv.field_size_limit(sys.maxsize)
tsv_reader = csv.reader(args.dist_tsv, delimiter='\t')
for row in tsv_reader:
    if not row: continue
    
    tmp = [x for x in row[1].split(', ')]
    distances[row[0]] = np.array(tmp, dtype='int32')

if args.verbose:
    print ('[{}] {} distances read.'.format(time.asctime(), len(distances)),
           file=sys.stderr)

# selective exclusion based on which genomic features are selected
include_no_info = False
include_genic = args.restrict in ['genic', 'genic_and_intergenic']
include_exonic = args.restrict == 'exonic'
include_intronic = args.restrict == 'intronic'
include_intergenic = args.restrict in ['intergenic', 'genic_and_intergenic']

# read in (0-based) locations of summits, and produce observed summit
# distributions (distances to gene starts/ends) for all scaffolds:
#   obs_summit_distrib[scaf] = [56, 789, 1012, ...]
# to generate expected summit distribution, also parse distances of
# all applicable positions to similar gene start/ends:
#   all_summit_distrib[scaf] = [34, 56, 778, 789, 1010, 1012, ...]
obs_summit_distrib = {}
all_summit_distrib = {}

tsv_reader = csv.reader(args.bismark_cov, delimiter='\t')
for row in tsv_reader:
    if not row: continue
    
    # selective exclusion based on which genomic features are selected
    feat = row[6]
    if feat == 'no_info' and not include_no_info: continue
    if feat == 'intergenic' and not include_intergenic: continue
    if 'gene' in feat.lower() and not include_genic: continue
    
    # only do this check if 'exonic' or 'intronic' is chosen
    if include_exonic or include_intronic:
        e_or_i = row[10]
        if e_or_i == 'no_info' and not include_no_info: continue
        if e_or_i[:4] == 'Exon' and not include_exonic: continue
        if e_or_i[:6] == 'Intron' and not include_intronic: continue
    
    scaf = row[0]
    coord = int(row[1]) - 1
    
    if scaf not in obs_summit_distrib:
        obs_summit_distrib[scaf] = []
    
    obs_summit_distrib[scaf].append(distances[scaf][coord])

if args.verbose:
    print ('[{}] Bismark cov processed.'.format(time.asctime()),
           file=sys.stderr)

# all_summit_distrib is simply distances < WINDOW bp in the scaffold
for scaf in obs_summit_distrib:
    all_summit_distrib[scaf] = distances[scaf]
    
    # remove "-1"s from all/obs_summit_distrib, as these values signify
    # summits at locations that are not considered in this analysis.
    # also filter out stuff that isn't within WINDOW bp.
    all_summit_distrib[scaf] = [x for x in all_summit_distrib[scaf]
                              if 0 <= x < WINDOW]
    obs_summit_distrib[scaf] = [x for x in obs_summit_distrib[scaf]
                              if 0 <= x < WINDOW]

# calculate mean distances!
if args.verbose:
    print ('[{}] Get observed and expected distances...'.format(time.asctime()), 
           file=sys.stderr)

obs_distances = np.array([y for x in obs_summit_distrib.values() for y in x],
                         dtype='int32')
exp_distances = []

if args.verbose:
    print ('[{}] Observed distances: {}'.format(
           time.asctime(), len(obs_distances)), file=sys.stderr)

for scaf in obs_summit_distrib:
    # pick m items from all_summit_distrib[scaf] without replacement
    m = len(obs_summit_distrib[scaf])
    
    # and do this n times
    for n in range(10000):
        m_items = random.sample(all_summit_distrib[scaf], m)
        exp_distances += m_items

exp_distances = np.array(exp_distances, dtype='int32')

if args.verbose:
    print ('[{}] Monte Carlo simulations done.'.format(time.asctime()), 
           file=sys.stderr)

# calculate chi-squared P value between observed and expected distributions;
# exp_matrix contain scaled-down values for the chisq test (in a way that 
# both matrices are equal in size).
obs_matrix = []
exp_matrix = []

increment = int(WINDOW / 40)
for n in range(0, WINDOW, increment):
    obs_matrix.append(len(obs_distances[obs_distances >= n]) - \
                      len(obs_distances[obs_distances >= n + increment]))
    exp_matrix.append(len(exp_distances[exp_distances >= n]) - \
                      len(exp_distances[exp_distances >= n + increment]))

multiplier = len(exp_distances) / len(obs_distances)
exp_matrix = [x/multiplier for x in exp_matrix]

chisq = scipy.stats.chisquare(obs_matrix, exp_matrix)
print ('SciPy chi-squared test')
print ('chisq test statistic = {}; p-value = {}'.format(*chisq))

# graph plotting with seaborn
sns.set_style('white')
sns.set_style('ticks')
plt.figure(figsize=(8,5))
sns.distplot(obs_distances, 
             bins=np.linspace(0, WINDOW, 41), norm_hist=True, color='#d8b365',
             hist_kws={'alpha': 0.3}, kde_kws={'bw': increment}, label='Observed')
sns.distplot(exp_distances, 
             bins=np.linspace(0, WINDOW, 41), norm_hist=True, color='#5ab4ac',
             hist_kws={'alpha': 0.3}, kde_kws={'bw': increment}, label='Expected')

if args.dist == 'upstream_to_start':
    sns.plt.xlim(int(WINDOW * 1.01), int(WINDOW * -0.01))
    sns.axlabel("Distance to gene 5' end (bp)", 'Density')
elif args.dist == 'downstream_from_start':
    sns.plt.xlim(int(WINDOW * -0.01), int(WINDOW * 1.01))
    sns.axlabel("Distance from gene 5' end (bp)", 'Density')
elif args.dist == 'upstream_to_end':
    sns.plt.xlim(int(WINDOW * 1.01), int(WINDOW * -0.01))
    sns.axlabel("Distance to gene 3' end (bp)", 'Density')
elif args.dist == 'downstream_from_end':
    sns.plt.xlim(int(WINDOW * -0.01), int(WINDOW * 1.01))
    sns.axlabel("Distance from gene 3' end (bp)", 'Density')

plt.legend(loc=9, ncol=2)
sns.despine(offset=10, trim=True)

# save plot
fig = plt.gcf()
fname = args.savefig if args.savefig else 'meth_distrib.pdf'
fig.savefig(fname, bbox_inches='tight')
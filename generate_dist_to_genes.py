#!/usr/bin/env python3

"""
> generate_dist_to_genes.py <

Given information about genic coordinates in scaffolds, generate a distance
matrix for each position's upstream/downstream distance to the closest 
genic start/stop.

This should be run before plot_feature_distrib.py, as it produces an output
file for the latter script to run on.
"""
import argparse
import csv
import sys
import time

import numpy as np

import natural_sort
import parse_fasta
import parse_gff3

parser = argparse.ArgumentParser(description="""
Given information about methylation positions across a stretch of sequence, 
calculate where these positions are relative to the start/ends of genes.""")

parser.add_argument('genome_fasta', metavar='fasta_file',
                    type=argparse.FileType('r'),
                    help='genome FASTA file.')
parser.add_argument('genome_gff3', metavar='gff3_file',
                    type=argparse.FileType('r'),
                    help='genome GFF3 annotation.')
parser.add_argument('--dist', metavar='type_of_distance',
                    type=str, default='upstream_to_start',
                    help='''calculates distance upstream/downstream to
                            start/ends of genes. Valid values are 
                            ['upstream_to_start', 'downstream_from_start',
                             'upstream_to_end', 'downstream_from_end'].''')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints diagnostic stuff to stderr.')
args = parser.parse_args()

# sanity check
assert args.dist in ['upstream_to_start', 'downstream_from_start',
    'upstream_to_end', 'downstream_from_end'], 'invalid --dist keyword used!'

def get_upstream_distances(gene_pos, scaf_len):
    """
    For consistency (a la annotate_bismark_cov.py), upstream distances can
    never be zero. Minimum upstream distance to the desired location is one.
    
    -1 denotes an invalid position (i.e. cannot calculate an upstream
    distance). Remember to exclude negative numbers subsequently!
    """
    upstream_dist_to_watson_pos = np.full(scaf_len, -1, dtype='int32')
    upstream_dist_to_crick_pos = np.full(scaf_len, -1, dtype='int32')
    upstream_dist_to_pos = np.full(scaf_len, -1, dtype='int32')
    
    watson_pos = [0] + [x for x in gene_pos if x > 0]
    paired_watsons = [x for x in zip(watson_pos[:-1], watson_pos[1:])]
    crick_pos = [-x for x in gene_pos if x < 0] + [scaf_len]
    paired_cricks = [x for x in zip(crick_pos[:-1], crick_pos[1:])]
    
    for p in paired_watsons:
        for n in range(p[0], p[1]):
            upstream_dist_to_watson_pos[n] = p[1] - n
    
    for p in paired_cricks:
        for n in range(p[0], p[1]):
            upstream_dist_to_crick_pos[n] = n - p[0] + 1
    
    # on a per-position basis, take the shortest distance to gene starts that 
    # lie either on Watson or Crick strands.
    for n in range(scaf_len):
        if upstream_dist_to_watson_pos[n] >= 0 and \
           upstream_dist_to_crick_pos[n] >= 0:
            upstream_dist_to_pos[n] = min(upstream_dist_to_watson_pos[n], 
                                          upstream_dist_to_crick_pos[n])
        else:
            # either or both distances are -1. if there's a non-negative number,
            # pick that; if both are -1, store a -1. DO VOODOO MATH!
            upstream_dist_to_pos[n] = -(upstream_dist_to_watson_pos[n] * \
                                        upstream_dist_to_crick_pos[n])
    
    return (upstream_dist_to_pos)

def get_downstream_distances(gene_pos, scaf_len):
    """
    For consistency (a la annotate_bismark_cov.py), downstream distances start
    at zero.
    
    -1 denotes an invalid position (i.e. cannot calculate a downstream
    distance). Remember to exclude negative numbers subsequently!
    """
    downstream_dist_to_watson_pos = np.full(scaf_len, -1, dtype='int32')
    downstream_dist_to_crick_pos = np.full(scaf_len, -1, dtype='int32')
    downstream_dist_to_pos = np.full(scaf_len, -1, dtype='int32')
    
    watson_pos = [x for x in gene_pos if x > 0] + [scaf_len]
    paired_watsons = [x for x in zip(watson_pos[:-1], watson_pos[1:])]
    crick_pos = [0] + [-x for x in gene_pos if x < 0]
    paired_cricks = [x for x in zip(crick_pos[:-1], crick_pos[1:])]
    
    for p in paired_watsons:
        for n in range(p[0], p[1]):
            downstream_dist_to_watson_pos[n] = n - p[0] 
    
    for p in paired_cricks:
        for n in range(p[0], p[1]):
            downstream_dist_to_crick_pos[n] = p[1] - (n + 1)
    
    # on a per-position basis, take the shortest distance to gene starts that 
    # lie either on Watson or Crick strands.
    for n in range(scaf_len):
        if downstream_dist_to_watson_pos[n] >= 0 and \
           downstream_dist_to_crick_pos[n] >= 0:
            downstream_dist_to_pos[n] = min(downstream_dist_to_watson_pos[n], 
                                            downstream_dist_to_crick_pos[n])
        else:
            # either or both distances are -1. if there's a non-negative number,
            # pick that; if both are -1, store a -1. DO VOODOO MATH!
            downstream_dist_to_pos[n] = -(downstream_dist_to_watson_pos[n] * \
                                          downstream_dist_to_crick_pos[n])
    
    return (downstream_dist_to_pos)

if args.verbose:
    print ('[', time.asctime(), '] Start.',
           sep='', file=sys.stderr)

# grab genome sequence to identify all CpG sites
genome_sequences = parse_fasta.get_all_sequences(args.genome_fasta, 'fasta')

# read coordinates of genes and exons from .gff3 file.
scaffold_gff3 = parse_gff3.parse_gff3(args.genome_gff3, 'gene')

for scaf in natural_sort.natural_sort(genome_sequences):
    # only allow scaffolds with annotated genic regions to pass the check
    if scaf not in scaffold_gff3: continue
    
    if args.verbose:
        print ('[', time.asctime(), '] Calculating distances for ', scaf, '...',
               sep='', file=sys.stderr)
    
    scaf_len = len(genome_sequences[scaf])
    scaf_gene_coords = []
    
    for gene in natural_sort.natural_sort(scaffold_gff3[scaf]):
        gene_coords = scaffold_gff3[scaf][gene].coords
        scaf_gene_coords.append(gene_coords)
    
    # based on the --dist keyword, calculate upstream or downstream distances
    # to gene starts or ends
    gene_starts = [x if x < y else -x for (x, y) in scaf_gene_coords]
    gene_ends = [y if x < y else -y for (x, y) in scaf_gene_coords]
    
    # THESE LINES MAKE SURE SCRIPT CAN DEAL WITH UNSORTED GFF3!
    gene_starts = sorted(gene_starts, key=abs)
    gene_ends = sorted(gene_ends, key=abs)
    
    if args.dist == 'upstream_to_start':
        distances = get_upstream_distances(gene_starts, scaf_len)
    elif args.dist == 'downstream_from_start':
        distances = get_downstream_distances(gene_starts, scaf_len)
    elif args.dist == 'upstream_to_end':
        distances = get_upstream_distances(gene_ends, scaf_len)
    elif args.dist == 'downstream_from_end':
        distances = get_downstream_distances(gene_ends, scaf_len)

    # write the distances into a file for downstream calculations/plotting
    #   Spis.scaffold... \t dist1, dist2, dist3, ...
    print (scaf, ', '.join([str(x) for x in distances]), sep='\t')

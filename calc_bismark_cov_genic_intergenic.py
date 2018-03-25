#!/usr/bin/env python3

"""
> calc_bismark_cov_genic_intergenic.py <

Script returns the number of methylated positions in intergenic/intronic/exonic
regions. This is theoretically also produced by annotate_bismark_cov.py,
but that script takes much longer as it covers more bases.
"""
import argparse
import collections
import csv
import time

import numpy as np

import natural_sort
import parse_fasta
import parse_gff3

parser = argparse.ArgumentParser(description="""
Script calculates the distribution of methylation over entire genic/intergenic
regions, and prints out lotsa details.""")

parser.add_argument('genome_fasta', metavar='fasta_file',
                    type=argparse.FileType('r'),
                    help='genome FASTA file.')
parser.add_argument('genome_gff3', metavar='gff3_file',
                    type=argparse.FileType('r'),
                    help='genome GFF3 annotation.')
parser.add_argument('bismark_cov', metavar="cov_filename",
                    type=argparse.FileType('r'),
                    help="Bismark .cov filename.")
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints diagnostic stuff to stderr.')
args = parser.parse_args()

def is_numpy_array_filled(np_array, coords):
    return np.any(np_array[min(coords):max(coords)])

def assign_value_to_numpy_array(np_array, coords, assigned_value):
    """
    Coords are in a tuple (coords1, coords2). Only +ve values are assigned
    to the region.
    """
    np_array[min(coords):max(coords)] = assigned_value
    
    return np_array

# read sequences
sequence_lengths = parse_fasta.get_all_sequences(args.genome_fasta, 
        'fasta', lengths_only=True)
if args.verbose:
    print ('[{}] Lengths for {} sequences parsed.'.format(
           time.asctime(), len(sequence_lengths)), file=sys.stderr)

# read coordinates of genes and exons from .gff3 file.
scaffold_gff3 = parse_gff3.parse_gff3(args.genome_gff3, 'exon')

# as genes might contain overlapping isoforms, the longest isoform is chosen,
# if multiples exist.
scaffold_gff3 = parse_gff3.pick_longest_mRNA(scaffold_gff3)

# make sure features in all mRNAs are sorted properly (for exon numbering).
scaffold_gff3 = parse_gff3.sort_features(scaffold_gff3)

# genic regions are denoted in a NumPy array as follows:
#    0: intergenic region
#    1: genic region
#
# exonic data is stored in a NumPy array:
#    0: intergenic region
#    1: exon
#    2: intron
#
# NOTE: if genes overlap, the first one that gets annotated (sorted 
# alphabetically) has priority.
gene_info = {}
exon_info = {}

for s in natural_sort.natural_sort(scaffold_gff3):
    gene_info[s] = np.zeros(sequence_lengths[s], np.int8)
    exon_info[s] = np.zeros(sequence_lengths[s], np.int8)
    
    for gene_id in natural_sort.natural_sort(scaffold_gff3[s]):
        # jump out of loop if the assigned region is not zeroes ([0, 0, ...])
        gene_coords = scaffold_gff3[s][gene_id].coords
        if is_numpy_array_filled(gene_info[s], gene_coords): continue
        
        # gene annot:
        gene_info[s] = assign_value_to_numpy_array(gene_info[s], gene_coords, 1)
        
        # exon annot:
        # sole_mRNA is a mRNA object, with the assumption that *.mRNAs
        # returns a dict containing {mRNA_name: mRNA_object}
        sole_mRNA = list(scaffold_gff3[s][gene_id].mRNAs.values())[0]
        
        exon_coords = sole_mRNA.details['exon']
        for e in exon_coords:
            exon_info[s] = assign_value_to_numpy_array(exon_info[s], e, 1)
        
        intron_coords = [(x[1], y[0]) for x, y in
                         zip(exon_coords[:-1], exon_coords[1:])]
        for i in intron_coords:
            exon_info[s] = assign_value_to_numpy_array(exon_info[s], i, 2)

# create empty list to store exon_info of methylated positions
pos_of_interest_info = []
counter = 0

# read per-position # of methylated and unmethylated bases
tsv_reader = csv.reader(args.bismark_cov, delimiter='\t')
for row in tsv_reader:
    # skip empty rows
    if not row: continue
    
    # skip lines with scaffold names not in the gff3 file
    scaf = row[0]
    pos = int(row[1]) - 1
    if scaf not in exon_info:
        pos_of_interest_info.append(0)
    else:
        pos_of_interest_info.append(exon_info[scaf][pos])

# print output
pos_of_interest_info = collections.Counter(pos_of_interest_info)
print ('total', 'intergenic', 'exonic', 'intronic', sep='\t')
print (sum(pos_of_interest_info.values()),
       *[pos_of_interest_info[x] for x in range(3)], sep='\t')

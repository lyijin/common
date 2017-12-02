#!/usr/bin/env python3

"""
> annotate_samtools_depth.py <

Script takes in the scaffold/position information of a genome, and prints out
whether the position is in a gene, and which exon it is in.

Note that script discards all intronic and intergenic info!
"""
import argparse
import csv
import gzip
import re
import sys
import time

import numpy as np

import natural_sort
import parse_fasta
import parse_gff3

parser = argparse.ArgumentParser(description="""
Script takes in the scaffold/position information of a genome, and prints out
whether the position is in a gene, and which exon it is in.""")

parser.add_argument('genome_fasta', metavar='fasta_file',
                    type=argparse.FileType('r'),
                    help='genome FASTA file.')
parser.add_argument('genome_gff3', metavar='gff3_file',
                    type=argparse.FileType('r'),
                    help='genome GFF3 annotation.')
parser.add_argument('s_depth_tsv', metavar="tsv_gz_filename",
                    type=argparse.FileType('r'), nargs='?',
                    default=sys.stdin, help="gzipped samtools depth file.")
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints diagnostic stuff to stderr.')
args = parser.parse_args()

def is_numpy_array_filled(np_array, coords):
    return np.any(np_array[min(coords):max(coords)])

def assign_value_to_numpy_array(np_array, coords, assigned_value):
    """
    Coords are in a tuple (coords1, coords2). +ve and -ve values are assigned
    respectively, depending on coords1 < coords2 or vice versa.
    """
    if coords[1] > coords[0]:
        np_array[coords[0]:coords[1]] = assigned_value
    elif coords[0] > coords[1]:
        np_array[coords[1]:coords[0]] = -assigned_value

    return np_array

def translate_exon_info(ei_integer, gene_id):
    """
    Change integers stored in NumPy array into strings, e.g.:
        1 --> Exon_1
        6 --> Intron_3
        
    Also return info in the negative form (i.e. Exon_-1 is the last exon,
    Exon_-2 is penultimate exon...).
    """
    # in rare cases, the longest transcript in a gene is shorter than the 
    # gene itself. Thus, there's a possibility that
    #     gene_info != AND exon_info == 0 !!
    if not ei_integer: return 'no_info'
    
    ei_integer = abs(int(ei_integer))   # NumPy ints aren't... really ints.    
    ei = 'Exon' if ei_integer % 2 else 'Intron'
    
    # round() is dangerous. Banker's rounding.
    ei_number = (ei_integer + 1) // 2       
    ei_inverse = ei_number - exon_count[gene_id] - (ei_integer % 2)
    return '_'.join([ei, str(ei_number), str(ei_inverse)])

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
#  +ve: gene in the + strand
#  -ve: gene in the - strand
#
# gene_names is a dict that converts ints in gene_info back to gene names.
#
# NOTE: if genes overlap, the first one that gets annotated (sorted 
# alphabetically) has priority.
gene_info = {}
gene_names = {}
gene_counter = 0

# similar to genic data, exonic data is stored in a NumPy array. However,
# exons/introns in all genes have the following numbering scheme:
#    0: intergenic region
# 1, 2: exon 1, intron 1
# 3, 4: exon 2, intron 2, ...
#  +ve: gene in the + strand
#  -ve: gene in the - strand
#
# exon_count stores # of exons in that gene.
exon_info = {}
exon_count = {}

for s in natural_sort.natural_sort(scaffold_gff3):
    gene_info[s] = np.zeros(sequence_lengths[s], np.int32)
    exon_info[s] = np.zeros(sequence_lengths[s], np.int8)
    for gene_id in natural_sort.natural_sort(scaffold_gff3[s]):
        # gene detail annotation
        gene_coords = scaffold_gff3[s][gene_id].coords
        
        # jump out of loop if the assigned region is not zeroes ([0, 0, ...]).
        if is_numpy_array_filled(gene_info[s], gene_coords): continue
        
        gene_counter += 1
        gene_info[s] = assign_value_to_numpy_array(gene_info[s], gene_coords, 
                                                   gene_counter)
        gene_names[gene_counter] = gene_id
        
        # exon detail annotation
        #        
        # sole_mRNA is a mRNA object, with the assumption that *.mRNAs
        # returns a dict containing {mRNA_name: mRNA_object}
        sole_mRNA = list(scaffold_gff3[s][gene_id].mRNAs.values())[0]
        
        exon_coords = sole_mRNA.details['exon']
        exon_count[gene_id] = len(exon_coords)
        unpacked_coords = [y for x in exon_coords for y in x]
        exon_intron_coords = [(x, y) for x, y in zip(unpacked_coords[:-1],
                                                     unpacked_coords[1:])]
        
        for n, ei in enumerate(exon_intron_coords):
            exon_info[s] = assign_value_to_numpy_array(exon_info[s], ei, n+1)
        
if args.verbose:
    print ('[{}] {} sequences with gene annotations parsed.'.format(
           time.asctime(), len(scaffold_gff3)), file=sys.stderr)
    print ('[{}] {} genes parsed in total.'.format(
           time.asctime(), gene_counter), file=sys.stderr)

# read per-position # of methylated and unmethylated bases
tsv_reader = csv.reader(gzip.open(args.s_depth_tsv.name, 'rt'), delimiter='\t')
counter = 0
for row in tsv_reader:
    # skip empty rows
    if not row: continue
    
    # to make script work on annotated samtools depth files (i.e. overwrite 
    # previous annotation), retain only first 3 columns
    row = row[:3]
    
    # skip lines with 0 coverage
    cov = int(row[2])
    if not cov: continue
    
    # skip lines with scaffold names not in the gff3 file
    scaf = row[0]
    if scaf not in gene_info: continue
    
    # all systems go!
    pos = int(row[1]) - 1
    
    # add gene info first
    extra_info = []
    
    # discard non-genic positions
    if gene_info[scaf][pos]:
        gene_id = gene_names[abs(gene_info[scaf][pos])]
        ei = translate_exon_info(exon_info[scaf][pos], gene_id)
        
        # only print exonic positions
        if 'Exon' in ei:
            print (*row, gene_id, ei, sep='\t')
    
    if args.verbose:
        counter += 1
        if counter % 10000 == 0:
            print ('[{}] {} methylated positions processed...'.format(
                   time.asctime(), counter), file=sys.stderr)

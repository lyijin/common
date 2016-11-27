#!/usr/bin/env python3

"""
> extract_upstream_regions.py <

Script extracts a user-defined region upstream of all gene models (well,
within the boundaries of the scaffold), preserving capitalisation of the
nucleotides.

This script makes no guarantee that all of the sequences extracted are
intergenic in nature. If the window is big enough, it'll encompass the entire
scaffold *shrug*.

Why this design decision? It's to deal with alternative promoters. If two
overlapping nested genes have different 5' genes, the shorter gene might
jolly well have a promoter region that includes genic region of the longer one.

Too bad. In practice though, if the window is small enough, MOST of the
extracted sequence WILL be from intergenic regions.
"""
import argparse

import natural_sort
import parse_fasta
import parse_gff3

parser = argparse.ArgumentParser(description="""
Script extracts a user-defined region upstream of all gene models (well,
within the boundaries of the scaffold), preserving capitalisation of the
nucleotides.""")

parser.add_argument('genome_fasta', metavar='fasta_file',
                    type=argparse.FileType('r'),
                    help='genome FASTA file.')
parser.add_argument('genome_gff3', metavar='gff3_file',
                    type=argparse.FileType('r'),
                    help='genome GFF3 annotation.')
parser.add_argument('--window', '-w', type=int, default=4000,
                    help='Define window size (default=4000).')

args = parser.parse_args()

def reverse_complement(seq):
    seq = seq.replace('U', 'T')

    translation_from = 'AaTtGgCcYyRrSsWwKkMmBbDdHhVvNn'
    translation_to   = 'TtAaCcGgRrYySsWwMmKkVvHhDdBbNn'
    translation_table = str.maketrans(translation_from, translation_to)

    seq = seq[::-1].translate(translation_table)

    return seq

# create dictionary to store upstream sequences
upstream_seq = {}

# read window size and sequences
WINDOW = args.window
genome_seq = parse_fasta.get_all_sequences(args.genome_fasta, 'fasta')
scaffold_gff3 = parse_gff3.parse_gff3(args.genome_gff3, 'gene')

for s in scaffold_gff3:
    for gene_id in scaffold_gff3[s]:
        gene_start, gene_end = scaffold_gff3[s][gene_id].coords
        
        if gene_start < gene_end:
            # gene is on watson strand
            upstream_start = max(gene_start - WINDOW, 0)
            upstream_end = gene_start
            upstream_seq[gene_id] = genome_seq[s][upstream_start:upstream_end]
        else:
            # gene is on crick strand
            upstream_start = min(gene_start + WINDOW, len(genome_seq[s]))
            upstream_end = gene_start
            upstream_seq[gene_id] = \
                reverse_complement(genome_seq[s][upstream_end:upstream_start])

# output
for gene_id in natural_sort.natural_sort(upstream_seq):
    print ('>{}_{}_upstream'.format(gene_id, WINDOW))
    print (upstream_seq[gene_id])

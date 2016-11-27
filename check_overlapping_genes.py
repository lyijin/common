#!/usr/bin/env python3

"""
> check_overlapping_genes.py <

Given a gff3 file, check whether genes overlap extensively--idea is to remove
isoforms, but allow overlapping genes that code for different protein products
to survive the cull.

To determine this, code checks whether genes are overlapping, and if they do,
check whether they share at least an exon. If they don't share a single exon,
it means that they code for different protein products!

Outputs overlapping genes, lengths and coordinates for downstream processing.

NOTE: it's probably easier to use the output to choose genes for EXCLUSION from
the set of gene models (i.e. pick_reads.py --exclude). In this case, use the 
data in the fourth column of the output as the exclusion list, as the shorter
gene name is printed in that column.
"""
import argparse

import natural_sort
import parse_gff3

parser = argparse.ArgumentParser(description="""
Given a gff3 file, check whether genes overlap (i.e. bits of them occupy same 
loci). Outputs which genes overlap, and their lengths for downstream processing.
""")

parser.add_argument('scaffold_gff3', metavar="gff3_file",
                    type=argparse.FileType('r'), 
                    help="gff3 file of a genome.")

args = parser.parse_args()

def calc_seq_length(coords):
    """
    Given a bunch of coords [(x1, y1), (x2, y2), ...], calculate the length
    of the sequence.
    """
    cumulative_len = 0
    for c in coords:
        cumulative_len += max(c) - min(c) + 1
    
    return cumulative_len

# read gff3
scaffold_gff3 = parse_gff3.parse_gff3(args.scaffold_gff3, 'exon')
scaffold_gff3 = parse_gff3.pick_longest_mRNA(scaffold_gff3)
for scaf in natural_sort.natural_sort(scaffold_gff3):
    genes_in_scaf = natural_sort.natural_sort(scaffold_gff3[scaf])
    genes_in_scaf_remainder = genes_in_scaf[:]
    for gene in genes_in_scaf:
        gene_coords = scaffold_gff3[scaf][gene].coords
        gene_len = calc_seq_length([gene_coords])
        
        tx = list(scaffold_gff3[scaf][gene].mRNAs.keys())[0]
        mrna_coords = scaffold_gff3[scaf][gene].mRNAs[tx].details['exon']
        mrna_len = calc_seq_length(mrna_coords)
        
        genes_in_scaf_remainder = genes_in_scaf_remainder[1:]
        # check whether gene overlaps with other genes on the same scaffold
        for gene2 in genes_in_scaf_remainder:
            gene2_coords = scaffold_gff3[scaf][gene2].coords
            
            if min(gene2_coords) <= max(gene_coords) <= max(gene2_coords) or \
               min(gene_coords) <= max(gene2_coords) <= max(gene_coords):
                # overlap between gene and gene2!
                gene2_len = calc_seq_length([gene2_coords])
                
                tx2 = list(scaffold_gff3[scaf][gene2].mRNAs.keys())[0]
                mrna2_coords = scaffold_gff3[scaf][gene2].mRNAs[tx2].details['exon']
                mrna2_len = calc_seq_length(mrna2_coords)
                
                # check whether there's at least a common exon: if not, abort!
                if not set(mrna_coords) & set(mrna2_coords): continue
                
                # always print shorter mRNA on the left and longer on the right.
                # alphabetical order is used to break ties:
                #   i.e. if Gene123 and Gene124 are equal, Gene123 is "longer"
                if mrna2_len > mrna_len:
                    print (gene2, mrna2_len, gene2_len, gene, mrna_len, gene_len,
                           scaf, mrna2_coords, gene2_coords, 
                           mrna_coords, gene_coords, sep='\t')
                else:
                    print (gene, mrna_len, gene_len, gene2, mrna2_len, gene2_len,
                           scaf, mrna_coords, gene_coords, 
                           mrna2_coords, gene2_coords, sep='\t')

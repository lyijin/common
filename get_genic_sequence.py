#!/usr/bin/env python3

"""
> get_genic_sequence.py <

CDS files do not contain intronic regions - this script reads the gff3 file
containing start and end coordinates for genes, and extracts the genic 
(i.e. exonic + intronic) sequences for each gene.
"""
import argparse

import parse_fasta
import parse_gff3

def reverse_complement(seq):
    seq = seq.replace('U', 'T')

    translation_from = 'AaTtGgCcYyRrSsWwKkMmBbDdHhVvNn'
    translation_to   = 'TtAaCcGgRrYySsWwMmKkVvHhDdBbNn'
    translation_table = str.maketrans(translation_from, translation_to)

    seq = seq[::-1].translate(translation_table)

    return seq

parser = argparse.ArgumentParser(description="""
CDS files do not contain intronic regions - this script reads the gff3 file
containing start and end coordinates for genes, and extracts the genic 
(i.e. exonic + intronic) sequences for each gene.""")

parser.add_argument('genome_fasta', metavar="fasta_file",
                    type=argparse.FileType('r'), 
                    help="FASTA file of the genome.")
parser.add_argument('scaffold_gff3', metavar="gff3_file",
                    type=argparse.FileType('r'), 
                    help="corresponding gff3 file of the genome.")

args = parser.parse_args()

# read genome details into memory
genome_fasta = parse_fasta.get_all_sequences(args.genome_fasta, 'fasta')
scaffold_gff3 = parse_gff3.parse_gff3(args.scaffold_gff3, 'gene')

# read the positions from the cov file
for scaf in scaffold_gff3:
    for gene in scaffold_gff3[scaf]:
        gene_coords = scaffold_gff3[scaf][gene].coords
        on_crick = gene_coords[0] > gene_coords[1]
        
        genic_seq = genome_fasta[scaf][min(gene_coords):max(gene_coords)]
        if on_crick:
            genic_seq = reverse_complement(genic_seq)
            
        print ('>' + gene)
        print (genic_seq)

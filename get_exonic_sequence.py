#!/usr/bin/env python3

"""
> get_exonic_sequence.py <

Based on input genome fasta and gff3 files, parse _exonic_ sequences of each
gene.
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
Based on input genome fasta and gff3 files, parse _exonic_ sequences of each
gene.""")

parser.add_argument('genome_fasta', metavar="fasta_file",
                    type=argparse.FileType('r'), 
                    help="FASTA file of the genome.")
parser.add_argument('scaffold_gff3', metavar="gff3_file",
                    type=argparse.FileType('r'), 
                    help="corresponding gff3 file of the genome.")

args = parser.parse_args()

# read genome details into memory
genome_fasta = parse_fasta.get_all_sequences(args.genome_fasta, 'fasta')
scaffold_gff3 = parse_gff3.parse_gff3(args.scaffold_gff3, 'exon')

# pick longest transcript in the gff3 file
scaffold_gff3 = parse_gff3.pick_longest_mRNA(scaffold_gff3)

# read the positions from the cov file
for scaf in scaffold_gff3:
    for gene in scaffold_gff3[scaf]:
        gene_coords = scaffold_gff3[scaf][gene].coords
        gene_on_crick = gene_coords[0] > gene_coords[1]

        tx = list(scaffold_gff3[scaf][gene].mRNAs.keys())[0]
        mrna_coords = scaffold_gff3[scaf][gene].mRNAs[tx].details['exon']
        
        exon_seq = ''
        for i in mrna_coords:
            temp = genome_fasta[scaf][min(i):max(i)]
            if gene_on_crick: temp = reverse_complement(temp)
            
            exon_seq += temp
        
        print ('>' + gene)
        print (exon_seq)

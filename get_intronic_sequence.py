#!/usr/bin/env python3

"""
> get_intronic_sequence.py <

Based on input genome fasta and gff3 files, parse _intronic_ sequences of each
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
        
        # get intron coords, then construct the intronic sequence
        exon_starts = [x for x, y in mrna_coords]
        exon_ends = [y for x, y in mrna_coords]
        
        if len(exon_starts) < 2:
            # there's only 1 exon, i.e. no introns are present
            continue
        
        intron_coords = [(x, y) for x, y in zip(exon_ends[:-1], exon_starts[1:])]
        
        intron_seq = ''
        for i in intron_coords:
            temp = genome_fasta[scaf][min(i):max(i)]
            if gene_on_crick: temp = reverse_complement(temp)
            
            intron_seq += temp
        
        print ('>' + gene)
        print (intron_seq)

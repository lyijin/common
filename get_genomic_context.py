#!/usr/bin/env python3

"""
> get_genomic_context.py <

Parses the genome for contextual sequences surrounding positions of interest.
Padding Ns are added if the positions are incredibly close to the start/end
of scaffolds.

Input file should contain scaffold name \t 1-based (not 0-based) position in 
the first two columns. Other columns would simply be ignored.

Default window sizes are -100 to +100 of the site (i.e. 201 bp sequences).

NOTE THAT THE OUTPUT, IN FASTA FORMAT, WILL ALWAYS BE ON THE WATSON STRAND.
"""
import argparse
import csv

import parse_fasta

def get_genomic_context(genome_fasta_file, pos_of_interest_file, window):
    """
    Function that does the heavy lifting. Returns a dictionary of sequences.
    """
    window_seqs = {}
    
    # read genome sequences
    genome_fasta = parse_fasta.get_all_sequences(genome_fasta_file, 'fasta')
    
    # read positions of interest
    tsv_reader = csv.reader(pos_of_interest_file, delimiter='\t')
    for row in tsv_reader:
        if len(row) < 2: continue
        
        scaf = row[0]
        pos = int(row[1]) - 1       # converting to 0-based numbering
        
        min_loc = max(0, pos - window)
        max_loc = min(len(genome_fasta[scaf]), pos + window + 1)
        window_sequence = genome_fasta[scaf][min_loc:max_loc]
        
        # if position is too close to the start/end of the scaffold, add padding
        # Ns in order to produce a sequence of (2 x WINDOW + 1) in length.
        front_n_needed = max(0, window - pos)
        back_n_needed = max(0, 
            window - (len(genome_fasta[scaf]) - (pos + 1)))
        
        # build the sequence of the window
        if front_n_needed:
            window_sequence = 'N' * front_n_needed + window_sequence
        
        if back_n_needed:
            window_sequence += 'N' * back_n_needed
        
        annot = scaf + '_' + str(pos+1)
        window_seqs[annot] = window_sequence
    
    return window_seqs

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    Parses the genome for contextual sequences surrounding positions of interest.
    Padding Ns are added if the positions are incredibly close to the start/end
    of scaffolds.""")

    parser.add_argument('genome_fasta', metavar="fasta_file",
                        type=argparse.FileType('r'), 
                        help="FASTA file of the genome.")
    parser.add_argument('pos_of_interest', metavar="tsv_file",
                        type=argparse.FileType('r'), 
                        help="tab-separated file containing 1-based coords.")
    parser.add_argument('--window', '-w', metavar="n", type=int, default=100,
                        help='''selects window of (-n, n) bases around position 
                                of interest.''')
    args = parser.parse_args()
    
    window_sequences = get_genomic_context(args.genome_fasta, 
                                           args.pos_of_interest, args.window)
    for annot in window_sequences:
        print ('>', annot, sep='')
        print (window_sequences[annot])

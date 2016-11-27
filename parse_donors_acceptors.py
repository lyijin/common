#!/usr/bin/env python3

"""
> parse_donors_acceptors.py <

When given a genome and gff3 file, parse the genome for *unique* exon/intron
splice sites, and parse -10 bp to +10 bp around the boundaries.
"""
import argparse
import collections
import random
import sys

import parse_fasta
import parse_gff3

def reverse_complement(seq):
    seq = seq.replace('U', 'T')

    translation_from = 'AaTtGgCcYyRrSsWwKkMmBbDdHhVvNn'
    translation_to   = 'TtAaCcGgRrYySsWwMmKkVvHhDdBbNn'
    translation_table = str.maketrans(translation_from, translation_to)

    seq = seq[::-1].translate(translation_table)

    return seq

def slice_window(seq, loc, window_size=10):
    """
    In a zero-based genome numbering system, loc is defined to be:
        E E E E I I I        # E = exon, I = intron
        0 1 2 3 4 5 6        # 0-based numbering
       0 1 2 3 4 5 6 7       # loc = 4
            <--|-->          # window_size = 2
    
    Use a -ve loc value to get the reverse-complement sequence of the window.
    """
    if loc == 0:
        print ('Sequence windows cannot be centered around location 0.',
               file=sys.stderr)
        return ''
    
    if loc < 0:
        loc = abs(loc)
        rev_comp_flag = True
    else:
        rev_comp_flag = False
    
    start = max(loc - window_size, 0)
    end = min(loc + window_size, len(seq))
    
    # add Ns as padding to produce sequence lengths of window_size
    upstream_window = 'N' * (window_size - loc) + seq[start:loc]
    downstream_window = seq[loc:end] + 'N' * (loc + window_size - len(seq))
    
    whole_window = upstream_window + downstream_window
    if rev_comp_flag: whole_window = reverse_complement(whole_window)
    
    return whole_window

def get_donor_acceptor_sequences(genome_fasta_file, scaffold_gff3_file, 
                                 window_size=10):
    """
    When given a FASTA file and GFF3 file corresponding to the same genome,
    parse out and return (in a giant list of tuples) all the donor and acceptor 
    sequences of all introns in the genome.
    """
    genome_fasta = parse_fasta.get_all_sequences(genome_fasta_file, 'fasta')
    scaffold_gff3 = parse_gff3.parse_gff3(scaffold_gff3_file, 'exon')
    
    # create dictionary to store exon/intron coordinates
    #   donor_acceptor_locs['scaffold_name'] = [tuple of ints]
    donor_acceptor_locs = {}
    for scaf in scaffold_gff3:
        donor_acceptor_locs[scaf] = []
        
    # get donor and acceptor locations
    for scaf in scaffold_gff3:
        for gene in scaffold_gff3[scaf]:
            for tx in scaffold_gff3[scaf][gene].mRNAs:
                # get all exon coordinates within the transcripts
                tx_coords = scaffold_gff3[scaf][gene].mRNAs[tx].details['exon']
                
                # check whether the transcript is reverse complemented
                rev_comp_flag = tx_coords[0][0] > tx_coords[0][1]
                
                if rev_comp_flag:
                    donor_locs = [-y for x, y in tx_coords][:-1]
                    acceptor_locs = [-x for x, y in tx_coords][1:]
                else:
                    donor_locs = [y for x, y in tx_coords][:-1]
                    acceptor_locs = [x for x, y in tx_coords][1:]
                    
                for x in zip(donor_locs, acceptor_locs):
                    donor_acceptor_locs[scaf].append(x)
        
        donor_acceptor_locs[scaf] = sorted(list(set(donor_acceptor_locs[scaf])))

    # parse the locations into sequences
    donor_acceptor_sequences = []
    for scaf in scaffold_gff3:
        for da in donor_acceptor_locs[scaf]:
            donor_seq = slice_window(genome_fasta[scaf], da[0])
            acceptor_seq = slice_window(genome_fasta[scaf], da[1])
            if donor_seq and acceptor_seq:
                donor_acceptor_sequences.append((donor_seq, acceptor_seq))
                 
    return donor_acceptor_sequences

def generate_null_distrib(donor_acceptor_sequences):
    """
    Given exons and introns of a certain base distribution at every position,
    is the donor-acceptor overlap (see next def) actually... real?
    
    Function shuffles all sequences while preserving per-base composition.
    """
    # create a dict to store randomly permuted sequences
    shuffled_base = {}
    
    concat_seqs = [d + a for d, a in donor_acceptor_sequences]
    d_length = len(donor_acceptor_sequences[0][0])
    da_length = len(concat_seqs[0])
    
    for n in range(da_length):
        bases_at_n = [concat_seqs[x][n] for x in range(len(concat_seqs))]
        random.shuffle(bases_at_n)
        
        shuffled_base[n] = bases_at_n
    
    shuffled_sequences = []
    for n in range(len(bases_at_n)):
        shuffled_donor = ''.join([shuffled_base[x][n] for x in range(d_length)])
        shuffled_acceptor = ''.join([shuffled_base[x][n] 
                                     for x in range(d_length, da_length)])
        shuffled_sequences.append((shuffled_donor, shuffled_acceptor))
    
    return shuffled_sequences

def check_donor_acceptor_overlap(donor_acceptor_sequences):
    """
    Function inspired by Mendez et al., 2015. They noticed that the
    dinoflagellate Crypthecodinium cohnii has splice sites that behave this way:
    
    Donor     cacacatcttgccagGCTGCCATCGTCGAG
                           |||||||             matching region = 7 bp
    Acceptor  CCTGCGACCATGAAGgctgcgctcccttcc
    
    Sequences around donor and acceptor sites are eerily similar! Maybe this is
    how dinoflagellates recognise and cut out introns.
    
    Function takes in tuples of sequences, returns # of matching bases around
    the splice site. This function assumes that the splice sites are dead
    middle of each sequence.
    """
    match_sizes = []
    for da in donor_acceptor_sequences:
        seq_len = len(da[0])
        half_seq_len = int(len(da[0]) / 2)
        d_match_a = [da[0][n] == da[1][n] for n in range(seq_len)]
        
        match_size = 0
        # start at len(seq)/2 - 1 -th base, and check upstream
        for n in range(half_seq_len - 1, -1, -1):
            if d_match_a[n]:
                match_size += 1
            else:
                # once a non-matching base is encountered (i.e. False), get out
                # of the for loop
                break
    
        # start at len(seq)/2 -th base, and check upstream
        for n in range(half_seq_len, seq_len):
            if d_match_a[n]:
                match_size += 1
            else:
                # once a non-matching base is encountered (i.e. False), get out
                # of the for loop
                break
    
        match_sizes.append(match_size)

    return match_sizes

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    When given a genome and gff3 file, parse the genome for *unique* exon/intron
    splice sites, and parse -10 bp to +10 bp around the boundaries.""")

    parser.add_argument('genome_fasta', metavar="fasta_file",
                        type=argparse.FileType('r'), 
                        help="FASTA file of the genome.")
    parser.add_argument('scaffold_gff3', metavar="gff3_file",
                        type=argparse.FileType('r'), 
                        help="corresponding gff3 file of the genome.")
    parser.add_argument('--window', '-w', metavar="n", type=int, default=10,
                        help="select window of (-n, n) bases around splice site.")
    
    args = parser.parse_args()

    splice_site_sequences = get_donor_acceptor_sequences(
        args.genome_fasta, args.scaffold_gff3, args.window)
    overlaps = check_donor_acceptor_overlap(splice_site_sequences)
    
    shuffled_ss_sequences = generate_null_distrib(splice_site_sequences)
    overlaps_shuffled = check_donor_acceptor_overlap(shuffled_ss_sequences)
    
    import collections
    x = collections.Counter(overlaps_shuffled)
    
    # for n in sorted(x):
        # print (n, x[n], sep='\t')
    
    #raise SystemExit()
    # remove empty strings
    # donor_sequences = filter(None, donor_sequences)
    # acceptor_sequences = filter(None, acceptor_sequences)

    # # print table of results
    # print ('Donors:')
    # all_possible_donor_bases = set(''.join([x for x in donor_sequences]))
    # for base in sorted(all_possible_donor_bases):
        # output = [base]
        # for pos in range(0, len(donor_sequences[0])):
            # bases_at_curr_pos = collections.Counter([x[pos] for x in donor_sequences])
            # output.append(bases_at_curr_pos[base] / sum(bases_at_curr_pos.values()))
        
        # output = [str(x) for x in output]
        # print ('\t'.join(output))

    # print ('Acceptors:')
    # all_possible_acceptor_bases = set(''.join([x for x in acceptor_sequences]))
    # for base in sorted(all_possible_acceptor_bases):
        # output = [base]
        # for pos in range(0, len(acceptor_sequences[0])):
            # bases_at_curr_pos = collections.Counter([x[pos] for x in acceptor_sequences])
            # output.append(bases_at_curr_pos[base] / sum(bases_at_curr_pos.values()))
        
        # output = [str(x) for x in output]
        # print ('\t'.join(output))
        
    # output all donor / acceptor sequences as separate fasta files for WebLogo.
    with open('donors.fa', 'w') as f:
        for n, seq in enumerate(splice_site_sequences):
            print ('>' + str(n), file=f)
            print (seq[0], file=f)

    with open('acceptors.fa', 'w') as f:
        for n, seq in enumerate(splice_site_sequences):
            print ('>' + str(n), file=f)
            print (seq[1], file=f)

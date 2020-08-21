#!/usr/bin/env python3

"""
> calc_primer3_stats.py <

When given a bunch of primer sequences as input, use PyPI package primer3-py
to calculate Tm and deltaG for the input sequences.

By default, assumes every provided primer is independent of everything else;
use --paired to pair first primer with the second, third with the fourth...,
use --all to pair all forwards ("_F") against all reverses ("_R"), and 
use --allall to pair everything against everything else, regardless of
annotation.

Primers are designed with qPCR conditions in mind.
  1. dv_conc = 4 (divalent cation concentration, mM; package default = 0)
  2. dna_conc = 200 (oligo concentration, nM; package default = 50)

otherwise, package default values are used.
  3. mv_conc = 50 (monovalent cation concentration, mM)
  4. dntp_conc = 0.8 (dNTP concentration, mM)

Output units
  1. Tm (in celsius)
  2. dG (in kcal / mol)

Use --noheader if there is no header in input file.

Input requirements: one primer per line, with format
    <primer_name> \t <sequence> \t <anything_else> ...

Primer3 does not normally handle degenerate IUPAC bases (e.g. R/Y/N), but script
generates all combinations of the primers and passes it onto Primer3.

WARNING: Primer3 does not like primers > 36 bp. Single and paired mode will
print those primers out without Tms, but all and allall mode will silently
ignore those primers.
"""
import argparse
import csv
import itertools
from pathlib import Path

import primer3

def fmt_tm(tm):
    if min(tm) == max(tm):
        # if min == max, return a single temp value
        return f'{min(tm):.2f}'
    else:
        return f'{min(tm):.2f}-{max(tm):.2f}'

def fmt_dg(dg):
    # just return lowest possible dG (most stable structure)
    try:
        # handles lists
        return f'{min(dg)/1000:.2f}'
    except TypeError:
        # handles floats/ints
        return f'{dg/1000:.2f}'

def min_dg(arg1, *args):
    if not args:
        return fmt_dg(min(arg1))
    else:
        return fmt_dg(min(arg1, *args))

def get_all_degenerate_possibilities(iupac_sequence):
    """
    Builds a list of sequences from a single IUPAC inputs e.g.
      'ACGY' --> ['ACGC', 'ACGT']
    
    Note that N is converted to 'ACGT'; it does not stay as 'N'.
    """
    iupac_nt = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
                'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
                'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
                'N': 'ACGT'}
    
    # this converts input string into list of lists
    #   e.g. 'ACGY' --> [['A'], ['C'], ['G'], ['C', 'T']]
    # list() is the secret sauce that breaks 'CT' into ['C', 'T']
    iupac_sequence = [list(iupac_nt[x]) for x in iupac_sequence]
    
    # itertools.product then joins them up (Cartesian product)
    all_sequences = [''.join(x) for x in itertools.product(*iupac_sequence)]
    
    return all_sequences

def calc_endstab(left_seq, right_seq, mv, dv, d, n):
    """
    End stability calculations ARE position-sensitive, i.e. 
        f(left_seq, right_seq) != f(right_seq, left_seq)
    
    Hence the permutation, saving results into a dict, and returning the
    more stable deltaG of f(seq1, seq1) vs. f(seq1, seq2).
    """
    # deconvolute sequence, as it might contain degenerate bases
    all_left_seq = get_all_degenerate_possibilities(left_seq)
    all_right_seq = get_all_degenerate_possibilities(right_seq)
    
    temp = {}
    for seq1 in all_left_seq + all_right_seq:
        temp[seq1] = {'max_tm': -999.99, 'min_dg': 9999999.99}
        for seq2 in all_left_seq + all_right_seq:
            e = primer3.bindings.calcEndStability(
                seq1, seq2, mv_conc=mv, dv_conc=dv, dntp_conc=n, dna_conc=d)
            temp[seq1]['max_tm'] = max(temp[seq1]['max_tm'], e.tm)
            temp[seq1]['min_dg'] = min(temp[seq1]['min_dg'], e.dg)
    
    # calculate per-primer max Tm and min dG
    endstab = {}
    endstab[left_seq] = {'tm': max(temp[x]['max_tm'] for x in all_left_seq),
                         'dg': min(temp[x]['min_dg'] for x in all_left_seq)}
    endstab[right_seq] = {'tm': max(temp[x]['max_tm'] for x in all_right_seq),
                          'dg': min(temp[x]['min_dg'] for x in all_right_seq)}
    
    return endstab

def calc_primer_stats(left_seq, right_seq=None, oneliner=False, 
                      mv=50, dv=4, d=200, n=0.8):
    """
    When given a single or a pair of primer sequences, compute annealing Tm and
    hairpin/homodimer/3'-end stability dGs.
    
    For paired primers, also calculate heterodimer dG.
    """
    # deconvolute sequence, as it might contain degenerate bases
    all_left_seq = get_all_degenerate_possibilities(left_seq)
    
    left_tm = [primer3.calcTm(x, mv_conc=mv, dv_conc=dv, 
                              dntp_conc=n, dna_conc=d)
               for x in all_left_seq]
    left_hairpin = [primer3.calcHairpin(x, mv_conc=mv, dv_conc=dv,
                                        dntp_conc=n, dna_conc=d).dg
                    for x in all_left_seq]
    left_homodimer = [primer3.calcHomodimer(x, mv_conc=mv, dv_conc=dv,
                                            dntp_conc=n, dna_conc=d).dg
                      for x in all_left_seq]
    
    if right_seq:
        # deconvolute sequence, as it might contain degenerate bases
        all_right_seq = get_all_degenerate_possibilities(right_seq)
        
        right_tm = [primer3.calcTm(x, mv_conc=mv, dv_conc=dv, 
                                  dntp_conc=n, dna_conc=d)
                   for x in all_right_seq]
        right_hairpin = [primer3.calcHairpin(x, mv_conc=mv, dv_conc=dv,
                                            dntp_conc=n, dna_conc=d).dg
                        for x in all_right_seq]
        right_homodimer = [primer3.calcHomodimer(x, mv_conc=mv, dv_conc=dv,
                                                dntp_conc=n, dna_conc=d).dg
                          for x in all_right_seq]
        
        heterodimer = [primer3.calcHeterodimer(x, y, mv_conc=mv, dv_conc=dv,
                                               dntp_conc=n, dna_conc=d).dg
                       for y in all_right_seq for x in all_left_seq]
        endstab = calc_endstab(left_seq, right_seq, mv, dv, d, n)
    else:
        endstab = calc_endstab(left_seq, left_seq, mv, dv, d, n)
    
    # returns section
    if right_seq:
        # two primers
        if oneliner:
            # return single tuple
            return (left_seq, 
                    right_seq, 
                    fmt_tm(left_tm), fmt_tm(right_tm),
                    fmt_dg(left_hairpin + right_hairpin),
                    fmt_dg(left_homodimer + right_homodimer),
                    fmt_dg(heterodimer),
                    fmt_dg([endstab[x]['dg'] for x in endstab]))
        else:
            # return list of two tuples
            return [(left_seq, fmt_tm(left_tm),
                     fmt_dg(left_hairpin), fmt_dg(left_homodimer),
                     fmt_dg(heterodimer), fmt_dg(endstab[left_seq]['dg'])),
                    (right_seq, fmt_tm(right_tm),
                     fmt_dg(right_hairpin), fmt_dg(right_homodimer),
                     fmt_dg(heterodimer), fmt_dg(endstab[right_seq]['dg']))]
    else:
        # single primer
        return (left_seq, fmt_tm(left_tm), fmt_dg(left_hairpin),
                fmt_dg(left_homodimer), fmt_dg(endstab[left_seq]['dg']))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    When given a bunch of primer sequences as input, use PyPI package primer3-py
    to calculate Tm and deltaG for the input sequences.
    
    By default, assumes every provided primer is independent of everything else;
    use --paired to pair first primer with the second, third with the fourth...,
    use --all to pair all forwards ("_F") against all reverses ("_R"), and 
    use --allall to pair everything against everything else, regardless of
    annotation.
    """)
    
    parser.add_argument('amplicon_file', metavar='amplicon_file', type=Path,
                        help='file with one amplicon per line, "<primer_name> \\t <seq>"')
    parser.set_defaults(pair_mode='single')
    pair_opt = parser.add_mutually_exclusive_group(required=False)
    pair_opt.add_argument('--paired', action='store_const',
                          dest='pair_mode', const='paired',
                          help='pair 1st with 2nd, 3rd with 4th, ...')
    pair_opt.add_argument('--all', action='store_const',
                          dest='pair_mode', const='all',
                          help='pair all forwards against all reverses')
    pair_opt.add_argument('--allall', action='store_const', 
                          dest='pair_mode', const='allall',
                          help='pair everything against each other')
    parser.add_argument('-o', '--oneliner', action='store_true',
                        help='print --paired output as single line, not two')
    parser.add_argument('--mv', type=float, default=50,
                        help='monovalent cation concentration, mM (default: 50)')
    parser.add_argument('--dv', type=float, default=4,
                        help='divalent cation concentration, mM (default: 4)')
    parser.add_argument('-d', type=float, default=200,
                        help='oligo concentration, nM (default: 200)')
    parser.add_argument('-n', type=float, default=0.8,
                        help='dNTP concentration, mM (default: 0.8)')
    parser.add_argument('--noheader', action='store_true',
                        help='input file has no header')
    
    args = parser.parse_args()
    
    tsv_reader = csv.reader(open(args.amplicon_file), delimiter='\t')
    # if there's a header, consume the first line
    if not args.noheader:
        header = next(tsv_reader)
    
    if args.pair_mode == 'single':
        # print header
        print ('annot', 'seq', 'Tm', 'hairpin_dG', 'homodimer_dG', 
               'endstab_dG', sep='\t')
        
        for row in tsv_reader:
            annot, seq = row[:2]
            
            # Primer3 refuses to generate Tms for primers > 36 bp. if that happens,
            # print original line out with no additional information
            if len(seq) > 36:
                print (annot, seq, sep='\t')
                continue
            
            # checks passed, grab primer stats
            output = calc_primer_stats(
                seq, mv=args.mv, dv=args.dv, d=args.d, n=args.n)
            print (annot, *output, sep='\t')
    
    elif args.pair_mode == 'paired':
        # print header
        if args.oneliner:
            print ('left_annot', 'left_seq', 'right_annot', 'right_seq',
                   'left_Tm', 'right_Tm', 'min_hairpin_dG', 'min_homodimer_dG',
                   'heterodimer_dG', 'min_endstab_dG', sep='\t')
        else:
            print ('annot', 'seq', 'Tm', 'hairpin_dG', 'homodimer_dG',
                   'heterodimer_dG', 'endstab_dG', sep='\t')
        
        # store contents of the first of each pair in a variable
        first_line = []
        for row in tsv_reader:
            if not first_line:
                first_line = row[:2]
                continue
            
            second_line = row[:2]
            left_annot, left_seq = first_line
            right_annot, right_seq = second_line
            
            # Primer3 refuses to generate Tms for primers > 36 bp. if that happens
            # for either primer, print them out with no additional information
            if len(left_seq) > 36 or len(right_seq) > 36:
                if args.oneliner:
                    print (*first_line, *second_line, sep='\t')
                else:
                    print (*first_line, sep='\t')
                    print (*second_line, sep='\t')
            else:
                # checks passed, grab primer stats
                output = calc_primer_stats(
                    left_seq, right_seq, oneliner=args.oneliner,
                    mv=args.mv, dv=args.dv, d=args.d, n=args.n)
                
                if args.oneliner:
                    print (left_annot, output[0], right_annot, *output[1:], sep='\t')
                else:
                    print (left_annot, *output[0], sep='\t')
                    print (right_annot, *output[1], sep='\t')
            
            # reset first line details
            first_line = []
    
    else:
        # handle --all and --allall similarly
        
        # print fixed header
        print ('left_annot', 'left_seq', 'right_annot', 'right_seq',
               'left_Tm', 'right_Tm', 'min_hairpin_dG', 'min_homodimer_dG',
               'heterodimer_dG', 'min_endstab_dG', sep='\t')
        
        # start by populating the empty dicts
        left_primers = {}
        right_primers = {}
        for row in tsv_reader:
            annot, seq = row[:2]
            
            # silently ignore primers > 36 bp
            if len(seq) > 36: continue
            
            # depending on mode, populate dicts differently
            if args.pair_mode == 'all':
                if 'F' in annot.split('_'):
                    left_primers[annot] = seq
                elif 'R' in annot.split('_'):
                    right_primers[annot] = seq
            
            else:
                left_primers[annot] = seq
                right_primers[annot] = seq
        
        for lp in left_primers:
            for rp in right_primers:
                # don't compare a primer to itself
                if lp == rp: continue
                
                # checks passed, grab primer stats
                output = calc_primer_stats(
                    left_primers[lp], right_primers[rp], oneliner=True,
                    mv=args.mv, dv=args.dv, d=args.d, n=args.n)
                print (lp, output[0], rp, *output[1:], sep='\t')

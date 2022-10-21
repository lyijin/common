#!/usr/bin/env python3

docstring = """
> convert_primers_normal_to_meth.py <

Converts primer sequences that work on normal, unconverted DNA to work on
bisulphite-/enzymatic-converted DNA.

Forward primers are C-to-T converted; reverse primers are G-to-A converted.
Assumes that forward primers have "_F" suffix; reverse have "_R" suffix.

By default, designs primers for post-conversion-top-strand version of the DNA;
use "-b" to override this.

Use the "-f" and "-r" flags to override the auto-detection, treating all primers
as forwards and reverses respectively.

Assumes that C can only be methylated in the CpG context. Tolerates wobble bases
in the input. If CpGs are present in input primers, output primers that could
amplify both the methylated and unmethylated DNA post-conversion, e.g.
  CG --> YG  (Y = C or T)

Remember that the converted primer sequences would not have the same properties
(Tm, GC%) as the unconverted primers!
""".strip()

import argparse
import csv
import itertools
from pathlib import Path

nt_to_iupac = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
               'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
               'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
               'N': 'ACGT'}
iupac_to_nt = {y:x for x,y in nt_to_iupac.items()}

def convert_c_to_t(seq):
    """
    Given a possibly degenerate FORWARD primer that targets unconverted DNA,
    return a degenerate primer that amplifies converted DNA.
    """
    all_degen_seq = get_all_degenerate_possibilities(seq)
    
    # to do C-to-T conversion, while being aware that C in CpG can be methylated
    # (C --> C) or unmethylated (C --> T)
    converted_seq = [x.replace('CG', 'XX').replace('C', 'T').replace('XX', 'YG')
                     for x in all_degen_seq]
    condensed_seq = merge_all_degenerate_possibilities(converted_seq)
    
    return condensed_seq

def convert_g_to_a(seq):
    """
    Given a possibly degenerate REVERSE primer that targets unconverted DNA,
    return a degenerate primer that amplifies converted DNA.
    """
    all_degen_seq = get_all_degenerate_possibilities(seq)
    
    # to do G-to-A conversion, while being aware that C in CpG can be methylated
    # (G --> G) or unmethylated (G --> A)
    converted_seq = [x.replace('CG', 'XX').replace('G', 'A').replace('XX', 'CR')
                     for x in all_degen_seq]
    condensed_seq = merge_all_degenerate_possibilities(converted_seq)
    
    return condensed_seq

def get_all_degenerate_possibilities(seq):
    """
    Builds a list of sequences from a single IUPAC-conforming input e.g.
      'ACGY' --> ['ACGC', 'ACGT']
    
    Note that N is converted to 'ACGT'; it does not stay as 'N'.
    """
    # this converts input string into list of lists
    #   e.g. 'ACGY' --> [['A'], ['C'], ['G'], ['C', 'T']]
    # list() is the secret sauce that breaks 'CT' into ['C', 'T']
    seq = [list(nt_to_iupac[x]) for x in seq]
    
    # itertools.product then joins them up (Cartesian product)
    all_sequences = [''.join(x) for x in itertools.product(*seq)]

    return all_sequences

def merge_all_degenerate_possibilities(sequences):
    """
    Basically a reverse of `get_all_degenerate_possibilities()`.
    
    Condenses a list of sequences into a single IUPAC-conforming input e.g.
      ['ACGC', 'ACGT'] --> 'ACGY'
    """
    # make sure all sequences in array have the same length
    assert isinstance(sequences, list), \
        f'ERROR: {sequences} does not appear to be a list!'
    assert len(set(len(x) for x in sequences)) == 1, \
        f'ERROR: some elements of {sequences} are of a different length!'
    
    condensed_seq = ''
    for n in range(len(sequences[0])):
        possible_chars = [x[n] for x in sequences]
        
        # there might be degenerate bases in "possible_chars", e.g. "Y"
        possible_chars = [nt_to_iupac[x] for x in possible_chars]
        
        # condense down to unique items
        possible_chars = ''.join(sorted(list(set(''.join(possible_chars)))))
        condensed_seq += iupac_to_nt[possible_chars]
    
    return condensed_seq

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=docstring, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('amplicon_file', metavar='amplicon_file', type=Path,
        help='file with one amplicon per line, "<primer_name> \\t <seq>"')
    parser.add_argument('-b', '--bottom_strand', action='store_true',
        help='produce primers that target post-conversion-bottom-strand DNA')
    pair_opt = parser.add_mutually_exclusive_group(required=False)
    pair_opt.add_argument('-f', '--all_forward', action='store_true',
        help='treat all primers as forwards')
    pair_opt.add_argument('-r', '--all_reverse', action='store_true',
        help='treat all primers as reverses')
    parser.add_argument('--noheader', action='store_true',
        help='input file has no header')

    args = parser.parse_args()
    
    tsv_reader = csv.reader(open(args.amplicon_file), delimiter='\t')
    # if there's a header, consume the first line
    if not args.noheader:
        header = next(tsv_reader)

    # start by populating the empty dicts
    left_primers = {}
    right_primers = {}
    for row in tsv_reader:
        annot, seq = row[:2]
        seq = seq.replace(' ', '')      # remove spaces in seq, if present
        
        # depending on mode, populate dicts differently
        if not args.all_forward and not args.all_reverse:
            # this is "normal mode", autodetects primer type based on last
            # letter in primer. "F" = forward; "R" = reverse
            
            # to create bottom-strand primers, just simply flip initial "_F"
            # and "_R" annotations around; then C-to-T conversion logic applies
            # to the new F while G-to-A conversion logic applies to new R
            if args.bottom_strand:
                if annot[-1] == 'F':
                    annot = annot[:-1] + 'R'
                elif annot[-1] == 'R':
                    annot = annot[:-1] + 'F'
            
            if annot[-1] == 'F':
                print (annot, seq, convert_c_to_t(seq), sep='\t')
            elif annot[-1] == 'R':
                print (annot, seq, convert_g_to_a(seq), sep='\t')
            else:
                print (annot, seq, 'unsure F or R', '', sep='\t')
        
        elif args.all_forward:
            print (annot, seq, convert_c_to_t(seq), sep='\t')
        
        elif args.all_reverse:
            print (annot, seq, convert_g_to_a(seq), sep='\t')

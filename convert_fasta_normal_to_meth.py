#!/usr/bin/env python3

docstring = """
> convert_fasta_normal_to_meth.py <

Converts four-base FASTA sequences to three-base converted FASTA sequences.
Degenerate bases (non-ACGT) are automatically converted to N, too much of a
hassle to properly convert degenerate bases. Takes FAR too long to produce
conversions that are 100% correct.

For each sequence, two converted sequences are produced:
  OT  : original top (C-to-T)
  CTOB: complementary to original bottom (G-to-A)

I prefer OT/CTOB as they adhere to the same chr:pos values as the original
input sequences.

By default, script assumes conversion of fully methylated sequences (i.e.,
all CGs stay CGs, whilst CHs are converted to THs), but this can be changed
to fully unmethylated (all Cs are converted to Ts).

Is gzip-aware, but output is always plaintext. Pipe into gzip or pigz if
compressed output is desired.
""".strip()

import argparse
import csv
import itertools
from pathlib import Path

import parse_fasta

# define string-insensitive translation tables
degen_to_N_from  = 'AaCcGgTtRrYySsWwKkMmBbDdHhVvNn'
degen_to_N_to    = 'AaCcGgTtNnNnNnNnNnNnNnNnNnNnNn'
degen_to_N_table = str.maketrans(degen_to_N_from, degen_to_N_to)

revcomp_from = 'AaTtGgCcYyRrSsWwKkMmBbDdHhVvNn'
revcomp_to   = 'TtAaCcGgRrYySsWwMmKkVvHhDdBbNn'
revcomp_table = str.maketrans(revcomp_from, revcomp_to)

unmeth_1base_from = 'AaCcGgTt'
unmeth_1base_to   = 'AaTtGgTt'
unmeth_1base_table = str.maketrans(unmeth_1base_from, unmeth_1base_to)

def convert_c_to_t(seq, seq_type='fully_meth'):
    """
    Given a potentially degenerate four-base sequence, C-to-T convert it to a 
    non-degenerate three-base equivalent.
    """
    assert seq_type in ['fully_meth', 'fully_unmeth'], \
        "`seq_type` argument can be one of ['fully_meth', 'fully_unmeth']!"
    
    # crush degeneracy out of sequence (all non-ACGT --> N)
    seq = seq.translate(degen_to_N_table)
    
    if seq_type == 'fully_meth':
        # CG stays CG, CH converted to TH
        conv_seq = seq.replace('CG', 'XX').replace('Cg', 'Xx').replace('cG', 'xX').replace('cg', 'xx')\
                      .translate(unmeth_1base_table)\
                      .replace('XX', 'CG').replace('Xx', 'Cg').replace('xX', 'cG').replace('xx', 'cg')
    elif seq_type == 'fully_unmeth':
        # C-to-T, wholesale conversion
        conv_seq = seq.translate(unmeth_1base_table)
    
    return conv_seq

def reverse_complement(seq):
    """
    Does not deal with U, assume input is DNA.
    """
    return seq[::-1].translate(revcomp_table)

def wrap_sequence(seq, width):
    return '\n'.join(seq[i:i+width] for i in range(0, len(seq), width))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=docstring, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('fasta_file', metavar='fasta_file', type=Path,
        help='(compressed) FASTA file')
    parser.add_argument('-u', '--unmeth', action='store_true',
        help='assume sequences are fully unmethylated')
    parser.add_argument('-w', '--wrap', metavar='n', type=int,
        help='wrap sequences with line width of n')
    args = parser.parse_args()
    
    # read sequences
    sequences = parse_fasta.get_all_sequences(args.fasta_file, 'fasta')
    
    for annot in sequences:
        seq = sequences[annot]
        # OT conversion, the top strand getting C-to-T converted
        ot_conv = convert_c_to_t(seq, 'fully_unmeth' if args.unmeth else 'fully_meth')
        
        # CTOB conversion. logic is, the reverse complement is getting C-to-T
        # converted (this is NOT complementary to OT). the re-reverse complement
        # produces a sequence with the same coordinate system as original & OT.
        ctob_conv = reverse_complement(
            convert_c_to_t(reverse_complement(seq), 'fully_unmeth' if args.unmeth else 'fully_meth'))
        
        # output
        print (f'>{annot}_OT')
        print (wrap_sequence(ot_conv, args.wrap) if args.wrap else ot_conv)
        
        print (f'>{annot}_CTOB')
        print (wrap_sequence(ctob_conv, args.wrap) if args.wrap else ctob_conv)

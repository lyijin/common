#!/usr/bin/env python3

docstring = """
> get_fasta_seqlen.py <

Script takes in a FASTA file, gets the lengths of the sequences, and outputs
this in a one-sequence-per-line format to stdout.

   annotation_1 \t seqlen_1
   annotation_2 \t seqlen_2
   ...
""".strip()
import argparse
from pathlib import Path

import parse_fasta

parser = argparse.ArgumentParser(
    description=docstring, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('genome_fasta', metavar='fasta_file', type=Path,
                    help='genome FASTA file.')

args = parser.parse_args()

# read sequences
sequence_lengths = parse_fasta.get_all_sequences(
    args.genome_fasta, 'fasta', lengths_only=True)

# output
for seq in sequence_lengths:
    print (seq, sequence_lengths[seq], sep='\t')

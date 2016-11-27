#!/usr/bin/env python3

"""
> collapse_reads.py <

Script is inspired by 'collapse_reads_md.pl' provided by the miRDeep2 package.

Takes in FASTA/FASTQ file and collapses the reads. Output is as follows:
>sid_0_x[y]
ATCGATCGATCG
>sid_y_x[z]
AAGTCAGCTAGCTGACT
>sid_y+z_x[a]  ...

where sid = three-letter species identifier, y/z/a.. are reads, sorted in reverse
of their abundances.
"""

import argparse

parser = argparse.ArgumentParser(description="""
Script takes in FASTA/FASTQ file and collapses the reads.""")

parser.add_argument('reads_file', metavar='reads_file', type=argparse.FileType('r'),
                    nargs=1, help='FASTA/FASTQ read file')
parser.add_argument('--species', metavar='species_identifier', type=str,
                    nargs=1, help='three-letter species identifier',
                    required=True)
fasta_opt = parser.add_mutually_exclusive_group(required=True)
fasta_opt.add_argument('--fasta', action='store_const', dest='file_format',
                       const='fasta', help='input file is in FASTA format.')
fasta_opt.add_argument('--fastq', action='store_const', dest='file_format',
                       const='fastq', help='input file is in FASTQ format.')

args = parser.parse_args()

# grab only the sequences - ignore all annotations
import parse_fasta
sequences = parse_fasta.get_all_sequences(args.reads_file[0], args.file_format,
                                          sequences_only=True)

# discard annotations
import collections
sequences = collections.Counter(sequences)

reads_counter = 0       # the y/z/a counter mentioned in script description
species_identifier = args.species[0][:3]
for m in sequences.most_common():
    # m[0] is the sequence; m[1] is the frequency
    print ('>{}_{}_x{}'.format(species_identifier, reads_counter, m[1]))
    print (m[0])
    
    reads_counter += m[1]

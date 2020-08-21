#!/usr/bin/env python3

"""
> reverse_complement.py <

Does what it says on the tin: reverse-complements the input sequence.
Autodetects input file formats of

  1. One-sequence-per-line
  2. Tab-separated: annot \t seq
  3. FASTA: >annot \n seq

Use flags to override detected file format. Files must be in plaintext.
"""

import argparse
import itertools
from pathlib import Path
import sys

def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

def rc_single(seq):
    seq = seq.replace('U', 'T')

    translation_from = 'AaTtGgCcYyRrSsWwKkMmBbDdHhVvNn'
    translation_to   = 'TtAaCcGgRrYySsWwMmKkVvHhDdBbNn'
    translation_table = str.maketrans(translation_from, translation_to)
    
    # if seq contains characters not in `translation_from`, then DO NOT
    # attempt to reverse complement the sequence, just return itself back
    if set(seq) <= set(translation_from):
        return seq[::-1].translate(translation_table)
    else:
        return seq

def rc_oneperline(txt_input):
    output = []
    for line in txt_input.splitlines():
        seq = line.strip()
        
        output.append(rc_single(seq))
    
    return '\n'.join(output)

def rc_tabseparated(tsv_input):
    output = []
    for line in tsv_input.splitlines():
        row = line.strip().split('\t')
        
        if len(row) > 1:
            row[1] = rc_single(row[1])
            output.append('\t'.join(row))
        else:
            # if there isn't a second column, skip conversion for the line
            output.append('\t'.join(row))
    
    return '\n'.join(output)

def rc_fasta(fasta_input):
    all_seqs = {}
    
    annot = ''
    for line, next_line in pairwise(fasta_input.splitlines()):
        if line[:1] == '>':
            annot = line[1:].strip()
            seq = next_line.strip()
        else:
            # ignore comment lines before the first '>'
            if annot:
                # produce output if next line = start of next sequence
                if next_line[:1] == '>':
                    all_seqs[annot] = seq
                else:
                    seq += next_line.strip()
    
    # handle the last sequence in file
    all_seqs[annot] = seq
    
    output = []
    for annot in all_seqs:
        output.append(f'>{annot}')
        output.append(rc_single(all_seqs[annot]))
    
    return '\n'.join(output)

def rc_autodetect(input_string):
    """
    Make an intelligent guess on what the input file format is, and call the
    corresponding function to interpret the string.
    """
    if '\t' in input_string:
        return rc_tabseparated(input_string)
    elif input_string[0] == '>':
        return rc_fasta(input_string)
    else:
        return rc_oneperline(input_string)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    Python script translates FASTA file containing DNA/RNA sequences into
    protein sequences, appending _+1/2/3 or _-1/-2/-3 based on the frame being
    translated. Standard genetic code is assumed.""")

    parser.add_argument('input_file', metavar='plaintext_file',
                        type=Path, default=sys.stdin.fileno(), nargs='?',
                        help='nucleotide sequence-containing plain text file.')
    filetype_opt = parser.add_mutually_exclusive_group(required=False)
    filetype_opt.add_argument('--oneperline', action='store_const',
                          dest='filetype', const='oneperline',
                          help='file contains a single sequence per line')
    filetype_opt.add_argument('--tabseparated', action='store_const',
                          dest='filetype', const='tabseparated',
                          help='file is tab-separated, "annot \\t seq"')
    filetype_opt.add_argument('--fasta', action='store_const',
                          dest='filetype', const='fasta',
                          help='file is FASTA-formatted')
    
    args = parser.parse_args()
    
    input_string = open(args.input_file).read()
    
    if args.filetype == 'tabseparated':
        rc = rc_tabseparated(input_string)
    elif args.filetype == 'fasta':
        rc = rc_fasta(input_string)
    elif args.filetype == 'oneperline':
        rc = rc_oneperline(input_string)
    else:
        rc = rc_autodetect(input_string)
    
    print (rc)

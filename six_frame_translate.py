#!/usr/bin/env python3

"""
> six_frame_translate.py <

Translate FASTA file containing DNA/RNA sequences into protein sequences,
appending _+1/2/3 or _-1/-2/-3 based on the frame being translated.

Assumes standard genetic code.

OPTIONS:
  --longest to get the longest ORF (Met to Stop / end of sequences) in the
  nucleotide sequences.

  --relaxed to relax the requirement of Met to start of sequence /
  non-Stop codon.

  --nt to get nucleotide output instead of protein sequences.
"""

import argparse
import re
import sys

import parse_fasta

# code taken from http://www.petercollingridge.co.uk/blog/python-codon-table
bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

def get_aa(codon):
    # resolves mixed-case scenario
    codon = codon.upper()
    if codon in codon_table:
        return codon_table[codon]
    else:
        return 'X'

def reverse_complement(seq):
    seq = seq.replace('U', 'T')

    translation_from = 'AaTtGgCcYyRrSsWwKkMmBbDdHhVvNn'
    translation_to   = 'TtAaCcGgRrYySsWwMmKkVvHhDdBbNn'
    translation_table = str.maketrans(translation_from, translation_to)

    seq = seq[::-1].translate(translation_table)

    return seq

def peptide_sequence(seq, frame=1):
    # note: 'frame' is in the biological sense. "+1" = no frameshift.
    if frame < 0:
        seq = reverse_complement(seq)
        frame = abs(frame)

    return ''.join([get_aa(seq[start:start + 3])
                    for start in range(frame - 1, len(seq), 3)
                    if start + 3 <= len(seq)])

def nucleotide_sequence(seq, frame=1):
    # note: 'frame' is in the biological sense. "+1" = no frameshift.
    # similar to peptide_sequence, but prints the equivalent, untranslated, seq
    if frame < 0:
        seq = reverse_complement(seq)
        frame = abs(frame)

    return ''.join([seq[start:start + 3]
                    for start in range(frame - 1, len(seq), 3)
                    if start + 3 <= len(seq)])

def six_frame_translate(annot, seq, get_nt_seq=False, display_length=False):
    """Prints out protein FASTA sequences of all six frames."""
    output = ''
    
    for fr in [1, 2, 3, -1, -2, -3]:
        new_annot = '>{}_{}'.format(annot, str(fr))
        if get_nt_seq:
            output_seq = nucleotide_sequence(seq, frame=fr)
        else:
            output_seq = peptide_sequence(seq, frame=fr)

        if display_length:
            new_annot += '_Length_{}'.format(len(output_seq))
        
        # save converted sequences into "output"
        if output:
            output = '\n'.join([output, new_annot, output_seq])
        else:
            output = '\n'.join([new_annot, output_seq])
    
    return output

def find_longest_orf(annot, seq, relaxed=False, display_length=False):
    """Prints out the protein FASTA sequence of the longest ORF only."""    
    longest_in_fr = {}
    for fr in [1, 2, 3, -1, -2, -3]:
        if relaxed:
            all_orfs = re.findall('[^\*]+', peptide_sequence(seq, fr))
        else:
            all_orfs = re.findall('M[^\*]+', peptide_sequence(seq, fr))

        if all_orfs:
            longest_in_fr[fr] = max(all_orfs, key=len)
    
    try:
        longest_overall = max(longest_in_fr, key=lambda x: len(longest_in_fr[x]))
    except ValueError:
        longest_overall = 0     # no ORF found in all frames
    
    if longest_overall != 0:     # if no ORFs found, skip the printing
        new_annot = '>{}_{}'.format(annot, longest_overall)
        if display_length:
            new_annot += '_Length_{}'.format(len(longest_in_fr[longest_overall]))

        output = '{}\n{}'.format(new_annot, longest_in_fr[longest_overall])
    
        return output
    else:
        return '>{}_NO_ORFs_FOUND'.format(annot)

def natural_sort(input_list):
    tryint = lambda x: int(x) if x.isdigit() else x
    chunked_text = lambda x: [tryint(y) for y in re.split('([0-9]+)', x)]
    sorted_list = sorted(input_list, key=chunked_text)

    return sorted_list

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    Python script translates FASTA file containing DNA/RNA sequences into
    protein sequences, appending _+1/2/3 or _-1/-2/-3 based on the frame being
    translated. Standard genetic code is assumed.""")

    parser.add_argument('infile', metavar="nucleotide_fasta_file",
                        type=argparse.FileType('r'), nargs='?',
                        default=sys.stdin, help="nucleotide FASTA file.")
    parser.add_argument('--longest', action='store_true',
                        help='find the longest ORF among the six frames.')
    parser.add_argument('--relaxed', action='store_true',
                        help='ORFs need not start with Met.')
    parser.add_argument('--nt', action='store_true',
                        help='print equivalent nucleotide sequences.')
    parser.add_argument('--print_length', action='store_true',
                        help='include the length of the sequence in the annot.')
    parser.add_argument('--nosort', action='store_true',
                        help='disable natural sorting on output.')

    args = parser.parse_args()

    fasta_seqs = parse_fasta.get_all_sequences(args.infile, 'fasta')

    if args.nosort:
        sorted_seqs = fasta_seqs
    else:
        sorted_seqs = natural_sort(fasta_seqs)

    if args.longest:
        for s in sorted_seqs:
            print (find_longest_orf(s, fasta_seqs[s], relaxed=args.relaxed, display_length=args.print_length))
    else:
        for s in sorted_seqs:
            print (six_frame_translate(s, fasta_seqs[s], get_nt_seq=args.nt, display_length=args.print_length))

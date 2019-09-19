#!/usr/bin/env python3

"""
> generate_expected_base_freq.amplicon_seq.py <

Given a flat file containing all possible amplicons sequences (one per line),
script generates an expected per-base frequency table (which should correspond
to the frequencies detected by the next gen sequencer). Assumes distinct
amplicons would be present in equal proportions.

Input accepts all IUPAC bases, e.g. B, D, H, V, N.
"""
import argparse

parser = argparse.ArgumentParser(description="""
Given a flat file containing all possible amplicons sequences (one per line),
script generates an expected per-base frequency table (which should correspond
to the frequencies detected by the next gen sequencer). Assumes distinct
amplicons would be present in equal proportions.""")

parser.add_argument('amplicon_seqs', metavar='seq_file',
                    type=argparse.FileType('r'),
                    help='sequence of all amplicons (one per line).')
parser.add_argument('-t', '--transpose', action='store_true',
                    help='transpose output (prints horizontally, instead '
                         'of vertically.')
args = parser.parse_args()

IUPAC_BASES = {'-': '', '.': '',
               'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
               'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
               'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
               'N': 'ACGT'}

def calc_base_proportions(bases):
    # convert every position as a 12-mer and concatenate them, e.g.
    #   'A' --> 'AAAAAAAAAAAA'
    #   'N' --> 'ACGTACGTACGT'
    b = [IUPAC_BASES[x] * int(12/len(IUPAC_BASES[x])) if x not in '-.' else ''
         for x in bases]
    b = ''.join(b)
    
    return [f'{b.count(x) / len(b):.2f}' for x in 'ACGT']

amplicon_seqs = []
with args.amplicon_seqs as f:
    for line in f:
        amplicon_seqs.append(line.strip())

# sort amplicons, from shortest to longest
amplicon_seqs = sorted(amplicon_seqs, key=len)

# generate and print output
print ('\t' * len(amplicon_seqs), 'A\tC\tG\tT')
for n in range(len(amplicon_seqs[-1])):
    pos_n_bases = [x[n] if n < len(x) else '-' for x in amplicon_seqs]
    print (*pos_n_bases, *calc_base_proportions(pos_n_bases), sep='\t')

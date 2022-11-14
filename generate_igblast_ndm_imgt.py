#!/usr/bin/env python3

docstring = """
> generate_igblast_ndm_imgt.py <

To run IgBLAST (1.20+) on custom organisms, there are two important flat files
that needs to be generated. This is so that IgBLAST can correctly delineate
where CDR1/2/3 starts and ends.

The first file, for IGHV sequences, is for CDR1/2 start and ends.
The second file, for IGHJ sequences, is for CDR3 ends.

This script generates the first file for custom organisms. To see what an
actual "first file" looks like on a known species, go to
  /path/you/installed/ncbi-igblast-1.20.0/internal_data/human/human.ndm.imgt

This script deduces the coordinates of FWR1 start/end, CDR1 start/end, ...
based on a sequence file downloaded from
  https://www.imgt.org/vquest/refseqh.html

Download the IGHV "F+ORF+in-frame P with IMGT gaps" nucleotides file for your
organism, and feed it through this script. To reiterate:
  1. IGHV (not D/J/C/whatever, and IG, not TR)
  2. **WITH IMGT GAPS** (these files have dots "..." in them), not ungapped
  3. NUCLEOTIDES, not PROTEINS

(which is why this script ends up producing text in the `*.ndm.imgt` format.)

Additional IMGT nomenclature information can be found at
https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
""".strip()

import argparse
from pathlib import Path
import re

import parse_fasta

def imgt_annot_to_igblast(orig_annot):
    # IMGT has annots in the form of
    # >M99641|IGHV1-18*01|Homo sapiens|F|V-REGION|188..483|296 nt|1| | | | |296+24=320| | |
    # IgBLAST only wants "IGHV1-18*01"
    return orig_annot.split('|')[1]
    
parser = argparse.ArgumentParser(
    description=docstring, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('imgt_ighv_file', metavar='imgt_ighv_file', type=Path,
    help='IMGT IGHV nucleotide file with gaps.')

args = parser.parse_args()

# reads all sequences from both files into memory
imgt_seqs = parse_fasta.get_all_sequences(args.imgt_ighv_file, 'fasta')

# print standard header plagiarised from `human.ndm.imgt`
print ('#gene/allele name, FWR1 start, FWR1 stop, CDR1 start, CDR1 stop, '
       'FWR2 start, FWR2 stop, CDR2 start, CDR2 stop, FWR3 start, FWR3 stop, '
       'chain type, coding frame start.')
print ('#FWR/CDR positions are 1-based while the coding frame start positions '
       'are 0-based')

for imgt_annot in imgt_seqs:
    seq = imgt_seqs[imgt_annot]
    # looking at `human.ndm.imgt`, the file does not contain coordinates
    # for sequences shorter than 312 bp (including gaps). this behaviour
    # is mimicked by skipping short sequences here
    if len(seq) < 312: continue
    
    igblast_annot = imgt_annot_to_igblast(imgt_annot)
    
    # this is WHY the gapped file is used--sequences in this file are already
    # aligned as per the IMGT numerotation concept:
    #   FWR1 starts from aa position  1 to  26, nt position   1 to  78
    #   CDR1 starts from aa position 27 to  38, nt position  79 to 114
    #   FWR2 starts from aa position 39 to  55, nt position 115 to 165
    #   CDR2 starts from aa position 56 to  65, nt position 166 to 195
    #   FWR3 starts from aa position 66 to 104, nt position 196 to 312
    #   (there's no CDR3 in IGHV)
    #
    # to figure out corresponding coordinates for the UNGAPPED file (which is
    # what IgBLAST wants), we basically just need to figure out how many .
    # are in each section, and subtract accordingly
    
    fwr1_start =   1
    fwr1_end   =  78 - seq[: 78].count('.')
    cdr1_start =  79 - seq[: 79].count('.')
    cdr1_end   = 114 - seq[:114].count('.')
    
    fwr2_start = 115 - seq[:115].count('.')
    fwr2_end   = 165 - seq[:165].count('.')
    cdr2_start = 166 - seq[:166].count('.')
    cdr2_end   = 195 - seq[:195].count('.')
    
    fwr3_start = 196 - seq[:196].count('.')
    fwr3_end   = 312 - seq[:312].count('.')
    
    chain_type = 'VH'
    coding_fs  = 0
    
    # some sequences have a lot of '.' at the start/end. this is best judged 
    # using 'fwr3_end'--in the `human.ndm.imgt` file, fwr3_end has a minimum 
    # value of 285. use this as a gauge for "real" sequences
    if fwr3_end < 280: continue
    
    print (igblast_annot, fwr1_start, fwr1_end, cdr1_start, cdr1_end,
           fwr2_start, fwr2_end, cdr2_start, cdr2_end,
           fwr3_start, fwr3_end, chain_type, coding_fs, sep='\t')

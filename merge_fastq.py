#!/usr/bin/env python3

"""
> merge_fastq.py <

Python script merges FASTQ files to get rid of duplicate sequences, but 
makes the HUGE assumption that sequences with the same annotation contain
the same sequence + quality.

The script is reasonably fast based on the fact that dictionary keys are unique.

It also assumes that the FASTQ file conforms to expectations (4 lines per
sequence, no weird newlines anywhere).

To make sure merged R1s and merged R2s are correctly paired, the keys are
sorted before printing it out (natural sorting is much slower).
"""
import argparse
import sys

parser = argparse.ArgumentParser(description="""
Python script merges FASTQ files to get rid of duplicate sequences, but 
makes the HUGE assumption that sequences with the same annotation contain
the same sequence + quality.""")
parser.add_argument('fastq_files', metavar='fastq_files',
                    type=argparse.FileType('r'), nargs='+',
                    help='input FASTQ files.')

args = parser.parse_args()

all_sequences = {}
for fq in args.fastq_files:
    for line in fq:
        if not line: continue
        
        # check for annotation line
        if line[0] == '@':
            annot = line.strip()
            seq = next(fq).strip()
            qual_annot = next(fq).strip()
            qual = next(fq).strip()
            
            all_sequences[annot] = (seq, qual)

# print all_sequences out as a FASTQ file
for annot in sorted(all_sequences):
    print (annot)
    print (all_sequences[annot][0])
    print ('+')
    print (all_sequences[annot][1])

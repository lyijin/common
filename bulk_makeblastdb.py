#!/usr/bin/env python3

"""
> bulk_makeblastdb.py <

Script to automate the makeblastdb process of (re-)making BLAST databases.
Useful for:

1. compiling new databases
2. changing descriptions of existing databases

To control which databases get recompiled, there's an external hardcoded file
that maps filename --> database type & description. Filename passed into
argparse MUST match these filenames identically.

Passing 'all' as the filename will process ALL files contained in the hardcoded
file. Use with caution.

NOTE: to prevent SequenceServer keeling over, FASTA sequences will have trailing
stop codons ('*') removed and internal stop codons replaced with 'X'.
"""
import argparse
import csv
import os
import re
import subprocess
import tempfile

parser = argparse.ArgumentParser(description='''
Script to automate the makeblastdb process of (re-)making BLAST databases.''')

parser.add_argument('fasta_file', metavar='fasta_file',
                    type=str, nargs='+',
                    help='FASTA files for makeblastdb.')

args = parser.parse_args()

# read hardcoded file
blast_dict = {}
script_abspath = os.path.abspath(__file__)
hardcoded_tsv = script_abspath[:-2] + 'tsv'

tsv_reader = csv.reader(open(hardcoded_tsv), delimiter='\t')
for row in tsv_reader:
    if not row: continue
    if row[0][0] == '#': continue       # ignore commented lines
    
    blast_dict[row[0]] = (row[1], row[2])

if args.fasta_file == ['all']:
    filenames = [x for x in blast_dict]
else:
    filenames = args.fasta_file

for f in filenames:
    if f not in blast_dict:
        print (f, 'not found in hardcoded input file!')
    else:
        # replace/remove stop codons appropriately
        # note: while these steps can be easily accomplished with sed (via
        #       subprocess), there's a chance that malicious filenames might
        #       compromise the system
        sequences = open(f).read()
        sequences = re.sub('\*\n>', '\n>', sequences)
        sequences = re.sub('\*', 'X', sequences)
        
        with tempfile.NamedTemporaryFile() as tf:
            tf.write(sequences.encode())
            tf.seek(0)
            
            subprocess.run(['makeblastdb', '-in', tf.name,
                                           '-out', f,
                                           '-dbtype', blast_dict[f][0],
                                           '-title', blast_dict[f][1], 
                                           '-parse_seqids'])

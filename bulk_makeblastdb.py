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
import gzip
import os
import re
import subprocess
import tempfile

parser = argparse.ArgumentParser(description='''
Script to automate the makeblastdb process of (re-)making BLAST databases.''')

parser.add_argument('database_name', metavar='database_name',
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
    
    blast_dict[row[0]] = (row[1], row[2], row[3])

if args.database_name == ['all']:
    filenames = [x for x in blast_dict]
else:
    filenames = args.database_name

for f in filenames:
    if f not in blast_dict:
        print (f, 'not found in hardcoded input file!')
    else:
        # replace/remove stop codons appropriately
        # note: while these steps can be easily accomplished with sed (via
        #       subprocess), there's a chance that malicious filenames might
        #       compromise the system
        db_filename = blast_dict[f][2]
        
        # crude way to determine whether file is gzipped or not
        if db_filename[-3:] == '.gz':
            sequences = gzip.open(blast_dict[f][2], 'rt').read()
        else:
            sequences = open(blast_dict[f][2]).read()
        
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
        
        print ('SUCCESS: database {} created, with type {} and title "{}"'.format(
            f, blast_dict[f][0], blast_dict[f][1]))

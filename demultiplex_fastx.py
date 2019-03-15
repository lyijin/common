#!/usr/bin/env python3

"""
> demultiplex_fastx.py <

GATK enforces the idea of "read groups" in many of its tools--it's loosely
defined as reads from the same biological sample, sequenced in the same
flowcell.

This script demultiplexes FASTA/FASTQ files that were merged previously by
checking the annotation of the read, then dumping the read in files
corresponding to its read group.

Works on compressed files.
"""
import argparse
import gzip
import os.path
import re
import sys
import time

parser = argparse.ArgumentParser(description="""
GATK enforces the idea of "read groups" in many of its tools--it's loosely
defined as reads from the same biological sample, sequenced in the same
flowcell.""")

parser.add_argument('genome_fastx', metavar='fastx_file',
                    type=argparse.FileType('r'),
                    help='genome FASTX file.')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints diagnostic stuff to stderr.')
args = parser.parse_args()

# script will deduce filetype/compression based on filename! script won't
# work properly if it's misled!
# supported extensions: .fa, .fq, .fasta, .fastq; and all of them + .gz
file_prefix, file_ext = os.path.splitext(args.genome_fastx.name)
if file_ext == '.gz':
    args.genome_fastx = gzip.open(args.genome_fastx.name, 'rt')
    file_prefix, file_ext = os.path.splitext(file_prefix)

if file_ext == '.fasta' or file_ext == '.fa':
    lines_per_seq = 2
elif file_ext == '.fastq' or file_ext == '.fq':
    lines_per_seq = 4

if args.verbose:
    print (f'[{time.asctime()}] Detected file extension: {file_ext}',
           file=sys.stderr)

# important bits from an annotation line (assumes reads from Illumina):
#   @HWI-ST1101:94:C1P8GACXX:1:1101:1233:2142 1:N:0:ACAGTG
#                  =====---- _
#   = : flowcell barcode
#   - : sequencer
#   _ : lane
#
# output files will have the barcode and lane in the filename
all_output_files = {}
output_file = ''
counter = 0
with args.genome_fastx as f:
    for line in f:
        # this if block captures annotation lines. decide which output file
        # the line should go to based on the flowcell/lane information
        if counter % lines_per_seq == 0:
            # figure out flowcell barcode and lane
            temp = line
            if ' ' in temp:
                temp = temp.split(' ')[0]
            
            # parsing goes from the back instead of the front as the back bits
            # are less variable than the front bits
            flowcell = temp.split(':')[-5][:-4]
            lane = temp.split(':')[-4]
            
            # always write gzip-compressed output
            output_file = f'{file_prefix}.{flowcell}_L{lane}{file_ext}.gz'
            if output_file not in all_output_files:
                all_output_files[output_file] = gzip.open(output_file, 'wt')
        
        print (line, end='', file=all_output_files[output_file])
        counter += 1
        
        if args.verbose and counter % 10000000 == 0:
            print (f'[{time.asctime()}] {counter} sequences binned into',
                   f'{len(all_output_files)} separate files...', file=sys.stderr)

# close all open files
for f in all_output_files:
    all_output_files[f].close()

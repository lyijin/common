#!/usr/bin/env python3

docstring = """
> demultiplex_undetermined_fastq.min_hamming.py <

Sometimes, our lab Illumina benchtop sequencer produces poor quality readouts
for the index regions (i5/i7). These reads will not be properly demultiplexed
with the manufacturer-provided `bcl2fastq` tool, as it can--at most--have
two mismatches to what's provided in SampleSheet.csv; these reads will then
end up in the dreaded "Undetermined" FASTQ file.

This script salvages reads from the dumpster--give it the Undetermined files,
and it will attempt to further demultiplex the file based on the (proper)
SampleSheet.csv file. It's designed to tolerate massive failures in the index
sequences. There is no maximum number of mismatches--for every read, this
script checks which index sequence on the sample sheet has the lowest Hamming
distance to it. If two indices on the sample sheet is equally distant from the
sequenced index, then the read is discarded.

Of course, this works better if there were fewer libraries pooled together;
beyond, say, 16 libraries, salvaging would likely not work, as it gets
far more difficult to have a minimum Hamming distance to a single library.

To use this script, run `bcl2fastq` with an empty "dummy" SampleSheet.csv,
so that all reads get pooled into a single Undetermined file for each read
direction. Feed those Undetermined files into this script--then hope for
the best. Hopefully saves you from a re-run!

Works on compressed files.
""".strip()

import argparse
import csv
import gzip
import os.path
import re
import sys
import time

parser = argparse.ArgumentParser(
    description=docstring, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('sample_sheet', metavar='samplesheet_csv_file',
                    type=argparse.FileType('r'),
                    help='the SampleSheet.csv file used in the run.')
parser.add_argument('fastq_files', metavar='fastq_files',
                    type=argparse.FileType('r'), nargs='+',
                    help='(paired) FASTQ files.')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints diagnostic stuff to stderr.')
args = parser.parse_args()

def deduce_most_likely_sample(i7_and_i5, indices_to_samples):
    """
    Given an i7 and i5 pair, calculate the number of matching bases, and
    assign it to the sample whose barcode has the highest number of matching
    bases. Discard ties.
    """
    num_matching_bases = {}
    for i in indices_to_samples:
        num_matching_bases[i] = sum([x == y for x, y in zip(i7_and_i5, i)])
    
    # sort num_matching_bases in reverse order
    most_likely_sample = sorted(num_matching_bases, key=num_matching_bases.get, 
                                reverse=True)
    # discard ties
    if num_matching_bases[most_likely_sample[0]] == \
        num_matching_bases[most_likely_sample[1]]: return 'unknown'
    
    # else, return the sample ID of the best matching sample
    return indices_to_samples[most_likely_sample[0]]

# read in sample sheet, and create two dicts that associate i7/i5 sequences
# to sample IDs
csv_reader = csv.reader(args.sample_sheet)
sample_id_encountered = False
indices_to_samples = {}
for row in csv_reader:
    # skip lines until it hits "Sample_ID", the header row for "[Data]"
    if row[0] != 'Sample_ID' and sample_id_encountered == False: continue
    
    if row[0] == 'Sample_ID':
        sample_id_encountered = True
        continue
    
    # this section processes index sequences and sample IDs
    sample_id = row[0]
    i7_index = row[5]
    i5_index = row[7]
    
    indices_to_samples[f'{i7_index}+{i5_index}'] = sample_id

# script will deduce filetype/compression based on filename! script won't
# work properly if it's misled!
# supported extensions: .fq, .fastq; and either of them + .gz
for fq in args.fastq_files:
    if args.verbose:
        print (f'[{time.asctime()}] Processing {fq.name}...',
               file=sys.stderr)
               
    file_prefix, file_ext = os.path.splitext(fq.name)
    if file_ext == '.gz':
        fq = gzip.open(fq.name, 'rt')
        file_prefix, file_ext = os.path.splitext(file_prefix)
    
    lines_per_seq = 4
    
    # important bits from an annotation line (assumes reads from Illumina):
    #   @MN00524:15:000H2T77K:1:11101:18067:1014 1:N:0:GNCCTNCT+ANANCNCN
    #                                                  ======== --------
    #   = : i7 barcode
    #   - : i5 barcode
    #
    # demultiplexed output files will have the sample name in the filename
    all_output_files = {}
    output_file = ''
    counter = 0
    for line in fq:
        # this if block captures annotation lines. decide which output file
        # the line should go to based on the flowcell/lane information
        if counter % lines_per_seq == 0:
            # figure out i7 and i5 barcodes
            temp = line.split(' ')[1]
            i7_and_i5 = re.search(r':(\w+\+\w+)$', temp).group(1)
            mls = deduce_most_likely_sample(i7_and_i5, indices_to_samples)
            
        # write uncompressed output into separate files
        output_file = f'{file_prefix}.{mls}.fastq'
        if output_file not in all_output_files:
            all_output_files[output_file] = open(output_file, 'w')
        
        print (line, end='', file=all_output_files[output_file])
        counter += 1
        
        if args.verbose and counter % 1000000 == 0:
            print (f'[{time.asctime()}] {counter} sequences binned into',
                   f'{len(all_output_files)} separate files...', file=sys.stderr)
    
    # close all open files
    for f in all_output_files:
        all_output_files[f].close()

#!/usr/bin/env python3

docstring = """
> filter_array-specific_CpGs.py <

There are far more deconvolution algorithms trained on HM450/EPIC data than
whole genome methylomes (e.g., WGBS/EM-seq). Deconvolution can still be done
on the latter by selectively filtering for array-specific CpGs and outputting
the subset of positions into a separate table.

This script expects the user to use the same genome version for everything
(if Bismark is run on hg38, then use the hg38 manifest). Mismatch in genome
version will cause nonsensical results.

Wanding Zhou maintains an excellent repo of array annots--massive props to him.
https://zwdzwd.github.io/InfiniumAnnotation

This script expects ONE-BASED COORDINATES for beta value tables (a la Bismark)
and ZERO-BASED COORDINATES for Infinium annotations (a la Wanding's tables).
The beta value tables should also have chr/start/end as the first 3 columns.
Basically, this is how the files are generated by default, so if you're unsure,
then just trust the script to do the right thing (TM).

Input: a TSV table with headers and beta values, e.g.,
  chr    start    end      S1    S2    S3    [...]
  chr1   123456   123457   0.5   0.6   0.8
  chr2   234567   234568   0.4   0.6   0.9
  [...]

Output: a CSV table with headers, replaces chr/start/end with probe ID if
there's an overlap, e.g. assume first row did not match but second did,
  CpGs,S1,S2,S3,[...]
  cg2345,0.4,0.6,0.9
  [...]

Autodetects gzip compression.
""".strip()

import argparse
import csv
from pathlib import Path
import glob
import gzip
import sys
import time

import numpy as np

def diag_print(*args):
    """
    Helper function to print diagnostic-level verbose output to stderr.
    """
    print (f'[{time.asctime()}]', *args, file=sys.stderr)


parser = argparse.ArgumentParser(
    description=docstring, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('array_manifest', metavar='tsv_file', type=Path,
                    help='Illumina manifest, in tsv format.')
parser.add_argument('beta_tsv', metavar='tsv_files', type=Path, nargs='+',
                    help='Whole genome methylation data, in tsv format.')
parser.add_argument('--verbose', '-v', action='store_true',
                    help='prints diagnostic stuff to stderr.')
args = parser.parse_args()

# start by reading the array manifest into a NumPy array. NumPy string
# arrays (storing the probe ID) assumes strings are of fixed lengths, but the
# manifest have IDs ranging from 8-36 chars. use ints for the array, and have
# another dict to map the ints back to string labels
counter = 0
counter_to_probe_id = {}
array_annot = {}
tsv_reader = csv.reader(open('grch38p13_lambda_puc.seqlen.tsv'), delimiter='\t')
for row in tsv_reader:
    scaf = row[0]
    scaf_size = int(row[1]) + 100    # embiggen scaf slightly, just to be safe
    array_annot[scaf] = np.zeros(scaf_size, dtype=np.uint32)  # overflows at 4 mil

if args.array_manifest.suffix == '.gz':
    # assume gzip-compressed file
    tsv_reader = csv.reader(gzip.open(args.array_manifest, 'rt'), delimiter='\t')
else:
    # assume it's plaintext
    tsv_reader = csv.reader(open(args.array_manifest, 'r'), delimiter='\t')

if args.verbose:
    diag_print(f'Reading {args.array_manifest.name}...')

# skip header
header = next(tsv_reader)
for row in tsv_reader:
    scaf = row[0]
    # HM450 has control probes with NA for scaf/start/end
    if scaf == 'NA': continue
    
    # the manifest does contain probes that map to ALT chromosomes--discard
    if scaf not in array_annot: continue
    
    start_pos = int(row[1])    # ZERO-BASED!!
    end_pos = int(row[2])
    probe_id = row[8]
    
    # store data into the array
    counter += 1
    array_annot[scaf][start_pos:end_pos] = counter
    counter_to_probe_id[counter] = probe_id

# okay, time to parse the input file(s) for array-specific CpGs. print these
# lines out as a separate file. output filenames are fixed: e.g.,
#   blabla.tsv.gz --> blabla.HM450.csv.gz (name of manifest precedes suffixes)
manifest_type = args.array_manifest.name.split('.')[0]  # "HM450"/"EPIC"
for bt in args.beta_tsv:
    bt_split_name = bt.name.split('.')
    if bt.suffix == '.gz':
        # assume gzip-compressed file
        tsv_reader = csv.reader(gzip.open(bt, 'rt'), delimiter='\t')
        output_filename = '.'.join([*bt_split_name[:-2], manifest_type, 'csv.gz'])
        output_file = gzip.open(output_filename, 'wt')
    else:
        # assume it's plaintext
        tsv_reader = csv.reader(open(bt, 'r'), delimiter='\t')
        output_filename = '.'.join([*bt_split_name[:-1], manifest_type, 'csv'])
        output_file = open(output_filename, 'w')
    
    if args.verbose:
        diag_print(f'Reading {bt.name}, writing array-specific CpGs to {output_filename}...')
    
    # assume printing out into file format that substitutes scaf/start/end with
    # probe ID. assume header present in the input file
    header = next(tsv_reader)
    print ('CpGs', *header[3:], sep=',', file=output_file)
    for row in tsv_reader:
        # assume scaf/start/end are the first three columns
        scaf = row[0]
        try:
            start_pos = int(row[1]) - 1   # ONE-BASED!! hence the -1
            end_pos = int(row[2])
        except ValueError:
            continue
        
        # sanity checks for valid scaf/pos values
        if scaf not in array_annot: continue
        if end_pos < 0: continue
        if start_pos > len(array_annot[scaf]): continue
        
        # lift annots from the array
        counter_ids = array_annot[scaf][start_pos:end_pos]    # multiple annots possible
        
        if len(set(counter_ids)) > 1: continue    # if the CpG has two IDs
                                                  # assigned, discard the CpG
        
        counter_id = counter_ids[0]
        if not counter_id: continue
        
        # finally--print this valid line out
        probe_id = counter_to_probe_id[counter_id]
        print (probe_id, *row[3:], sep=',', file=output_file)
    
    output_file.close()

if args.verbose:
    diag_print('Done.')

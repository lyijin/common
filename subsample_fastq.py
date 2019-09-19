#!/usr/bin/env python3

"""
> subsample_fastq.py <

Subsamples all provided FASTQ files to the same depth (to the LOWEST amongst
the input files).

Files with "[.|_]R1[.|_]" in the filename are mandatory; similarly named files
with "R2" replacing "R1" are assumed to be the paired counterpart, and will be
subsampled in a similar fashion.
"""
import argparse
import csv
import re
import sys
import time

import numpy as np

parser = argparse.ArgumentParser(description="""
Subsamples all provided FASTQ files to the same depth (to the LOWEST amongst
the input files).

Files with "[.|_]R1[.|_]" in the filename are mandatory; similarly named files
with "R2" replacing "R1" are assumed to be the paired counterpart, and will be
subsampled in a similar fashion.""")

parser.add_argument('fastqs', metavar='fastq_files',
                    type=argparse.FileType('r'), nargs='+',
                    help='FASTQ files (accepts all R1 and R2 files).')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints diagnostic stuff to stderr.')
args = parser.parse_args()

def get_line_count(fname):
    with open(fname, 'r') as f:
        for n, line in enumerate(f):
            pass
    
    return n + 1

# classify input files into three bins: R1, R2 and unpaired
args_fnames = [x.name for x in args.fastqs]
r1_fnames = [x for x in args_fnames if re.search(r'[._]R1[._]', x)]

# drop filenames that contain "subsamp" in R1--these are output files
# produced by this script
r1_fnames = [x for x in r1_fnames if '.subsamp.f' not in x]

# store R2 filenames as a dict, with the R1 filename as the key. R1 files 
# without a R2 counterpart would not be a key in this dict
r2_fnames = {x:re.sub(r'(\.|_)R1(\.|_)', r'\1R2\2', x) for x in r1_fnames}
r2_fnames = {x:y for x, y in r2_fnames.items() if y in args_fnames}

# unpaired files are those that are neither R1 nor R2, excluding "subsamp"
unpaired_fnames = [x for x in args_fnames if x not in r1_fnames and x not in r2_fnames]
unpaired_fnames = [x for x in unpaired_fnames if '.subsamp.f' not in x]

if args.verbose:
    print (f'[{time.asctime()}] R1 files detected: {r1_fnames}', file=sys.stderr)
    print (f'[{time.asctime()}] R2 files detected: {r2_fnames}', file=sys.stderr)
    print (f'[{time.asctime()}] Unpaired files detected: {unpaired_fnames}',
           file=sys.stderr)

# for downstream purposes, unpaired files can be treated as orphaned R1s
r1_fnames += unpaired_fnames

# parse files to check what their original depths are
r1_line_counts = {x:get_line_count(x) for x in r1_fnames}
assert all([x % 4 == 0 for x in r1_line_counts.values()]), \
    f'One or more FASTQ files do not have line counts with multiple of 4: {r1_line_counts}'
r1_read_counts = {x:int(y/4) for x, y in r1_line_counts.items()}

# input files all good. subsample i things down to the lowest depth j
j = min(r1_read_counts.values())

if args.verbose:
    print (f'[{time.asctime()}] Minimum depth: {j:,} reads', file=sys.stderr)

for r1 in r1_fnames:
    i = r1_read_counts[r1]
    
    output_fname = r1.split('.')
    output_fname.insert(-1, 'subsamp')
    output_fname = '.'.join(output_fname)
    
    if args.verbose:
        print (f'[{time.asctime()}] Subsampling and writing to {output_fname}...',
               file=sys.stderr)
    
    with open(output_fname, 'w') as o:
        # first, initialise a numpy array of length i containing (i-j) 0s and
        # j 1s (where 1 = include the read, 0 = discard the read)
        inclusion_array = np.zeros(i, np.int8)
        inclusion_array[:j] = 1
        
        # second, shuffle the array once. seed is set to a fixed value for
        # replicability
        np.random.seed(1337)
        np.random.shuffle(inclusion_array)
        
        # based on the inclusion array, keep or toss reads away
        with open(r1, 'r') as f:
            for n in inclusion_array:
                if n:
                    # include read, i.e. 4 lines
                    print (f.readline().strip(), file=o)
                    print (f.readline().strip(), file=o)
                    print (f.readline().strip(), file=o)
                    print (f.readline().strip(), file=o)
                else:
                    # discard read
                    f.readline()
                    f.readline()
                    f.readline()
                    f.readline()
        
    # do the same for R2, if it exists
    if r1 not in r2_fnames: continue
    
    # R2 exists! run code below
    output_fname = r2_fnames[r1].split('.')
    output_fname.insert(-1, 'subsamp')
    output_fname = '.'.join(output_fname)
    
    if args.verbose:
        print (f'[{time.asctime()}] Subsampling and writing to {output_fname}...',
               file=sys.stderr)

    with open(output_fname, 'w') as o:
        with open(r2_fnames[r1], 'r') as f:
            for n in inclusion_array:
                if n:
                    # include read, i.e. 4 lines
                    print (f.readline().strip(), file=o)
                    print (f.readline().strip(), file=o)
                    print (f.readline().strip(), file=o)
                    print (f.readline().strip(), file=o)
                else:
                    # discard read
                    f.readline()
                    f.readline()
                    f.readline()
                    f.readline()

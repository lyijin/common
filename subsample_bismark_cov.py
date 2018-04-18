#!/usr/bin/env python3

"""
> subsample_bismark_covs.py <

Subsamples all input covs to the same depth (the LOWEST among the input covs).
"""
import argparse
import csv

import numpy as np

parser = argparse.ArgumentParser(description="""
Subsamples all input covs to the same depth (the LOWEST among the input covs.""")

parser.add_argument('covs', metavar='cov_files',
                    type=argparse.FileType('r'), nargs='+',
                    help='Bismark covs containing meth/unmeth information.')
args = parser.parse_args()

# parse files to check what their original depths are
overall_depths = {}

cov_filenames = sorted([x.name for x in args.covs])
for c in cov_filenames:
    tsv_reader = csv.reader(open(c), delimiter='\t')
    cov_counter = 0
    for row in tsv_reader:
        cov_counter += int(row[4]) + int(row[5])

    overall_depths[c] = cov_counter

# subsample i things down to the lowest depth j
j = min(overall_depths.values())
for c in cov_filenames:
    i = overall_depths[c]
    
    output_file = c[:-4] + '.subsamp.cov'
    with open(output_file, 'w') as o:
        # first, initialise a numpy array of length i containing (i-j) 0s and
        # j 1s (where 1 = include the read, 0 = discard the read)
        inclusion_array = np.zeros(i, np.int8)
        inclusion_array[:j] = 1
        
        # second, shuffle the array once. seed is set to a fixed value for
        # replicability
        np.random.seed(0)
        np.random.shuffle(inclusion_array)
        
        # based on the inclusion array, keep or toss reads away
        start_pos = 0
        end_pos = 0
        
        tsv_reader = csv.reader(open(c), delimiter='\t')
        for row in tsv_reader:
            meth = int(row[4])
            unmeth = int(row[5])
            
            end_pos = start_pos + meth + unmeth
            
            sub_incl_array = inclusion_array[start_pos:end_pos]
            # if sub_incl_array is entirely 0, toss the line
            if not np.any(sub_incl_array):
                start_pos = end_pos
                continue
            # if sub_incl_array is entirely 1, keep entire line
            elif np.all(sub_incl_array):
                print (*row, sep='\t', file=o)
                start_pos = end_pos
                continue
            else:
                # overwrite old meth and unmeth values with the sub_incl_array
                new_meth = sum(sub_incl_array[:meth])
                new_unmeth = sum(sub_incl_array) - new_meth
                
                row[3] = round(new_meth / (new_meth + new_unmeth) * 100, 4)
                row[4] = new_meth
                row[5] = new_unmeth
                
                print (*row, sep='\t', file=o)
                start_pos = end_pos
                continue

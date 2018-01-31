#!/usr/bin/env python3

"""
> tally_cov_stats.py <

Script that produces numbers that succinctly describes the coverages from 
cov files produced by bismark_methylation_extractor.

NumPy is used because x[x < 5] is super convenient for list-slicing.
"""
import argparse
import csv
import gzip

import numpy as np

parser = argparse.ArgumentParser(description="""
Script that produces numbers that succinctly describes the coverages from 
cov files produced by bismark_methylation_extractor.""")

parser.add_argument('cov_files', metavar='cov_files',
                    type=argparse.FileType('r'), nargs='+',
                    help='cov files that will be tallied.')
parser.add_argument('--gzip', '-g', action='store_true', default=False,
                    help='Tallied files are gzip-compressed.')

args = parser.parse_args()

# print header
print ('filename', '>= 1 cov', '>= 5 cov', '>= 10 cov', '>= 100 cov',
       'mean cov', 'median cov',
       '>= 1 meth', '>= 5 meth', '>= 10 meth', '>= 100 meth', sep='\t')

for f in args.cov_files:
    all_covs = []
    all_meths = []
    
    # yeah, very crude way to detect whether file is gzip-compressed
    if f.name[-3:] == '.gz':
        tsv_reader = csv.reader(gzip.open(f.name, 'rt'), delimiter='\t')
    else:
        tsv_reader = csv.reader(f, delimiter='\t')
    
    for row in tsv_reader:
        if not row: continue
        
        meth = int(row[4])
        unmeth = int(row[5])
        cov = meth + unmeth
        
        all_covs.append(cov)
        all_meths.append(meth)
    
    all_covs = np.array(all_covs, np.int32)
    all_meths = np.array(all_meths, np.int32)
    
    # print results
    print (f.name, len(all_covs[all_covs >= 1]), len(all_covs[all_covs >= 5]), 
           len(all_covs[all_covs >= 10]), len(all_covs[all_covs >= 100]),
           np.mean(all_covs), np.median(all_covs),
           len(all_meths[all_meths >= 1]), len(all_meths[all_meths >= 5]),
           len(all_meths[all_meths >= 10]), len(all_meths[all_meths >= 100]),
           sep='\t')

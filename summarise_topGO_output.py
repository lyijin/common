#!/usr/bin/env python3

"""
> summarise_topGO_output.py <

Script checks folder for output files from topGO pipeline (bp_*.txt, cc_*.txt
and mf_*.txt), then summarises everything into a summary file.

This script used to exist as a shell script that parsed the same files using
awk, but I discovered that awk might treat floats in scientific notation
as strings, therefore incorrectly concluding that '2.3e-09' is not smaller than
0.05. This weird behaviour was first diagnosed for the awk distributed in 
Debian on WSL--it didn't happen on Cygwin/Ubuntu on WSL.

Regardless, as awk is much less predictable than Python, this script was
written to give me peace of mind.
"""
import csv
import glob

bp_files = sorted(glob.glob('bp_*.txt'))

for bp_file in bp_files:
    # the '1' is to replace bp with the replacement string only once
    cc_file = bp_file.replace('bp', 'cc', 1)
    mf_file = bp_file.replace('bp', 'mf', 1)
    output_file = bp_file.replace('bp', 'summary', 1)
    
    with open(output_file, 'w') as opf:
        for f in [bp_file, cc_file, mf_file]:
            tsv_reader = csv.reader(open(f), delimiter='\t')
            
            # print header
            print ('-- {} --'.format(f), file=opf)
            header = next(tsv_reader)
            print ('', *header, sep='\t', file=opf)
            
            for row in tsv_reader:
                if not row: continue
                
                # two criteria to pass:
                # 1. At least 5 terms in the universe (row[3])
                if int(row[3]) < 5: continue
                
                # 2. P value < 0.05 OR if it contains the character '<' 
                #    (because R outputs stuff like '< 1e-30') (row[6])
                if row[6][0] == '<':
                    print ('\t'.join(row), file=opf)
                elif float(row[6]) < 0.05:
                    print ('\t'.join(row), file=opf)
        
            print (file=opf)

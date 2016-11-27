#!/usr/bin/env python3

"""
> parse_vcf.py <

Python script to help parse the vcf files produced by SnpEff.
"""
import argparse
import csv
import re
import sys

def parse_vcf(filename):
    """
    Expects a properly-formatted VCF file. Returns a dict containing 
    positions grouped by editing type.
    """
    all_effects = {}
    
    tsv_reader = csv.reader(filename, delimiter='\t')
    for row in tsv_reader:
        if not row: continue
        if row[0][0] == '#': continue      # skip comment lines
        
        pos = '_'.join(row[:2])
        pos_effects = re.findall('EFF=(\w+)\(', row[7]) + \
                      re.findall(',(\w+)\(', row[7])
        
        # pos_effects sometimes will have duplicates (due to multiple
        # transcripts). suppress duplicates with set().
        for p in set(pos_effects):
            if p not in all_effects:
                all_effects[p] = []
        
            all_effects[p].append(pos)
        
    return all_effects
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    Python script to help parse the vcf files produced by SnpEff.""")
    parser.add_argument('vcf_file', metavar='vcf_filename',
                        type=argparse.FileType('r'), nargs='?',
                        default=sys.stdin,
                        help='VCF file, or piped from stdin.')
    
    args = parser.parse_args()
    
    overall_effects = parse_vcf(args.vcf_file)
    
    for e in sorted(overall_effects):
        print (e, len(overall_effects[e]), sep='\t')

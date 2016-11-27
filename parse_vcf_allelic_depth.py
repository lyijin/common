#!/usr/bin/env python3

"""
> parse_vcf_allelic_depth.py <

Python script scans through a VCF file and, on a per-position basis, parses
the number of ACGTs in that position.
"""
import argparse
import csv

parser = argparse.ArgumentParser(description="""
Python script scans through a VCF file and, on a per-position basis, parses
the number of ACGTs in that position.""")

parser.add_argument('vcf', metavar='vcf_file',
                    type=argparse.FileType('r'),
                    help='*.vcf file that has to be parsed.')

args = parser.parse_args()

bases = ['A', 'C', 'G', 'T']

tsv_reader = csv.reader(args.vcf, delimiter='\t')
for row in tsv_reader:
    if not row: continue
    if row[0][0] == '#': continue       # skip comment lines
    
    chrom = row[0]
    pos = row[1]
    ref = row[3]
    alt = row[4]                        # 'G,T,<*>'
    ad = row[7].split(';')[3]           # 'AD=29,3,1,1'
    
    # match allelic depths to ref and alt
    ref_alt = [ref] + alt.split(',')
    ad = [int(x) for x in ad[3:].split(',')]
    
    # construct a dict to store per-base frequencies
    ad_dict = {x:y for x,y in zip(ref_alt, ad)}
    for b in bases:
        if b not in ad_dict: ad_dict[b] = 0
    
    # print parsed stuff out
    print (chrom, pos, ref, *[ad_dict[x] for x in bases], sep='\t')

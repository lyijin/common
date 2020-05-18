#!/usr/bin/env python3

"""
> filter_nonconverted_reads.py <

Script takes in BAM files produced by Bismark, and filters out reads that
appear nonconverted. Reads with non-CpG methylation above a specified rate are
tossed, i.e. (methylated CHG + methylated CHH) / (all CHG + all CHH) > n.

Like other Bismark scripts, this script outputs filtered files in the directory
where the command was called. This behaviour can be altered with the -o flag.

Requires pysam to parse BAM files.
"""

import argparse
import collections
from pathlib import Path
import sys

import pysam

parser = argparse.ArgumentParser(description="""
Script takes in BAM files produced by Bismark, and filters out reads that
appear nonconverted (non-CpG methylation in the read above a specified rate).
""")
parser.add_argument('bam_files', metavar='bam_file',
                    type=Path, nargs='+',
                    help='bam file(s) produced by Bismark.')
parser.add_argument('-m', '--meth_rate', metavar='n', type=float, default=0.2,
                    help='remove reads with non-CpG methylation rate of >= n.')
parser.add_argument('-o', '--output_directory', metavar='folder_name',
                    type=Path, default='./',
                    help='set output directory.')

args = parser.parse_args()

# read bamfiles
for bamfile in args.bam_files:
    prev_read = {'name': ''}
    bf = pysam.AlignmentFile(bamfile, 'rb')
    ofile = args.output_directory / bamfile.name.replace('.bam', '.filt.bam')
    of = pysam.AlignmentFile(ofile, 'wb', header=bf.header)
    
    # define a few variables for logging
    n_reads = 0
    n_pass = 0
    logfile = args.output_directory / bamfile.name.replace('.bam', '.filt.txt')
    lf = logfile.open('w')
    
    for read in bf:
        # handle paired-end reads
        if read.is_paired:
            read_qname = read.qname
            prev_read_qname = prev_read['name']
            
            if read_qname == prev_read_qname:
                # `read`/`prev_read` are pairs of the same read!
                n_reads += 1
                
                # start parsing the 'XM' tag
                read_xm_tag = read.tags[2][1]
                prev_read_xm_tag = prev_read['tags'][2].replace('XM:Z:', '')
                combined_xm_tag = read_xm_tag + prev_read_xm_tag
                
                # tally the number of 'zZxXhH' characters in the XM tag
                #   z: unmethylated C in CpG context
                #   Z: methylated C in CpG context
                #   x: unmethylated C in CHG context
                #   X: methylated C in CHG context
                #   h: unmethylated C in CHH context
                #   H: methylated C in CHH context
                c = collections.Counter(combined_xm_tag)
                
                # set meth proportion to 0 if either meth/unmeth version doesn't
                # exist in the read
                cpg_meth = 0 if not c['Z'] else c['Z'] / (c['Z'] + c['z'])
                chg_meth = 0 if not c['X'] else c['X'] / (c['X'] + c['x'])
                chh_meth = 0 if not c['H'] else c['H'] / (c['H'] + c['h'])
                chn_meth = 0 if not chg_meth + chh_meth else \
                    (c['X'] + c['H']) / (c['X'] + c['x'] + c['H'] + c['h'])
                
                # selectively print reads below specified nonconversion rate
                if chn_meth < args.meth_rate:
                    n_pass += 1
                    of.write(pysam.AlignedSegment.from_dict(prev_read, bf.header))
                    of.write(read)
            
            prev_read = read.to_dict()
        
        # handle single-end stuff
        else:
            n_reads += 1
            
            # parse 'XM' tag, similar logic like paired-end reads
            read_xm_tag = read.tags[2][1]
            c = collections.Counter(read_xm_tag)
            
            cpg_meth = 0 if not c['Z'] else c['Z'] / (c['Z'] + c['z'])
            chg_meth = 0 if not c['X'] else c['X'] / (c['X'] + c['x'])
            chh_meth = 0 if not c['H'] else c['H'] / (c['H'] + c['h'])
            chn_meth = 0 if not chg_meth + chh_meth else \
                (c['X'] + c['H']) / (c['X'] + c['x'] + c['H'] + c['h'])
            
            if chn_meth < args.meth_rate:
                n_pass += 1
                of.write(read)
        
    bf.close()
    of.close()

    # print logged variables
    if n_reads:
        logstr = f'{bamfile.name}: retained {n_pass:,} of {n_reads:,} reads ' + \
                 f'({n_pass / n_reads * 100:.1f}%).'
    else:
        logstr = f'{bamfile.name}: ERROR; has no valid reads in it.'
    
    print (logstr, file=lf)
    lf.close()
    
    print (logstr, file=sys.stderr)

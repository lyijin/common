#!/usr/bin/env python3

"""
> tabulate_gcfrac_meth.py <

Script reads in deduplicated BAM files produced by Bismark, and on a per-read
basis, print

  1. read annotation
  2. read length
  3. position of 5' end
  4. read direction
  5. GC fraction
  6/7. meth C in CpG context; all C in CpG context
  8/9. meth C in CHG context; all C in CHG context
  10/11. meth C in CHH context; all C in CHH context
  12/13. meth C in CHN context; all C in CHN context

These values can then be used by downstream files for further analysis/plotting.

Requires pysam to parse BAM files.
"""

import argparse
import collections
import gzip
from pathlib import Path
import sys
import time

import pysam
import numpy as np

import parse_fasta

def benchmark_print(message):
    """
    Prints a timestamped message to stderr, to aid visual assessment of
    chokepoints in code. 
    """
    print (f'[{time.asctime()}] {message}', file=sys.stderr)

def calc_gc_frac(val_array, start_pos, end_pos):
    '''
    `val_array`: genomic sequence transformed into np array
    `start_pos`, `end_pos`: genomic coordinates of the start and end of the 
                            mapped read
    GC fraction can be thus calculated as:
      GC = 0.5 + 0.5 * sum(val_array[start_pos:end_pos]) /
                       count_nonzero(val_array[start_pos:end_pos])
    this is because the sum() / count_nonzero() portion ranges from [-1, 1],
    while GC fraction is [0, 1]. the 0.5 parts is just a rescale function.
    '''
    s = sum(val_array[start_pos:end_pos])
    cnz = np.count_nonzero(val_array[start_pos:end_pos])
    
    return 0.5 + 0.5 * s/cnz

parser = argparse.ArgumentParser(description="""
Script reads in deduplicated BAM files produced by Bismark, and on a per-read
basis, check for GC fraction and methylation levels.""")

parser.add_argument(
    'genome_fasta', metavar='genome_fasta', type=Path,
    help='genome FASTA file used in mapping BAMs.')
parser.add_argument(
    'bam_files', metavar='bam_file', type=Path, nargs='+',
    help='bam file(s) produced by Bismark.')
parser.add_argument(
    '-o', '--output_directory', metavar='folder_name', type=Path, default='./',
    help='set output directory.')
parser.add_argument('-v', '--verbose', action='store_true',
    help='prints diagnostic stuff to stderr.')
args = parser.parse_args()

# parse the genome. there are three pieces of information that are important:
#   1. G/C/S [+1]
#   2. A/T/W [-1]
#   3. other unmentioned bases e.g. N, other degenerate bases [0]
# square brackets denote the value assigned to the base in the int8 array.
if args.verbose:
    benchmark_print(f'Reading {args.genome_fasta.name}...')

genome_fasta = parse_fasta.get_all_sequences(args.genome_fasta, 'fasta')
genome_vals = {}

for seq in genome_fasta:
    genome_vals[seq] = np.zeros(len(genome_fasta[seq]), dtype=np.int8)
    
    for n in range(len(genome_fasta[seq])):
        base = genome_fasta[seq][n].upper()
        if base == 'C' or base == 'G' or base == 'S':
            genome_vals[seq][n] = 1
        elif base == 'A' or base == 'T' or base == 'W':
            genome_vals[seq][n] = -1

if args.verbose:
    benchmark_print(f'Finished parsing {args.genome_fasta.name}.')

# read bamfiles
for bamfile in args.bam_files:
    if args.verbose:
        benchmark_print(f'Reading {bamfile.name}...')
    
    prev_read = pysam.AlignedSegment()    # empty object
    bf = pysam.AlignmentFile(bamfile, 'rb')
    ofile = args.output_directory / bamfile.name.replace('.bam', '.gcfrac.tsv.gz')
    of = gzip.open(ofile.name, 'wt')
    
    # define a few variables for logging
    n_reads = 0
    logfile = args.output_directory / bamfile.name.replace('.bam', '.gcfrac.txt')
    lf = logfile.open('w')
    # disagree = agree = 0  # for debugging
    
    for read in bf:
        # handle paired-end reads
        if read.is_paired:
            if read.qname == prev_read.qname:
                # only use paired-end reads that map to the same reference seq
                if not read.reference_name == prev_read.reference_name: continue
                
                # `read`/`prev_read` are pairs of the same read!
                n_reads += 1
                
                # use union of covered positions across both reads to
                # determine start:end of the entire template
                covered_pos = read.positions + prev_read.positions
                start_pos = min(covered_pos) 
                end_pos = max(covered_pos) + 1    # x.positions produces a list
                                                  # of 0-based positions
                
                # use R1 to define 5' position and directionality of read
                r1 = read if read.is_read1 else prev_read
                # note that +1 is to transform 0-based to 1-based positions
                fivep_pos = max(r1.positions) + 1 if r1.is_reverse else min(r1.positions) + 1
                direction = '-' if r1.is_reverse else '+'
                                
                # calculate GC fraction
                read_gc_frac = calc_gc_frac(genome_vals[read.reference_name], start_pos, end_pos)
                
                # parse the 'XM' tag to get meth levels
                read_xm_tag = read.tags[2][1]
                prev_read_xm_tag = prev_read.tags[2][1]
                
                # by default, if R1 and R2 overlaps, for the overlapping bit,
                # `bismark` keeps the calls on R1 and discards the corresponding
                # calls on R2 (rationale: ends of R1 have higher read quals
                # than R2)
                #
                # here, a dict is used to store key-value pairs of
                # "pos:XM_character" for R2 first, then stores R1. in doing so,
                # if R1 has a "pos" call that overlaps with R2, the former
                # will overwrite the latter!
                if prev_read.is_read2:
                    combined_xm_tag = {x:y for x, y in zip(prev_read.positions, prev_read_xm_tag)}
                    for x in zip(read.positions, read_xm_tag):
                        combined_xm_tag[x[0]] = x[1]
                else:
                    combined_xm_tag = {x:y for x, y in zip(read.positions, read_xm_tag)}
                    for x in zip(prev_read.positions, prev_read_xm_tag):
                        # DEBUG: disagree / (agree + disagree) should be ~0.02
                        # if x[0] in combined_xm_tag:
                        #     if combined_xm_tag[x[0]] != x[1]:
                        #         disagree += 1
                        #     else:
                        #         agree += 1
                        combined_xm_tag[x[0]] = x[1]
                
                # sort the dict by pos, then join the values together to form
                # an actual "merged" xm_tag!
                combined_xm_tag = ''.join(combined_xm_tag[x] for x in sorted(combined_xm_tag))
                
                # tally the number of 'zZxXhH' characters in the XM tag
                #   z: unmethylated C in CpG context
                #   Z: methylated C in CpG context
                #   x: unmethylated C in CHG context
                #   X: methylated C in CHG context
                #   h: unmethylated C in CHH context
                #   H: methylated C in CHH context
                c = collections.Counter(combined_xm_tag)
                
                # print stats for each read
                # read annot, length, 5' start pos, direction;
                # GC fraction;
                # meth C in CpG, all C in CpG;
                # meth C in CHG, all C in CHG;
                # meth C in CHH, all C in CHH;
                # meth C in CHN, all C in CHN
                print (read.qname, end_pos - start_pos, fivep_pos, direction,
                       f'{read_gc_frac:.5f}',
                       c['Z'], c['Z'] + c['z'],
                       c['X'], c['X'] + c['x'],
                       c['H'], c['H'] + c['h'],
                       c['X'] + c['H'], c['X'] + c['x'] + c['H'] + c['h'],
                       sep='\t', file=of)
            
            prev_read = read
        
        # handle single-end stuff
        else:
            n_reads += 1
            
            # use covered positions to determine start:end of the entire template
            covered_pos = read.positions
            start_pos = min(covered_pos) 
            end_pos = max(covered_pos) + 1    # x.positions produces a list
                                              # of 0-based positions
            
            # note that +1 is to transform 0-based to 1-based positions
            fivep_pos = max(covered_pos) + 1 if read.is_reverse else min(covered_pos) + 1
            direction = '-' if read.is_reverse else '+'
            
            # calculate GC fraction
            read_gc_frac = calc_gc_frac(genome_vals[read.reference_name], start_pos, end_pos)
            
            # parse 'XM' tag, similar logic like paired-end reads
            read_xm_tag = read.tags[2][1]
            c = collections.Counter(read_xm_tag)
            
            print (read.qname, end_pos - start_pos, fivep_pos, direction,
                   f'{read_gc_frac:.5f}',
                   c['Z'], c['Z'] + c['z'],
                   c['X'], c['X'] + c['x'],
                   c['H'], c['H'] + c['h'],
                   c['X'] + c['H'], c['X'] + c['x'] + c['H'] + c['h'],
                   sep='\t', file=of)
        
        # verbose logging
        if args.verbose and n_reads % 1_000_000 == 0:
            benchmark_print(f'Processed {n_reads:,} reads...')
        
    bf.close()
    of.close()

    # print logged variables
    if n_reads:
        logstr = f'{bamfile.name}: processed {n_reads:,} reads.'
    else:
        logstr = f'{bamfile.name}: ERROR; has no valid reads in it.'
    
    print (logstr, file=lf)
    lf.close()
    
    print (logstr, file=sys.stderr)

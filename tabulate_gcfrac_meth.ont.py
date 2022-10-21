#!/usr/bin/env python3

docstring = """
> tabulate_gcfrac_meth.ont.py <

Script reads in "mod_mappings.bam" files produced by megalodon, and on a
per-read basis, print
  
  1. read ID
  2. length of read that was mapped to genome (not length of full read)
  3. position of 5' end
  4. read direction
  5. GC fraction
  6. # 5mC in CpGs
  7. # all CpGs

These values can then be used by downstream files for further analysis/plotting.

Allows for the specification of a "region of interest"--CpG information outside
this window would not be considered.

Methylated positions are filtered on a per-read, per-site basis. In megalodon,
methylation state is a continuum (i.e., non-binary), ranging from 0% confidence
to 100% confidence. Based on their recommended cutoffs, confidence values are
treated as
   0 -  20% confidence: unmethylated C
  20 -  80% confidence: C with unknown state
  80 - 100% confidence: methylated C

Reads below 100 bp are discarded.

Requires pysam to parse BAM files.
""".strip()

import argparse
import gzip
from pathlib import Path
import re
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

parser = argparse.ArgumentParser(
    description=docstring, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument(
    'genome_fasta', metavar='genome_fasta', type=Path,
    help='genome FASTA file used in mapping BAMs.')
parser.add_argument(
    'bam_files', metavar='bam_file', type=Path, nargs='+',
    help='bam file(s) produced by megalodon.')
parser.add_argument(
    '-o', '--output_directory', metavar='folder_name', type=Path, default='./',
    help='set output directory.')
parser.add_argument(
    '-r', '--region_of_interest', metavar='roi', type=int, nargs=2,
    help='set region of interest; use 1-based coords for start/end.')
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

# do some calcs if "region of interest" mode is specified
if args.region_of_interest:
    start_roi = min(args.region_of_interest) - 1    # convert this to 0-based
    end_roi = max(args.region_of_interest)

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
    
    for read in bf:
        # reads are all single-end, no need to bother about handling paired- vs.
        # single-end reads
        n_reads += 1
        
        # use covered positions to determine start:end of the entire template
        covered_pos = read.positions
        start_pos = min(covered_pos) 
        end_pos = max(covered_pos) + 1    # x.positions produces a list
                                          # of 0-based positions
        
        # grab 'Mm' and 'Ml' tag, too complicated to describe in this script;
        # check this out https://samtools.github.io/hts-specs/SAMtags.pdf
        # in the context of this script, Ml quantifies likelihood that C is 5mC.
        # there are THREE possible scenarios:
        #   1. C is deemed methylated (> 80% confidence; Ml score > 204): "M"
        #   2. C is deemed unmethylated (< 20% confidence; Ml score < 51): "U"
        #   3. Low confidence call (20% <= confidence <= 80%; 51 <= Ml <= 204)
        #      Discard these Cs, so that meth% calcs are M / (U + M)
        # see https://github.com/nanoporetech/megalodon/issues/206
        if len(read.tags) < 3: continue     # some rows do not have Ml tag
        
        read_mm_tag = read.tags[1][1]
        read_ml_tag = read.tags[2][1]
        
        bool_partial_overlap = False
        if args.region_of_interest:
            # check whether read has any overlap with region of interest; if
            # no overlap, skip read completely
            if start_pos >= end_roi or end_pos <= start_roi: continue
            
            # check partial overlap
            if max(end_pos, end_roi) - min(start_pos, start_roi) > end_roi - start_roi:
                bool_partial_overlap = True
        
        # do computationally expensive parsing when there's a partial overlap
        if bool_partial_overlap:
            # disregard really short reads post-overlap
            readlen = min(end_pos, end_roi) - max(start_pos, start_roi)
            if readlen < 100: continue
            fivep_pos = min(end_pos, end_roi) if read.is_reverse else max(start_pos, start_roi) + 1
            
            # if there is some overlap, then force GC% calculation to the
            # overlapped region
            read_gc_frac = calc_gc_frac(genome_vals[read.reference_name],
                                        max(start_pos, start_roi), 
                                        min(end_pos, end_roi))
            
            # get 5mC positions (if read on Watson, then "C" in "CG"; if read
            # on Crick, then "G" in "CG")
            if read.is_reverse:
                c_pos = [x.start() + start_pos + 1 for x in re.finditer(r'G', read.seq)]
                c_pos = c_pos[::-1]  # meth confidences in Ml is always in 5'-to-3'
            else:
                c_pos = [x.start() + start_pos + 1 for x in re.finditer(r'C', read.seq)]
            
            # define an array to store not-CpG / unmeth CpG / meth CpG statuses
            # this can be figured out by 'Mm' and 'Ml' tags
            #   0: C/G not in CpG, or caller not confident 
            #   1: unmeth CpG: "U"
            #   2: meth CpG: "M"
            read_mm_tag = read_mm_tag.replace('C+m,', '').replace(';', '').split(',')
            read_mm_tag = [int(x) for x in read_mm_tag]
            status_array = []
            for mm in read_mm_tag:
                status_array += [0] * mm
                
                ml = read_ml_tag.pop(0)
                if ml > 204:
                    status_array.append(2)
                elif ml < 51:
                    status_array.append(1)
                else:
                    status_array.append(0)
            
            # then pad out status array with 0 so that lengths of "c_pos" and 
            # "status_array" is identical
            status_array += [0] * (len(c_pos) - len(status_array))
            
            # go through "c_pos" and "status_array" position-by-position,
            # and filter for pos within region-of-interest, 
            filt_status_array = []
            for n in range(len(status_array)):
                if start_roi < c_pos[n] <= end_roi:
                    filt_status_array.append(status_array.pop(0))
                else:
                    status_array.pop(0)
            
            n_5mc = sum(1 for x in filt_status_array if x > 1)
            n_cpg = sum(1 for x in filt_status_array if x > 0)
        
        else:
            # disregard really short reads
            readlen = end_pos - start_pos
            if readlen < 100: continue
            fivep_pos = end_pos if read.is_reverse else start_pos + 1
            
            # calculate GC fraction normally
            read_gc_frac = calc_gc_frac(genome_vals[read.reference_name], start_pos, end_pos)
            
            # calculate number of 5mC and CpG just from the 'Ml' tag
            n_5mc = sum(1 for ml in read_ml_tag if ml > 204)  # "> 204" == "> 80%"
            n_cpg = len(read_ml_tag)
        
        print (read.qname, readlen, fivep_pos, '-' if read.is_reverse else '+',
               f'{read_gc_frac:.5f}', n_5mc, n_cpg, sep='\t', file=of)
        
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

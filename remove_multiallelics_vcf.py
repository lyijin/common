#!/usr/bin/env python3

docstring = """
> remove_multiallelics_vcf.py <

This script was written specifically to deal with gnomAD VCF files. I wanted
a file containing only *biallelic* SNPs, but existing tools seem to stop at
picking SNPs out... but not removing multiallelics (or removing all but one).

There are many ways to skin the pick-SNPs-out cat. I used `gatk SelectVariants`,
with `--exclude-filtered true` to pick PASSes and `--restrict-alleles-to
BIALLELIC --select-type-to-include SNP` to get "biallelic" SNPs.

The biallelic flag doesn't work on gnomAD files because their VCFs have one
line per variant, e.g.,
  chr1	14402	rs1345225645	G	A	331.31	PASS	...
  chr1	14402	rs1345225645	G	C	331.31	PASS	...

It should work if it was
  chr1	14402	rs1345225645	G	A,C	331.31	PASS	...

But alas, there's a million file formats in bioinformatics and none seem to
ever be ideal.

Am I reinventing the wheel, other tools surely would be able to do this no?
I don't think so... from what I've found:

It was a feature request for `bcftools`, but that flag seemed to have
disappeared/changed its functionality over time.
https://github.com/samtools/bcftools/issues/974

`bcftools norm --rm-dup` seem to always pick the first of duplicated positions.
https://samtools.github.io/bcftools/bcftools.html#norm

"Oh you can run `bcftools --multiallelics +snps` to join biallelic records into
multiallelics and then cull it"--ugh gnomAD files are huge man, can't a simpler
script deal with this.

There's an `awk` way of doing things but it reads the entire file into the
damn memory.
  `awk -F '\t' '{chrpos[$1,$2]++; text[$1,$2]=$0} END 
   {for (cp in chrpos) if (chrpos[cp] == 1) print text[cp]}'`

Fine for smaller files, not fine for 30+ GB .vcf.gz. But bonus is that this
`awk` way doesn't assume VCF file is sorted.

This script ASSUMES
  1. VCF input file is sorted
  2. VCF input file only has SNPs (too lazy to write SNP detection code)

Then it works like a charm, culls all multiallelic SNPs in the biallelic VCF
format, prints plain-text to standard out. Super memory-efficient because
it only stores the previous line in memory.
"""

import argparse
import csv
import gzip
from pathlib import Path

parser = argparse.ArgumentParser(
    description=docstring, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('vcf_file', metavar='vcf_file', type=Path,
                    help='Sorted VCF file in biallelic format. '
                         'Can be gzip-compressed.')
args = parser.parse_args()

prev_line = []
duplicate_chrpos = ['', '']

# yeah, very crude way to detect whether file is gzip-compressed
if args.vcf_file.name[-3:] == '.gz':
    tsv_reader = csv.reader(gzip.open(args.vcf_file.name, 'rt'), delimiter='\t')
else:
    tsv_reader = csv.reader(open(args.vcf_file), delimiter='\t')

for line in tsv_reader:
    # print all comment lines
    if line[0][0] == '#':
        print (line)
        continue

    # for the first non-duplicated line
    if not prev_line:
        prev_line = line
        continue

    # start to filter lines. if the first two columns matches a duplicated
    # chr/pos pairing, suppress printing
    if line[:2] == duplicate_chrpos:
        continue

    # check whether duplication happened
    if line[:2] == prev_line[:2]:
        duplicate_chrpos = prev_line[:2]
        continue

    # ok this line's chr/pos doesn't match previous line's
    if prev_line[:2] != duplicate_chrpos:
        print (*prev_line, sep='\t')
    prev_line = line
    duplicate_chrpos = ['', '']

# handle last line
if line[:2] != duplicate_chrpos:
    print (*line, sep='\t')

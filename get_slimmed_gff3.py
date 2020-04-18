#!/usr/bin/env python3

"""
> get_slimmed_gff3.py <

Based on input gff3 files, pick the gene model that produces the longest
cDNA/protein sequence, then produces a gff3 file that only contains one
gene model per gene.
"""
import argparse
import re

import parse_gff3

parser = argparse.ArgumentParser(description="""
Based on input gff3 files, pick the gene model that produces the longest
cDNA/protein sequence, then produces a gff3 file that only contains one
gene model per gene.""")

parser.add_argument('scaffold_gff3', metavar="gff3_file",
                    type=argparse.FileType('r'), 
                    help="corresponding gff3 file of the genome.")

args = parser.parse_args()

# read genome details into memory
scaffold_gff3 = parse_gff3.parse_gff3(args.scaffold_gff3, 'exon')

# pick longest transcript in the gff3 file
scaffold_gff3 = parse_gff3.pick_longest_mRNA(scaffold_gff3)

# store which transcript is the longest in a dict
longest_mRNA = {}
for scaf in scaffold_gff3:
    for gene in scaffold_gff3[scaf]:
        # double check that there is only one gene model per gene; crash 
        # noisily if not
        tx = list(scaffold_gff3[scaf][gene].mRNAs.keys())
        assert len(tx) == 1, f'{gene} does not have one gene model per gene!'

        longest_mRNA[gene] = tx[0]

# dump all longest mRNAs into a list
all_longest_mRNAs = list(longest_mRNA.values())

# rewind the gff3 file now, then depending on what's in the line, different
# strategies are adopted
#   1. comment lines: retain only comment lines for longest gene model
#   2. content lines with no "Parent=": these are genes. retain everything
#   3. content lines with "Parent=": retain if parent is the longest gene model
args.scaffold_gff3.seek(0)
for line in args.scaffold_gff3:
    line = line.strip()
    
    # if line is blank, just print it out, then go onto next iteration
    if not line:
        print (line)
        continue
    
    if line[0] == '#':
        # comment line has two types: those that contain protein sequence, and 
        # those that do not
        if 'PROT' in line:
            genem = line.split(' ')[1]
        else:
            genem = line.split(' ')[2]
            # remove commas in gene model name, if they exist
            if genem[-1] == ',':
                genem = genem[:-1]
        
        if genem in all_longest_mRNAs:
            print (line)
    
    else:
        # these are content lines
        if 'Parent=' in line:
            # check whether this line is an "mRNA" line. if it is, the
            # transcript name is under 'Name='
            if '\tmRNA\t' in line:
                genem = re.search(r'Name=(.*?)$', line)[1]
            else:
                genem = re.search(r'Parent=(.*?)(;|$)', line)[1]
            
            if genem in all_longest_mRNAs:
                print (line)
        else:
            print (line)

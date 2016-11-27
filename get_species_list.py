#!/usr/bin/env python3

"""
> get_species_list.py <

When provided with a .tsv file produced by parse_blast_xml.py, parse and list
the species present in the table.
"""

import argparse
import csv
import re
import sys

parser = argparse.ArgumentParser(description="""
When provided with a .tsv file produced by parse_blast_xml.py, parse and list
the species present in the table.""")

parser.add_argument('pbx_tsv', metavar="tsv_filename",
                    type=argparse.FileType('r'), nargs='?',
                    default=sys.stdin, 
                    help="name of file produced by parse_blast_xml.py.")
args = parser.parse_args()

species_list = []

tsv_reader = csv.reader(args.pbx_tsv, delimiter='\t')

# skip header row
next(tsv_reader)
for row in tsv_reader:
    # skip blank rows
    if not row: continue
    
    s = re.findall('\[(.*?)\]', row[2])
    
    if s:
        # there can be multiple "[xxxx]" in the annotation. assume the last one
        # is the correct one.
        species = s[-1]
        
        # exclude stuff in brackets e.g. "Xenopus (Silurana) tropicalis"
        species = re.sub('\(.*?\) ', '', species)
        
        # some poorly-formatted annotations produces random ' ' and '[' at
        # the start of strings
        species = species.replace('[', '')
        species = species.strip()
        
        # only record genus and species
        species = ' '.join(species.split(' ')[:2])
        
        species_list.append(species)
    else:
        # this makes sure that even though no species exist in the line,
        # a blank line is created so that the XML file and output from this
        # script can be merged sideways
        species_list.append('')
    
# printing happens in a chunk to avoid broken pipe errors if this output
# gets piped to other scripts
output = ''
for s in species_list:
    output += s + '\n'

output = output.strip()
print (output)

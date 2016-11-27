#!/usr/bin/env python3

"""
> tally_species.py <

Based on the output of resolve_species.py, tally the number of species 
that falls under the superkingdom / kingdom levels of:
- Bacteria
- Archaea
- Viruses
- Eukaryota

Eukaryota is further subdivided into:
- Viridiplantae / Plantae
- Metazoa
- Fungi
- (for all other things) Protista

Program logic is:
- Check whether "kingdom" field is filled. If it's filled, it's eukaryotic,
  and one of plant/metazoa/fungi.
- If it's not filled, check for "bacteria" or "archaea". If both aren't found,
  then the species is a protist.
"""

import argparse
import collections
import csv
import re

parser = argparse.ArgumentParser(description="""
Based on the output of resolve_species.py, tally the number of species 
that falls under the superkingdom / kingdom levels.""")

parser.add_argument('tsv_file', metavar="resolved_species_tsv",
                    type=argparse.FileType('r'),
                    help="Resolved species file in .tsv format.")
parser.add_argument('-c', '--check_protists', action='store_true',
                    help='Display all protists for visual assessment.')
args = parser.parse_args()

tsv_reader = csv.reader(args.tsv_file, delimiter='\t')

kingdom_tally = {'Archaea': [], 'Bacteria': [], 'Viruses': [], 
                 'Viridiplantae': [], 'Plantae': [], 'Metazoa': [],
                 'Fungi': [], 'Chromista': [], 'Protista': []}

for row in tsv_reader:
    if not row: continue
    
    species = row[0]
    authority = row[1]
    kingdom = row[2]
    full_desc = row[8]
    organism_type = full_desc.split('>')[0].strip()
    
    # jump out if we hit any of these keywords
    if organism_type == 'Unclassified sequences': continue
    
    # manually fix some exceptions
    if species == 'Bacteria':
        kingdom_tally['Bacteria'].append(species)
        continue
    elif authority == 'not found':
        if 'bacteri' in species.lower():
            kingdom_tally['Bacteria'].append(species)
        elif 'phage' in species.lower():
            kingdom_tally['Viruses'].append(species)
        continue
    
    if kingdom in kingdom_tally:
        kingdom_tally[kingdom].append(species)
    else:
        # Virus check
        if organism_type == 'Viruses':
            kingdom_tally[organism_type].append(species)
            continue
        
        # Bacteria/Archaea check
        superk = full_desc.split('>')[1].strip()
        if superk in ['Bacteria', 'Archaea']:
            kingdom_tally[superk].append(species)
        elif superk in ['Artificial']:
            continue
        else:
            # bin everything else as Protista
            # print (species)
            kingdom_tally['Protista'].append(species)

# merge Plantae and Viridiplantae
kingdom_tally['Plantae'] += kingdom_tally['Viridiplantae']
del kingdom_tally['Viridiplantae']

# merge Chromista and Protista
kingdom_tally['Protista'] += kingdom_tally['Chromista']
del kingdom_tally['Chromista']

if args.check_protists:
    protists = collections.Counter(kingdom_tally['Protista'])
    for p in protists.most_common():
        print (p[0], p[1], sep='\t')
else:
    for k in sorted(kingdom_tally):
        print (k, len(kingdom_tally[k]), sep='\t')

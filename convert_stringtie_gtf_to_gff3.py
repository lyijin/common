#!/usr/bin/env python3

"""
> convert_gtf_to_gff3.py <

Given a StringTie-produced *.gtf file, convert it to a gff3 file that can be 
displayed properly in JBrowse. Nomenclature rules are elaborated below.

Expects 1+ *.gtf files, the output files replaces "gtf" with "gff3".

TODO: also produce a BigWig file containing the expression values of the 
StringTie exons.

TRANSCRIPT NAMES AND EXON NAMES ARE CONSISTENT ACROSS INPUT FILES!

Nomenclature rules:
1. If a StringTie transcript overlaps with only one gene model, then the 
   StringTie transcript inherits the name of the gene + a "Sxx" tag.
   This applies even if transcript is on opposite strand of the gene model!
     e.g. third StringTie transcript in Spis580 is "Spis580.S3".
   
   Transcripts that cut through multiple gene models would acquire the name of
   the gene that they overlap most extensively.
     e.g. if half of transcript is in Spis1337, a quarter in intergenic, 
          and last quarter in Spis1338, it would be called "Spis1337.Sxx".
   
   Transcripts that appear in intergenic regions are assigned a new gene name.
     e.g. SpisStringTie888

2. Exon naming: known exon names (as defined in the genome gff3) take precedence
   over new exons identified by StringTie.

3. For column 3 ("type"):
     "transcripts" and "exon" stays the same.
     For bin/flatfile-to-json.pl in JBrowse, use flag "--className transcript".
   For columns 4-8:
     Coordinates, directionality is preserved.
   For column 9 ("attributes"):
     transcript: ID=___;cov=___;FPKM=___
     exon: ID=___;Parent=___;cov=___;FPKM=___
     (note that cov and FPKM are non-standard fields)
"""
import argparse
import collections
import csv
import itertools
import re

import numpy as np

import natural_sort
import parse_fasta
import parse_gff3

parser = argparse.ArgumentParser(description="""
Given a StringTie-produced *.gtf file, convert it to a gff3 file that can be 
displayed PROPERLY in JBrowse. Nomenclature rules are elaborated in 
the source code.""")

parser.add_argument('species_abbreviation', metavar='species_abbreviation',
                    help='3/4-letter abbreviation of the species (e.g. "Spis")')
parser.add_argument('genome_fasta', metavar='fasta_filename',
                    type=argparse.FileType('r'),
                    help='Genome FASTA file.')
parser.add_argument('genome_gff3', metavar='gff3_filename',
                    type=argparse.FileType('r'),
                    help='Genome GFF3 annotation.')
parser.add_argument('stringtie_gtf', metavar='gtf_filename',
                    type=argparse.FileType('r'), nargs='+',
                    help='StringTie *.gtf file.')

args = parser.parse_args()

def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..., (sn-1, sn)"""
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

def get_single_transcript(gtf_file):
    """
    When given a gtf_file, this generator yields blocks containing one line of
    "transcript" and all the successive "exon" lines.
    
    gtf_file is a file object.
    """
    transcript_block = []
    for line, next_line in pairwise(gtf_file):
        if not line: continue
        
        transcript_block.append(line)
        if next_line:
            # produce output if next line = start of next sequence
            if next_line[2] == 'transcript':
                yield transcript_block
                transcript_block = []
    
    # add last line into block if it's not an empty row
    if next_line:
        transcript_block.append(next_line)
    yield transcript_block

def assign_value_to_numpy_array(np_array, coords, assigned_value):
    """
    Coords are in a tuple (coords1, coords2). +ve and -ve values are assigned
    respectively, depending on coords1 < coords2 or vice versa.
    
    Values are assigned in a first-come-first-serve mamner (this only
    concerns overlapping genes), i.e. alphabetically-first genes would have 
    entire sequence on the numpy array, but not the later ones.
    """
    for c in range(min(coords), max(coords)):
        if not np_array[c]:
            if coords[1] > coords[0]:
                np_array[c] = assigned_value
            elif coords[0] > coords[1]:
                np_array[c] = -assigned_value
    
    return np_array

def scrape_gtf_attribute(entire_attribute):
    """
    Given a *.gtf attribute,
      e.g. gene_id "spis.pH1.6"; exon_number "9"; cov "2.805369";
    scrape the data into a dictionary, and return it.
    """
    return dict(re.findall('(\S+) "(\S+)";', entire_attribute))

# read sequences
#'/home/liewy/kaust/dual_genome/raw_data/Spis.genome.scaffold.final.fa'
#'/home/liewy/kaust/dual_genome/raw_data/Spis.genome.annotation.gff3'
sequence_lengths = parse_fasta.get_all_sequences(args.genome_fasta, 'fasta', 
                                                 lengths_only=True)

# read coordinates of genes and exons from .gff3 file.
scaffold_gff3 = parse_gff3.parse_gff3(args.genome_gff3, 'exon')

# transcript_names stores the transcript names assigned to particular
# coordinates on the scaffold. if another transcript happens to share the exact
# same coordinates (all exons in transcript have the same coordinates),
# then the transcript will have the same name:
#   transcript_names['scaffold_name']['(start1, end1), (start2, end2), ...'] =\
#           'transcript_name'
transcript_names = {}

# exon_names functions similarly to transcript_names:
#   exon_names['scaffold_name'][(start, end)] = 'exon_name'
exon_names = {}

# genic regions are denoted in a NumPy array as follows:
#    0: intergenic region
#  +ve: gene in the + strand
#  -ve: gene in the - strand
#
# gene_names is a dict that converts ints in gene_info back to gene names.
#
# NOTE: if genes overlap, the first one that gets annotated (sorted 
# alphabetically) has priority.
gene_info = {}
gene_names = {}
gene_counter = 0

for scaf in sequence_lengths:
    gene_info[scaf] = np.zeros(sequence_lengths[scaf], np.int32)
    transcript_names[scaf] = {}
    exon_names[scaf] = {}

# pre-populate transcript_names and exon_names from the genome gff3
for s in natural_sort.natural_sort(scaffold_gff3):
    for gene in natural_sort.natural_sort(scaffold_gff3[s]):
        # gene detail annotation
        gene_coords = scaffold_gff3[s][gene].coords
        
        gene_counter += 1
        gene_info[s] = assign_value_to_numpy_array(gene_info[s], gene_coords, 
                                                   gene_counter)
        gene_names[gene_counter] = gene
        
        for tx in scaffold_gff3[s][gene].mRNAs:        
            # get all exon coordinates within the transcripts
            tx_coords = scaffold_gff3[s][gene].mRNAs[tx].details['exon']
            transcript_names[s][str(tx_coords)] = tx
        
            e_counter = 0
            for e in tx_coords:
                exon_start = min(e) + 1
                exon_end = max(e)
                e_counter += 1
                
                # if exon already has a name, do not annotate it again
                if (exon_start, exon_end) not in exon_names[s]:
                    exon_names[s][(exon_start, exon_end)] = '{}.exon{}'.format(
                        tx, e_counter)

transcript_counter = {}
exon_counter = {}

# read StringTie gtf. Make sure the first line is a "transcript" line!
for strg_gtf in args.stringtie_gtf:
    output_file = open(strg_gtf.name.replace('gtf', '') + 'gff3', 'w')
    
    tsv_reader = csv.reader(strg_gtf, delimiter='\t')
    for tx in get_single_transcript(tsv_reader):
        # sanity check: there are only lines with "transcript" and "exon"
        row_types = [x[2] for x in tx]
        assert set(row_types) == set(['transcript', 'exon']), \
            'rows contain more than just "transcript" and "exon"!'
    
        # tx is a transcript block: one line of "transcript" (tx[0]), 
        # followed by multiple lines of "exon" (tx[1:]).
        scaf = tx[0][0]
        tx_attrs = scrape_gtf_attribute(tx[0][8])
        # scaling_factor: the k in  FPKM = k . cov
        if float(tx_attrs['cov']) > 0:
            scaling_factor = float(tx_attrs['FPKM']) / float(tx_attrs['cov'])
        else:
            scaling_factor = 0
        
        if 'reference_id' in tx_attrs:
            # hallelujah, life is simpler.
            current_transcript = tx_attrs['reference_id']
        else:
            # check whether this transcript has existed in previous files
            tx_coords = []
            for row in tx[1:]:
                start = int(row[3])
                end = int(row[4])
                if row[6] == '+':
                    tx_coords.append((start - 1, end))
                else:
                    tx_coords.append((end, start - 1))
            
            tx_coords = str(tx_coords)
            if tx_coords in transcript_names[scaf]:
                current_transcript = transcript_names[scaf][tx_coords]
            else:
                # sigh. assign new name to the transcript.
                
                # check whether this transcript overlaps an existing gene
                start = int(tx[0][3]) - 1
                end = min(int(tx[0][4]), sequence_lengths[scaf])
                
                # if transcript_slice cuts across genes / intergenic regions,
                # it belongs to the gene that transcript_slice overlaps most.
                transcript_slice = collections.Counter(gene_info[scaf][start:end])
                
                # in determining overlaps, we only care about genic regions
                transcript_slice[0] = 0
                
                # check whether transcript overlaps at least part of a gene
                if transcript_slice.most_common(1)[0][1] > 0:
                    existing_gene = abs(transcript_slice.most_common(1)[0][0])
                    
                    # transcript is in an existing gene!
                    gene_id = gene_names[existing_gene]
                    
                    if gene_id not in transcript_counter:
                        transcript_counter[gene_id] = 0
                    
                    transcript_counter[gene_id] += 1
                    current_transcript = gene_id.replace('Gene', '') \
                        + '.S{}'.format(transcript_counter[gene_id])
                else:
                    if 'unknown' not in transcript_counter:
                        transcript_counter['unknown'] = 0
                    
                    transcript_counter['unknown'] += 1
                    current_transcript = '{}StringTie{}'.format(
                        args.species_abbreviation,
                        transcript_counter['unknown'])
                
                transcript_names[scaf][tx_coords] = current_transcript
        
        # print the "transcript" line out
        new_attr_line = 'ID={};cov={};FPKM={}'.format(
            current_transcript, tx_attrs['cov'], tx_attrs['FPKM'])
        print (tx[0][0], tx[0][1], 'transcript', '\t'.join(tx[0][3:8]), 
               new_attr_line, sep='\t', file=output_file)
        
        for row in tx[1:]:
            start = int(row[3])
            end = int(row[4])
            if (start, end) not in exon_names[scaf]:
                if current_transcript not in exon_counter:
                    exon_counter[current_transcript] = 0
                
                exon_counter[current_transcript] += 1
                exon_names[scaf][(start, end)] = '{}.exon{}'.format(
                    current_transcript, exon_counter[current_transcript])
                
            exon_id = exon_names[scaf][(start, end)]
            
            exon_attrs = scrape_gtf_attribute(row[8])
            new_attr_line = 'ID={};Parent={};cov={};FPKM={}'.format(
                exon_id, current_transcript, exon_attrs['cov'],
                round(scaling_factor * float(exon_attrs['cov']), 6))
            
            print (row[0], row[1], 'exon', '\t'.join(row[3:8]), 
                   new_attr_line, sep='\t', file=output_file)

    output_file.close()
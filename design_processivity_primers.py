#!/usr/bin/env python3

"""
> design_processivity_primers.py <

When given a gene name and the exon/intron number for which primers have to
be designed for, brute-force primers that fit all primer design criteria, as 
well as being are roughly spaced apart equally from each other.

Primer design is as follows:
1. 5' primer is in the exon of interest (i.e. exon N).
2. 3' primers are in the following intron (i.e. intron N). By default, 5
   primers are chosen in the intron and equally-spaced out, i.e. if intron N
   was of length 1, then the primers would be centered at
   [0.1, 0.3, 0.5, 0.7, 0.9].
   
As primers are ~20 bp in length, exons should be > 20 bp while introns should
be at least 20n bp, where n is number of primers designed within it. Less than
that, and primers will start to overlap. Too bad.
"""
import argparse
import csv

import numpy as np

import natural_sort
import parse_fasta
import parse_gff3
import primer3_api

parser = argparse.ArgumentParser(description="""
When given a gene name and the exon/intron number for which primers have to
be designed for, brute-force primers that fit all primer design criteria, as 
well as being are roughly spaced apart equally from each other.""")

parser.add_argument('genome_fasta', metavar='fasta_filename',
                    type=argparse.FileType('r'),
                    help='Genome FASTA file.')
parser.add_argument('genome_gff3', metavar='gff3_filename',
                    type=argparse.FileType('r'),
                    help='Genome GFF3 annotation.')
parser.add_argument('desired_amplicons', metavar='tsv_filename',
                    type=argparse.FileType('r'),
                    help='''Tab-delimited file containing gene and desired 
                            exon/intron (i.e. 'SpisGene10907    4'). One line 
                            per desired amplicon.''')
parser.add_argument('--intron_primers', '-i', type=int, default=5,
                    help='Number of primers designed within intron.')
parser.add_argument('--polyx', action='store_true', default=False,
                    help='Allow primers to have poly-X stretches of > 5.')
parser.add_argument('--wobble', type=int, default=20,
                    help='Change wobble window (default is 20).')

args = parser.parse_args()

# allow for primer to start +- WOBBLE_WINDOW base pairs of the intended site
WOBBLE_WINDOW = args.wobble

# SECTION DISABLED: PRIMER3 DOES NOT ACCEPT PRIMERS > 36nt. -__-"
# Illumina 16S Metagenomic Sequencing Kit sequences. The sequences are in:
#   5'--(overhang)--(amplicon-specific primer)--3'
# FWD_OVERHANG = 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
# REV_OVERHANG = 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'


def reverse_complement(seq):
    seq = seq.replace('U', 'T')

    translation_from = 'AaTtGgCcYyRrSsWwKkMmBbDdHhVvNn'
    translation_to   = 'TtAaCcGgRrYySsWwMmKkVvHhDdBbNn'
    translation_table = str.maketrans(translation_from, translation_to)

    seq = seq[::-1].translate(translation_table)

    return seq

def find_optimal_primer(scaffold_name, desired_location, primer_type,
                        gene_on_crick=False):
    '''
    desired_location is an int.
    primer_type is either of ['forward', 'reverse'].
    additional_overhang is a str, used in checking hairpin/dimers structures.
    
    This function checks for directionality, as hairpin structures have
    different melting points when reverse-complemented.
    
    Produces a list of 3 tuples, corresponding to top 3 most suitable primers
    close to desired_location. Data in tuple are:
        (primer_seq, start_location, length, Tm, GC%)
    '''
    assert primer_type in ['forward', 'reverse'], \
        'primer_type is either "forward" or "reverse".'
    
    strand = 'crick' if gene_on_crick else 'watson'
    possible_primers = []
    
    # grab all possible primers that fit criteria in WOBBLE_WINDOW around the
    # intended site
    for a in range(-WOBBLE_WINDOW, WOBBLE_WINDOW):
        # optimal primer lengths are 18-23 (primer3 defaults), but extended
        # to 36 when designing for bisulphite-converted DNA - GC% is low.
        max_primer_len = 37
        for p_len in range(18, max_primer_len):
            amp_start = a + desired_location
            
            scaffold_sequence = genome_sequences[scaffold_name]
            if (strand == 'watson' and primer_type == 'forward') \
                or (strand == 'crick' and primer_type == 'reverse'):
                p_temp = scaffold_sequence[amp_start:amp_start + p_len]
            else:
                p_temp = scaffold_sequence[amp_start - p_len:amp_start]
            
            # reverse-complement primer sequence if
            #   1. mRNA on Watson strand, reverse primer
            #   2. mRNA on Crick strand, forward primer
            if strand == 'watson' and primer_type == 'reverse':
                p_temp = reverse_complement(p_temp)
            elif strand == 'crick' and primer_type == 'forward':
                p_temp = reverse_complement(p_temp)
            
            # skip primers containing any degenerate bases (Y/R/...), primer3 
            # doesn't like it
            if p_temp.count('A') + p_temp.count('C') + p_temp.count('G') + \
                    p_temp.count('T') < len(p_temp):
                continue
            
            # check for consecutive Ns ("poly-X") (default: max of 5)
            if not args.polyx:
                if 'AAAAA' in p_temp or 'CCCCC' in p_temp \
                        or 'GGGGG' in p_temp or 'TTTTT' in p_temp:
                    continue
            
            p = primer3_api.check_individual_primers(p_temp)
            
            # secondary structure: exclude anything that evaluates to True
            bool_any_th = float(p[p_temp]['PRIMER_LEFT_0_SELF_ANY_TH']) > 45
            bool_end_th = float(p[p_temp]['PRIMER_LEFT_0_SELF_END_TH']) > 45
            bool_hairpin = float(p[p_temp]['PRIMER_LEFT_0_HAIRPIN_TH']) > 45
            
            # ignore secondary structure requirements for bis-converted primers
            if any([bool_any_th, bool_end_th, bool_hairpin]): continue
            
            # grab tm and gc
            tm = float(p[p_temp]['PRIMER_LEFT_0_TM'])
            gc = float(p[p_temp]['PRIMER_LEFT_0_GC_PERCENT'])
            if 57.5 < tm < 62.5 and 30 < gc < 70:
                possible_primers.append((p_temp, amp_start, p_len, tm, gc))

    # pick primer as close to the intended conditions:
    #   1. location: 0.2 penalty per base away from intended site (max WW * .2)
    #   2. length: 1 per base away from 20 (max 3) (primer3 default)
    #   3. tm: 1 per degree away from 60 (max 2.5) (primer3 default)
    #   4. gc: 0.01 per gc% away from 50 (max 0.2) (primer3 default=0)
    # retain top 3 forward primers that fit criteria
    penalty = lambda x: abs(x[1] - desired_location) * 0.2 + \
                        abs(x[2] - 20) * 1 + \
                        abs(x[3] - 60) * 1 + \
                        abs(x[4] - 50) * 0.01
    possible_primers = sorted(possible_primers, key=penalty)
    
    return possible_primers[:3]

def check_primer_pair(forward_primer, reverse_primer):
    '''
    Checks whether the pair of primers show high degree of complementarity.
    Returns a boolean True if all checks pass; False if not.
    '''
    p_pair = (forward_primer, reverse_primer)
    p = primer3_api.check_paired_primers((forward_primer, reverse_primer))
    
    # exclude anything that evaluates to True
    bool_any_th = float(p[p_pair]['PRIMER_PAIR_0_COMPL_ANY_TH']) > 45
    bool_end_th = float(p[p_pair]['PRIMER_PAIR_0_COMPL_END_TH']) > 45
    if any([bool_any_th, bool_end_th]):
        return False
    else:
        return True

def generate_relative_locations(n):
    '''
    Evenly split the range {0..1} depending on n (number of divisions), and 
    returns the midpoint of the sub-ranges.
    
    i.e. if n == 5, return [0.1, 0.3, 0.5, 0.7, 0.9]
    '''
    return [(x + 0.5)/n for x in range(n)]

# read sequences
genome_sequences = parse_fasta.get_all_sequences(args.genome_fasta, 'fasta')

# read coordinates of genes and exons from .gff3 file
scaffold_gff3 = parse_gff3.pick_longest_mRNA(
    parse_gff3.parse_gff3(args.genome_gff3, 'exon'))

# create dictionary to map genes to their respective scaffold (this is needed
# to obtain gene coords based solely on gene names)
gene_to_scaffold = {}
for s in scaffold_gff3:
    for g in scaffold_gff3[s]:
        gene_to_scaffold[g] = s
    
# print header row for results
print ('Gene', 'Intron relative location', 'Scaffold',
       'Desired region', 'Outer region', 
       'Outer amplicon length', 'Outer forward', 'Outer reverse', 
       'OF: loc | len | Tm | GC', 'OR: loc | len | Tm | GC',
       'Warnings', sep='\t')

# test gene + exon/intron number
# desired_amplicons = [('AIPGENE5756_gene', 6),   # watson
                     # ('AIPGENE21488_gene', 2)]  # crick

# read in desired amplicon
desired_amplicons = []
tsv_reader = csv.reader(args.desired_amplicons, delimiter='\t')
for row in tsv_reader:
    if not row: continue
    
    desired_amplicons.append((row[0], int(row[1])))

# create outer primers
for d in desired_amplicons:
    # get exon and intron coordinates
    gene = d[0]
    exon = d[1]
    
    scaf = gene_to_scaffold[gene]
    tx = list(scaffold_gff3[scaf][gene].mRNAs.keys())[0]
    
    gene_details = scaffold_gff3[scaf][gene].mRNAs[tx].details['exon']
    
    exon_coords = gene_details[exon - 1]
    intron_coords = (gene_details[exon - 1][1], gene_details[exon][0])
    
    gene_on_crick = exon_coords[0] > exon_coords[1]
    
    # set amplicon start (within exon) at 40 bp upstream of donor site
    if gene_on_crick:
        amplicon_start = exon_coords[1] + 40
    else:
        amplicon_start = exon_coords[1] - 40
    
    # set multiple amplicon ends within introns
    for relative_loc in generate_relative_locations(args.intron_primers):
        amplicon_end = int(intron_coords[0] + \
            relative_loc * (intron_coords[1] - intron_coords[0]))
        
        outer_starts = find_optimal_primer(scaf, amplicon_start, 
                                           'forward', gene_on_crick)
        outer_ends = find_optimal_primer(scaf, amplicon_end,
                                         'reverse', gene_on_crick)
        
        # if either of outer_starts or outer_ends does not exist, best to let user
        # re-pick desired locations, instead of over-engineering the script to 
        # widen search range etc.
        
        # proceed onwards if working primers are found
        desired_region = '..'.join([str(amplicon_start), str(amplicon_end)])
        error_msg = []
        if not outer_starts:
            error_msg.append('outer forward primer')
        if not outer_ends:
            error_msg.append('outer reverse primer')
        if error_msg:
            # dummy columns at the end makes it easier to merge files sideways
            print (d[0], exon + relative_loc, scaf, desired_region,
                   'ERROR: no valid ' + ' & '.join(error_msg), 
                   '', '', '', '', '', '', sep='\t')
            continue
        
        outer_fwd_primer = outer_starts[0]
        outer_rev_primer = outer_ends[0]
        outer_region = '{}..{}'.format(outer_fwd_primer[1], outer_rev_primer[1])
        
        warnings = []
        # check whether exonic primer lies outside of exon
        fwd_start = outer_fwd_primer[1]
        fwd_end = outer_fwd_primer[1] - outer_fwd_primer[2] if gene_on_crick else \
                  outer_fwd_primer[1] + outer_fwd_primer[2]
        if not min(exon_coords) < fwd_start < max(exon_coords):
            warnings.append('Exon primer starts outside desired exon')
        if not min(exon_coords) < fwd_end < max(exon_coords):
            warnings.append('Exon primer ends outside desired exon')
        
        # same for intronic primer
        rev_start = outer_rev_primer[1]
        rev_end = outer_rev_primer[1] + outer_rev_primer[2] if gene_on_crick else \
                  outer_rev_primer[1] - outer_rev_primer[2]
        if not min(intron_coords) < rev_start < max(intron_coords):
            warnings.append('Intron primer starts outside desired intron')
        if not min(intron_coords) < rev_end < max(intron_coords):
            warnings.append('Intron primer ends outside desired intron')
        
        # print everything out nicely
        outer_amplicon_len = abs(outer_fwd_primer[1] - outer_rev_primer[1])
        print (d[0], exon + relative_loc, scaf, desired_region, outer_region, 
               outer_amplicon_len, outer_fwd_primer[0], outer_rev_primer[0], 
               ' | '.join([str(x) for x in outer_fwd_primer[1:]]),
               ' | '.join([str(x) for x in outer_rev_primer[1:]]), 
               '; '.join(warnings), sep='\t')
#!/usr/bin/env python3

"""
> design_miseq_amplicon_primers.py <

When given a gene name and rough coordinates of the desired amplified region,
design primers that are suitable for Illumina's 16S Metagenomic Sequencing 
Library Preparation kit.

A nested approach to primer design is used:
1. Outer primer pair (~20 bp) amplifies a ~250 bp region.
2. Inner primer pair (~20 bp + general overhangs).

The amplicon amplified by the inner primer pair is the sequence that gets 
sequenced by the MiSeq.

The kit recommends the locus-specific portion of the primer to have a Tm of
60-65 deg C, while hairpin/dimer calculations should be based on the entire
primer (locus-specific + general overhangs).

When designing primers for bisulphite-converted DNA, remember that
non-methylated Cs on the primers (i.e. all of them) are converted to Ts!

Input requirements:
<gene_name> \t <scaffold_name> \t <start_coords> \t <end_coords>
(gene_name is required)
"""
import argparse
import csv

import numpy as np

import natural_sort
import parse_fasta
import parse_gff3
import primer3_api

parser = argparse.ArgumentParser(description="""
When given a gene name and rough coordinates of the desired amplified region,
design primers that are suitable for Illumina's 16S Metagenomic Sequencing 
Library Preparation kit.""")

parser.add_argument('genome_fasta', metavar='fasta_filename',
                    type=argparse.FileType('r'),
                    help='Genome FASTA file.')
parser.add_argument('genome_gff3', metavar='gff3_filename',
                    type=argparse.FileType('r'),
                    help='Genome GFF3 annotation.')
parser.add_argument('desired_amplicons', metavar='tsv_filename',
                    type=argparse.FileType('r'),
                    help='''Tab-delimited file containing gene and desired 
                            coordinates (in the format of "xxx..yyy"). One
                            line per desired amplicon.''')
parser.add_argument('--cov', metavar='cov_filename',
                    type=argparse.FileType('r'),
                    help='File with coordinates that primers cannot overlap.')
parser.add_argument('--bis', action='store_true', default=False,
                    help='Design primers for bisulphite-converted DNA.')
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

def find_optimal_primer(scaffold_name, desired_location, primer_type,
                        additional_overhang=''):
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
    
    if gene_info[scaffold_name][desired_location] > 0:
        strand = 'watson'
    elif gene_info[scaffold_name][desired_location] < 0:
        strand = 'crick'
    else:
        # sanity check
        print ('''Primers cannot be in intergenic regions, as it is impossible
                  to tell whether it's on the Watson or Crick strand.''')
        raise SystemExit()
    
    possible_primers = []
    
    # grab all possible primers that fit criteria in WOBBLE_WINDOW around the
    # intended site
    for a in range(-WOBBLE_WINDOW, WOBBLE_WINDOW):
        # optimal primer lengths are 18-23 (primer3 defaults), but extended
        # to 36 when designing for bisulphite-converted DNA - GC% is low.
        max_primer_len = 37 if args.bis else 24
        for p_len in range(18, max_primer_len):
            amp_start = a + desired_location
            
            scaffold_sequence = genome_sequences[scaffold_name]
            if (strand == 'watson' and primer_type == 'forward') \
                or (strand == 'crick' and primer_type == 'reverse'):
                gi = gene_info[scaffold_name][amp_start:amp_start + p_len]
                p_temp = scaffold_sequence[amp_start:amp_start + p_len]
            else:
                gi = gene_info[scaffold_name][amp_start - p_len:amp_start]
                p_temp = scaffold_sequence[amp_start - p_len:amp_start]
            
            # restrict p_temp to only contain ACGT; skip sequence otherwise
            if not all(c in 'ACGT' for c in p_temp): continue
            
            # check whether there's any +2 / -2 in gene_info in this window
            if any(gi == 2) or any(gi == -2):
                continue
            
            # reverse-complement primer sequence if
            #   1. mRNA on Watson strand, reverse primer
            #   2. mRNA on Crick strand, forward primer
            if strand == 'watson' and primer_type == 'reverse':
                p_temp = reverse_complement(p_temp)
            elif strand == 'crick' and primer_type == 'forward':
                p_temp = reverse_complement(p_temp)

            # Convert primers to work on bisulphite-converted DNA
            if args.bis and primer_type == 'forward':
                p_temp = p_temp.replace('C', 'T')
            if args.bis and primer_type == 'reverse':
                p_temp = p_temp.replace('G', 'A')
            
            # check for consecutive Ns ("poly-X") (default: max of 5)
            if not args.polyx:
                if 'AAAAA' in p_temp or 'CCCCC' in p_temp \
                    or 'GGGGG' in p_temp or 'TTTTT' in p_temp:
                    continue
            
            p = primer3_api.check_individual_primers(p_temp)
            
            if additional_overhang:
                p2 = primer3_api.check_individual_primers(additional_overhang \
                                                          + p_temp)
            else:
                p2 = p
            
            # secondary structure: exclude anything that evaluates to True
            bool_any_th = float(p2[p_temp]['PRIMER_LEFT_0_SELF_ANY_TH']) > 45
            bool_end_th = float(p2[p_temp]['PRIMER_LEFT_0_SELF_END_TH']) > 45
            bool_hairpin = float(p2[p_temp]['PRIMER_LEFT_0_HAIRPIN_TH']) > 45
            
            # ignore secondary structure requirements for bis-converted primers
            if not args.bis:
                if any([bool_any_th, bool_end_th, bool_hairpin]): continue
            
            # grab tm and gc
            tm = float(p[p_temp]['PRIMER_LEFT_0_TM'])
            gc = float(p[p_temp]['PRIMER_LEFT_0_GC_PERCENT'])
            # relax GC%-requirements when designing for bis-converted DNA
            if args.bis:
                if 55 < tm < 65:
                    possible_primers.append((p_temp, amp_start, p_len, tm, gc))
            else:
                if 60 < tm < 65 and 30 < gc < 70:
                    possible_primers.append((p_temp, amp_start, p_len, tm, gc))

    # pick primer as close to the intended conditions:
    #   1. location: 0.2 penalty per base away from intended site (max WW * .2)
    #   2. length: 1 per base away from 20 (max 3) (primer3 default)
    #   3. tm: 1 per degree away from 62.5 (max 2.5) (primer3 default)
    #   4. gc: 0.01 per gc% away from 50 (max 0.2) (primer3 default=0)
    # retain top 3 forward primers that fit criteria
    penalty = lambda x: abs(x[1] - desired_location) * 0.2 + \
                        abs(x[2] - 20) * 1 + \
                        abs(x[3] - 62.5) * 1 + \
                        abs(x[4] - 50) * 0.01
    possible_primers = sorted(possible_primers, key=penalty)
    
    return possible_primers[:3]

def check_primer_pair(forward_primer, reverse_primer, 
                      additional_overhangs=('', '')):
    '''
    Checks whether the pair of primers show high degree of complementarity.
    Returns a boolean True if all checks pass; False if not.
    '''
    p_pair = (forward_primer, reverse_primer)
    if additional_overhangs:
        p2 = primer3_api.check_paired_primers(
                (additional_overhangs[0] + forward_primer,
                 additional_overhangs[1] + reverse_primer))
    else:
        p2 = primer3_api.check_paired_primers((forward_primer, reverse_primer))
    
    # exclude anything that evaluates to True
    bool_any_th = float(p2[p_pair]['PRIMER_PAIR_0_COMPL_ANY_TH']) > 45
    bool_end_th = float(p2[p_pair]['PRIMER_PAIR_0_COMPL_END_TH']) > 45
    if any([bool_any_th, bool_end_th]):
        return False
    else:
        return True
    
# read sequences
genome_sequences = parse_fasta.get_all_sequences(args.genome_fasta, 'fasta')

# read coordinates of genes and exons from .gff3 file.
scaffold_gff3 = parse_gff3.parse_gff3(args.genome_gff3, 'gene')

# genic regions are denoted in a NumPy array as follows:
#    0: intergenic region
#   +1: gene in the + strand
#   -1: gene in the - strand
#
#   +2: gene in the + strand, POSITION EXCLUDED FROM PRIMERS
#   -2: gene in the - strand, POSITION EXCLUDED FROM PRIMERS
# NOTE: if genes overlap, the first one that gets annotated (sorted 
# alphabetically) has priority.
gene_info = {}

# gene_location['gene_name'] = 'scaffold_name'
gene_location = {}

for scaf in genome_sequences:
    gene_info[scaf] = np.zeros(len(genome_sequences[scaf]), np.int8)

# annotate the NumPy array with genes
for scaf in natural_sort.natural_sort(scaffold_gff3):
    for gene in natural_sort.natural_sort(scaffold_gff3[scaf]):
        gene_coords = scaffold_gff3[scaf][gene].coords
        
        gene_info[scaf] = assign_value_to_numpy_array(
            gene_info[scaf], gene_coords, 1)
            
        gene_location[gene] = scaf

# read in information from .cov file, disallow primers from overlapping
# methylated DNA / edited RNA bases (optional).
if args.cov:
    tsv_reader = csv.reader(args.cov, delimiter='\t')

    for row in tsv_reader:
        if not row: continue
        
        # mark location as non-primer-overlappable
        gene_info[row[0]][int(row[1]) - 1] *= 2

# print header row for results
print ('Gene', 'Scaffold', 'Desired region', 'Desired length',
       'Outer region', 'Inner region',
       'Outer amplicon length', 'Inner amplicon length',
       'Outer forward', 'Inner forward', 'Inner reverse', 'Outer reverse', 
       'OF: loc | len | Tm | GC', 'IF: loc | len | Tm | GC',
       'IR: loc | len | Tm | GC', 'OR: loc | len | Tm | GC', sep='\t')

# test gene/coordinate. coordinates are 1-based.
# test_data = [('SpisGene12521', '138325..138525'), # watson
             # ('SpisGene10907', '218900..219150')] # crick
# gene_info['Spis.scaffold224|size490686'][138541] = 2

# read in desired amplicon
desired_amplicons = []
tsv_reader = csv.reader(args.desired_amplicons, delimiter='\t')
for row in tsv_reader:
    if not row: continue
    
    desired_amplicons.append(row)

# create outer primers
for d in desired_amplicons:
    gene, scaf = d[0], d[1]
    amplicon_start, amplicon_end = int(d[2]), int(d[3])
    
    # convert 1-based coordinates to 0-based
    amplicon_start, amplicon_end = \
        sorted([amplicon_start, amplicon_end])
    amplicon_start -= 1
    
    amplicon_scaf = gene_location[d[0]]
    if gene_info[amplicon_scaf][amplicon_start] < 0:
        # gene is in Crick strand
        amplicon_start, amplicon_end = amplicon_end, amplicon_start
    
    # if either of outer_starts or outer_ends does not exist, best to let user
    # re-pick desired locations, instead of over-engineering the script to 
    # widen search range etc.
    outer_starts = find_optimal_primer(amplicon_scaf, amplicon_start, 
                                       'forward')
    outer_ends = find_optimal_primer(amplicon_scaf, amplicon_end,
                                     'reverse')
    
    desired_region = '{}..{}'.format(amplicon_start, amplicon_end)
    desired_length = abs(amplicon_end - amplicon_start)
    
    # print ('os:', outer_starts)
    # print ('oe:', outer_ends)
    
    # proceed onwards if working primers are found.
    error_msg = []
    if not outer_starts:
        error_msg.append('outer forward primer')
    if not outer_ends:
        error_msg.append('outer reverse primer')
    
    if error_msg:
        print (gene, scaf, desired_region, desired_length,
               'ERROR: no valid ' + ' & '.join(error_msg), sep='\t')
        continue
    
    outer_inner_ok = False
    for s in outer_starts:
        for e in outer_ends:
            # test the primer pair in primer3 for inter-primer annealing
            if check_primer_pair(s[0], e[0]):
                # create new inners based on the good outers. only accept
                # outers if they allow for inners to be designed.
                if s[1] < e[1]:
                    inner_starts = find_optimal_primer(
                        amplicon_scaf, s[1] + s[2], 'forward')
                    inner_ends = find_optimal_primer(
                        amplicon_scaf, e[1] - e[2], 'reverse')
                else:
                    inner_starts = find_optimal_primer(
                        amplicon_scaf, s[1] - s[2], 'forward')
                    inner_ends = find_optimal_primer(
                        amplicon_scaf, e[1] + e[2], 'reverse')
                
                # print ('is:', inner_starts)
                # print ('ie:', inner_ends)
                
                if inner_starts and inner_ends:
                    # test the primer pairs for inter-primer annealing
                    for ss in inner_starts:
                        for ee in inner_ends:
                            if check_primer_pair(ss[0], ee[0]):
                                outer_inner_ok = True
                                break
                        if outer_inner_ok: break
            if outer_inner_ok: break
        if outer_inner_ok: break
    
    outer_fwd_primer = s
    outer_rev_primer = e
    outer_region = '{}..{}'.format(outer_fwd_primer[1], outer_rev_primer[1])
    
    if not outer_inner_ok:
        # failed to find inner/outer pair
        print (gene, scaf, desired_region, desired_length, outer_region, 
               'ERROR: No valid inner primers.', sep='\t')
        continue
    
    inner_fwd_primer = ss
    inner_rev_primer = ee
    
    # print everything out nicely.
    outer_amplicon_len = abs(outer_fwd_primer[1] - outer_rev_primer[1])
    inner_amplicon_len = abs(inner_fwd_primer[1] - inner_rev_primer[1])
    inner_region = '{}..{}'.format(inner_fwd_primer[1], inner_rev_primer[1])    
    print (gene, scaf, desired_region, desired_length,
           outer_region, inner_region,
           outer_amplicon_len, inner_amplicon_len,
           outer_fwd_primer[0], inner_fwd_primer[0], 
           inner_rev_primer[0], outer_rev_primer[0], 
           ' | '.join([str(x) for x in outer_fwd_primer[1:]]),
           ' | '.join([str(x) for x in inner_fwd_primer[1:]]),
           ' | '.join([str(x) for x in inner_rev_primer[1:]]),
           ' | '.join([str(x) for x in outer_rev_primer[1:]]), sep='\t')
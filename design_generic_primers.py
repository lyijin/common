#!/usr/bin/env python3

"""
> design_generic_primers.py <

When given a scaffold name and rough coordinates of the desired amplicon region,
design primers that fit generic primer design criteria.

The script doesn't do sanity checks (e.g. testing that amplicons are of
reasonable size) so it's up to the user to vet the primers for their
suitability!

Input requirements:
<gene/amplicon_name> \t <scaffold_name> \t <start_coord> \t <end_coord>
"""
import argparse
import csv

import natural_sort
import parse_fasta
import primer3_api

parser = argparse.ArgumentParser(description="""
When given a scaffold name and rough coordinates of the desired amplicon region,
design primers that fit generic primer design criteria.""")

parser.add_argument('genome_fasta', metavar='fasta_filename',
                    type=argparse.FileType('r'),
                    help='Genome FASTA file.')
parser.add_argument('desired_amplicons', metavar='tsv_filename',
                    type=argparse.FileType('r'),
                    help='''Tab-delimited file with four columns:
                            <gene/amplicon_name> <scaffold_name>
                            <start_coord> <end_coord>.
                            One line per desired amplicon.''')
parser.add_argument('--polyx', action='store_true', default=False,
                    help='Allow primers to have poly-X stretches of > 5.')
parser.add_argument('--wobble', type=int, default=20,
                    help='Change wobble window (default is 20).')

args = parser.parse_args()

# allow for primer to start +- WOBBLE_WINDOW base pairs of the intended site
WOBBLE_WINDOW = args.wobble

MAX_PRIMER_LEN = 37

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
        for p_len in range(18, MAX_PRIMER_LEN):
            amp_start = a + desired_location
            
            scaffold_sequence = genome_sequences[scaffold_name]
            if (strand == 'watson' and primer_type == 'forward') \
                or (strand == 'crick' and primer_type == 'reverse'):
                p_temp = scaffold_sequence[amp_start:amp_start + p_len]
            else:
                p_temp = scaffold_sequence[amp_start - p_len:amp_start]
            
            # restrict p_temp to only contain ACGT; skip sequence otherwise
            if not all(c in 'ACGT' for c in p_temp): continue
            
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

# read sequences
genome_sequences = parse_fasta.get_all_sequences(args.genome_fasta, 'fasta')

# print header row for results
print ('Gene', 'Scaffold',
       'Desired region', 'Desired length', 
       'Outer region', 'Outer amplicon length',
       'Outer forward', 'Outer reverse', 
       'OF: loc | len | Tm | GC', 'OR: loc | len | Tm | GC',
       'Warnings', sep='\t')

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
    gene_on_crick = amplicon_start > amplicon_end
    if gene_on_crick:
        amplicon_end -= 1
    else:
        amplicon_start -= 1
    
    # if either of outer_starts or outer_ends does not exist, best to let user
    # re-pick desired locations, instead of over-engineering the script to 
    # widen search range etc.
    outer_starts = find_optimal_primer(scaf, amplicon_start, 
                                       'forward', gene_on_crick)
    outer_ends = find_optimal_primer(scaf, amplicon_end,
                                     'reverse', gene_on_crick)
    
    # proceed onwards if working primers are found
    desired_region = '{}..{}'.format(amplicon_start, amplicon_end)
    desired_length = abs(amplicon_end - amplicon_start)
    
    error_msg = []
    if not outer_starts:
        error_msg.append('outer forward primer')
    if not outer_ends:
        error_msg.append('outer reverse primer')
    
    if error_msg:
        # dummy columns at the end makes it easier to merge files sideways
        print (gene, scaf, desired_region, desired_length,
               'ERROR: no valid ' + ' & '.join(error_msg), 
               '', '', '', '', '', '', sep='\t')
        continue
    
    outer_fwd_primer = outer_starts[0]
    outer_rev_primer = outer_ends[0]
    outer_region = '{}..{}'.format(outer_fwd_primer[1], outer_rev_primer[1])
    
    warnings = []
    
    # print everything out nicely
    outer_amplicon_len = abs(outer_fwd_primer[1] - outer_rev_primer[1])
    print (gene, scaf, desired_region, desired_length, outer_region, 
           outer_amplicon_len, outer_fwd_primer[0], outer_rev_primer[0], 
           ' | '.join([str(x) for x in outer_fwd_primer[1:]]),
           ' | '.join([str(x) for x in outer_rev_primer[1:]]), 
           '; '.join(warnings), sep='\t')
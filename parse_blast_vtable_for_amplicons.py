#!/usr/bin/env python3

docstring = """
> parse_blast_vtable_for_amplicons.py <

After BLASTN-ing against the genome to check where primers could potentially
bind, parse the vtable output from `parse_blast_xml.py` and check for Fs and Rs
that bind within ~1 kb of each other (changeable option --max_amplicon_length).

Helps assuage concerns that primers might be amplifying off-target regions
(e.g., pseudogenes, or just purely by chance).

Example input
  CDKN2A_F  gnl|BL_ORD_ID|17  chr9_CTOB  [...]  21974802  21974781  [...]  1.31e-03
  CDKN2A_F  gnl|BL_ORD_ID|26  chr14_OT   [...]  97811438  97811456  [...]  8.06e-02
  CDKN2A_R  gnl|BL_ORD_ID|20  chr11_OT   [...]  66559302  66559284  [...]  4.03e-02
  CDKN2A_R  gnl|BL_ORD_ID|17  chr9_CTOB  [...]  21974733  21974751  [...]  4.03e-02

into
  'CDKN2A': [('chr9_CTOB', 21974802, 21974733), ('somewhere_else', 1234, 1235)]

Assume primers can be paired purely by changing FINAL char from F to R.
Multiple primers for the same gene can be accommodated, but script assumes the
second primer for CDKN2A would be either
  CDKN2A_2_F
  CDKN2A_2F
""".strip()

import argparse
import csv
import gzip
from pathlib import Path
import sys
import time

def diag_print(*args):
    """
    Helper function to print diagnostic-level verbose output to stderr.
    """
    print (f'[{time.asctime()}]', *args, file=sys.stderr)

def calculate_tuple_dists(f_coord_list, r_coord_list, max_size):
    """
    Calculates distances between all F primers and R primers.
    
    Returns list of F/R primers that create an amplicon below the predefined
    size.
    """
    potential_amplicons = []
    for f_coords in f_coord_list:
        for r_coords in r_coord_list:
            # the trick in figuring out amplicon size is that... it's just the
            # distance between the first values in both tuples, as exemplified
            # in the docstring of this script:
            #   F: (21974802, 21974781); R: (21974733, 21974751)
            #   i.e., amplicon: (21974802, 21974733)
            amplicon_size = abs(f_coords[0] - r_coords[0]) + 1
            
            if amplicon_size < max_size:
                potential_amplicons.append([
                    min(f_coords[0], r_coords[0]),  # start of amplicon
                    max(f_coords[0], r_coords[0]),  # end of amplicon
                    amplicon_size])
    
    return potential_amplicons


parser = argparse.ArgumentParser(
    description=docstring, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('vtable_file', metavar='tsv_file', type=Path,
                    help='output from `parse_blast_xml.py` --vtable format.')
parser.add_argument('--amplicon_size', '-a', metavar='bp', type=int, default=1000,
                    help='set maximum amplicon size (default 1k bp) for F/R pairing.')
parser.add_argument('--e_value', '-e', metavar='e_value', type=float,
                    help='set maximum E value.')
parser.add_argument('--coverage', '-c', metavar='coverage_pct', type=float,
                    help='set minimum coverage (in %%) for individual hits.')
parser.add_argument('--identity', '-i', metavar='identity_pct', type=float,
                    help='set minimum local identity (in %%).')
parser.add_argument('--bit_score', '-b', metavar='bit_score', type=float,
                    help='set minimum total bit score value.')
parser.add_argument('--verbose', '-v', action='store_true',
                    help='prints diagnostic stuff to stderr.')
args = parser.parse_args()

if args.vtable_file.suffix == '.gz':
    # assume gzip-compressed file
    tsv_reader = csv.reader(gzip.open(args.vtable_file, 'rt'), delimiter='\t')
else:
    # assume it's plaintext
    tsv_reader = csv.reader(open(args.vtable_file, 'r'), delimiter='\t')

# create dict to store vtable details, and for downstream pairing
query_dict = {}

if args.verbose:
    diag_print(f'Starting to parse {args.vtable_file}...')

# start reading tsv
for row in tsv_reader:
    # skip header row
    if row[0] == 'Query': continue
    
    # valid rows should have 15 cols, skip rows that are not valid
    if len(row) < 15: continue
    
    # filter based on specified cutoffs
    remove_flag = False
    if args.e_value:
        e_value = float(row[14])
        remove_flag = args.e_value < e_values
    
    if args.coverage:
        coverage_pct = float(row[13].replace('%', ''))
        remove_flag = args.coverage > coverage_pct
    
    if args.identity:
        identity_pct = float(row[12].replace('%', ''))
        remove_flag = args.identity > identity_pcts
    
    if args.bit_score:
        bit_score = float(row[10])
        remove_flag = args.bit_score > bit_score
    
    if remove_flag: continue
    
    # filters (if any) passed; check for F and R primers in close proximity
    query_annot = row[0]
    
    # raise a stink if query annot doesn't end with either F or R
    primer_f_or_r = query_annot[-1:].upper()
    assert primer_f_or_r in ['F', 'R'], \
        f'{query_annot} must end with either F or R.'
    
    # infer gene name by lopping off two chars if it ends with "_F" or "_R",
    # but only lop one char off if it's non-standard e.g., "_4F", "_4R"
    gene_name = query_annot[:-2] if query_annot[-2:] in ['_F', '_R', '_f', '_r'] \
                else query_annot[:-1]
    
    if gene_name not in query_dict:
        query_dict[gene_name] = {}
    
    hit_annot = row[2]
    if hit_annot not in query_dict[gene_name]:
        query_dict[gene_name][hit_annot] = {'F': [], 'R': []}
    
    hit_start = int(row[7])
    hit_end = int(row[8])
    
    query_dict[gene_name][hit_annot][primer_f_or_r].append((hit_start, hit_end))

if args.verbose:
    diag_print(f'Done parsing {args.vtable_file}.')

# check which Fs are close to which Rs
for gene_name in query_dict:
    for hit_annot in query_dict[gene_name]:
        potential_amplicons = calculate_tuple_dists(
            query_dict[gene_name][hit_annot]['F'],
            query_dict[gene_name][hit_annot]['R'],
            args.amplicon_size)
        
        if not potential_amplicons: continue
        
        # print the potential amplicon out
        for pa in potential_amplicons:
            print (gene_name, hit_annot, pa[0], pa[1], pa[2], sep='\t')

if args.verbose:
    diag_print(f'Done predicting potential amplicons.')

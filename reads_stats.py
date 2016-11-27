#!/usr/bin/env python3

"""
> reads_stats.py <

Python script aims to do two things:
- Calculates the length of the reads and returns the distribution (histogram) of
  those lengths.
- For every read length, calculates the nucleotide distribution at each position
  along the read. Useful to check whether adapter trimming has been too
  successful (eg cutadapt only requires a minimum match length of 3 nts, which
  might lead to the removal of non-adapter sequences from 5'/3' ends), or..
  reveal new patterns!
"""

import argparse

parser = argparse.ArgumentParser(description="""
Script calculates the length of the reads and returns the distribution
(histogram) of those lengths; for every read length, calculates the nucleotide
distribution at each position along the read; also, outputs the top 1000 most
abundant reads. """)

parser.add_argument('reads_file', metavar='reads_file', type=argparse.FileType('r'),
                    nargs=1, help='FASTA/FASTQ/collapsed FASTA file')
parser.add_argument('--include', metavar='included_reads', type=argparse.FileType('r'),
                    nargs='+', help='only calculate stats of reads in the given files')
parser.add_argument('--exclude', metavar='excluded_reads', type=argparse.FileType('r'),
                    nargs='+', help='exclude stats of reads in the given files')
parser.add_argument('--top_seqs', metavar='number_of_top_sequences', type=int,
                    help='print top_seqs sequences at the end of the file. \
                    Default is 1000. Use 0 to print all sequences.', default=1000)
fasta_opt = parser.add_mutually_exclusive_group(required=True)
fasta_opt.add_argument('--fasta', action='store_const', dest='file_format',
                       const='fasta', help='input file is in FASTA format.')
fasta_opt.add_argument('--fastq', action='store_const', dest='file_format',
                       const='fastq', help='input file is in FASTQ format.')
fasta_opt.add_argument('--cfasta', action='store_const',
                       dest='file_format', const='collapsed_fasta',
                       help='input file is in collapsed FASTA format.')

args = parser.parse_args()

# file exists, counting begins
read_lengths = {}       # read_lengths[length] = total number
nucleotide_stats = {}   # nucleotide_stats[length] = {1: {'A': m, 'T': n ...},
                        #                             2: {'A': m2, 'T': n2 ...}}
VALID_NUCLEOTIDES = ['A', 'T', 'U', 'G', 'C', 'N']

# parse input file, and convert the option "cfasta" to "fasta" because
# parse_fasta.py doesn't discriminate between those two formats
import parse_fasta
sequences = parse_fasta.get_all_sequences(args.reads_file[0], args.file_format[-5:])

# handle read inclusions/exclusions
if args.include == None:
    included_reads = set(sequences.keys())
else:
    included_reads = set()
    for i in args.include:
        included_reads |= set([x.strip() for x in i])
    
    # make sure that the reads exist in the original reads file
    included_reads &= set(sequences.keys())
                                        
if args.exclude != None:
    for e in args.exclude:
        included_reads -= set([x.strip() for x in e])

included_sequences = [] # needed to tally most-occuring sequences
# note: included_seq isn't a collections.Counter() object as appending to the
#       object is horribly slow.

read_freq = 1           # default for FASTA/FASTQ files: every annot = 1 read

for i in included_reads:
    seq = sequences[i]
    seq = seq.rstrip().upper()
    
    # for collapsed_fasta files, each annot will have different # reads.
    if args.file_format == 'collapsed_fasta':
        read_freq = int(i.split('x')[-1])
    
    included_sequences += [seq] * read_freq     # not Counter(), see above
    
    s = len(seq)
    
    # tallying sequence lengths
    if s not in read_lengths:
        read_lengths[s] = read_freq
        nucleotide_stats[s] = {}
        for l in range(1, s + 1):
            nucleotide_stats[s][l] = {}
            for v in VALID_NUCLEOTIDES:
                nucleotide_stats[s][l][v] = 0
    else:
        read_lengths[s] += read_freq

    # tallying nucleotide distribution
    nt_location = 0
    for nt in seq:
        nt_location += 1
        nucleotide_stats[s][nt_location][nt] += read_freq

list_of_lengths = list(read_lengths.keys())
list_of_lengths.sort()

# output
print ('Distribution of read lengths:')
print ('length', 'count', sep='\t')
for lol in list_of_lengths:
    print ('\t'.join([str(lol), str(read_lengths[lol])]))
print ('sum', sum([read_lengths[x] for x in read_lengths]), sep='\t')
    
print ('', end='\n\n')

print ('Nucleotide stats:', end='\n\n')

for lol in list_of_lengths:
    print ('Stats for sequences of length', lol)
    print ('', *range(1, lol + 1), sep='\t')
    for v in VALID_NUCLEOTIDES:
        output_list = []
        for l in range(1, lol + 1):
            output_list.append(nucleotide_stats[lol][l][v])
        
        # omit printing entire line if that nucleotide has 0 reads for that length
        if sum(output_list) > 0:
            print (v, *output_list, sep='\t')
    
    print ('')

print ('', end='\n\n')

# print top n reads as well
import collections
raw_sequences_counter = collections.Counter(included_sequences)

ts = args.top_seqs
if ts == 0:
    print ('All reads:', end='\n\n')
    for m in raw_sequences_counter.most_common():
        print ('\t'.join([str(x) for x in m]))
else:
    print ('Top', ts, 'reads:', end='\n\n')
    for m in raw_sequences_counter.most_common(ts):
        print ('\t'.join([str(x) for x in m]))

#!/usr/bin/env python3

"""
> pick_reads.py <

Script takes in FASTA/FASTQ file, picks reads that are in the included list (by
default, all reads are included) and removes those in the excluded list.

Note that --include operates on a FILE, while --contain operates on string(s).
"""

import re
import sys

import natural_sort

def reverse_complement(seq):
    seq = seq.replace('U', 'T')

    translation_from = 'AaTtGgCcYyRrSsWwKkMmBbDdHhVvNn'
    translation_to   = 'TtAaCcGgRrYySsWwMmKkVvHhDdBbNn'
    translation_table = str.maketrans(translation_from, translation_to)

    seq = seq[::-1].translate(translation_table)

    return seq

def slice_sequence(sequence, startpos, endpos):
    """
    Function takes in 1-based genomic coordinates, and slices the input
    sequence based on these coordinates (after sanity checks).
    
    If startpos is greater than endpos, the sequence is also reverse-
    complemented.
    """
    # handle edge cases
    if startpos < 1 and endpos < 1:
        return ''
    
    # set numerical boundaries to startpos and endpos
    startpos = max(1, startpos)
    startpos = min(startpos, len(sequence))
    
    endpos = max(1, endpos)
    endpos = min(endpos, len(sequence))
    
    # reverse complement sequence flag: True when startpos > endpos
    revcomp_flag = startpos > endpos
    if revcomp_flag:
        startpos, endpos = endpos, startpos
    
    sliced_sequence = sequence[startpos-1:endpos]
    if revcomp_flag: sliced_sequence = reverse_complement(sliced_sequence)
    
    return sliced_sequence

def wrap_sequence(sequence, width):
    return '\n'.join(sequence[i:i+width] for i in range(0, len(sequence), width))

if __name__ == '__main__':
    import argparse
    
    import parse_fasta
    
    parser = argparse.ArgumentParser(description="""
    Script takes in FASTA file, picks reads that are in the included list (by
    default, all reads are included) and removes those in the excluded list.""")
    
    parser.add_argument('reads_file', metavar='reads_file',
                        type=argparse.FileType('r'), nargs='?',
                        default=sys.stdin, help='FASTA/collapsed FASTA file')
    parser.add_argument('--fastq', action='store_true',
                        help='input is a FASTQ file, not FASTA')
    parser.add_argument('--include', metavar='reads_file',
                        type=argparse.FileType('r'), nargs='+',
                        help='include reads with annots in the given files')
    parser.add_argument('--exclude', metavar='reads_file',
                        type=argparse.FileType('r'), nargs='+',
                        help='exclude reads with annots in the given files')
    parser.add_argument('--contain', metavar='list_of_strings',
                        type=str, nargs='+',
                        help='pick reads that contain specified string')
    parser.add_argument('--donotcontain', metavar='list_of_strings',
                        type=str, nargs='+',
                        help='pick reads that do not contain specified string')
    parser.add_argument('--min', metavar='n', type=int,
                        help='exclude reads with length < n bp')
    parser.add_argument('--max', metavar='n', type=int,
                        help='exclude reads with length > n bp')
    parser.add_argument('--start', metavar='n', type=int,
                        help='get sequences from position n')
    parser.add_argument('--end', metavar='n', type=int,
                        help='get sequences to position n')
    parser.add_argument('--wrap', metavar='n', type=int,
                        help='wrap sequences with line width of n')
    parser.add_argument('--longest', metavar='regex_string',
                        help='pick longest read based on annotation')
    parser.add_argument('--alphasort', action='store_true',
                        help='sort sequences alphabetically (chr11, chr2)')
    parser.add_argument('--nosort', action='store_true',
                        help='disable sorting of read file')
    parser.add_argument('--order_include', action='store_true',
                        help='follow order of singular file in --include')

    
    args = parser.parse_args()
    
    # sanity checking
    if bool(args.start) != bool(args.end):
        raise ValueError('--start and --end has to be used in conjuction!')
    
    if args.order_include:
        # order_include can only be True when include is True
        if not args.include: args.order_include == False
        if len(args.include) > 1: args.order_include == False
    
    # start of script - get sequence data
    if args.fastq:
        sequences = parse_fasta.get_all_sequences(args.reads_file, 'fastq')
    else:
        sequences = parse_fasta.get_all_sequences(args.reads_file, 'fasta')
    
    # handle read inclusions/exclusions
    if not args.include:
        included_reads = set(sequences.keys())
    else:
        included_reads = set()
        for i in args.include:
            included_reads |= set([x.strip() for x in i])
        
        # make sure that the reads exist in the original reads file
        included_reads &= set(sequences.keys())

    if args.exclude:
        for e in args.exclude:
            included_reads -= set([x.strip() for x in e])

    # checks annotations and picks any read that contains any of args.contains
    if args.contain:
        included_reads = [x for x in included_reads if
                          any(c in x for c in args.contain)]
    
    if args.donotcontain:
        included_reads = [x for x in included_reads if
                          all(c not in x for c in args.donotcontain)]
    
    if args.min:
        included_reads = [x for x in included_reads if
                          len(sequences[x]) >= args.min]

    if args.max:
        included_reads = [x for x in included_reads if
                          len(sequences[x]) <= args.max]

    if args.longest:
        # create dictionary that creates abbreviations for the annotations
        # (and based on these groupings, pick the longest one in the group!)
        short_annot = {}
        for r in included_reads:
            # assuming that common annotations share a same numerical
            # identifier near the start of the sequence, e.g. 'maker-12345'
            s = re.search(args.longest, r).group(1)

            if s not in short_annot:
                short_annot[s] = []
            short_annot[s].append(r)

        included_reads = []
        for s in short_annot:
            longest_read = max(short_annot[s], key=lambda x: len(sequences[x]))
            included_reads.append(longest_read)

            # debug:
            # print (longest_read, len(sequences[longest_read]), [len(sequences[x]) for x in short_annot[s]])

    if args.alphasort:
        sorted_reads = sorted(included_reads)
    elif args.nosort:
        sorted_reads = [s for s in sequences if s in included_reads]
    elif args.order_include:
        # rewind to start of file
        args.include[0].seek(0)
        sorted_reads = [x.strip() for x in args.include[0]]
    else:
        sorted_reads = natural_sort.natural_sort(included_reads)
    
    for i in sorted_reads:
        print('>{}'.format(i))
        seq = sequences[i]
        if args.start and args.end:
            seq = slice_sequence(seq, args.start, args.end)
            
        if args.wrap:
            print (wrap_sequence(seq, args.wrap))
        else:
            print (seq)

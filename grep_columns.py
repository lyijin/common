#!/usr/bin/env python3

"""
> grep_columns.py <

Given a query file (tab-separated) and a target file (tab-separated), by
default, searches the specified column in target file for all lines that 
contains EXACT MATCHES to any values of the first column of the query file.

-c (zero-based) allows for multi-column matching;
-r (string) controls what the search string is in the target file.
"""
import argparse
import csv
import re
import sys

parser = argparse.ArgumentParser(description='''
Given a query file (tab-separated) and a target file (tab-separated), by
default, searches the specified column in target file for all lines that 
contains EXACT MATCHES to any values of the first column of the query file.''')

parser.add_argument('query_file', metavar='query_file',
                    type=argparse.FileType('r'), 
                    help='query file containing one query per line.')
parser.add_argument('target_file', metavar='target_file',
                    type=argparse.FileType('r'), 
                    help='tab-separated target file.')
parser.add_argument('--column', '-c', metavar='n', type=int,
                    nargs='+', default=[0],
                    help='''When no arguments supplied: print if target[0] == 
                            query[0]. If single n supplied: print if 
                            target[n1] == query[0]. If multiple n supplied, 
                            print if target[n1] == query[0],
                            target[n2] == query[1], ...''')
parser.add_argument('--numerical', '-n', action='store_true',
                    help='''assumes that query file contains lines in the
                            pattern query \\t x \\t y in each line. Check
                            whether target[n1] matches query, and whether
                            column n + 1 satisfies the inequality
                            x <= target[n1 + 1] <= y. Print if it does. 
                            NOTE THE "n1": THIS FLAG ONLY USES FIRST PROVIDED
                            VALUE OF -c!!''')
parser.add_argument('--regex', '-r', metavar='regex_string',
                    help='''use regex string to define search target in target 
                            file. Ignores -c and -n flag.''')

args = parser.parse_args()

# read queries
queries = []
tsv_reader = csv.reader(args.query_file, delimiter='\t')
for row in tsv_reader:
    if not row: continue
    
    if args.numerical:
        query, start, end = row
        start = int(start)
        end = int(end)
        queries.append([query, start, end])
    else:
        queries.append(row[:len(args.column)])

# read and print lines in target file
tsv_reader = csv.reader(args.target_file, delimiter='\t')
for row in tsv_reader:
    if not row: continue
    if len(row) < max(args.column): continue
    
    if args.regex:
        s = re.search(args.regex, '\t'.join(row))
        if not s: continue
        
        # do crude detection of whether s.group(0) or s.group(1) should be used,
        # which depends on whether brackets appear in the regex string.
        if '(' in args.regex and ')' in args.regex:
            if s.group(1) not in queries: continue
        else:
            if s.group(0) not in queries: continue
    else:
        # very inefficient, but oh well. computing power wasted for the
        # sake of convenience
        if args.numerical:
            valid_seq = False
            c = args.column[0]
            for q in queries:
                if row[c] == q[0] and q[1] <= int(row[c+1]) <= q[2]:
                    valid_seq = True
                    break
            
            if not valid_seq: continue
        else:
            query_string = [row[x] for x in args.column]
            if query_string not in queries: continue
    
    print ('\t'.join(row))

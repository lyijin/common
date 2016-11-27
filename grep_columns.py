#!/usr/bin/env python3

"""
> grep_columns.py <

Given a query file (one string per line) and a target file (tab-separated), 
searches the target file for all lines that contains EXACT MATCHES to the 
queries.

-c (zero-based) controls which column the search is performed (default is 0);
-r (string) controls what the search string is in the target file.
"""
import argparse
import csv
import re
import sys

parser = argparse.ArgumentParser(description='''
Given a query file (one string per line) and a target file (tab-separated), 
searches the target file for all lines that contains EXACT MATCHES to the 
queries.''')

parser.add_argument('query_file', metavar='query_file',
                    type=argparse.FileType('r'), 
                    help='query file containing one query per line.')
parser.add_argument('target_file', metavar='target_file',
                    type=argparse.FileType('r'), 
                    help='tab-separated target file.')
parser.add_argument('--column', '-c', metavar='n', type=int, default=0,
                    help='check column n for query in target file.')
parser.add_argument('--numerical', '-n', action='store_true',
                    help='''assumes that query file contains lines in the
                            pattern x..y in each line. Check whether value in 
                            column n satisfies the inequality x < n < y.
                            Print if it does.''')
parser.add_argument('--regex', '-r', metavar='regex_string',
                    help='''use regex string to define search target in target 
                            file. Ignores -c and -n flag.''')

args = parser.parse_args()

# read queries
queries = []
with args.query_file as f:
    for line in f:
        line = line.strip()
        if not line: continue
        
        if args.numerical:
            queries.append([float(x) for x in line.split('..')])
        else:
            queries.append(line)

# read and print lines in target file
c = args.column

tsv_reader = csv.reader(args.target_file, delimiter='\t')
for row in tsv_reader:
    if not row: continue
    if len(row) < c: continue
    
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
        if args.numerical:
            valid_seq = False
            rc = float(row[c])
            
            for q in queries:
                if min(q) < rc < max(q):
                    valid_seq = True
            
            if not valid_seq: continue
        else:
            if row[c] not in queries: continue
    
    print ('\t'.join(row))

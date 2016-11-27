#!/usr/bin/env python3

"""
> calc_max_min_avg.py <

Given one/more files containing a list of numbers, calculate the max, min,
mean, median and mode of the numbers.
"""
import argparse
import csv
import statistics

import natural_sort

parser = argparse.ArgumentParser(description="""
Given one/more files containing a list of numbers, calculate the max, min,
mean, median and mode of the numbers.""")

parser.add_argument('tsv_file', metavar='tsv_file',
                    type=argparse.FileType('r'), nargs='+',
                    help='one/more input lists.')
args = parser.parse_args()

def int_or_float(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

# read data
data = {}
for t in args.tsv_file:
    data[t.name] = []
    tsv_reader = csv.reader(t, delimiter='\t')
    for row in tsv_reader:
        if not row: continue
        
        data[t.name].append(int_or_float(row[0]))

# header row
print ('File', 'Max', 'Min', 'Mean', 'Mode', 'Median', sep='\t')
for d in natural_sort.natural_sort(data):
    # print results out
    print (d, max(data[d]), min(data[d]), round(statistics.mean(data[d]), 3),
           round(statistics.mode(data[d]), 3),
           round(statistics.median(data[d]), 3), sep='\t')

#!/usr/bin/env python3

"""
> convert_csv_to_tsv.py <

Does what it says on the tin - converts csv to tsv.
"""
import argparse
import csv

import natural_sort

parser = argparse.ArgumentParser(description="""
Does what it says on the tin - converts csv to tsv.""")

parser.add_argument('csv_file', metavar='csv_file',
                    type=argparse.FileType('r'),
                    help='a csv file.')
args = parser.parse_args()

# read data
csv_reader = csv.reader(args.csv_file)
for row in csv_reader:
    print ('\t'.join(row))

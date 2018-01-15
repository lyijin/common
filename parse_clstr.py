#!/usr/bin/env python3

"""
> parse_clstr.py <

Parses cd-hit-est's .clstr files, produces a dictionary containing clusters
in the following format:
  cluster[representative sequence] = [sequences_that_belong_to_the_cluster]

and prints it out to standard output.
"""
import argparse
import re

import natural_sort

def parse_clstr(clstr_file):
    # start parsing
    clstr = clstr_file.read()
    
    # splits clstr into clusters of clusters (ha!)
    clstr = clstr.split('>Cluster ')
    clstr.pop(0)
    
    # start populating 'cluster' dictionary
    cluster = {}
    
    
    for indiv_clusters in clstr:
        seqs_in_cluster = []
        representative_seq = ''
        for line in indiv_clusters.split('\n'):
            if '...' in line:
                r = re.search('>(.*?)\.\.\.', line).group(1)
                seqs_in_cluster.append(r)
                
                if '... *' in line:
                    representative_seq = r
            
        cluster[representative_seq] = seqs_in_cluster
    
    return cluster
    
parser = argparse.ArgumentParser(description="""
Parses cd-hit-est's .clstr files""")

parser.add_argument('clstr_file', metavar="clstr_filename",
                    type=argparse.FileType('r'), help="clstr filename.")

args = parser.parse_args()

clusters = parse_clstr(args.clstr_file)

for c in natural_sort.natural_sort(clusters):
    d = natural_sort.natural_sort(clusters[c])
    print (c, ','.join(d), sep='\t')

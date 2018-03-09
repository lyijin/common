#!/usr/bin/env python3

"""
> parse_clstr.py <

Contains a function to parse cd-hit-est's .clstr files and returns a dictionary
containing clusters in the following format:
  cluster[representative sequence] = [sequences_that_belong_to_the_cluster]

If the script is called in a standalone manner, it outputs a table containing
the relative frequencies of each cluster (in descending order).
"""
import argparse
import re

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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    Contains a function to parse cd-hit-est's .clstr files, and outputs a table
    containing the relative frequencies of each cluster (in descending order).""")
    
    parser.add_argument('clstr_file', metavar='clstr_filename',
                        type=argparse.FileType('r'),
                        help="name of clstr file.")

    args = parser.parse_args()
    
    clusters = parse_clstr(args.clstr_file)
    
    # get cluster values, so that we can find the total number of clusters
    cluster_values = [len(clusters[x]) for x in clusters]
    total = sum(cluster_values)
    
    # print results out
    for c in sorted(clusters, key=lambda x: len(clusters[x])):
        print (c, len(clusters[c]), round(len(clusters[c]) / total, 4), sep='\t')

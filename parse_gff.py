#!/usr/bin/env python

"""
> parse_gff.py <

Python script intends to be a helper function that can be called by other
scripts to handle gff annotations.

Script requires a _file object_ (not filename).

The output is in the form of:
    dict[seqname] = {five_prime_UTR: {ID001: (start_coord1, end_coord1), ...}
                     three_prime_UTR: {ID001: (start_coord2, end_coord2), ...}
                     gene: ...
                     exon: ...
                     CDS: ...
                     ...: ...}
"""

# from http://genome.ucsc.edu/FAQ/FAQformat.html#format3:
# Here is a brief description of the GFF fields:
#   seqname - The name of the sequence. Must be a chromosome or scaffold. 
#   source - The program that generated this feature. 
#   feature - The name of this type of feature. Some examples of standard
#             feature types are "CDS", "start_codon", "stop_codon", and "exon". 
#   start - The starting position of the feature in the sequence. The first base
#           is numbered 1. 
#   end - The ending position of the feature (inclusive). 
#   score - A score between 0 and 1000. If the track line useScore attribute is
#           set to 1 for this annotation data set, the score value will
#           determine the level of gray in which this feature is displayed
#           (higher numbers = darker gray). If there is no score value, enter ".". 
#   strand - Valid entries include '+', '-', or '.' (for don't know/don't care). 
#   frame - If the feature is a coding exon, frame should be a number between
#           0-2 that represents the reading frame of the first base. If the
#           feature is not a coding exon, the value should be '.'. 
#   group - All lines with the same group are linked together into a single item.

# sample line:
#   RNA-1;1.gff:1000131     maker   five_prime_UTR  1968    1999    .   
#   +       .       ID=maker-1000131-exonerate...

# IMPORTANT NOTE ABOUT startpos AND endpos: the biological definition of
# position (starting base numbered 1) is BAD at dealing with the edge case
# of startpos == endpos. If startpos == endpos, then there's no way to tell
# whether it's on the '+' strand or the '-' strand.
# Thus, startpos and endpos produced by the parse_gff def follows the Python
# convention (starting base is sequence[0:1]), which allows for the
# discrimination of '+' and '-' strands.
#   e.g. if startpos = 2 and endpos = 3: '+' strand    }  both refer to
#        if endpos = 2 and startpos = 3: '-' strand    }  the same base

import csv
import re

def get_attribute(gff3_row, attribute):
    attr_search = re.search('{}=(.*?)(;|$)'.format(attribute), gff3_row[8])
    if attr_search:
        return attr_search.group(1)
    else:
        raise AttributeError("'{}' does not contain '{}='.".format(
            '\t'.join(gff3_row), attribute))

def parse_gff(gff_file, select_feature='all', dict_key='seqname'):
    """
    'gff_file' refers to file object containing gff file.
    'select_feature' can be used to select for one or more features of interest
    in the gff file (e.g. "three_prime_UTR", ['three_prime_UTR', 'five_prime_UTR'])
    'dict_key" refers to which variable should be used as the key to the dict
    returned - allowed values are ['seqname', 'group_id', 'group_parent'].
    """

    # note: all variables are named according to convention!

    if dict_key not in ['seqname', 'group_id', 'group_parent', 'group_name']:
        raise ValueError("dict_key can only take one of four values: " +\
                         "'seqname', 'group_id', 'group_parent' or 'group_name'.")

    gff_details = {}
    tsv_reader = csv.reader(gff_file, delimiter='\t')
    for row in tsv_reader:
        # ignore blank lines and comments (lines starting with '#')
        if not row: continue
        if row[0].startswith('#'): continue
        
        feature = row[2]
        # ignore lines that do not correspond to the feature wanted:
        if select_feature == 'all' or feature in select_feature:
            # by default, seqname's value is "seqname" itself, but modify the
            # value based on which dictionary key is preferred!
            seqname = row[0]
            if dict_key == 'group_id':
                seqname = get_attribute(row[8], 'ID')
            elif dict_key == 'group_parent':
                seqname = get_attribute(row[8], 'Parent')
            elif dict_key == 'group_name':
                seqname = get_attribute(row[8], 'Name')

            if seqname not in gff_details:
                gff_details[seqname] = {}

            if feature not in gff_details[seqname]:
                gff_details[seqname][feature] = {}
            
            strand = row[6]
            
            # note that startpos/endpos follow programming convention,
            # not biological convention (see note before def)
            if strand == '+':
                startpos = int(row[3]) - 1
                endpos = int(row[4])
            elif strand == '-':
                startpos = int(row[4])
                endpos = int(row[3]) - 1
            else:
                print ('Error:', line.strip(), 'has no valid strand information.')
                raise SystemExit()
            
            feature_id = get_attribute(row, 'ID')
            gff_details[seqname][feature][feature_id] = (startpos, endpos)

    return gff_details
    
if __name__ == '__main__':
    print (parse_gff(open('test.gff3')))
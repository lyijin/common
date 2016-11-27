#!/usr/bin/env python3

"""
> parse_uniprot_dat.py <

Python script to parse through UniProt's *.dat files, and then output sequences
in chunks till '//' is encountered.

Script requires a _file object_ (not filename) and filetype.
"""

import re

def iter_sequences(dat_file):
    """
    A generator that handles UniProt *.dat that follows the pattern of:
      ID xxxx    (starting line)
      ...        (random stuff)
      //         (ending line)
    """    
    data = ''
    for line in dat_file:
        if line[:2] == '//':
            yield data
            data = ''
        else:
            data += line

def parse_uniprot_seq(seq):    
    id = re.search('^ID\s+(\w+_\w+)', seq, re.M).group(1)
    ac = re.search('^AC\s+(\w+);', seq, re.M).group(1)
    src = re.search('^DT.*?UniProtKB/([\w|-]+)', seq, re.M).group(1)
    sv = re.search('sequence version (\d+)', seq).group(1)
    de = re.search('[Rec|Sub]Name: Full=(.*?);', seq).group(1)
    try:
        gn = re.search('^GN\s+\w+=(\w+)', seq, re.M).group(1)
    except:
        gn = ''
    os = re.findall('^OS\s+(.*?)\n', seq, re.M)
    oc = re.findall('^OC\s+(.*?)\n', seq, re.M)
    pe = re.search('^PE\s+(\d+)', seq, re.M).group(1)
    kw = re.findall('^KW\s+(.*?)\n', seq, re.M)
    sq = re.search('^SQ.*;(.*)', seq, re.DOTALL + re.M).group(1)
    
    # postprocessing
    os = ' '.join(os)[:-1]      # removes the '.' at the end of the string
    oc = ' '.join(oc)[:-1]
    kw = ' '.join(kw)
    sq = sq.replace(' ', '').replace('\n', '')
    
    parsed_data = {'id': id, 'ac': ac, 'src': src, 'sv': sv, 'de': de,
                   'gn': gn, 'os': os, 'oc': oc, 'pe': pe, 'kw': kw, 'sq': sq}
    
    return parsed_data
#!/usr/bin/env python3

"""
> resolve_species.py <

When given a file containing species of interest, query eol.org and obtain 
more info e.g. kingdom, phylum for those species.

API usage is documented at http://eol.org/api.
"""

import argparse
import json
import os
import sqlite3
import sys
import time
import urllib.request

def query_eol_ping():
    """
    Checks that the service is up.
    """
    query_url = 'http://eol.org/api/ping/1.0.json'
    
    json_output = urllib.request.urlopen(query_url).read()
    json_output = json_output.decode('utf-8', 'ignore')
    
    json_dict = json.loads(json_output)
    
    if json_dict['response']['message'] == 'Success':
        return True
    else:
        return False
    
def query_eol_search(species):
    """
    Searches the EOL database for the species. Expects full species name,
    as it expects an EXACT match.
    """
    # remove 'sp.' in the species name, broadens search space a bit more
    #species = species.replace(' sp.', '')
    # line disabled because it caused 504 Timeout Error on "Bradyrhizobium sp."
    
    # queries require spaces to be replaced by '+'
    species = species.replace(' ', '+')
    
    query_url = 'http://eol.org/api/search/1.0.json?' + \
                'q={}&page=1&exact=true'.format(species)
    
    json_output = urllib.request.urlopen(query_url).read()
    json_output = json_output.decode('utf-8', 'ignore')
    
    json_dict = json.loads(json_output)
    
    return json_dict

def query_eol_pages(species_id):
    """
    Searches the EOL database for details associated with a given species ID.
    """    
    query_url = 'http://eol.org/api/pages/1.0/' + \
                '{}.json?'.format(species_id) + \
                'images=0&videos=0&text=0'
    
    json_output = urllib.request.urlopen(query_url).read()
    json_output = json_output.decode('utf-8', 'ignore')
    
    json_dict = json.loads(json_output)
    
    return json_dict

def query_eol_hierarchy_entries(hierarchy_id):
    """
    Searches the EOL database for details associated with a given species ID.
    """    
    query_url = 'http://eol.org/api/hierarchy_entries/1.0/' + \
                '{}.json?'.format(hierarchy_id) + \
                'common_names=false&synonyms=false'
    
    json_output = urllib.request.urlopen(query_url).read()
    json_output = json_output.decode('utf-8', 'ignore')
    
    json_dict = json.loads(json_output)
    
    return json_dict

def parse_eol_search(json_dict):
    """
    Expects a JSON dictionary produced by query_eol_search(), then parses 
    the dictionary for the species ID.
    """
    
    if not json_dict['results']:
        return 'not found'
    
    # if there's more than one match, query all and preferentially
    # select for the highest richness score
    if len(json_dict['results']) > 1:
        poss_ids = [x['id'] for x in json_dict['results']]            
        richness_scores = {}
        for p in poss_ids:
            richness_scores[p] = float(query_eol_pages(p)['richness_score'])

        species_id = max(richness_scores, key=richness_scores.get)
    else:
        species_id = json_dict['results'][0]['id']

    return species_id

def parse_eol_pages(json_dict):
    """
    Expects a JSON dictionary produced by query_eol_pages(), then parses 
    the dictionary for the best species hierarchy ID.
    
    Preferred sources (stored in 'nameAccordingTo') are as follows:
    1. NCBI Taxonomy
    2. Integrated Taxonomic Information System (ITIS)
    3. any of the others
    
    As sources might contain nothing in the 'ancestors' key (i.e.
    hierarchy_dict['ancestors'] = []), the for loop goes over the possible 
    hierarchy IDs and exits when the first one returns valid hierarchy
    information.
    """
    preferred_sources = ['NCBI Taxonomy',
                         'Integrated Taxonomic Information System (ITIS)']
    possible_ids = []
    valid_info = ''
    
    for p in preferred_sources + [x['nameAccordingTo'] for x in json_dict['taxonConcepts']]:
        if valid_info: break
        
        pref_details = [x for x in json_dict['taxonConcepts']
                        if x['nameAccordingTo'] == p]
        
        for pd in pref_details:
            possible_ids.append(pd['identifier'])
        
        possible_ids = sorted(possible_ids, reverse=True)
        
        for pi in possible_ids:
            if valid_info: break
            
            hierarchy_dict = query_eol_hierarchy_entries(pi)
            
            if hierarchy_dict['ancestors']:
                # note: this is [] for non-valid hierarchy_dicts
                # valid hierarchy information has been found!
                valid_info = hierarchy_dict

    return valid_info

def parse_eol_hierarchy_entries(json_dict):
    """
    Expects a JSON dictionary produced by query_eol_hierarchy_entries(), 
    then parses the dictionary for hierarchy information.
    
    At the moment, subphylum/subkingdom/sub-xxx stuff is ignored; only
    KPCOFGS is recorded.
    """
    controlled_vocab = ['kingdom', 'phylum', 'class', 'order', 'family',
                        'genus', 'species']
    kpcofgs = dict([(x, '') for x in controlled_vocab])
    
    # also store source of names
    kpcofgs['authority'] = json_dict['nameAccordingTo'][0]
    
    # store all 'scientificName' in a string
    all_sn = ' > '.join([x['scientificName'] for x in json_dict['ancestors']])
    kpcofgs['full_taxonomy'] = all_sn
    
    # store query time
    kpcofgs['query_time'] = int(time.time())
    
    for a in json_dict['ancestors']:
        if 'taxonRank' not in a: continue
        
        if a['taxonRank'] in controlled_vocab:
            kpcofgs[a['taxonRank']] = a['scientificName']
    
    return kpcofgs

def get_kpcofgs(species_of_interest, verbose=False):
    """
    One def to rule them all.
    """
    eol_search_json = query_eol_search(species_of_interest)
    species_id = parse_eol_search(eol_search_json)
    
    if species_id == 'not found':
        if verbose:
            print ('{} not found!'.format(species_of_interest), file=sys.stderr,
                   flush=True)
        return {'query_time': int(time.time()),
                'authority': 'not found'}
    
    eol_pages_json = query_eol_pages(species_id)
    valid_info = parse_eol_pages(eol_pages_json)
    
    # there's a chance that valid_info is empty - happens when all sources
    # fail to list a single ancestor for the species
    if not valid_info:
        if verbose:
            print ('{} not found!'.format(species_of_interest), file=sys.stderr,
                   flush=True)
        return {'query_time': int(time.time()),
                'authority': 'not found'}
                
    kpcofgs = parse_eol_hierarchy_entries(valid_info) 

    return kpcofgs

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    When given a file containing species of interest, query eol.org and obtain 
    more info e.g. kingdom, phylum for those species.""")

    parser.add_argument('species_file', metavar='species_name',
                        type=argparse.FileType('r'), nargs='?',
                        default=sys.stdin, 
                        help='file containing species names, one per line.')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='be verbose when querying EOL.')
    parser.add_argument('-f', '--force', action='store_true',
                        help='disable use of cache, force query EOL.')
    args = parser.parse_args()

    # re-query EOL when entry is older than 180 days
    MAX_QUERY_AGE = 15552000
    if args.force: MAX_QUERY_AGE = 0

    # if query_eol_ping() is False:
        # raise SystemExit('EOL website is down!')

    line_total = sum(1 for line in args.species_file)
    line_counter = 0

    args.species_file.seek(0)
    with args.species_file as f:
        script_abspath = os.path.abspath(__file__)
        sqdb_filename = script_abspath[:-2] + 'sqdb'
        con = sqlite3.connect(sqdb_filename)
        cur = con.cursor()
        
        # create a species_details table if it doesn't exist
        cur.execute('CREATE TABLE IF NOT EXISTS species_details \
                        (species TEXT PRIMARY KEY, \
                         query_time INTEGER, \
                         authority TEXT, \
                         tax_kingdom TEXT, \
                         tax_phylum TEXT, \
                         tax_class TEXT, \
                         tax_order TEXT, \
                         tax_family TEXT, \
                         tax_genus TEXT, \
                         full_taxonomy TEXT)')
        
        con.commit()
        con.close()
        
        for line in f:
            line_counter += 1
            species_of_interest = line.strip()
            
            if not species_of_interest:
                # no data in line, skip to next one
                print ()
                continue
            
            # connect to the sqdb every loop, so even if script crashes in the
            # middle, EOL doesn't need to be re-queried for entries prior to crash.
            con = sqlite3.connect(sqdb_filename)
            cur = con.cursor()
            
            # check whether species_of_interest exists in the sqlite3 table -- if
            # it does, then return the cached info!
            q = cur.execute('SELECT * FROM species_details WHERE \
                             species="{}"'.format(species_of_interest)).fetchone()
            
            if q:
                query_age = int(time.time()) - q[1]
                print ('[{}/{}] {} exists in cache ({} days old).'.format(
                        line_counter, line_total, species_of_interest,
                       round(query_age/86400,2)), file=sys.stderr, flush=True)
                
                # delete entry if older than MAX_QUERY_AGE
                if query_age > MAX_QUERY_AGE:
                    cur.execute('DELETE FROM species_details WHERE \
                                 species="{}"'.format(species_of_interest))
            
            # make query if entry doesn't exist anymore
            if not q or query_age > MAX_QUERY_AGE:
                print ('[{}/{}] Querying for {}...'.format(
                        line_counter, line_total, species_of_interest), 
                       file=sys.stderr, flush=True)
                kpcofgs = get_kpcofgs(species_of_interest, verbose=args.verbose)
                
                # cache all searches, whether failed or successful, avoid hammering
                # EOL with failed searches
                if kpcofgs['authority'] == 'not found':
                    cur.execute('INSERT INTO species_details VALUES \
                                 ("{}", {}, "{}", "{}", "{}", "{}", "{}", "{}", \
                                  "{}", "{}")'.format(species_of_interest, 
                                                      kpcofgs['query_time'], 
                                                      kpcofgs['authority'],
                                                      '-', '-', '-', '-', '-',
                                                      '-', '-'))
                else:
                    cur.execute('INSERT INTO species_details VALUES \
                                 ("{}", {}, "{}", "{}", "{}", "{}", "{}", "{}", \
                                  "{}", "{}")'.format(species_of_interest, 
                                                      kpcofgs['query_time'], 
                                                      kpcofgs['authority'], 
                                                      kpcofgs['kingdom'], 
                                                      kpcofgs['phylum'], 
                                                      kpcofgs['class'], 
                                                      kpcofgs['order'], 
                                                      kpcofgs['family'], 
                                                      kpcofgs['genus'], 
                                                      kpcofgs['full_taxonomy']))
                    
                q = cur.execute('SELECT * FROM species_details WHERE \
                                 species="{}"'.format(species_of_interest)).fetchone()
                
                # flooding the EOL server is bad!
                time.sleep(10)
            
            # second column is query_time -- do not print this column
            output = [x for n, x in enumerate(q) if n != 1]
            print ('\t'.join(output))
        
            con.commit()
            con.close()

#!/usr/bin/env python3

"""
> primer3_api.py <

Contains a collection of functions that make it much easier to use Primer3,
which is natively written in C.
"""
import subprocess
import tempfile

PRIMER3_PATH = '/home/liewy/tools/primer3-2.3.6/src/'
PRIMER3_CORE = PRIMER3_PATH + 'primer3_core'
PRIMER3_CONFIG = PRIMER3_PATH + 'primer3_config/'

def test_primer3_exists():
    """
    Tests whether primer3 exists on the system. If it does, return the
    version of primer3 used (using the -about flag in the program).
    
    If primer3 isn't at said path, script will crash (FileNotFoundError).
    """
    primer3_ver = subprocess.check_output([PRIMER3_CORE, '-about'], 
                                          universal_newlines=True)
    
    return primer3_ver.strip()

def parse_primer3_output(primer3_output):
    """
    A dict is created using SEQUENCE_PRIMER as keys, and other properties
    (e.g. PRIMER_LEFT_0_TM, PRIMER_LEFT_0_GC_PERCENT) as values.
    """
    primer3_data = {}
    for po in primer3_output.split('=\n'):
        if not po: continue
        
        temp_dict = {}
        for line in po.split('\n'):
            # each line should be in the pattern of
            #   PROPERTY=VALUE
            if '=' not in line: continue
            
            prop, val = line.split('=')
            temp_dict[prop] = val
        
        if 'SEQUENCE_PRIMER_REVCOMP' not in temp_dict:
            dict_key = temp_dict['SEQUENCE_PRIMER']
        else:
            dict_key = (temp_dict['SEQUENCE_PRIMER'],
                        temp_dict['SEQUENCE_PRIMER_REVCOMP'])
        primer3_data[dict_key] = temp_dict
    
    return primer3_data

def check_individual_primers(primer_sequences):
    """
    primer_sequences can be either a string (single primer) or a list
    (multiple primers). For each individual primer, it is run through primer3,
    then output is parsed into a dict.
    """
    if isinstance(primer_sequences, str):
        primer_sequences = [primer_sequences]
    
    primer3_input = ''
    for p in primer_sequences:
        primer3_input += 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH={}\n'.format(PRIMER3_CONFIG) + \
                         'SEQUENCE_PRIMER={}\n'.format(p) + \
                         'PRIMER_TASK=check_primers\n' + \
                         'PRIMER_PICK_ANYWAY=1\n=\n'
    
    t = tempfile.NamedTemporaryFile()
    t.write(primer3_input.strip().encode())
    t.seek(0)
    primer3_output = subprocess.check_output([PRIMER3_CORE, t.name],
                                             universal_newlines=True)
    t.close()
    
    parsed_data = parse_primer3_output(primer3_output)
    
    return parsed_data

def check_paired_primers(primer_sequences):
    """
    primer_sequences can be either a tuple (single pair) or a list of tuples
    (multiple primers). For each individual primer, it is run through primer3,
    then output is parsed into a dict.
    
    Tuple should contain sequences in the form (forward, reverse).
    """
    if isinstance(primer_sequences, tuple):
        primer_sequences = [primer_sequences]
    
    primer3_input = ''
    for p in primer_sequences:
        primer3_input += 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH={}\n'.format(PRIMER3_CONFIG) + \
                         'SEQUENCE_PRIMER={}\n'.format(p[0]) + \
                         'SEQUENCE_PRIMER_REVCOMP={}\n'.format(p[1]) + \
                         'PRIMER_TASK=check_primers\n' + \
                         'PRIMER_PICK_ANYWAY=1\n=\n'
    
    t = tempfile.NamedTemporaryFile()
    t.write(primer3_input.strip().encode())
    t.seek(0)
    primer3_output = subprocess.check_output([PRIMER3_CORE, t.name],
                                             universal_newlines=True)
    t.close()
    
    parsed_data = parse_primer3_output(primer3_output)
    
    return parsed_data
    
if __name__ == '__main__':
    # print (test_primer3_exists())
    
    print (check_individual_primers(['TACGATCGATGCTACGAGTAA', 'TACGATCGATGCTACGAGTAAA', 'TACGATCGATGCTACGAGTAAAA']))
    print (check_paired_primers(('TACTAGCGCTAGTCGACGTAC', 'TGACGAGCAGCTGTGTGTGA')))
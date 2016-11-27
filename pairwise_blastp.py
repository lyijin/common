#!/usr/bin/env python3

"""
> pairwise_blastp.py <

Script takes in two FASTA files of equal number of sequences, and does pairwise
BLASTP (#1 vs. #1, #2 vs. #2, ...). Prints results to stdout.

The FASTA file cannot contain duplicates, as this will introduce off-by-one
comparison errors.
"""
import argparse
import subprocess
import tempfile

import parse_fasta

parser = argparse.ArgumentParser(description="""
Script takes in two FASTA files of equal number of sequences, and does pairwise
BLASTP (#1 vs. #1, #2 vs. #2, ...). Prints results to stdout.""")
parser.add_argument('protein_fastas', metavar='fasta_file',
                    type=argparse.FileType('r'), nargs=2,
                    help='Pair of protein FASTA files.')

args = parser.parse_args()

first_fasta = parse_fasta.get_all_sequences(args.protein_fastas[0], 'fasta')
second_fasta = parse_fasta.get_all_sequences(args.protein_fastas[1], 'fasta')

assert len(first_fasta) == len(second_fasta), \
    'number of unique sequences in both files are different!'

# all systems go!
first_fasta = list(first_fasta.items())
second_fasta = list(second_fasta.items())

# header line
print ('Query', 'Hit accession', 'Hit description','Query length', 'Hit length', 
       'Query (start, end)', 'Hit (start, end)', 'Frame', 'Max bit score',
       'Total bit score', 'Identity', 'Identity %', 'Coverage %', 'Expect',
       sep='\t')

for n in range(len(first_fasta)):
    # create dummy files for NCBI BLASTP to operate on
    first_seq = '>{}\n{}'.format(*first_fasta[n])
    second_seq = '>{}\n{}'.format(*second_fasta[n])
    
    t = tempfile.NamedTemporaryFile()
    t.write(first_seq.encode())
    t.seek(0)
    
    u = tempfile.NamedTemporaryFile()
    u.write(second_seq.encode())
    u.seek(0)
    
    blastp_command = ['blastp', '-query', t.name, '-subject', u.name, '-outfmt', '5']
    v = tempfile.NamedTemporaryFile()
    v.write(subprocess.check_output(blastp_command))
    v.seek(0)
    
    pbx_command = ['parse_blast_xml.py', v.name, '--table', '--noheader']
    output = subprocess.check_output(pbx_command, universal_newlines=True)
    print (output.strip())
    
    t.close()
    u.close()
    v.close()

#!/bin/bash

# > deduplicate_fasta.sh <
#
# USAGE: deduplicate_fasta.sh <file1> <file2> ...
#
# Some FASTA files contain multiple sequences with identical annotations.
# To counter this, simply add a "_<incrementing_number>" as a suffix to
# all annotations!
#
# Use with care; it replaces the original files.

for fastafile in "$@"
do
    awk 'BEGIN {s=1; OFS=""} {if (/^>/) {print $0, "_", s; s+=1} else {print $0}}' $fastafile > tmp && mv -f tmp $fastafile
done
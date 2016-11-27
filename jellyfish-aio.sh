#!/bin/bash

# Jellyfish all-in-one script, simplifiying the entire Jellyfish pipeline:
# 1. jellyfish count
# 2. jellyfish merge
# 3. jellyfish histo
# 4. jellyfish stats
#
# Two arguments needed: first one is FASTA/FASTQ file; second one is k-mer size.
# (note that -C collapses any k-mer with its complementary, use this if the 
# library is non-directional)

echo Counting $1...
jellyfish count -m $2 -s 1G -t 15 -o ${1}.k${2} $1

if [ -f ${1}.k${2}_1 ];
then
    echo Merging variants of ${1}.k{2}...
    jellyfish merge ${1}.k${2}_* -o temp && rm -f ${1}.k${2}_* && mv temp ${1}.k${2}
else
    echo Skipping jellyfish merge, as ${1}.k${2} is not in multiple pieces.
fi

echo Creating histo and stats for ${1}.k${2}...
jellyfish stats ${1}.k${2} > ${1}.k${2}.stats
jellyfish histo -h 10000 -f -t 8 ${1}.k${2} > ${1}.k${2}.histo

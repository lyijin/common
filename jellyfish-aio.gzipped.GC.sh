#!/bin/bash

# Jellyfish all-in-one script, simplifiying the entire Jellyfish pipeline:
# 1. jellyfish count
# 2. jellyfish stats
# 3. jellyfish dump
#
# Two arguments needed: first one is k-mer size; second one is FASTA/FASTQ file.
# (note that -C collapses any k-mer with its complementary, use this if the 
# library is non-directional)
#
# NOTE: jellyfish merge not needed in Jellyfish v2 - merging is done whenever
# needed by jellyfish count.

set -e		# script exits upon a command producing nonzero exit value

echo Counting gzip-compressed $2...
zcat ${2} | jellyfish count -m $1 -s 1G -t 10 -o ${2}.k${1} /dev/fd/0

echo Creating dump and stats for ${2}.k${1}...
jellyfish stats ${2}.k${1} > ${2}.k${1}.stats.txt
jellyfish dump -c -t ${2}.k${1} > ${2}.k${1}.dump.tsv

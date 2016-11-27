#!/bin/bash

# only arguments '06', '07', '18' and '19' would work as ${1}

echo Interleaving for ${1}...
# pre-normalisation: interleave reads, as raw paired-end data is provided in _1 and _2.fastq files
#interleave-reads.py ../raw_data/pe${1}_1.paired.fastq ../raw_data/pe${1}_2.paired.fastq  > pe${1}_interleaved.fastq

echo Digital normalisation, k=20, C=20...
# step one: digital normalisation to k=20, C=20.
normalize-by-median.py -k 20 -C 20 -N 4 -x 4e9 -p -s pe${1}_norm_hash.kh pe${1}_interleaved.fastq > pe${1}.s1.out

echo Removing low-abundance k-mers...
# step two: remove low-abundance k-mers
filter-abund.py pe${1}_norm_hash.kh pe${1}_interleaved.fastq.keep > pe${1}.s2.out

echo Extracting paired reads...
# step 2.5: remove non-paired end reads
extract-paired-reads.py pe${1}_interleaved.fastq.keep.abundfilt > pe${1}.s25.out

echo Digital normalisation, k=20, C=5...
# step three: second round of normalisation to C=5
normalize-by-median.py -k 20 -C 5 -N 4 -x 4e9 -p pe${1}_interleaved.fastq.keep.abundfilt.pe > pe${1}.s3.out
#!/bin/bash

# cutadapt finally supports paired-end trimming!
#
# provide arguments as follows: 
#   cutadapt_16S_miseq.sh <8 bp i5 fwd tag> <8 bp i5 revcomp>
#                         <8 bp i7 rev tag> <8 bp i7 revcomp> <read1> <read2>
# 
# make sure read 1 contains "R1" and read 2 "R2".
#
# optional flags:
#   -q FRONT,BACK   quality (from front, from back)
#   -m LENGTH       retain reads of minimum length
#   --trim-n        removes flanking Ns (NNNATCGNNN --> ATCG)

outfile=`echo ${5} | sed 's/.fastq$//' | sed 's/.fq$//'`.trim.fastq
logfile=`echo ${5} | sed 's/.fastq$//' | sed 's/.fq$//'`.trim.log

# -q 20,20 changed to 5,5 to deal with low quality library
cutadapt -O 10 -q 5,5 -m 25 --trim-n  \
         -b AATGATACGGCGACCACCGAGATCTACAC${1}TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
         -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC${4}ATCTCGTATGCCGTCTTCTGCTTG \
         -B CTGTCTCTTATACACATCTGACGCTGCCGACGA${2}GTGTAGATCTCGGTGGTCGCCGTATCATT \
         -B CAAGCAGAAGACGGCATACGAGAT${3}GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
         ${5} ${6} -o ${outfile} -p ${outfile/R1/R2} > ${logfile}

#!/bin/bash

# cutadapt finally supports paired-end trimming!
#
# provide arguments as follows: 
#   cutadapt_16S_miseq.sh <8-char i5 fwd tag> <8 bp i7 rev tag> <read1> <read2>
# 
# make sure read 1 contains "R1" and read 2 "R2".
#
# optional flags:
#   -q FRONT,BACK   quality (from front, from back)
#   -m LENGTH       retain reads of minimum length
#   --trim-n        removes flanking Ns (NNNATCGNNN --> ATCG)

outfile=`echo ${3} | sed 's/.fastq$//' | sed 's/.fq$//'`.trim.fastq
logfile=`echo ${3} | sed 's/.fastq$//' | sed 's/.fq$//'`.trim.log

cutadapt -O 10 -q 20,20 -m 25 --trim-n -b AATGATACGGCGACCACCGAGATCTACAC${1}TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -b CAAGCAGAAGACGGCATACGAGAT${2}GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -B AATGATACGGCGACCACCGAGATCTACAC${1}TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -B CAAGCAGAAGACGGCATACGAGAT${2}GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG ${3} ${4} -o ${outfile} -p ${outfile/R1/R2} > ${logfile}

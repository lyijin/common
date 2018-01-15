#!/bin/bash

# provide arguments as follows: 
#   cutadapt_teloprime.sh <read1>
#
# optional flags:
#   -m LENGTH       retain reads of minimum length
#   --trim-n        removes flanking Ns (NNNATCGNNN --> ATCG)
#   -n              removes up to n adapters per read 

outfile=`echo ${1} | sed 's/.fa$//'`.trim.fa
logfile=`echo ${1} | sed 's/.fa$//'`.trim.log

cutadapt -O 15 -m 25 --trim-n -n 5 -g TGGATTGATATGTAATACGACTCACTATAG -a AAAAAAAAAAAAAAAAAACGCCTGAGA -b TCCGTAGCCATTTTGGCTCAAG ${1} -o ${outfile} > ${logfile}

#!/bin/bash

# > get_sra_from_ebi.sh <
#
# USAGE: get_sra_from_ebi.sh <acc_no1> <acc_no2> ...
#
# Relies on plain ol' `wget`. Gave up on getting `ascp` to work reliably.
#
# `wget` flags:
#   -c    Continue downloads
#   -nv   "no verbose", which is "quiet" + error messages
#   -r    Recursive
#   -l 1  One-level deep
#   -np   Do not ascend to parent directories
#   -nd   Do not create directories
for acc_no in "$@"
do
    if [ ${#acc_no} -gt 10 ]; then
        echo [`date`] Downloading ${acc_no}...
        wget -c -nv -r -l 1 -np -nd "ftp://ftp.ebi.ac.uk/ena/vol1/fastq/${acc_no:0:6}/0${acc_no: -2}/${acc_no}/"
    elif [ ${#acc_no} -eq 10 ]; then
        echo [`date`] Downloading ${acc_no}...
        wget -c -nv -r -l 1 -np -nd "ftp://ftp.ebi.ac.uk/ena/vol1/fastq/${acc_no:0:6}/00${acc_no: -1}/${acc_no}/"
    else
        echo [`date`] Downloading ${acc_no}...
        wget -c -nv -r -l 1 -np -nd "ftp://ftp.ebi.ac.uk/ena/vol1/fastq/${acc_no:0:6}/${acc_no}/"
    fi
    
    echo [`date`] Completed ${acc_no}.
    echo
done

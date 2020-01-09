#!/bin/bash

# > get_sra_from_ebi.sh <
#
# USAGE: get_sra_from_ebi.sh <acc_no1> <acc_no2> ...
#
# Depends on Aspera SCP. Script assumes sequencing files are from a paired-end
# run, then downloads <acc_no>_1.fastq.gz and <acc_no>_2.fastq.gz from EBI.
#
# ascp flags:
#   -pvT    Preserve timestamps, verbose, no encryption
#   -k 1    Resume downloads if the current and original attributes match
for acc_no in "$@"
do
    echo [`date`] Getting R1 for ${acc_no}...
    ~/.aspera/connect/bin/ascp -pvT -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -P 33001 -O 33001 -k 1 -l 50M --mode recv --overwrite older era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${acc_no:0:6}/00${acc_no: -1}/${acc_no}/${acc_no}_1.fastq.gz .
    
    echo [`date`] Getting R2 for ${acc_no}...
    ~/.aspera/connect/bin/ascp -pvT -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -P 33001 -O 33001 -k 1 -l 50M --mode recv --overwrite older era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${acc_no:0:6}/00${acc_no: -1}/${acc_no}/${acc_no}_2.fastq.gz .
    
    echo [`date`] Completed ${acc_no}.
    echo
done
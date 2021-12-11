#!/bin/bash

if [ "$#" == 0 ] || [ $1 == "-h" ]; then
  printf "Usage:\t tblout-to-fasta.sh hmmsearchResult RefSeq outputFasta\n\n"
  exit
fi

cut -d ' ' -f1 $1 > ids.txt
seqtk subseq $2 ids.txt > $3
rm ids.txt

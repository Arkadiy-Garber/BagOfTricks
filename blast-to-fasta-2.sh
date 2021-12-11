#!/bin/bash

if [ "$#" == 0 ] || [ $1 == "-h" ]; then
  printf "Usage:\t blast-to-fasta.sh blastResult QuerySeq outputFasta\n\n"
  exit
fi

cut -f1 $1 > ids.txt
seqtk subseq $2 ids.txt > $3
rm ids.txt

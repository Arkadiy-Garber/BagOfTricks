#!/bin/bash

if [ "$#" == 0 ] || [ $1 == "-h" ]; then
  printf "Usage:\t test.sh blastResult RefSeq outputFasta\n\n"
  exit
fi

echo “File name is ” $0
echo “First arg. is ” $1
echo “Second arg. is ” $2
echo “Second arg. is ” $3

cut -f2 $1 > ids.txt
seqtk subseq $2 ids.txt > $3
rm ids.txt

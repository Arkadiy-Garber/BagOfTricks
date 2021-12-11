#!/bin/bash


if [ "$#" == 0 ] || [ $1 == "-h" ]; then
  printf "Go from SAM files to sorted BAM files\n
Usage:\t samtools_looper.sh directory\n\n"
  exit
fi

cut -d ' ' -f1 $1 > ids.txt
seqtk subseq $2 ids.txt > $3
rm ids.txt


for i in $1/*sam; do samtools view -bS -o $1/${i%.*}.bam $1/$i; samtools sort --threads 16 -o $1/${i%.*}.bam.sorted $1/${i%.*}.bam; samtools index -b $1/${i%.*}.bam.sorted $1/${i%.*}.bam.sorted.bai; done

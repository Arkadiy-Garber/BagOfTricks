#!/bin/bash

if [ "$#" == 0 ] || [ $1 == "-h" ]; then
  printf "Usage:\t checkm-quality.sh fileExtension BinDirectory threads
Completion and redundancy results written to file: checkm_qaResults\n"
  exit
fi

checkm tree -x $1 $2 checkm-output/ -t $3
checkm tree_qa checkm-output/
checkm lineage_set checkm-output/ checkm-markers
checkm analyze -x $1 checkm-markers $2 checkm-output/ -t $3
checkm qa -o 1 -f checkm_qaResults checkm-markers checkm-output/ -t $3

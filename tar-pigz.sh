#!/bin/bash

if [ "$#" == 0 ] || [ $1 == "-h" ]; then
  printf "Usage:\t tar-pigz.sh directory\n\n"
  exit
fi


tar -cf $i.tar $i
pigz --best -p 16 $i.tar



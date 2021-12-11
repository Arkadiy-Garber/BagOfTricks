#!/bin/bash

if [ "$#" == 0 ] || [ $1 == "-h" ]; then
  printf "Usage:\t touch.sh dir\n"
  exit
fi

find $i -exec touch -c '{}' +

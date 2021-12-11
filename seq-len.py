#!/usr/bin/env python3
# !/bin/sh
# Author: Arkadiy Garber

from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys

parser = argparse.ArgumentParser(
    prog="seq-len.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    ************************************************************************
    ************************************************************************
    
    Script for printing the length of each sequence in a provided FASTA file
    Developed by Arkadiy Garber: agarber4@asu.edu
    ************************************************************************
    ************************************************************************
    '''))


parser.add_argument('-fasta', help='fasta file')

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]


def fasta(fasta_file):
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                seq = ''
            else:
                header = i[1:]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    return Dict


file = open(args.fasta, "r")
file = fasta(file)
ls = []
Dict = defaultdict(list)
for i in file.keys():
    # print(i + "\t" + str(len(file[i])))
    Dict[len(file[i])].append(i.split(" ")[0])


for i in sorted(Dict.keys(), reverse=True):
    for j in Dict[i]:
        print(str(j) + "\t" + str(i))






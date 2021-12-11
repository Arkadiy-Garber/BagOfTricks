#!/usr/bin/env python3
# !/bin/sh
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys


parser = argparse.ArgumentParser(
    prog="SeqLengthSep.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************
    *******************************************************
    *******************************************************
    
    Script for segregating a file into two by read length
    Developed by Arkadiy Garber: agarber4@asu.edu 

    *******************************************************
    *******************************************************
    *******************************************************
    '''))

parser.add_argument('-reads', type=str, help='File with sequences in FASTA format')

parser.add_argument('-long', type=str, help="output file for long reads")

parser.add_argument('-short', type=str, help="output file for short reads")

parser.add_argument('-len', type=str, help="max read length for short-reads file (default=600)", default=600)

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


file = open(args.reads, "r")
file = fasta(file)

long = open(args.long, "w")
short = open(args.short, "w")
for i in file.keys():
    if len(file[i]) > int(args.len):
        long.write(">" + i + "\n")
        long.write(file[i] + "\n")
    else:
        short.write(">" + i + "\n")
        short.write(file[i] + "\n")
long.close()
short.close()
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
    prog="GC-calc.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    
    Script for calculating GC content from FASTA files
    
    Developed by Arkadiy Garber: agarber4@asu.edu
    '''))

parser.add_argument('-g', help='FASTA file input (must be nucleotide genes or contigs)')

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]


def lastItem(string, delim):
    x = ''
    ls = string.split(delim)
    for i in ls:
        x = i
    return x


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
                header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    return Dict

genome = open(args.g, "r")
genome = fasta(genome)
GC = 0

total = 0
for i in genome.keys():
    seq = genome[i]
    total += len(seq)
    gc = 0
    for bp in seq:
        if bp == "C" or bp == "G":
            GC += 1
            gc += 1
    print(i + "\t" + str(round(float(GC/total), 3)))
print("------------------------------------")
print(lastItem(args.g, "/") + "\t" + str(round(float(GC/total), 3)))
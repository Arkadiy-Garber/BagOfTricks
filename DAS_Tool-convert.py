#!/usr/bin/env python3
# !/bin/sh
# Author: Arkadiy Garber

import textwrap
import argparse
import os
import re
from collections import defaultdict
import sys

parser = argparse.ArgumentParser(
    prog="DAS_Toool-convert.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    Program for creating TSV-formatted bin files for DAS_Tool;
    Developed by Arkadiy Garber; agarber4@asu.edu
    '''))
parser.add_argument('-dir', help='Directory with bins. Please provide full path.')
parser.add_argument('-name', help='Basename of output file.')
parser.add_argument('-id', help='suffix or prefix of bins')

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
    Dict[header]["seq"] = seq
    return Dict

cwd = os.getcwd()


outfile = open(cwd + "/" + args.name, "w")
dir = os.listdir(args.dir)
for i in dir:
    if args.id == i[len(i) - len(args.id): len(i)]:
        binX = open(args.directory + "/" + i, "r")
        binX = fasta(binX)
        for j in binX.keys():
            outfile.write(j + "\t" + i + "\n")
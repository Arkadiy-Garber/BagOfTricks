#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys


parser = argparse.ArgumentParser(
    prog="blast-to-fasta.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    ************************************************************************

    Script for generating a FASTA file from BLAST results

    Developed by Arkadiy Garber; agarber4@asu.edu
    ************************************************************************
    '''))

parser.add_argument('-b', type=str, help="blast result", default="NA")
parser.add_argument('-h', type=str, help="column to use as the header (must be a number; 1 means first column, 2 means second, and so forth)", default=13)
parser.add_argument('-s', type=str, help="column to use as the seq (default = 14)", default=14)
parser.add_argument('-o', type=str, help="output fasta file", default="NA")


if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]

redunDict = defaultdict(list)
file = open(args.b)
out = open(args.o, "w")
for i in file:
    ls = i.rstrip().split("\t")
    if ls[12] not in redunDict.keys():
        out.write(">" + ls[args.h-1] + "\n")
        out.write(ls[args.s-1] + "\n")
        redunDict[args.s-1].append(ls[1])
        redunDict[args.h-1].append(ls[1])
out.close()





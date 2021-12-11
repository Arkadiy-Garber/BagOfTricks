#!/usr/bin/env python3
# !/bin/sh
# Author: Arkadiy Garber
from collections import defaultdict
import statistics
import os
import re
import textwrap
import argparse
import sys


parser = argparse.ArgumentParser(
    prog="filter-seq-len.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    
    Script to filter out sequences below a certain length
    
    Developed by Arkadiy Garber: agarber4@asu.edu
    '''))

parser.add_argument('-f', help='input FASTA sequence file')
parser.add_argument('-m', help="minimum sequence length")
parser.add_argument('-o', help="name output file")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]

file = open(args.f, "r")
seq = ''
header = ''
if args.o != args.f:
    out = open(args.output, "w")
    for i in file:
        i = i.rstrip()
        if re.match(r'^>', i):
            if len(seq) >= int(args.m):
                print(header)
                out.write(">" + header + "\n")
                out.write(seq + "\n")
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
        else:
            seq += i

    file.close()
    out.close()
else:
    print("Output file name same as the fasta file. You want your assembly overwritten!?")

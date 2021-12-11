#!/usr/bin/env python3
# !/bin/sh
# Author: Arkadiy Garber


from collections import defaultdict
import re, os
import argparse
import textwrap
import statistics
import sys


parser = argparse.ArgumentParser(
    prog="genecall-mod.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    ************************************************************************
    ************************************************************************
    Script for adding a prefix to FASTA headers.
    Developed by Arkadiy Garber: agarber4@asu.edu
    ************************************************************************
    ************************************************************************
    '''))


parser.add_argument('-fasta', help='fasta file')

parser.add_argument('-prefix', help="what to add after carrot, but before the header")

parser.add_argument('-out', help="outfile")

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
                header = args.prefix + i[1:]
                seq = ''
            else:
                header = args.prefix + i[1:]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    return Dict


file = open(args.fasta, "r")
file = fasta(file)
out = open(args.out, "w")
for i in file.keys():
    out.write(">" + i + "\n")
    out.write(file[i] + "\n")
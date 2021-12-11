#!/usr/bin/env python3
from collections import defaultdict
import re, os
import textwrap
import argparse
import sys


parser = argparse.ArgumentParser(
    prog="binning-results.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    
    Script for summarizing into a tab-delimited file binning results
    
    Developed by Arkadiy Garber: agarber4@asu.edu
    
    '''))

parser.add_argument('-dir', help='directory of FASTA files')

parser.add_argument('-out', help="output table")

parser.add_argument('-x', help="extension for fasta file")

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


dir = os.listdir(args.dir)
out = open(args.out, "w")
for i in dir:
    if re.findall(args.x, i):
        file = open(args.dir + "/" + i, "r")
        file = fasta(file)
        for j in file.keys():
            out.write(j + "\t" + i + "\n")















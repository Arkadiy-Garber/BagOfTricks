#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys


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
    # print(count)
    return Dict


parser = argparse.ArgumentParser(
    prog="seq-pull.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************

    Script for pulling out sequences based on a character or string that should be present in the header
    Developed by Arkadiy Garber: agarber4@asu.edu

    *******************************************************
    '''))

parser.add_argument('-fasta', type=str, help="input FASTA file")
parser.add_argument('-id', type=str, help="string or character to use in regex for header identification")
parser.add_argument('-out', type=str, help="name of output FASTA file")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]

file = open(args.fasta)
file = fasta(file)
out = open(args.out, "w")
for i in file.keys():
    if re.findall(args.id, i):
        out.write(">" + i + '\n')
        out.write(file[i] + "\n")

out.close()
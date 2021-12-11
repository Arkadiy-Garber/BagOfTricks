#!/usr/bin/env python3
# !/bin/sh
# Author: Arkadiy Garber
from collections import defaultdict
import re
import os
import textwrap
import argparse
import urllib.request
import ssl
from urllib.error import HTTPError
import sys


parser = argparse.ArgumentParser(
    prog="de-align.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    Program for removing the dashes \'(-)\' and dots \'(.)\' generated from sequence aligners.
    
    Developed by Arkadiy Garber: agarber4@asu.edu
    '''))

parser.add_argument('-i', help='input aligned fasta file')

parser.add_argument('-o', help="output unaligned fasta file")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]


def remove(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


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
            seq += remove(i, ["-", "."])
    Dict[header] = seq
    return Dict

IN = open(args.i, "r")
IN = fasta(IN)
OUT = open(args.o, "w")
for i in IN.keys():
    OUT.write(">" + i + "\n")
    OUT.write(IN[i] + "\n")
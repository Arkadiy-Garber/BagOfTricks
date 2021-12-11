#!/usr/bin/env python3
# !/bin/sh
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
    return Dict


parser = argparse.ArgumentParser(
    prog="masker.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    
    Script for masking FASTA sequences - removes positions were there are too many gaps 
    (how many is too many is decided by the user)
    Developed by Arkadiy Garber: agarber4@asu.edu
    '''))

parser.add_argument('-i', help='alignment file to mask (in FASTA format)', default="NA")
parser.add_argument('-o', help="output file (also in FASTA format)", default="NA")
parser.add_argument('-m', help="maximum proportion of sequences that can contain a gap for the position to be kept in the alignment. Should be a number between 0 and 1 (default = 0.2)", default=0.2)

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]


infile = open(args.i)
infile = fasta(infile)


gapDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
firstHeader = (list(infile.keys())[0])
firstSeq = infile[firstHeader]
for i in range(0, len(firstSeq)):
    count = 0
    for j in infile.keys():
        char = (infile[j][i])
        if char == "-":
            count += 1
    gapDict[i] = count/len(infile.keys())

count = 0
out = open(args.o, "w")
infile = open(args.i)
infile = fasta(infile)
for i in infile.keys():
    out.write(">" + i + "\n")
    for j in gapDict.keys():
        if gapDict[j] < float(args.m):
            out.write(infile[i][j])
        else:
            count += 1
    out.write("\n")

print(str(int(count/len(infile.keys()))) + " positions were removed")
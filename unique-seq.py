#!/usr/bin/env python3
from collections import defaultdict
import re, os
import textwrap
import argparse
import sys


def fasta(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


parser = argparse.ArgumentParser(
    prog="unique-seq.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    
    Script for removing redundant/identical entries (based on header name) from FASTA file
    Developed by Arkadiy Garber: agarber4@asu.edu

    '''))

parser.add_argument('-i', help='input fasta file')

parser.add_argument('-o', help="output fasta file")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]


if args.i == args.o:
    print("Input and output file names identical. Please choose a different name for output file")
    raise SystemExit


file = open(args.i)
file = fasta(file)
out = open(args.o, "w")
for i in file.keys():
    out.write(">" + i + "\n")
    out.write(file[i] + "\n")
out.close()





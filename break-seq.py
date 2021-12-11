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


def allButTheFirst(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(1, length):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)]


parser = argparse.ArgumentParser(
    prog="break-seq.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    
    Script for splitting contig sequences into 2000 bp chunks
    
    Developed by Arkadiy Garber: agarber4@asu.edu
    
    '''))

parser.add_argument('-i', help='input fasta file')

parser.add_argument('-o', help="output fasta file")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]


file = open(args.i)
file = fasta(file)
out = open(args.o, "w")

for i in file.keys():
    print(i)
    seq = file[i]
    contig = i
    firstItem = contig.split(" ")[0]
    rest = allButTheFirst(contig, " ")
    count = 0
    counter = 0
    for j in range(2000, len(seq), 2000):
        start = counter
        end = j
        segment = seq[start:end]
        out.write(">" + firstItem + "_" + str(count) + " " + rest + "\n")
        out.write(segment + "\n")
        counter += 2000
        count += 1















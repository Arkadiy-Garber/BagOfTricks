#!/usr/bin/env python3
# !/bin/sh
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys

parser = argparse.ArgumentParser(
    prog="Seqs.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************
    
    Script for printing out a summary of provided FASTA file
    Developed by Arkadiy Garber: agarber4@asu.edu

    *******************************************************
    '''))

parser.add_argument('-fasta', type=str, help='FASTA file')

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

file = open(args.fasta, "r")
file = fasta(file)
count = 0
counter = 0
seq300 = 0
seq1000 = 0
seq5000 = 0
seq10000 = 0
seq25000 = 0
seq50000 = 0
seq75000 = 0
seq100000 = 0
for i in file.keys():
    count += 1
    counter += len(file[i])
    if len(file[i]) >= 300:
        seq300 += 1
    if len(file[i]) >= 1000:
        seq1000 += 1
    if len(file[i]) >= 5000:
        seq5000 += 1
    if len(file[i]) >= 10000:
        seq10000 += 1
    if len(file[i]) >= 25000:
        seq25000 += 1
    if len(file[i]) >= 50000:
        seq50000 += 1
    if len(file[i]) >= 75000:
        seq75000 += 1
    if len(file[i]) >= 100000:
        seq100000 += 1
average = counter/count
numSeqs = count

print("Total length in bp: " + str(counter))
print("Number of Sequences: " + str(numSeqs))
print("Average Sequence Length: " + str(average))
print("Sequences over 300 bp: " + str(seq300) + " (" + str(round((seq300/numSeqs)*100, 3)) + "%)")
print("Sequences over 1000 bp: " + str(seq1000) + " (" + str(round((seq1000/numSeqs)*100, 3)) + "%)")
print("Sequences over 5000 bp: " + str(seq5000) + " (" + str(round((seq5000/numSeqs)*100, 3)) + "%)")
print("Sequences over 10000 bp: " + str(seq10000) + " (" + str(round((seq10000/numSeqs)*100, 3)) + "%)")
print("Sequences over 25000 bp: " + str(seq25000) + " (" + str(round((seq25000/numSeqs)*100, 3)) + "%)")
print("Sequences over 50000 bp: " + str(seq50000) + " (" + str(round((seq50000/numSeqs)*100, 3)) + "%)")
print("Sequences over 75000 bp: " + str(seq75000) + " (" + str(round((seq75000/numSeqs)*100, 3)) + "%)")
print("Sequences over 100000 bp: " + str(seq100000) + " (" + str(round((seq100000/numSeqs)*100, 3)) + "%)")

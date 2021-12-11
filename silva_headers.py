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
    prog="silva-headers.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    ************************************************************************
    ************************************************************************
    Developed by Arkadiy Garber; University of Delaware, Geological Sciences
    Please send comments and inquiries to arkg@udel.edu
    ************************************************************************
    ************************************************************************
    '''))


parser.add_argument('-silva_db', help='location of SILVA database file (in fasta format).')


parser.add_argument('-output', help="Name of output file (i.e. the file containing sampled sequences from SILVA db",
                    default="silva_sample.fa")

parser.add_argument('-SILVA_depth', help="SAMtools-generated depth file for silva database")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]


def sum(list):
    sum = 0
    for i in list:
        i = int(i)
        sum += i
    return sum


print("\nstarting analysis...")
print("reading in SILVA database file")
print("...")
dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
file = open(args.silva_db, "r")
for i in file:
    i = i.rstrip()
    if re.search(r'>', i):
        ls = i.split(" ")
        id = ls[0]
        id = id[1:]
        dict[id] = i

print("analyzing depths file")
print("...")
dictA = defaultdict(list)
dictB = defaultdict(list)
out = open(args.output, "w")
depth = open(args.SILVA_depth, "r")
for i in depth:
    i = i.rstrip()
    ls = i.split("\t")
    if int(ls[2]) == 0 and int(ls[3]) == 0:
        pass
    else:
        id = ls[0]
        header = (dict[id][1:])
        headerLS = header.split(" ")
        try:
            Class = headerLS[1].split(";")[2]
        except IndexError:
            Class = headerLS[1].split(";")[1]
        # print(id)
        # print(header)
        # print("")
        if int(ls[1]) > 250 and int(ls[1]) < 1250:
            dictA[Class].append(int(ls[2]))
            dictB[Class].append(int(ls[3]))
        # out.write(header + "\t" + str(ls[2]) + "\t" + str(ls[3]) + "\n")

print("computing coverage...")
for i in dictA.keys():
    depthsA = dictA[i]
    depthsB = dictB[i]
    meanA = statistics.mean(depthsA)
    meanB = statistics.mean(depthsB)
    # stdev = statistics.stdev(depths)
    print(i)
    print("meanA: " + str(meanA))
    print("meanB: " + str(meanB))
    print("")
    out.write(i + "\t" + str(meanA) + "\t" + str(meanB) + "\n")

print("\nAll done.")
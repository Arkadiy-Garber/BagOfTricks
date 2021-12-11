#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys


parser = argparse.ArgumentParser(
    prog="count-to-rpkm.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    ************************************************************************

    Script for generating RPKM values from an htseq-generated counts table

    Developed by Arkadiy Garber; agarber4@asu.edu
    ************************************************************************
    '''))


parser.add_argument('-gff', type=str, help="GFF file used to generate the counts table", default="NA")
parser.add_argument('-id', type=str, help="attribute used from the GFF file (default = ID)", default="ID")
parser.add_argument('-counts', type=str, help="counts table output from htseq-count", default="NA")

parser.add_argument('-delim', type=str, help="delimiter for counts table (tab or comma; default = tab)", default="\t")
parser.add_argument('-out', type=str, help="name output file (default is the name of your counts file appended with \'.rpkm\')", default="NA")

parser.add_argument('-readsMap', type=int, help="Tab-delimited two column file. First column is name of htseq output counts file; "
                                                "second column has the corresponding read totals for each sample", default="NA")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]


if args.gff == "NA":
    print("GFF file not provided")
    raise SystemExit

if args.counts == "NA":
    print("Counts file not provided")
    raise SystemExit

if args.reads == "NA":
    print("total number of reads not provided")
    raise SystemExit

counts = open(args.counts)
gff = open(args.gff)
if args.out == "NA":
    outfilename = args.counts
    outfilename = outfilename + ".rpkm"
    out = open(outfilename, "w")
else:
    out = open(args.out, "w")


gffDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
for i in gff:
    if not re.match(r'#', i):
        ls = i.rstrip().split("\t")
        start = float(ls[3])
        end = float(ls[4])
        length = end - start
        attributes = ls[8].split(";")
        for j in attributes:
            attribute = args.id + "="
            if re.findall(attribute, j):
                attributeList = j.split("=")
                ID = attributeList[1]
                gffDict[ID] = length
                break

readDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
readmap = open(args.readsMap)
for i in readmap:
    ls = i.rstrip().split("\t")
    readDict[ls[0]] = ls[1]


DELIM = ''
if args.delim == "comma":
    DELIM = ","
elif args.delim == "tab":
    DELIM = "\t"
else:
    DELIM = "\t"


totalRPK = 0
T = 0
RPM = int(readDict[args.counts])/1000000
for i in counts:
    ls = i.rstrip().split("\t")
    if not re.match(r'#', i) and not re.match(r'__', i) and len(ls) > 1:
        ID = ls[0]
        count = float(ls[1])
        length = gffDict[ID]

        rpm = count / RPM
        rpkm = rpm / (length/1000)
        out.write(ID + DELIM + str(rpkm) + "\n")

out.close()

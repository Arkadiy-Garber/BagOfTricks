#!/usr/bin/env python3
#!/bin/sh
from collections import defaultdict
import re, os
import textwrap
import argparse
import sys


parser = argparse.ArgumentParser(
    prog="reformat-fastas",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    Program for reformatting multiple fasta files alla Anvi'o
    '''))

parser.add_argument('-i', help='input directory with fasta files')
parser.add_argument('-x', help='extension of fasta files (include the preceeding period)')
parser.add_argument('-min_length', help="minimum sequence length for fasta file (default = 1)", default=1)

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]


def stabilityCounter(int):
    if len(str(int)) == 1:
        string = (str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 2:
        string = (str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 3:
        string = (str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 4:
        string = (str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 5:
        string = (str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 6:
        string = (str(0) + str(0) + str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 7:
        string = (str(0) + str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 8:
        string = (str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 9:
        string = (str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 10:
        string = (str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 11:
        string = (str(0) + str(int))
        return (string)


dir = os.listdir(args.i)


def fasta(dir):
    count = 0
    for file in dir:
        if re.findall(args.x, file):
            ext = len(args.x)
            infile = open(args.i + "/" + file, "r")
            outfile = open(args.i + "/" + file[0:len(file)-ext] + "-fixed" + args.x, "w")
            seq = ''
            header = ''
            for i in infile:
                i = i.rstrip()
                if re.match(r'^>', i):
                    if len(seq) >= int(args.min_length):
                        outfile.write(">" + header + "\n")
                        outfile.write(seq + "\n")
                        count += 1
                        headerNum = stabilityCounter(count)
                        header = ("c_" + str(headerNum))
                        seq = ''
                    else:
                        count += 1
                        headerNum = stabilityCounter(count)
                        header = ("c_" + str(headerNum))
                        seq = ''
                else:
                    seq += i
            outfile.write(">" + header + "\n")
            outfile.write(seq + "\n")
            infile.close()
            outfile.close()

fasta(dir)

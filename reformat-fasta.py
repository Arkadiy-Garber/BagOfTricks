from collections import defaultdict
import re, os
import textwrap
import argparse
import sys


parser = argparse.ArgumentParser(
    prog="reformat-fasta",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    Program for reformatting fasta file alla Anvi'o
    '''))

parser.add_argument('-i', help='input fasta file')
parser.add_argument('-o', help="output re-formatted fasta file")
parser.add_argument('-min_length', help="minimum sequence length for fasta file")

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


def fasta(fasta_file):
    outfile = open(args.o, "w")
    seq = ''
    header = ''
    count = 0
    for i in fasta_file:
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

file = open(args.i, "r")
fasta(file)
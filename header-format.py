from collections import defaultdict
import re, os
import textwrap
import argparse
import sys


parser = argparse.ArgumentParser(
    prog="fasta-reformat.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    
    Script for removing everything after the first space in a FASTA file
    Developed by Arkadiy Garber: agarber4@asu.edu
    
    '''))

parser.add_argument('-file', help='input fasta file')

parser.add_argument('-out', help="output re-formatted fasta file")

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
                header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    return Dict


file = open(args.file, "r")
file = fasta(file)
out = open(args.out, "w")
for i in file.keys():
    out.write(">" + i + "\n")
    out.write(file[i] + "\n")
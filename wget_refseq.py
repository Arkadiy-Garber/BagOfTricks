#!/usr/bin/env python3
# !/bin/sh
from collections import defaultdict
import re
import os
import textwrap
import argparse
import ssl
import sys

gcontext = ssl.SSLContext(ssl.PROTOCOL_TLSv1)


parser = argparse.ArgumentParser(
    prog="FindMeHemes.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************
    *******************************************************
    *******************************************************

    Developed by Arkadiy Garber^1; 
    ^1University of Delaware, Geological Sciences

    Please send comments and inquiries to arkg@udel.edu

    *******************************************************
    *******************************************************
    *******************************************************
    '''))

parser.add_argument('-f', type=str, help='input file', default="NA")


if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]


def replace(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
    outString = "".join(emptyList)
    return outString


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
            seq += i
            seq += "\n"
    Dict[header] = seq
    return Dict


def lastItem(ls):
    x = ''
    for i in ls:
        if i != "":
            x = i
    return x

total = 0
file = open(args.f)
for i in file:
    if re.findall(r'genomic.fna', i):
        ls = i.rstrip().split(" ")
        size = lastItem(ls)
        size = remove(size, ["M", "K"])
        size = float(size)
        total += size
        sys.stdout.write(": %d%% Mb   \r" % (total))
        sys.stdout.flush()
        os.system("wget https://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/%s" % ls[0])











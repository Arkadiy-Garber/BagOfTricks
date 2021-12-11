#!/usr/bin/env python3
# !/bin/sh
from collections import defaultdict
import re
import os
import textwrap
import argparse
import urllib.request
import ssl

gcontext = ssl.SSLContext(ssl.PROTOCOL_TLSv1)


parser = argparse.ArgumentParser(
    prog="Seqs.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************

    Developed by Arkadiy Garber^1 and Gustavo Ramiréz^2; 
    1^University of Delaware, Geological Sciences
    2^University of Rhode Island, Graduate School of Oceanography

    Please send comments and inquiries to arkg@udel.edu

    *******************************************************
    '''))

parser.add_argument('-fasta', type=str, help='FASTA file')

args = parser.parse_args()

total = 0
file = open(args.fasta, "r")
for i in file:
    if not re.match(r'>', i):
        total += len(i)
print(total)
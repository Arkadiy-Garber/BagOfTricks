#!/usr/bin/env python3
# !/bin/sh
from collections import defaultdict
import re
import os
import textwrap
import argparse
import urllib.request
import ssl
import sys

gcontext = ssl.SSLContext(ssl.PROTOCOL_TLSv1)


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


def mass(seq):
    dal = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    dal["A"] = 89
    dal["R"] = 174
    dal["N"] = 132
    dal["D"] = 133
    dal["C"] = 121
    dal["Q"] = 146
    dal["E"] = 147
    dal["G"] = 75
    dal["H"] = 155
    dal["I"] = 131
    dal["L"] = 131
    dal["K"] = 146
    dal["M"] = 149
    dal["F"] = 165
    dal["P"] = 115
    dal["S"] = 105
    dal["T"] = 119
    dal["W"] = 204
    dal["Y"] = 181
    dal["V"] = 117
    dal["*"] = 0
    mw = 0
    for i in seq:
        mw += int(dal[i])
    return mw/1000


parser = argparse.ArgumentParser(
    prog="daltons.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    ************************************************************************

    Script for predicting the molecular mass in kilo daltons (kDa) of proteins

    Developed by Arkadiy Garber; agarber4@asu.edu
    ************************************************************************
    '''))


parser.add_argument('-f', type=str, help="input protein seqs in FASTA format", default="NA")


if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]

prots = open(args.f)
prots = fasta(prots)
for i in prots.keys():
    print(i.split(" ")[0] + '\t' + str(mass(prots[i])) + " kDa")



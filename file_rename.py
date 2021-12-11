#!/usr/bin/env python3
# !/bin/sh
# Author: Arkadiy Garber


from collections import defaultdict
import re, os
import argparse
import textwrap
import statistics


def replace(stringOrlist, ilItem, item):
    emptyList = []
    for i in stringOrlist:
        if i != ilItem:
            emptyList.append(i)
        else:
            emptyList.append(item)
    outString = "".join(emptyList)
    return outString


parser = argparse.ArgumentParser(
    prog="file_rename.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    ************************************************************************
    ************************************************************************
    Developed by Arkadiy Garber; University of Delaware, Geological Sciences
    Please send comments and inquiries to arkg@udel.edu
    ************************************************************************
    ************************************************************************
    '''))


parser.add_argument('-folder', help='folder with fasta files')


parser.add_argument('-illegal_char', help="illiegal character to replace")

parser.add_argument('-replace', help="character to replace illegal character")

args = parser.parse_args()

dir = os.listdir(args.folder)
for i in dir:
    os.system("mv " + args.folder + "/" + i + " " + args.folder + "/" + replace(i, args.illegal_char, args.replace))
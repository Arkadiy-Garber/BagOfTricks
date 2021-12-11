#!/usr/bin/env python3
from collections import defaultdict
import time
import re
import os
import sys
import textwrap
import argparse
import urllib.request
import ssl
from urllib.error import HTTPError
gcontext = ssl.SSLContext(ssl.PROTOCOL_TLSv1)


parser = argparse.ArgumentParser(
    prog="silva-survey.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    Developed by Arkadiy Garber; University of Delaware, Geological Sciences
    Please send comments and inquiries to arkg@udel.edu
    '''))

parser.add_argument('-taxa', help='taxanomic affiliation of organism that you would like to investigate')
parser.add_argument('-silva_DB', help="location of SILVA database")
parser.add_argument('-out_dir', help="directory to which the output file will be written")

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
    Dict[header] = seq
    return Dict


silva = open(args.silva_DB, "r")
silva =fasta(silva)
outfile = open(args.out_dir + "/" + args.taxa + "_SILVA.csv", "w")
outfile.write(args.taxa + "_strain" + "," + "Source" + "," + "Reference" + "," + "16S_seq" + "\n")

for i in silva.keys():
    if re.findall(args.taxa, i):
        print(i)
        time.sleep(2)
        outfile.write(i + ",")
        id = i.split(" ")[0]
        id = id.split(".")[0]
        fp = urllib.request.urlopen("https://www.ebi.ac.uk/ena/data/view/%s&display=text" % id, context=gcontext)
        mybytes = fp.read()
        mystr = mybytes.decode("utf8")
        fp.close()
        lines = mystr.split('\n')
        string = ''
        count = 0
        ref = ''
        for j in lines:
            if re.findall(r'isolation_source', j) or re.findall(r'tissue_type', j) or \
                    re.findall(r'Source of Isolate', j):
                count += 1
                try:
                    string = j.split("\"")[1]
                    string = replace(string, [","], ";")
                except IndexError:
                    string = j
                    string = replace(j, [","], ";")
                count += 1
            if re.match(r'^RT', j) and len(j) > 6:
                ref += j
        if count == 0:
            outfile.write("No_data" + ",")
        else:
            outfile.write(replace(string, [","], ";") + ",")
        ref = (replace(ref, [","], ";"))
        ref = remove(ref, ["RT"])
        outfile.write(ref + ",")
        outfile.write(silva[i] + "\n")

#!/usr/bin/env python3
# !/bin/sh
# Author: Arkadiy Garber


from collections import defaultdict
import re, os
import argparse
import textwrap
import sys


parser = argparse.ArgumentParser(
    prog="silva-sampler.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    Program for sampling sequences from the SILVA database, based on taxanomical resolution;
    Developed by Arkadiy Garber; University of Delaware, Geological Sciences
    Please send comments and inquiries to arkg@udel.edu
    '''))


parser.add_argument('-silva_db', help='location of SILVA database file (in fasta format).')


parser.add_argument('-level', help="Options: phylum (default), class, family, order, genus",
                    default="phylum")

parser.add_argument('-output', help="Name of output file (i.e. the file containing sampled sequences from SILVA db",
                    default="silva_sample.fa")

parser.add_argument('-RNA_to_DNA', help="convert uracil nucleotides to thymine? (y/n). Default = n",
                    default="n")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]

def main1():
    def fastaPhylum(fasta_file):
        seq = ''
        header = ''
        Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in fasta_file:
            i = i.rstrip()
            if re.match(r'^>', i):
                if len(seq) > 0:
                    Dict[head]["seq"] = seq
                    Dict[head]["header"] = header
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    head = headerList[1]
                    seq = ''
                else:
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    head = headerList[1]
                    seq = ''
            else:
                seq += i
        Dict[head]["seq"] = seq
        Dict[head]["header"] = header
        outfile = open(args.output, "w")
        for j in Dict.keys():
            outfile.write(">" + Dict[j]["header"] + "\n")
            outfile.write(Dict[j]["seq"] + "\n")


    def fastaClass(fasta_file):
        seq = ''
        header = ''
        Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in fasta_file:
            i = i.rstrip()
            if re.match(r'^>', i):
                if len(seq) > 0:
                    Dict[head]["seq"] = seq
                    Dict[head]["header"] = header
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    head = headerList[2]
                    seq = ''
                else:
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    head = headerList[2]
                    seq = ''
            else:
                seq += i
        Dict[head]["seq"] = seq
        Dict[head]["header"] = header
        outfile = open(args.output, "w")
        for j in Dict.keys():
            outfile.write(">" + Dict[j]["header"] + "\n")
            outfile.write(Dict[j]["seq"] + "\n")


    def fastaOrder(fasta_file):
        seq = ''
        header = ''
        Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in fasta_file:
            i = i.rstrip()
            if re.match(r'^>', i):
                if len(seq) > 0:
                    Dict[head]["seq"] = seq
                    Dict[head]["header"] = header
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    head = headerList[3]
                    seq = ''
                else:
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    head = headerList[3]
                    seq = ''
            else:
                seq += i
        Dict[head]["seq"] = seq
        Dict[head]["header"] = header
        outfile = open(args.output, "w")
        for j in Dict.keys():
            outfile.write(">" + Dict[j]["header"] + "\n")
            outfile.write(Dict[j]["seq"] + "\n")


    def fastaFamily(fasta_file):
        seq = ''
        header = ''
        Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in fasta_file:
            i = i.rstrip()
            if re.match(r'^>', i):
                if len(seq) > 0:
                    Dict[head]["seq"] = seq
                    Dict[head]["header"] = header
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    head = headerList[4]
                    seq = ''
                else:
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    head = headerList[4]
                    seq = ''
            else:
                seq += i
        Dict[head]["seq"] = seq
        Dict[head]["header"] = header
        outfile = open(args.output, "w")
        for j in Dict.keys():
            outfile.write(">" + Dict[j]["header"] + "\n")
            outfile.write(Dict[j]["seq"] + "\n")


    def fastaGenus(fasta_file):
        seq = ''
        header = ''
        Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in fasta_file:
            i = i.rstrip()
            if re.match(r'^>', i):
                if len(seq) > 0:
                    Dict[head]["seq"] = seq
                    Dict[head]["header"] = header
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    head = headerList[5]
                    seq = ''
                else:
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    head = headerList[5]
                    seq = ''
            else:
                seq += i
        Dict[head]["seq"] = seq
        Dict[head]["header"] = header
        outfile = open(args.output, "w")
        for j in Dict.keys():
            outfile.write(">" + Dict[j]["header"] + "\n")
            outfile.write(Dict[j]["seq"] + "\n")


    silva = open(args.silva_db, "r")
    if args.level == "phylum":
        fastaPhylum(silva)
    if args.level == "class":
        fastaClass(silva)
    if args.level == "order":
        fastaOrder(silva)
    if args.level == "family":
        fastaFamily(silva)
    if args.level == "genus":
        fastaGenus(silva)

    print("Finished sampling from SILVA database.")
    print("Thanks for using the silva-sampler program!")


# *************************************************************************************************
# *************************************************************************************************
# *************************************************************************************************
# *************************************************************************************************
# *************************************************************************************************


def main2():
    def replace(stringOrlist, list, item):
        emptyList = []
        for i in stringOrlist:
            if i not in list:
                emptyList.append(i)
            else:
                emptyList.append(item)
        outString = "".join(emptyList)
        return outString


    def fastaPhylum(fasta_file):
        seq = ''
        header = ''
        Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in fasta_file:
            i = i.rstrip()
            if re.match(r'^>', i):
                if len(seq) > 0:
                    Dict[head]["seq"] = replace(seq, ["U"], "T")
                    Dict[head]["header"] = header
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    head = headerList[1]
                    seq = ''
                else:
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    head = headerList[1]
                    seq = ''
            else:
                seq += i
        Dict[head]["seq"] = seq
        Dict[head]["header"] = header
        outfile = open(args.output, "w")
        for j in Dict.keys():
            outfile.write(">" + Dict[j]["header"] + "\n")
            outfile.write(Dict[j]["seq"] + "\n")

    def fastaClass(fasta_file):
        seq = ''
        header = ''
        Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in fasta_file:
            i = i.rstrip()
            if re.match(r'^>', i):
                if len(seq) > 0:
                    Dict[head]["seq"] = replace(seq, ["U"], "T")
                    Dict[head]["header"] = header
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    try:
                        head = headerList[2]
                    except IndexError:
                        pass
                    seq = ''
                else:
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    try:
                        head = headerList[2]
                    except IndexError:
                        pass
                    seq = ''
            else:
                seq += i
        Dict[head]["seq"] = seq
        Dict[head]["header"] = header
        outfile = open(args.output, "w")
        for j in Dict.keys():
            outfile.write(">" + Dict[j]["header"] + "\n")
            outfile.write(Dict[j]["seq"] + "\n")

    def fastaOrder(fasta_file):
        seq = ''
        header = ''
        Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in fasta_file:
            i = i.rstrip()
            if re.match(r'^>', i):
                if len(seq) > 0:
                    Dict[head]["seq"] = replace(seq, ["U"], "T")
                    Dict[head]["header"] = header
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    try:
                        head = headerList[3]
                    except IndexError:
                        try:
                            head = headerList[2]
                        except IndexError:
                            try:
                                head = headerList[1]
                            except IndexError:
                                pass
                    seq = ''
                else:
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    try:
                        head = headerList[3]
                    except IndexError:
                        try:
                            head = headerList[2]
                        except IndexError:
                            try:
                                head = headerList[1]
                            except IndexError:
                                pass
                    seq = ''
            else:
                seq += i
        Dict[head]["seq"] = seq
        Dict[head]["header"] = header
        outfile = open(args.output, "w")
        for j in Dict.keys():
            outfile.write(">" + Dict[j]["header"] + "\n")
            outfile.write(Dict[j]["seq"] + "\n")

    def fastaFamily(fasta_file):
        seq = ''
        header = ''
        Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in fasta_file:
            i = i.rstrip()
            if re.match(r'^>', i):
                if len(seq) > 0:
                    Dict[head]["seq"] = replace(seq, ["U"], "T")
                    Dict[head]["header"] = header
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    try:
                        head = headerList[4]
                    except IndexError:
                        try:
                            head = headerList[3]
                        except IndexError:
                            try:
                                head = headerList[2]
                            except IndexError:
                                try:
                                    head = headerList[1]
                                except IndexError:
                                    pass
                    seq = ''
                else:
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    try:
                        head = headerList[4]
                    except IndexError:
                        try:
                            head = headerList[3]
                        except IndexError:
                            try:
                                head = headerList[2]
                            except IndexError:
                                try:
                                    head = headerList[1]
                                except IndexError:
                                    pass
                    seq = ''
            else:
                seq += i
        Dict[head]["seq"] = seq
        Dict[head]["header"] = header
        outfile = open(args.output, "w")
        for j in Dict.keys():
            outfile.write(">" + Dict[j]["header"] + "\n")
            outfile.write(Dict[j]["seq"] + "\n")

    def fastaGenus(fasta_file):
        seq = ''
        header = ''
        Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for i in fasta_file:
            i = i.rstrip()
            if re.match(r'^>', i):
                if len(seq) > 0:
                    Dict[head]["seq"] = replace(seq, ["U"], "T")
                    Dict[head]["header"] = header
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    try:
                        head = headerList[5]
                    except IndexError:
                        try:
                            head = headerList[4]
                        except IndexError:
                            try:
                                head = headerList[3]
                            except IndexError:
                                try:
                                    head = headerList[2]
                                except IndexError:
                                    try:
                                        head = headerList[1]
                                    except IndexError:
                                        pass
                    seq = ''
                else:
                    header = i[1:]
                    headerLS = header.split(" ")
                    headerList = headerLS[1].split(";")
                    try:
                        head = headerList[5]
                    except IndexError:
                        try:
                            head = headerList[4]
                        except IndexError:
                            try:
                                head = headerList[3]
                            except IndexError:
                                try:
                                    head = headerList[2]
                                except IndexError:
                                    try:
                                        head = headerList[1]
                                    except IndexError:
                                        pass
                    seq = ''
            else:
                seq += i
        Dict[head]["seq"] = seq
        Dict[head]["header"] = header
        outfile = open(args.output, "w")
        for j in Dict.keys():
            outfile.write(">" + Dict[j]["header"] + "\n")
            outfile.write(Dict[j]["seq"] + "\n")

    silva = open(args.silva_db, "r")
    if args.level == "phylum":
        fastaPhylum(silva)
    if args.level == "class":
        fastaClass(silva)
    if args.level == "order":
        fastaOrder(silva)
    if args.level == "family":
        fastaFamily(silva)
    if args.level == "genus":
        fastaGenus(silva)

    print("Finished sampling from SILVA database.")
    print("Thanks for using the silva-sampler program!")


if args.RNA_to_DNA == "n":
    main1()
if args.RNA_to_DNA == "y":
    main2()
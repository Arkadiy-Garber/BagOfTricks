#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys


parser = argparse.ArgumentParser(
    prog="mmseqs-helper.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************
    '''))

parser.add_argument('-gff', type=str, help="GFF file used to generate the counts table", default="NA")
parser.add_argument('-id', type=str, help="attribute used from the GFF file (default = ID)", default="ID")
parser.add_argument('-counts', type=str, help="counts table output from htseq-count", default="NA")
parser.add_argument('-reverse', type=str, help="counts table output from htseq-count for reverse "
                                               "strand-mapping reads (optional)", default="NA")
parser.add_argument('-delim', type=str, help="delimiter for counts table (default = TAB)", default="\t")
parser.add_argument('-len', type=str, help="average read length", default="NA")
parser.add_argument('-out', type=str, help="name output file (default is the name of your counts file appended with \'.tpm\')", default="NA")
parser.add_argument('-out2', type=str, help="name output file for tpm calculated from reverse strand-mapping counts table "
                                            "(default is the name of your reverse count file appended with \'.tpm\')", default="NA")

parser.add_argument('--wagner', type=str, help="include this flag to use the equation from Wagner et al., 2012 paper. Must include read length via the -len argument with thos flag", const=True, nargs="?")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_args()

if args.gff == "NA":
    print("GFF file not provided")
    raise SystemExit

if args.counts == "NA":
    print("Counts file not provided")
    raise SystemExit

if args.wagner:
    if args.len == "NA":
        print("provided the --wagner flag, but no overage read length. Using 150 bp")
        reads = 150
    else:
        reads = float(args.len)

counts = open(args.counts)
gff = open(args.gff)
if args.out == "NA":
    outfilename = args.counts
    outfilename = outfilename + ".tpm"
    out = open(outfilename, "w")
else:
    out = open(args.out, "w")


if args.wagner:
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

    totalRPK = 0
    T = 0
    for i in counts:
        ls = i.rstrip().split("\t")
        if not re.match(r'#', i) and not re.match(r'__', i):
            ID = ls[0]
            count = float(ls[1])
            length = gffDict[ID]

            # WAGNER ET AL 2012
            t = (count * float(reads)) / length
            T += t

            # ORIGINAL
            # pk = length / 1000
            # rpk = count / pk
            # totalRPK += rpk


            ############################################################################
    # SECOND COUNTS TABLE
    if args.reverse != "NA":
        reverse = open(args.reverse)

        for i in reverse:
            ls = i.rstrip().split("\t")
            if not re.match(r'#', i) and not re.match(r'__', i):
                ID = ls[0]
                count = float(ls[1])
                length = gffDict[ID]

                # WAGNER ET AL 2012
                t = (count * float(reads)) / length
                T += t

                # ORIGINAL
                # pk = length / 1000
                # rpk = count / pk
                # totalRPK += rpk
    ############################################################################

    # WAGNER ET AL 2012
    # print(T)

    # ORIGINAL
    # denom = totalRPK / 1000000

    totalTPM = 0
    counts = open(args.counts)
    for i in counts:
        ls = i.rstrip().split("\t")
        if not re.match(r'#', i) and not re.match(r'__', i):
            ID = ls[0]
            count = float(ls[1])
            length = gffDict[ID]

            # WAGNER ET AL 2012
            tpm = (float(count) * float(reads) * 1000000) / (length * T)

            # ORIGINAL
            # pk = length / 1000
            # rpk = count / pk
            # tpm = rpk / denom

            out.write(ID + "\t" + str(tpm) + "\n")
            totalTPM += tpm

    ############################################################################
    if args.reverse != "NA":
        if args.out2 == "NA":
            outfilename = args.reverse
            outfilename = outfilename + ".tpm"
            out2 = open(outfilename, "w")
        else:
            out2 = open(args.out2, "w")

        reverse = open(args.reverse)
        for i in reverse:
            ls = i.rstrip().split("\t")
            if not re.match(r'#', i) and not re.match(r'__', i):
                ID = ls[0]
                count = float(ls[1])
                length = gffDict[ID]

                # WAGNER ET AL 2012
                tpm = (float(count) * float(reads) * 1000000) / (length * T)

                # ORIGINAL
                # pk = length / 1000
                # rpk = count / pk
                # tpm = rpk / denom

                out2.write(ID + "\t" + str(tpm) + "\n")
                totalTPM += tpm
    ############################################################################

    out.close()
    print(totalTPM)

else:
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

    totalRPK = 0
    T = 0
    for i in counts:
        ls = i.rstrip().split("\t")
        if not re.match(r'#', i) and not re.match(r'__', i):
            ID = ls[0]
            count = float(ls[1])
            length = gffDict[ID]

            # WAGNER ET AL 2012
            # t = (count * float(args.len)) / length
            # T += t

            # ORIGINAL
            pk = length/1000
            rpk = count/pk
            totalRPK += rpk


            ############################################################################
    # SECOND COUNTS TABLE
    if args.reverse != "NA":
        reverse = open(args.reverse)

        for i in reverse:
            ls = i.rstrip().split("\t")
            if not re.match(r'#', i) and not re.match(r'__', i):
                ID = ls[0]
                count = float(ls[1])
                length = gffDict[ID]

                # WAGNER ET AL 2012
                # t = (count * float(args.len)) / length
                # T += t

                # ORIGINAL
                pk = length / 1000
                rpk = count / pk
                totalRPK += rpk
    ############################################################################

    # WAGNER ET AL 2012
    # print(T)

    # ORIGINAL
    denom = totalRPK/1000000


    totalTPM = 0
    counts = open(args.counts)
    for i in counts:
        ls = i.rstrip().split("\t")
        if not re.match(r'#', i) and not re.match(r'__', i):
            ID = ls[0]
            count = float(ls[1])
            length = gffDict[ID]

            # WAGNER ET AL 2012
            # tpm = (float(count) * float(args.len) * 1000000) / (length * T)

            # ORIGINAL
            pk = length/1000
            rpk = count/pk
            tpm = rpk/denom

            out.write(ID + "\t" + str(tpm) + "\n")
            totalTPM += tpm

    ############################################################################
    if args.reverse != "NA":
        if args.out2 == "NA":
            outfilename = args.reverse
            outfilename = outfilename + ".tpm"
            out2 = open(outfilename, "w")
        else:
            out2 = open(args.out2, "w")

        reverse = open(args.reverse)
        for i in reverse:
            ls = i.rstrip().split("\t")
            if not re.match(r'#', i) and not re.match(r'__', i):
                ID = ls[0]
                count = float(ls[1])
                length = gffDict[ID]

                # WAGNER ET AL 2012
                # tpm = (float(count) * float(args.len) * 1000000) / (length * T)

                # ORIGINAL
                pk = length / 1000
                rpk = count / pk
                tpm = rpk / denom

                out2.write(ID + "\t" + str(tpm) + "\n")
                totalTPM += tpm
    ############################################################################

    out.close()
    print(totalTPM)













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


# TODO: in this version, I made pfam cross-validation optional


def ribosome(seq):
    NTs = ['T', 'C', 'A', 'G']
    stopCodons = ['TAA', 'TAG', 'TGA']
    Codons = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                codon = NTs[i] + NTs[j] + NTs[k]
                # if not codon in stopCodons:
                Codons.append(codon)

    CodonTable = {}
    AAz = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    AAs = list(AAz)
    k = 0
    for base1 in NTs:
        for base2 in NTs:
            for base3 in NTs:
                codon = base1 + base2 + base3
                CodonTable[codon] = AAs[k]
                k += 1

    prot = []
    for j in range(0, len(seq), 3):
        codon = seq[j:j + 3]
        try:
            prot.append(CodonTable[codon])
        except KeyError:
            prot.append("X")
    protein = ("".join(prot))
    return protein


def reverseComplement(seq):
    out = []
    for i in range(len(seq)-1, -1, -1):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        #->>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        #->>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    # data = sorted(data)
    data.sort(key=int)
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups


def lastItem(ls):
    x = ''
    for i in ls:
        x = i
    return x


def RemoveDuplicates(ls):
    empLS = []
    counter = 0
    for i in ls:
        if i not in empLS:
            empLS.append(i)
        else:
            pass
    return empLS


def allButTheLast(iterable, delim):
    x = []
    count = 0
    length = len(iterable.split(delim))
    for i in range(0, length):
        count += 1
        x.append(iterable.split(delim)[i])
        if count == length:
            x.append(delim)
    return delim.join(x[0:len(x)-2])


def secondToLastItem(ls):
    x = ''
    for i in ls[0:len(ls)-1]:
        x = i
    return x


def pull(item, one, two):
    ls = []
    counter = 0
    for i in item:
        if counter == 0:
            if i != one:
                pass
            else:
                counter += 1
                ls.append(i)
        else:
            if i != two:
                ls.append(i)
            else:
                ls.append(i)
                counter = 0
    outstr = "".join(ls)
    return outstr


def stabilityCounter(int):
    if len(str(int)) == 1:
        string = (str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 2:
        string = (str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 3:
        string = (str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 4:
        string = (str(0) + str(int))
        return (string)


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


def removeLS(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    return emptyList


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
    # print(count)
    return Dict


def filter(list, items):
    outLS = []
    for i in list:
        if i not in items:
            outLS.append(i)
    return outLS


def delim(line):
    ls = []
    string = ''
    for i in line:
        if i != " ":
            string += i
        else:
            ls.append(string)
            string = ''
    ls = filter(ls, [""])
    return ls


parser = argparse.ArgumentParser(
    prog="ORFS-to-BED.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************
    Developed by Arkadiy Garber^1;

    ^1University of Delaware, Geological Sciences
    Please send comments and inquiries to arkg@udel.edu
    *******************************************************
    '''))

parser.add_argument('-gff', type=str, help='contigs in FASTA format')

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]


# gffDict = defaultdict(list)
# totalLength = 0
# totalCoding = 0
# gff = open(args.gff)
# for i in gff:
#     if re.match(r'#', i):
#         if re.match(r'##FASTA', i):
#             break
#         elif re.match(r'##sequence-region', i):
#             ls = i.rstrip().split(" ")
#             length = int(ls[3])
#             totalLength += length
#         else:
#             pass
#     else:
#         ls = i.rstrip().split("\t")
#         if len(ls) > 1:
#             if ls[2] == "gene":
#                 gffDict[ls[0]].append([ls[3], ls[4]])
#                 geneLength = int(ls[4]) - int(ls[3])
#                 totalCoding += geneLength
#                 print(ls)
#                 print(geneLength)
#                 print("")
#
#
# print(totalLength)
# print(totalCoding)
#
# print(args.gff + "\t" + str(round(totalCoding/totalLength, 3)))


totalLength = 0
gffDict = defaultdict(list)

gffFile = args.gff
base = allButTheLast(gffFile, ".")
gbkFile = base + ".gb"



gff = open(args.gff)
for i in gff:
    if re.findall(r'FASTA', i):
        break
    else:
        if not re.match(r'#', i):
            ls = i.rstrip().split("\t")
            if len(ls) > 1:
                if ls[2] == "gene":
                    gffDict[ls[0]].append([int(ls[3]), int(ls[4])])
                    # print(ls)
                    # print(ls[8].split("ID=")[1].split(";")[0])
                    # print("")
        else:
            if re.match(r'##sequence-region', i):
                ls = i.rstrip().split(" ")
                length = int(ls[3])
                totalLength += length

totalIG = 0
for i in gffDict.keys():
    if len(gffDict[i]) > 2:
        end = 0
        start = 0
        # print(i)
        for j in sorted(gffDict[i]):
            start = str(j[0])
            if end > 0:
                IG = int(start) - int(end)
                if IG > 0:
                    totalIG += IG
            end = int(j[1])

coding = (totalLength - totalIG)
cd = coding/totalLength
print(args.gff + "\t" + str(cd))















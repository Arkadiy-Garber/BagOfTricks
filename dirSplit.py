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
    prog="dirSplit.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************
    
    Script for splitting a directory into multiple directories, each containing a subset of files
    
    Developed by Arkadiy Garber: agarber4@asu.edu

    *******************************************************
    '''))

parser.add_argument('-d', type=str, help='directory to split')

parser.add_argument('-b', type=str, help='number of batches to split into', default=2)

parser.add_argument('-f', type=str, help='comma-separated list of filename extensions to look for')

parser.add_argument('--gz', type=str,
                    help="files are gzipped", const=True, nargs="?")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]


DIR = os.listdir(args.d)

extensions = args.f

extensionsList = extensions.split(",")


if args.gz:
    count = 0
    for i in DIR:
        file = allButTheLast(i, ".")
        if lastItem(file.split(".")) == extensionsList[0]:
            count += 1

    print("total files with the %s extension: %s" % (extensionsList[0], str(count)))
    batches = count / int(args.b)
    batchSize = int(round(batches, 0))
    os.system("sleep 1")
    print("splitting into %s batches, each with about %s files" % (args.b, str(batchSize)))
    os.system("sleep 1\n")
    count = 1
    DIRdict = defaultdict(list)
    for i in DIR:
        file = allButTheLast(i, ".")
        extension = lastItem(file.split("."))
        if len(DIRdict[count]) == batchSize:
            print("batch_%s filled" % str(count))
            count += 1
            if extension == extensionsList[0]:
                DIRdict[count].append(i)
        else:
            if extension == extensionsList[0]:
                DIRdict[count].append(i)

    print("")
    for i in DIRdict:
        os.system("mkdir batch_%s" % str(i))
        print(i)
        for j in DIRdict[i]:
            file = allButTheLast(j, ".")
            basename = allButTheLast(file, ".")
            for k in extensionsList:
                outfile = basename + "." + k + ".gz"
                os.system("mv %s/%s batch_%s/%s > /dev/null 2>&1" % (args.d, outfile, i, outfile))

else:
    count = 0
    for i in DIR:
        if lastItem(i.split(".")) == extensionsList[0]:
            count += 1

    print("total files with the %s extension: %s" % (extensionsList[0], str(count)))
    batches = count/int(args.b)
    batchSize = int(round(batches, 0))
    os.system("sleep 1")
    print("splitting into %s batches, each with about %s files" % (args.b, str(batchSize)))
    os.system("sleep 1\n")
    count = 1
    DIRdict = defaultdict(list)
    for i in DIR:
        extension = lastItem(i.split("."))
        if len(DIRdict[count]) == batchSize:
            print("batch_%s filled" % str(count))
            count += 1
            if extension == extensionsList[0]:
                DIRdict[count].append(i)
        else:
            if extension == extensionsList[0]:
                DIRdict[count].append(i)

    print("")
    for i in DIRdict:
        os.system("mkdir batch_%s" % str(i))
        print(i)
        for j in DIRdict[i]:
            basename = allButTheLast(j, ".")
            for k in extensionsList:
                file = basename + "." + k
                os.system("mv %s/%s batch_%s/%s > /dev/null 2>&1" % (args.d, file, i, file))





























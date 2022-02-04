#!/usr/bin/env python3
from collections import defaultdict
import re
import statistics
import numpy as np
import os
import textwrap
import argparse

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


def Complement(seq):
    out = []
    for i in range(0, len(seq)):
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


def ribosome(seq):
    Dict = defaultdict(lambda: defaultdict(list))
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


def codonTable(seq):
    Dict = defaultdict(lambda: defaultdict(list))
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
            Dict[CodonTable[codon]][codon].append(codon)
            prot.append(CodonTable[codon])
        except KeyError:
            prot.append("X")
    protein = ("".join(prot))
    return Dict


def SeqCoord(seq, start, end):
    return seq[start:end]


def howMany(ls, exclude):
    counter = 0
    for i in ls:
        if i != exclude:
            counter += 1
    return counter


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
    if len(str(int)) > 4:
        string = str(int)
        return (string)


def sum(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count


def ave(ls):
    count = 0
    for i in ls:
        try:
            count += float(i)
        except ValueError:
            pass
    return count/len(ls)


def derep(ls):
    outLS = []
    for i in ls:
        if i not in outLS:
            outLS.append(i)
    return outLS


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


def GCcalc(seq):
    count = 0
    for i in seq:
        if i == "G" or i == "C":
            count += 1
    return count/len(seq)


def reject_outliers(data):
    m = 2
    u = np.mean(data)
    s = np.std(data)
    filtered = [e for e in data if (u - 2 * s < e < u + 2 * s)]
    return filtered


def lastItem(ls):
    x = ''
    for i in ls:
        if i != "":
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
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length-1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)-1]


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
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def fasta2(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

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


def allButTheFirst(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(1, length):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)]


def filter(list, items):
    outLS = []
    for i in list:
        if i not in items:
            outLS.append(i)
    return outLS


def filterRe(list, regex):
    ls1 = []
    ls2 = []
    for i in list:
        if re.findall(regex, i):
            ls1.append(i)
        else:
            ls2.append(i)
    return ls1, ls2


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
    prog="codonUsage",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    '''))

parser.add_argument('-g', help='genes sequences (nucleotides)')

parser.add_argument('-o', help='output')

parser.add_argument('--starts', type=str,
                    help="only look at start codons",
                    const=True, nargs="?")

args = parser.parse_args()


if args.starts:
    codonDict = defaultdict(list)
    hs = open(args.g)
    hs = fasta2(hs)
    AAs = {
        "M": "methionine",
        "F": "phenylalanine",
        "Y": "tyrosine",
        "P": "proline",
        "D": "aspartate",
        "H": "histidine",
        "V": "valine",
        "I": "isoleucine",
        "G": "glycine",
        "A": "alanine",
        "T": "threonine",
        "E": "glutamate",
        "S": "serine",
        "R": "arginine",
        "C": "cystine",
        "Q": "glutamine",
        "L": "leucine",
        "N": "asparagine",
        "K": "lysine",
        "W": "tryptophan",
        "*": "stop",

    }

    totalStarts = len(hs.keys())
    for i in hs.keys():
        seq = hs[i]
        start = seq[0:3]
        codonDict[start].append(start)

    out = open(args.o, "w")
    out.write("start_codon,amino_acid_letter,amino_acid_name,totalCount,perc_of_all_start_codons\n")
    for i in codonDict.keys():
        codons = len(codonDict[i])
        prop = len(codonDict[i])/totalStarts
        # prop = round(prop, 6)
        # print(i + "\t" + ribosome(i) + "\t" + AAs[ribosome(i)] + "\t" + str(codons) + "\t" + str(prop))
        out.write(i + "," + ribosome(i) + "," + AAs[ribosome(i)] + "," + str(codons) + "," + str(prop) + "\n")


else:
    codonDict = defaultdict(lambda: defaultdict(list))
    hs = open(args.g)
    hs = fasta2(hs)
    for i in hs.keys():
        seq = (hs[i])
        Dict = (codonTable(seq))
        for j in Dict.keys():
            for k in Dict[j]:
                codonDict[j][k].append(len(Dict[j][k]))

    totalCodons = 0
    for i in codonDict.keys():
        for j in codonDict[i]:
            totalCodons += sum(codonDict[i][j])

    AAs = {
        "M":"methionine",
        "F": "phenylalanine",
        "Y": "tyrosine",
        "P": "proline",
        "D": "aspartate",
        "H": "histidine",
        "V": "valine",
        "I": "isoleucine",
        "G": "glycine",
        "A": "alanine",
        "T": "threonine",
        "E": "glutamate",
        "S": "serine",
        "R": "arginine",
        "C": "cystine",
        "Q": "glutamine",
        "L": "leucine",
        "N": "asparagine",
        "K": "lysine",
        "W": "tryptophan",
        "*": "stop",

    }

    out = open(args.o, "w")
    out.write("amino_acid_letter,amino_acid_name,codon,totalCount,ratio,perc_of_all_codons\n")
    for i in sorted(codonDict.keys()):
        total = 0
        for j in sorted(codonDict[i]):
            Codons = sum(codonDict[i][j])
            total += Codons

        for j in sorted(codonDict[i]):
            Codons = sum(codonDict[i][j])
            print(i + "\t\t" + AAs[i] + "\t\t" + j + "\t" + str(Codons) + "\t\t" + str(Codons/total) + "\t\t" + str((Codons/totalCodons)*100 ))
            out.write(i + "," + AAs[i] + "," + j + "," + str(Codons) + "," + str(Codons / total) + "," + str((Codons / totalCodons) * 100) + "\n")
        out.write("#############################################\n")
        # print("")
    out.close()








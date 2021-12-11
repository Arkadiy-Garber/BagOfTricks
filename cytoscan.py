#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys
import time


# TODO: ADD CYTOCHROME 579 HMM
# TODO: ADD COLUMN WITH ORF STRAND


def SUM(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count

def firstNum(string):
    outputNum = []
    for i in string:
        try:
            int(i)
            outputNum.append(i)
        except ValueError:
            break
    Num = "".join(outputNum)
    return Num


def Num(ls):
    outputNum = "0"
    for i in ls:
        try:
            int(i)
            outputNum = i
        except ValueError:
            pass
    return outputNum


def Strip(ls):
    outList = []
    for i in ls:
        gene = i.split("|")[0]
        outList.append(gene)
    return outList

def unique(ls, ls2):
    unqlist = []
    for i in ls:
        if i not in unqlist and i in ls2:
            unqlist.append(i)
    return len(unqlist)

def Unique(ls):
    unqList = []
    for i in ls:
        if i not in unqList:
            unqList.append(i)
    return unqList

def Unique2(ls):
    unqList = []
    for i in ls:
        hmm = i.split("|")[0]
        if hmm not in unqList:
            unqList.append(hmm)
    return unqList


def checkDFE1(ls):
    count = 0
    uniqueLS = []
    for i in ls:
        hmm = i.split("|")[0]
        if hmm not in uniqueLS:
            uniqueLS.append(hmm)
            if hmm in ["DFE_0461", "DFE_0462", "DFE_0463", "DFE_0464", "DFE_0465"]:
                count += 1
    return count

def checkDFE2(ls):
    count = 0
    uniqueLS = []
    for i in ls:
        hmm = i.split("|")[0]
        if hmm not in uniqueLS:
            uniqueLS.append(hmm)
            if hmm in ["DFE_0448", "DFE_0449", "DFE_0450", "DFE_0451"]:
                count += 1
    return count

def checkGACE(ls):
    count = 0
    uniqueLS = []
    for i in ls:
        hmm = i.split("|")[0]
        if hmm not in uniqueLS:
            uniqueLS.append(hmm)
            if hmm in ["GACE_1843", "GACE_1844", "GACE_1845", "GACE_1846", "GACE_1847"]:
                count += 1
    return count

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

def lastItem(ls):
    x = ''
    for i in ls:
        x = i
    return x

def RemoveDuplicates(ls):
    empLS = []
    for i in ls:
        if i not in empLS:
            empLS.append(i)
        else:
            pass
    return empLS


def allButTheFirst(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(1, length):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)-1]


def allButTheLast(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length - 1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x) - 1]

def secondToLastItem(ls):
    x = ''
    for i in ls[0:len(ls) - 1]:
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

def remove2(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    # outString = "".join(emptyList)
    return emptyList

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
                seq = ''
            else:
                header = i[1:]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def fasta2(fasta_file):
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


def fastaRename(fasta_file):
    counter = 0
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
                counter += 1
                header = header + "_" + str(counter)
                seq = ''
            else:
                header = i[1:]
                header = header.split(" ")[0]
                counter += 1
                header = header + "_" + str(counter)
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
    prog="cytoscan.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    ************************************************************************

    Script for predicting periplasmic, outer-membrane, and extracellular c-type cytochromes

    Developed by Arkadiy Garber; agarber4@asu.edu
    ************************************************************************
    '''))


parser.add_argument('-p', type=str, help="input phobius output", default="NA")

parser.add_argument('-f', type=str, help="input protein seqs", default="NA")

parser.add_argument('-o', type=str, help="output CSV", default="NA")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]

print("Analyzing: " + args.p.split(".")[0])
Dict = defaultdict(list)
hemeDict = defaultdict(list)
protDict = defaultdict(lambda: 'EMPTY')
# prots = open(args.f)
prots = open(args.f)
prots = fasta2(prots)
for i in prots.keys():

    if re.findall(r'C(..)CH', prots[i]) or re.findall(r'C(...)CH', prots[i]):
        hemes = len(re.findall(r'C(..)CH', prots[i])) + len(re.findall(r'C(...)CH', prots[i]))
        Dict[i.split(" ")[0]].append("hemes")
        hemeDict[i.split(" ")[0]].append(hemes)

    protDict[i.split(" ")[0]] = allButTheFirst(i, " ")

# phobius = open(args.p)
phobius = open(args.p)
for i in phobius:
    ls = delim(i.rstrip())
    if ls[1] == "0" and ls[2] == "Y":
        Dict[ls[0]].append("SP")
        hemeDict[ls[0]].append("SP")

clusterDict = defaultdict(list)
for i in Dict.keys():
    numOrf = (int(lastItem(i.split("_"))))
    contig = allButTheLast(i, "_")
    clusterDict[contig].append(numOrf)


cluDict = defaultdict(list)
cluDict2 = defaultdict(list)
cluDict3 = defaultdict(list)
for i in clusterDict:
    clu = cluster(clusterDict[i], 1)
    for j in clu:
        rep = i + "_" + str(j[0])
        for k in j:
            locus = i + "_" + str(k)
            cluDict[rep].append(locus)
            cluDict2[rep].append(Dict[locus])
            cluDict3[rep].append(hemeDict[locus])

if args.o == "NA":
    outfile = args.p + "." + args.f  + ".cytoscan.csv"
else:
    outfile = args.o

out = open(outfile, "w")
out.write("locus,periplasmic_or_outer_membrane,hemes,seq\n")
for i in cluDict.keys():
    if len(cluDict[i]) > 0:
        countDict = defaultdict(list)
        counter = 0
        for j in cluDict2[i]:
            if len(list(j)) > 1:
                counter += 1
            for k in list(j):
                countDict[k].append(k)
        # print("-----------")
        # for j in countDict.keys():
        #     print(j + "\t" + str(len(countDict[j])))

        if len(countDict.keys()) > 1 and counter > 0:
            # print(cluDict[i])
            # print(cluDict3[i])
            for j in countDict.keys():
                # print(j + ":\t" + str(len(countDict[j])))
                pass

            for j in range(len(cluDict[i])):
                # print(cluDict[i][j] + "\t" + str(cluDict3[i][j]))
                if "SP" in cluDict3[i][j]:
                    SP = "Y"
                else:
                    SP = "N"

                hemes = Num(cluDict3[i][j])
                out.write(cluDict[i][j] + "," + SP + "," + str(hemes) + "," + prots[cluDict[i][j]] + "\n")

            out.write("#####################################################\n")
            # print("\n")
out.close()















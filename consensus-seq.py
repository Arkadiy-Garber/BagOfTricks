#!/usr/bin/env python3
from collections import defaultdict
import re
import textwrap
import argparse


def remove(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


def mostCommon(Dict):
    for i in sorted(Dict.keys()):
        highest = 0
        for j in Dict:
            aa = j
            count = Dict[j]
            if count > highest:
                highest = count
                consensus = aa
            if count == highest:
                if aa != "-":
                    consensus = aa
        return consensus


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


parser = argparse.ArgumentParser(
    prog="consensus-seq.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    Script for generating a representative consensus sequence from an alignment
    
    contact: agarber4@asu.edu
    '''))

parser.add_argument('-i', help='alignment file')


args = parser.parse_args()

aln = open(args.i)
aln = fasta(aln)
firstheader = list(aln.keys())[0]
seqLength = len(aln[firstheader])
clusterLength = len(aln.keys())

checkDict = defaultdict(list)
for i in aln.keys():
    checkDict[len(aln[i])].append(i)

if len(checkDict.keys()) > 1:
    print("Alignment file not formatted correctly. All sequences must be the same length.")
    raise SystemExit

Dict = defaultdict(lambda: defaultdict(list))
for j in aln.keys():
    seq = aln[j]

    for k in range(0, seqLength):
        Dict[k][seq[k]].append(1)

Dict2 = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
for j in Dict.keys():
    for k in Dict[j]:
        Dict2[j][k] = sum(Dict[j][k])

consensus = ""
for j in range(0, seqLength):

    residue = mostCommon(Dict2[j])
    consensus += residue
print(consensus)















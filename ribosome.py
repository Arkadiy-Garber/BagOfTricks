#!/usr/bin/env python3
from collections import defaultdict
import argparse
import textwrap
import sys
import re


def capitalizeCodon(codon):
    codonOut = ''
    for i in codon:
        if i == "a":
            codonOut += "A"
        elif i == "g":
            codonOut += "G"
        elif i == "c":
            codonOut += "C"
        elif i == "t":
            codonOut += "T"
        elif i == "u":
            codonOut += "U"
        else:
            codonOut += i
    return codonOut


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
    prog="riobosome.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************

    Script for translating gene sequences (nucleic acid format) into proteins (amino acid format)
    
    Developed by Arkadiy Garber: agarber4@asu.edu

    *******************************************************
    '''))

parser.add_argument('-i', type=str, help="input DNA sequence file in FASTA format")
parser.add_argument('-o', type=str, help="name output file to which protein sequences will be written")
parser.add_argument('-x', type=str, help="include the asterist denoting the stop codon? (y/n)", default="y")


if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]

print("reading in FASTA file")
Seqs = open(args.i)
Seqs = fasta(Seqs)


genestring = ''
NTs= [ 'T', 'C', 'A', 'G' ]
stopCodons = [ 'TAA', 'TAG', 'TGA' ]
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
            codon = base1+base2+base3
            CodonTable[codon] = AAs[k]
            k += 1

if args.o == args.i:
    print("Name of input file is same as output. Exiting...")
    raise SystemExit

count = 0
out = open(args.o, "w")
for i in Seqs.keys():
    prot = []
    count += 1
    perc = (count / len(Seqs.keys())) * 100
    sys.stdout.write("completed: %d%%   \r" % (perc))
    sys.stdout.flush()
    for j in range(0, len(Seqs[i]), 3):
        codon = Seqs[i][j:j + 3]
        codon = (capitalizeCodon(codon))
        try:
            if args.x == "y":
                prot.append(CodonTable[codon])
            else:
                if CodonTable[codon] != "*":
                    prot.append(CodonTable[codon])
                else:
                    pass
        except KeyError:
            prot.append("X")
    protein = ("".join(prot))
    out.write(">" + i + "\n")
    out.write(protein + "\n")

print("")
print("Done!")
out.close()
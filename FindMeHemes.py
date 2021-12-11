#!/usr/bin/env python3
# !/bin/sh
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys

parser = argparse.ArgumentParser(
    prog="FindMeHemes.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************
    *******************************************************
    *******************************************************

    Script for identifying proteins with heme c-binding motifs

    Developed by Arkadiy Garber: agarber4@asu.edu

    *******************************************************
    *******************************************************
    *******************************************************
    '''))

parser.add_argument('-contigs', type=str, help='File with contigs in FASTA format. if you provide this, the program '
                                               'will run Prodigal to identify candidate open reading frames, then use'
                                               'those proteins to search for candidate heme-binding '
                                               'sites.', default="NA")

parser.add_argument('-proteins', type=str, help="File with gene-calls in FASTA amino acid format. "
                                                "These proteins will be used to search for candidate heme-binding "
                                                "sites", default="NA")

# parser.add_argument('-outDir', type=str, help="name of directory to which output files will be written")

parser.add_argument('-mode', type=str, help="Is the provided FASTA file a single genome or a metagenome? "
                                                 "(default=genome). This is only relevant if you provide contigs, instead of proteins.", default="genome")

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
            seq += "\n"
    Dict[header] = seq
    return Dict


if args.contigs != "NA":
    print("Contigs file provided. Starting ORF prediction using Prodigal")
    if args.mode == "genome":
        os.system("prodigal -q -i " + args.contigs + " -a " + args.contigs + ".proteins.faa -o progidal.out")
        print("finished ORF prediction. Beginning identification of heme-binding motifs")
        prots = open(args.contigs + ".proteins.faa", "r")
        prots = fasta(prots)
        outFASTA = open(args.contigs + ".hemeProts.fasta", "w")
        outCSV = open(args.contigs + ".hemeProts.csv", "w")
        outCSV.write("ORF" + "," + "number_of_hemes" + "\n")
        for i in prots.keys():
            if re.findall(r'C(..)CH', prots[i]) or re.findall(r'C(...)CH', prots[i]) or re.findall(r'C(....)CH', prots[i]) or re.findall(r'C(..............)CH', prots[i]) or re.findall(r'C(...............)CH', prots[i]):
                outCSV.write(replace(i, [","], ";") + "," + str(
                    len(re.findall(r'C(..)CH', prots[i])) + len(re.findall(r'C(...)CH', prots[i])) + len(re.findall(r'C(....)CH', prots[i])) +
                    len(re.findall(r'C(..............)CH', prots[i])) + len(re.findall(r'C(...............)CH', prots[i]))) + "\n")
                outFASTA.write(">" + i + "\n")
                outFASTA.write(prots[i] + "\n")
        outFASTA.close()
        outCSV.close()

    if args.mode == "metagenome":
        os.system("prodigal -q -i " + args.contigs + " -a " + args.contigs + ".proteins.faa -p meta -o progidal.out")
        print("finished ORF prediction. Beginning identification of heme-binding motifs")
        prots = open(args.contigs + ".proteins.faa", "r")
        prots = fasta(prots)
        outFASTA = open(args.contigs + ".hemeProts.fasta", "w")
        outCSV = open(args.contigs + ".hemeProts.csv", "w")
        outCSV.write("ORF" + "," + "number_of_hemes" + "\n")
        for i in prots.keys():
            if re.findall(r'C(..)CH', prots[i]) or re.findall(r'C(...)CH', prots[i]) or re.findall(r'C(....)CH', prots[i]) or re.findall(r'C(..............)CH', prots[i]) or re.findall(r'C(...............)CH', prots[i]):
                outCSV.write(replace(i, [","], ";") + "," + str(
                    len(re.findall(r'C(..)CH', prots[i])) + len(re.findall(r'C(...)CH', prots[i])) + len(
                        re.findall(r'C(....)CH', prots[i])) +
                    len(re.findall(r'C(..............)CH', prots[i])) + len(re.findall(r'C(...............)CH', prots[i]))) + "\n")
                outFASTA.write(">" + i + "\n")
                outFASTA.write(prots[i] + "\n")
        outFASTA.close()
        outCSV.close()
    else:
        print("Please provide the program with some information about the type of file that is being processed. Is "
              "it a single genome, or a set of assembled contigs from a metagenome or multiple set of genomes? "
              "This can be answered via the -mode flag")


if args.proteins != "NA":
    print("Amino acid FASTA file is provided. Skipping Prodigal, and moving on to identification of heme-binding motifs")
    prots = open(args.proteins, "r")
    prots = fasta(prots)

    outFASTA = open(args.proteins + ".hemeProts.fasta", "w")
    outCSV = open(args.proteins + ".hemeProts.csv", "w")
    outCSV.write("ORF" + "," + "number_of_hemes" + "\n")
    for i in prots.keys():
        if re.findall(r'C(..)CH', prots[i]) or re.findall(r'C(...)CH', prots[i]) or re.findall(r'C(....)CH', prots[i]) or re.findall(r'C(..............)CH', prots[i]) or re.findall(r'C(...............)CH', prots[i]):
            outCSV.write(replace(i, [","], ";") + "," + str(
                len(re.findall(r'C(..)CH', prots[i])) + len(re.findall(r'C(...)CH', prots[i])) + len(
                    re.findall(r'C(....)CH', prots[i])) +
                len(re.findall(r'C(..............)CH', prots[i])) + len(re.findall(r'C(...............)CH', prots[i]))) + "\n")
            outFASTA.write(">" + i + "\n")
            outFASTA.write(prots[i] + "\n")
    outFASTA.close()
    outCSV.close()

print("All done. Thank you for using FindMeHemes!")
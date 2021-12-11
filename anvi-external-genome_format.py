#!/usr/bin/env python3
# !/bin/sh
# Author: Arkadiy Garber


from collections import defaultdict
import re, os
import argparse
import textwrap
import statistics


parser = argparse.ArgumentParser(
    prog="anvi-external-genome_format.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    ************************************************************************
    
    Script for making a tab-delimited file for ANVIO. 
    This file is required for inputting external genomes into the metagenomic workflow
    
    Developed by Arkadiy Garber; agarber4@asu.edu
    ************************************************************************
    '''))


parser.add_argument('-folder', help='folder with fasta files')


parser.add_argument('-output', help="Name of output TSV file")

parser.add_argument('-fasta', help="suffix for fasta file")

args = parser.parse_known_args()[0]

folder = os.listdir(args.folder)
out = open(args.output, "w")
out.write("name" + "\t" + "contigs_db_path" + "\n")
for i in folder:
    if re.search(args.fasta, i):
        length = len(args.fasta)
        out.write(i[0:len(i)-length] + "\t" + i[0:len(i)-length] + ".db" + "\n")
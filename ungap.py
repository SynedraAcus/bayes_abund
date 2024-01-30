#! /usr/bin/env python3

from argparse import ArgumentParser

from Bio import SeqIO

parser = ArgumentParser("Un-gap SILVA alignment for use as reference")
parser.add_argument("-f", type=str, help="Input SILVA alignment")
parser.add_argument("-o", type=str, help="Output file")
args = parser.parse_args()

with open(args.o, mode="w") as outfile:
    for record in SeqIO.parse(args.f, "fasta"):
        ungapped = record.seq.replace(".", "").replace("-", "")
        print(f">{record.description}\n{ungapped}", file=outfile)

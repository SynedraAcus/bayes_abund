#! /usr/bin/env python3

from argparse import ArgumentParser

from Bio import SeqIO

from src.utils import extract_taxon

parser = ArgumentParser("Calculate expected abundances")
parser.add_argument("-d", type=str, help="SILVA reference alignment")
parser.add_argument("-r", type=str, help="Reads file")
parser.add_argument("-o", help="Phylotype output TSV")
parser.add_argument(
    "--rank",
    type=str,
    default="Genus",
    help="Taxonomic rank for classification. Defaults to genus",
)
parser.add_argument(
    "--pseudocount",
    type=float,
    default=0.1,
    help="Pseudocount for k-mers absent in training sequences (between 0 and 1, defaults to 0.1)",  # noqa: E501
)
args = parser.parse_args()


valid_ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
try:
    rank = valid_ranks.index(args.rank.lower())
except ValueError:
    print(
        "Invalid rank {args.rank}. Please use one of kingdom, phylum, class, order, family, genus or species"  # noqa: E501
    )
# Parse reference DB
for record in SeqIO.parse(args.d, "fasta"):
    # TODO: process eukaryotes' taxonomies properly
    if extract_taxon(record.description, 0) == "Eukaryota":
        # Skipping eukaryotes for now, they use a different taxonomy scheme
        continue
    taxon = extract_taxon(record.description, rank)
    # TODO: prepare an un-aligned SILVA instead of doing this on the fly
    ungapped = record.seq.replace(".", "").replace("-", "")
    print(ungapped)
    quit()

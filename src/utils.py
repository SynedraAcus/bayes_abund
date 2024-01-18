"""
Various utility functions for classifier
"""


def extract_taxon(header: str, position: int = 0) -> list[str]:
    """
    Extract the taxon of a required rank from SILVA-formatted FASTA header.

    If this element is empty, return the smallest non-empty taxon prefixed with
    'unclassified '.
    """
    s = header.split("\t")[2]
    taxa = s.split(";")
    if position <= len(taxa) - 1 and taxa[position] != "":
        return taxa[position]
    else:
        taxa.remove("")
        return "Unclassified " + taxa[-1]

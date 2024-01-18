"""
Tests for utilities
"""
from src.utils import extract_taxon


def test_taxa_extraction():
    header = "AY214344.RvmSpec4\t100\tBacteria;Proteobacteria;Alphaproteobacteria;Rhodobacterales;Rhodobacteraceae;Oceanicella;"  # noqa: E501
    taxa = [
        "Bacteria",
        "Proteobacteria",
        "Alphaproteobacteria",
        "Rhodobacterales",
        "Rhodobacteraceae",
        "Oceanicella",
        "Unclassified Oceanicella",
    ]
    assert [extract_taxon(header, x) for x in range(7)] == taxa

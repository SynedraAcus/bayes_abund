"""
Tests
"""
import pytest

from src.tree import KmerTree, KmerTreeException


def test_tree_gen():
    t = KmerTree(3, "ACTG", "ACACACC")
    assert t.get_count("ACA") == 2
    assert t.get_count("CAC") == 2
    assert t.get_count("ACC") == 1
    assert t.get_count("GTG") == 0


def test_freq():
    t = KmerTree(3, "ACTG", "ACACACC")
    assert t.get_freq("ACA") == 0.4
    assert t.get_freq("CAC") == 0.4
    assert t.get_freq("ACC") == 0.2
    assert t.get_freq("GTG") == 0


def test_seq_probs():
    """
    Just a few sanity checks, nothing particularly deep
    :return:
    """
    t1 = KmerTree(3, "ACTG", "ACACACC")
    t2 = KmerTree(3, "ACTG", "ACACACT")
    # No info to differentiate these two
    assert t1.get_seq_prob("GGACAC") == t2.get_seq_prob("GGACAC")
    # Definitely more like t1
    assert t1.get_seq_prob("CACCGG") > t2.get_seq_prob("CACCGG")


def test_pseudocount_effect():
    t = KmerTree(3, "ACTG", "ACACACC")
    assert t.get_seq_prob("CACCGG", pseudocount=0.2) > t.get_seq_prob(
        "CACCGG", pseudocount=0.1
    )


def test_init_exceptions():
    with pytest.raises(KmerTreeException):
        KmerTree(3, "ACTG", "ACTGACTGU\n")


def test_method_exceptions():
    t = KmerTree(3, "ACTG", "ACTGACTGACTG")
    with pytest.raises(KmerTreeException):
        # Wrong k
        t.get_count("ACTG")
    with pytest.raises(KmerTreeException):
        t.get_seq_prob("ACUGACUUUC")

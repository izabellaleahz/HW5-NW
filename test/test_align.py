# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    Test that the 3 alignment matrices are filled correctly
    for MYQR vs MQR using BLOSUM62, gap open -10, gap extend -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    score, alignA, alignB = nw.align(seq1, seq2)

    # matrices should be (len(seq1)+1) x (len(seq2)+1)
    assert nw._align_matrix.shape == (5, 4)
    assert nw._gapA_matrix.shape == (5, 4)
    assert nw._gapB_matrix.shape == (5, 4)

    # check score — optimal alignment is MYQR / M-QR
    # M-M(5) + gap(-10 + -1) + Q-Q(5) + R-R(5) = 4
    assert score == 4

    # spot-check align matrix along the optimal path
    assert nw._align_matrix[0][0] == 0
    assert nw._align_matrix[1][1] == 5    # M-M match
    assert nw._align_matrix[2][2] == 4    # Y-Q via diagonal
    assert nw._align_matrix[4][3] == 4    # R-R final

    # check gap boundaries are initialized correctly
    assert nw._gapB_matrix[1][0] == -11
    assert nw._gapB_matrix[2][0] == -12
    assert nw._gapA_matrix[0][1] == -11
    assert nw._gapA_matrix[0][2] == -12

    # gapB entry used in optimal path (Y aligned to gap)
    assert nw._gapB_matrix[2][1] == -6


def test_nw_backtrace():
    """
    Test backtrace on MAVHQLIRRP vs MQLIRHP using BLOSUM62.
    Expected score: 17, expected alignment: MAVHQLIRRP / M---QLIRHP
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    score, alignA, alignB = nw.align(seq3, seq4)

    assert score == 17
    assert alignA == "MAVHQLIRRP"
    assert alignB == "M---QLIRHP"


def test_nw_identical():
    """Aligning a sequence to itself should produce no gaps."""
    seq1, _ = read_fasta("./data/test_seq1.fa")

    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    score, alignA, alignB = nw.align(seq1, seq1)

    assert "-" not in alignA
    assert "-" not in alignB
    assert alignA == alignB
    assert score > 0
#!/usr/bin/env python3
"""Counts the number of contigs in a given FASTA file."""

import sys

from commons import number_of_contigs

if __name__ == "__main__":
    print(number_of_contigs(sys.argv[1]))

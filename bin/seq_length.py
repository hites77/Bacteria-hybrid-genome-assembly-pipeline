#!/usr/bin/env python3

"""Check the total length of all the contigs in a fasta file and print it to stdout."""

import sys

from Bio import SeqIO


def fasta_length(fasta_file):
    records = SeqIO.parse(fasta_file, "fasta")
    return sum(len(rec) for rec in records)


if __name__ == "__main__":
    print(fasta_length(sys.argv[1]))

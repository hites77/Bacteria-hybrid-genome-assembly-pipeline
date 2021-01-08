#!/usr/bin/env python3

from Bio import SeqIO


def get_contig_names(fasta_file):
    records = SeqIO.parse(fasta_file, "fasta")
    return [rec.id for rec in records]

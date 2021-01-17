#!/usr/bin/env python3

import re
from Bio import SeqIO


def make_flag(name):
    return "--" + name


def bool_to_str(boolean):
    if not isinstance(boolean, bool):
        raise ValueError("boolean must be a bool")
    return "yes" if boolean else "no"


def string_to_bases(bases_str):
    """Converts a string such as '10k', '10m', '10g', or '10' into the number of bases."""

    SUFFIX_DICT = {"k": 1000, "m": 10 ** 6, "g": 10 ** 9}

    if not isinstance(bases_str, str):
        raise ValueError

    if not bases_str[-1].isdigit():
        suffix = bases_str[-1].lower()
        number = int(bases_str[:-1])
        try:
            multiplier = SUFFIX_DICT[suffix]
        except KeyError:
            raise ValueError("Unrecognised suffix. Must be one of k,m,g.")
    else:
        number = int(bases_str)
        multiplier = 1

    return number * multiplier


# TODO: port to biopython ?
def number_of_contigs(assembly_file):
    """Counts the number of contigs in a fasta file."""
    text = open(assembly_file).read()
    matches = re.findall("^>", text, flags=re.MULTILINE)
    return len(matches)


def contig_lengths(assembly_fasta):
    records = SeqIO.parse(assembly_fasta, "fasta")
    return {rec.id: len(rec) for rec in records}

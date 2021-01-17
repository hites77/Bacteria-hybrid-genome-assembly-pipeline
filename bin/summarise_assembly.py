#!/usr/bin/env python3

import json
import sys

from commons import contig_lengths, number_of_contigs


def make_assembly_summary(assembly_fasta, circularity_summary):
    """
    Args:
        assembly_fasta: Path to assembly fasta file.
        circularity_summary: Path to circularity summary, which should be a JSON file of
          contig name: circularity. eg. { "contig 1" : "circular" }
    """

    return {
        "number of contigs": number_of_contigs(assembly_fasta),
        "contig lengths": contig_lengths(assembly_fasta),
        "contig circularity": json.load(open(circularity_summary)),
    }


if __name__ == "__main__":
    summary = make_assembly_summary(sys.argv[1], sys.argv[2])
    print(json.dumps(summary, indent=4))

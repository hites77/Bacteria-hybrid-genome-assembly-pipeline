#!/usr/bin/env python3

import json
import sys
from pathlib import Path

from commons import contig_lengths, number_of_contigs


def flye_may_have_found_plasmids(flye_dir):
    """Returns true if Flye *may* have found plasmids.

    As Flye does not identify which contigs are plasmids in its final assembly, we can
    only check the contents of Flye's intermediate files created when it tries to identify
    possible plasmids. It is possible that these plasmids initially identified by Flye get
    merged into other contigs in later stages.
    """
    raw_plasmid_file = Path(flye_dir) / "22-plasmids/plasmids_raw.fasta"
    return raw_plasmid_file.stat().st_size > 0


def make_assembly_summary(assembly_fasta, circularity_summary, flye_dir):
    """
    Args:
        assembly_fasta (str): Path to assembly fasta file.
        circularity_summary (str): Path to circularity summary, which should be a JSON file of
          contig name: circularity. eg. { "contig 1" : "circular" }
        flye_dir (str): Path to Flye output directory.
    """

    return {
        "number of contigs": number_of_contigs(assembly_fasta),
        "contig lengths": contig_lengths(assembly_fasta),
        "contig circularity": json.load(open(circularity_summary)),
        "flye found plasmids": "maybe" if flye_may_have_found_plasmids(flye_dir) else "no",
    }


if __name__ == "__main__":
    summary = make_assembly_summary(
        assembly_fasta=sys.argv[1], circularity_summary=sys.argv[2], flye_dir=sys.argv[3]
    )
    print(json.dumps(summary, indent=4))

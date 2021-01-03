#!/usr/bin/env python3

"""Checks whether 2 genome assemblies are different.

See the documentation for the main function.
"""

import re
import sys

import pandas as pd
from Bio import SeqIO

### CONSTANTS ###

YES = "yes"
NO = "no"
# columns
COL_DIFFERENT = "different"
COL_CONTIGS_CHANGE = "contigsChange"
COL_BEFORE_CONTIGS = "beforeContigs"
COL_AFTER_CONTIGS = "afterContigs"

### FUNCTIONS ###


def sequences_different(assembly_file1, assembly_file2):
    """Returns whether 2 assemblies are different.

    Precondition: number of contigs are the same, and 1st contigs correspond, 2nd contigs
    correspond etc

    Args:
        assembly_file1 (str): Path to 1st assembly, in FASTA format.
        assembly_file2 (str): Path to 2nd assembly, in FASTA format.

    Returns:
        True if the assemblies are different, False if they are the same.
    """
    a1_record = SeqIO.parse(assembly_file1, "fasta")
    a2_record = SeqIO.parse(assembly_file2, "fasta")

    for record1, record2 in zip(a1_record, a2_record):
        if record1.seq != record2.seq:
            return True
    return False


# TODO: port to biopython ?
def number_of_contigs(assembly_file):
    text = open(assembly_file).read()
    matches = re.findall("^>", text, flags=re.MULTILINE)
    return len(matches)


def main(assembly1, assembly2):
    """Analyses the difference between 2 genome assemblies and outputs the results in the tsv
    format to stdout.

    Assumptions:
    - When checking the equality of 2 assemblies with the same number of contigs, it is
      assumed that the 1st contigs correspond, 2nd contigs correspond etc.
    - Sequences are treated as linear.

    Args:
        assembly1 (str): Path to file of original assembly (FASTA format).
        assembly2 (str): Path to file of new assembly (FASTA format).

    Returns:
        Dataframe of differences between the 2 assemblies. The dataframe has a single row,
        and has the following columns:
        - different: whether the assemblies are different
        - contigsChange: whether the number of contigs differs
        - beforeContigs: number of contigs in assembly 1
        - afterContigs: number of contigs in assembly 2
    """
    before_contigs = number_of_contigs(assembly1)
    after_contigs = number_of_contigs(assembly2)
    if before_contigs != after_contigs:
        different = YES
        contigs_change = YES
    elif sequences_different(assembly1, assembly2):
        different = YES
        contigs_change = NO
    else:
        different = NO
        contigs_change = NO

    df = pd.DataFrame(
        {
            COL_DIFFERENT: [different],
            COL_CONTIGS_CHANGE: [contigs_change],
            COL_BEFORE_CONTIGS: [before_contigs],
            COL_AFTER_CONTIGS: [after_contigs],
        }
    )

    df.to_csv(sys.stdout, index=False)

    return df


def test_main():
    main(
        "/home/chloe/Documents/NUS/UROPS/server-data/test-pipeline/assembly/pilon/final_assembly.fasta",
        "/home/chloe/Documents/NUS/UROPS/server-data/test-pipeline/assembly/pilon/final_assembly.fasta",
    )

    main(
        "/home/chloe/Documents/NUS/UROPS/server-data/test-pipeline/assembly/pilon/final_assembly.fasta",
        "/home/chloe/Documents/NUS/UROPS/server-data/test-pipeline/assembly/pilon/final_assembly_doubled.fasta",
    )

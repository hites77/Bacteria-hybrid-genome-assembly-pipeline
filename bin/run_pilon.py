#!/usr/bin/env python3

import argparse
import re
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from commons import make_flag

### CONSTANTS ###

PILON_PREFIX = "pilon"
PILON_INFO_FILE = "pilon_info.tsv"

### FUNCTIONS ###


def run(command, check=True):
    return subprocess.run(command, shell=True, check=check)


def preprocess(assembly_name, illumina1_fq, illumina2_fq, threads):
    bam_file = Path(assembly_name).name + ".bam"  # remove directory, if any
    run(f"bwa index {assembly_name}")
    run(
        f"bwa mem -t {threads} {assembly_name} {illumina1_fq} {illumina2_fq} | samtools sort --threads {threads} > {bam_file}"
    )
    run(f"samtools index {bam_file}")
    return bam_file


def run_pilon(assembly_name, bam_file, index, pilon_args):
    run(
        f"java -Xmx13G -jar $PILONJAR --genome {assembly_name} --frags {bam_file} --changes --output {PILON_PREFIX}{index} {pilon_args}"
    )


def remove_temp_files():
    run("rm *.bai", check=False)
    run("rm *.bwt", check=False)
    run("rm *.amb", check=False)
    run("rm *.ann", check=False)
    run("rm *.pac", check=False)
    run("rm *.sa", check=False)


def no_changes(changes_file):
    return changes_file.read_text().strip() == ""


SUFFIX = re.compile("(_pilon)+$")


def strip_suffix(name):
    return SUFFIX.sub("", name)


def strip_pilon_suffix(input_fasta, output_fasta):
    records = SeqIO.parse(input_fasta, "fasta")
    new_records = (SeqRecord(rec.seq, id=strip_suffix(rec.id), description="") for rec in records)
    SeqIO.write(new_records, output_fasta, "fasta")


def main(
    assembly_name, output_assembly, illumina1_fq, illumina2_fq, threads, max_iters, pilon_args
):
    if max_iters < 1:
        raise ValueError(f"{max_iters} should be at least 1")

    changes_file = None
    i = 0
    for i in range(1, max_iters + 1):
        bam_file = preprocess(assembly_name, illumina1_fq, illumina2_fq, threads)
        run_pilon(assembly_name, bam_file, i, pilon_args)
        changes_file = Path(f"pilon{i}.changes")
        assembly_name = f"pilon{i}.fasta"
        if no_changes(changes_file):
            break

    if no_changes(changes_file):  ## pilon converged
        final_assembly_name = f"pilon{i-1}.fasta"
        num_iters = i - 1
        converged = True

        Path(assembly_name).unlink()  # remove since no changes
    else:  ## pilon did not converge
        # print warning
        warning = (
            f"WARNING: pilon failed to converge after {max_iters} iterations."
            + f"\nYou may want to run pilon for more iterations."
            + f"\nSee {PILON_PREFIX}X.changes for the changes made in each iteration."
        )
        print(warning, file=sys.stderr)

        final_assembly_name = assembly_name
        num_iters = i
        converged = False

    open(PILON_INFO_FILE, "w").write("iterations\tconverged\n" + f"{num_iters}\t{converged}")

    strip_pilon_suffix(final_assembly_name, output_assembly)

    remove_temp_files()


### PARSER ####

# constants
PARSER_ASSEMBLY = "assembly"
PARSER_READS_1 = "reads1"
PARSER_READS_2 = "reads2"
PARSER_THREADS = "threads"
PARSER_MAX_ITERS = "maxiters"
PARSER_OUT_ASSEMBLY = "out"
PARSER_PILON_ARGS = "args"


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(make_flag(PARSER_ASSEMBLY), required=True, help="Path to assembly")
    parser.add_argument(
        make_flag(PARSER_READS_1),
        required=True,
        help="Path to one of the short read files",
    )
    parser.add_argument(
        make_flag(PARSER_READS_2),
        required=True,
        help="Path to the other short read file",
    )
    parser.add_argument(
        make_flag(PARSER_THREADS),
        type=int,
        required=True,
        help="Number of threads to use",
    )
    parser.add_argument(
        make_flag(PARSER_MAX_ITERS),
        type=int,
        required=True,
        help="Maximum number of iterations to run pilon",
    )
    parser.add_argument(
        make_flag(PARSER_OUT_ASSEMBLY), required=True, help="Path to output assembly"
    )
    parser.add_argument(
        make_flag(PARSER_PILON_ARGS),
        required=False,
        default="",
        help="Command line args other than inputs, outputs and --changes to pass to pilon",
    )
    return parser


### MAIN ###

if __name__ == "__main__":
    parser = make_parser()
    args = vars(parser.parse_args())
    main(
        assembly_name=args[PARSER_ASSEMBLY],
        illumina1_fq=args[PARSER_READS_1],
        illumina2_fq=args[PARSER_READS_2],
        threads=args[PARSER_THREADS],
        max_iters=args[PARSER_MAX_ITERS],
        output_assembly=args[PARSER_OUT_ASSEMBLY],
        pilon_args=args[PARSER_PILON_ARGS],
    )

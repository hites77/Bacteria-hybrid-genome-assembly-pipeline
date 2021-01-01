#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path

### CONSTANTS ###

PILON_PREFIX = "pilon"
FINAL_ASSEMBLY_SYMLINK = "final_assembly.fasta"

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


def run_pilon(assembly_name, bam_file, index):
    run(
        f"java -Xmx13G -jar $PILONJAR --genome {assembly_name} --frags {bam_file} --changes --output {PILON_PREFIX}{index}"
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


def main(assembly_name, illumina1_fq, illumina2_fq, threads, max_iters):
    if max_iters < 1:
        raise ValueError(f"{max_iters} should be at least 1")

    # this file is present while pilon has yet to / failed to converge
    incomplete_file = Path("incomplete")

    incomplete_file.touch()

    changes_file = None
    i = 0
    for i in range(1, max_iters + 1):
        bam_file = preprocess(assembly_name, illumina1_fq, illumina2_fq, threads)
        run_pilon(assembly_name, bam_file, i)
        changes_file = Path(f"pilon{i}.changes")
        assembly_name = f"pilon{i}.fasta"
        if no_changes(changes_file):
            break

    # run pilon
    # bamFile = preprocess(assemblyFasta, illumina1Fastq, illumina2Fastq, threads)
    # run_pilon(assemblyFasta, bamFile, 1)
    # changes_file = None
    # for i in range(2, maxIters + 1):
    #     changes_file = Path(f"pilon{i-1}.changes")
    #     assemblyFasta = f"pilon{i-1}.fasta"
    #     if no_changes(changes_file):
    #         break
    #     bamFile = preprocess(assemblyFasta, illumina1Fastq, illumina2Fastq, threads)
    #     run_pilon(assemblyFasta, bamFile, i)

    if no_changes(changes_file):  ## pilon converged
        incomplete_file.unlink()

        final_assembly_name = f"pilon{i-1}.fasta"

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

    remove_temp_files()

    Path(FINAL_ASSEMBLY_SYMLINK).symlink_to(final_assembly_name)


### PARSER ####

# constants
PARSER_ASSEMBLY = "assembly"
PARSER_READS_1 = "reads1"
PARSER_READS_2 = "reads2"
PARSER_THREADS = "threads"
PARSER_MAX_ITERS = "maxiters"


def make_flag(name):
    return "--" + name


parser = argparse.ArgumentParser()
parser.add_argument(
    make_flag(PARSER_ASSEMBLY), nargs=1, required=True, help="Path to assembly"
)
parser.add_argument(
    make_flag(PARSER_READS_1),
    nargs=1,
    required=True,
    help="Path to one of the short read files",
)
parser.add_argument(
    make_flag(PARSER_READS_2),
    nargs=1,
    required=True,
    help="Path to the other short read file",
)
parser.add_argument(
    make_flag(PARSER_THREADS),
    nargs=1,
    type=int,
    required=True,
    help="Number of threads to use",
)
parser.add_argument(
    make_flag(PARSER_MAX_ITERS),
    nargs=1,
    type=int,
    required=True,
    help="Maximum number of iterations to run pilon",
)

### EXECUTION ###

if __name__ == "__main__":
    args = vars(parser.parse_args())
    print(args)
    main(
        assembly_name=args[PARSER_ASSEMBLY][0],
        illumina1_fq=args[PARSER_READS_1][0],
        illumina2_fq=args[PARSER_READS_2][0],
        threads=args[PARSER_THREADS][0],
        max_iters=args[PARSER_MAX_ITERS][0],
    )

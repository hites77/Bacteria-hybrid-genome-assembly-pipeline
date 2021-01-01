#!/usr/bin/env python3

import subprocess
import sys
from pathlib import Path

PILON_PREFIX = "pilon"
FINAL_ASSEMBLY_SYMLINK = "final_assembly.fasta"


def run(command, check=True):
    return subprocess.run(command, shell=True, check=check)


def preprocess(assembly_name, illumina1_fq, illumina2_fq, threads):
    bam_file = assembly_name + ".bam"
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

        final_assembly_name = f"pilon{i-2}.fasta"

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

    Path(final_assembly_name).symlink_to(FINAL_ASSEMBLY_SYMLINK)


# run('$assemblyFasta', '$illumina1Fq', '$illumina2Fq', $params.threads, 6)

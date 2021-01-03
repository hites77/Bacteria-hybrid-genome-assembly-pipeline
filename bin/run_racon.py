#!/usr/bin/env python3

"""Repeatedly runs racon for a certain number of times or until no more changes are made
to the assembly, whichever comes first.

See the documentation for the main function for details.
"""

import argparse
import subprocess
from pathlib import Path

import assembly_diff

### CONSTANTS ###

# prefix used for intermediate files
RACON_PREFIX = "racon"

COL_RACON_ITERATIONS = "raconIterations"


### FUNCTIONS ###


def run(command, check=True):
    return subprocess.run(command, shell=True, check=check)


def run_racon(assembly, pacbio, output_assembly, log, threads, racon_args):
    samfile = Path(assembly).name + ".sam"  # current directory only
    run(f"minimap2 -t {threads} -ax map-pb {assembly} {pacbio} > {samfile}")
    run(f"racon {racon_args} -t {threads} {pacbio} {samfile} {assembly} > {output_assembly}")
    log_df = assembly_diff.main(assembly, output_assembly)
    log_df.to_csv(log, index=False)
    return log_df


def has_difference(log_df):
    return log_df[assembly_diff.COL_DIFFERENT][0] == assembly_diff.YES


def main(assembly, pacbio, output_prefix, threads, max_iters, racon_args):
    """Run racons up to max_iters number of times or until no more changes are made to the
    assembly, whichever comes first.

    Output files created:
    - `raconX.fa` (where X = iteration number): Assembly after iteration X of racon.
    - `racon_log_X.tsv`: Changes to the assembly after iteration X of racon, generated using
      assembly_diff.py.
    - `<output_prefix>_assembly.fa`: Final assembly. This is a symlink to the relevant raconX.fa file.
    - `<output_prefix>_log.tsv`: Difference between initial and final assembly, generated using
      assembly_diff.py. There is also an additional column containing the number of iterations
      racon was run for.
    - (note that `<output_prefix>` is substituted with the value of the output_prefix parameter)
    - various sam files: intermediate files.

    Args:
        assembly (str): Path to genome assembly file in FASTA format.
        pacbio (str): Path to pacbio reads file.
        output_prefix (str): prefix to use for final assembly file and final log file.
        threads (int): Number of threads.
        max_iters (int): Maximum number of times to run racon. Must be at least 1.
        racon_args (str): Additional argument to pass to racon.
    """
    if max_iters < 1:
        raise ValueError("max_iters must be at least 1")

    current_output = ""
    num_iters = 0
    for i in range(1, max_iters + 1):
        num_iters = i
        current_input = f"{RACON_PREFIX}{i-1}.fa" if i > 1 else assembly
        current_output = f"{RACON_PREFIX}{i}.fa"
        log_df = run_racon(
            assembly=current_input,
            pacbio=pacbio,
            output_assembly=current_output,
            log=f"{RACON_PREFIX}_log_{i}.tsv",
            threads=threads,
            racon_args=racon_args,
        )

        if not has_difference(log_df):
            Path(current_output).unlink()  # remove assembly with no changes
            num_iters = i - 1
            if i > 1:
                final_assembly = current_input
            break

    final_assembly = current_output

    # make final log
    # append column: number of racon iterations
    log_df = assembly_diff.main(assembly, final_assembly)
    log_df[COL_RACON_ITERATIONS] = num_iters
    log_df.to_csv(f"{output_prefix}_log.tsv", sep="\t", index=False)

    # make symlink pointing to final assembly
    Path(f"{output_prefix}_assembly.fa").symlink_to(final_assembly)


def test_main():
    main(
        "/home/chloe/Documents/NUS/UROPS/server-data/test-pipeline/assembly/pilon/final_assembly.fasta",
        "/home/chloe/Documents/NUS/UROPS/server-data/S8E_3_1/reads/cleaned/pacbio.fq",
        "final_racon",
        8,
        2,
    )


### PARSER ###

# constants
FLAG_IN_ASSEMBLY = "in_assembly"
FLAG_IN_PACBIO = "in_pacbio"
FLAG_OUT_PREFIX = "out_prefix"
FLAG_THREADS = "threads"
FLAG_MAXITERS = "maxiters"
FLAG_RACON_ARGS = "args"


def make_flag(name):
    return "--" + name


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(make_flag(FLAG_IN_ASSEMBLY), required=True)
    parser.add_argument(make_flag(FLAG_IN_PACBIO), required=True)
    parser.add_argument(make_flag(FLAG_OUT_PREFIX), required=True)
    parser.add_argument(make_flag(FLAG_THREADS), required=True, type=int)
    parser.add_argument(make_flag(FLAG_MAXITERS), required=True, type=int)
    parser.add_argument(make_flag(FLAG_RACON_ARGS), required=False, default="")
    return parser


### MAIN ###

if __name__ == "__main__":
    parser = make_parser()
    args = vars(parser.parse_args())
    main(
        assembly=args[FLAG_IN_ASSEMBLY],
        pacbio=args[FLAG_IN_PACBIO],
        output_prefix=args[FLAG_OUT_PREFIX],
        threads=args[FLAG_THREADS],
        max_iters=args[FLAG_MAXITERS],
        racon_args=args[FLAG_RACON_ARGS],
    )

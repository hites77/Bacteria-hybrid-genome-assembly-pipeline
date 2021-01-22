#!/usr/bin/env python3

import subprocess
import sys
from pathlib import Path


def run_fastani(query_genome, reference_genomes_list, output_tsv):
    subprocess.run(
        f"fastANI -q {query_genome} --rl {reference_genomes_list} -o {output_tsv}",
        shell=True,
        check=True,
    )


def make_reference_genomes_list_from_dir(reference_genomes_dir, reference_list_name):
    reference_list = Path(reference_list_name)
    reference_genomes_dir = Path(reference_genomes_dir)
    with open(reference_list, "w") as f:
        for genome_file in reference_genomes_dir.iterdir():
            f.write(str(genome_file.absolute()) + "\n")


def make_reference_genomes_list_from_file(
    accessions_file, reference_genomes_dir, reference_list_name
):
    accessions = open(accessions_file).readlines()
    accessions = list(acc.strip() for acc in accessions)
    reference_list = Path(reference_list_name)
    reference_genomes_dir = Path(reference_genomes_dir).absolute()

    with open(reference_list, "w") as f:
        for acc in accessions:
            genome_file = reference_genomes_dir / (acc + ".fa")
            if not genome_file.exists():
                raise ValueError(f"The file {genome_file} does not exist")
            f.write(str(genome_file) + "\n")


if __name__ == "__main__":
    query_genome = sys.argv[1]
    ref_genome_accessions = sys.argv[2]
    ref_genomes_dir = sys.argv[3]

    ref_list_name = "fastani_ref_genome_paths.txt"
    make_reference_genomes_list_from_file(ref_genome_accessions, ref_genomes_dir, ref_list_name)
    run_fastani(query_genome, ref_list_name, "fastani_results.tsv")

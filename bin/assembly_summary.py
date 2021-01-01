#!/usr/bin/env python3

# TODO: add support for plasmids

import json
import sys
from pathlib import Path

from alignment_stats import average_coverage, reads_mapped
from cds import cds
from checkm_stats import checkm_stats
from circlator_genome_topology import is_circular
from quast_stats import quast_stats


def make_assigner(results_dict, error_list):
    """Helper function for make_summary"""
    if not isinstance(results_dict, dict):
        raise ValueError("results_dict is needs to be a dict.")
    if not isinstance(error_list, list):
        raise ValueError("error_list is needs to be a list.")

    def try_assign(key, get_value_function, update=False):
        try:
            if update:
                results_dict.update(get_value_function())
            else:
                results_dict[key] = get_value_function()
        except Exception as e:
            error_list.append(f"[{key}] {str(e)}")

    return try_assign


def make_summary(
    short_reads_coverage_dir,
    long_reads_coverage_dir,
    quast_dir,
    prokka_dir,
    circlator_dir,
    checkm_dir,
):
    """
    Naming assumptions:
    - short_reads_coverage_dir contains stats.txt and bbmap_stderr.txt
    - long_reads_coverage_dir contains stats.txt and pileup_stderr.txt
    - prokka was used with the flag \"--prefix prokka\"

    Returns:
        A dictionary containing useful summary values about the genome assembly.
        Any errors encoutered are stored under they key \"error\".
    """
    results = {}
    error_list = []

    assign_to_dict = make_assigner(results, error_list)

    short_reads_coverage_dir = Path(short_reads_coverage_dir)
    long_reads_coverage_dir = Path(long_reads_coverage_dir)
    prokka_dir = Path(prokka_dir)

    assign_to_dict(
        "avg short reads coverage",
        lambda: average_coverage(short_reads_coverage_dir / "stats.txt"),
    )
    assign_to_dict(
        "short reads mapped",
        lambda: reads_mapped(short_reads_coverage_dir / "bbmap_stderr.txt"),
    )
    assign_to_dict(
        "avg long reads coverage",
        lambda: average_coverage(long_reads_coverage_dir / "stats.txt"),
    )
    assign_to_dict(
        "long reads mapped",
        lambda: reads_mapped(long_reads_coverage_dir / "pileup_stderr.txt"),
    )
    assign_to_dict("quast", lambda: quast_stats(quast_dir), update=True)
    assign_to_dict("CDS", lambda: cds(prokka_dir / "prokka.txt"))
    assign_to_dict("is circular", lambda: is_circular(circlator_dir))
    assign_to_dict("checkm", lambda: checkm_stats(checkm_dir), update=True)

    results["errors"] = error_list

    if error_list:
        print("Errors:\n", error_list, file=sys.stderr)

    return results


def make_summary_and_save(
    save_to,
    short_reads_coverage_dir,
    long_reads_coverage_dir,
    quast_dir,
    prokka_dir,
    circlator_dir,
    checkm_dir,
):
    """Saves summary statistsics as a pretty-printed JSON file.

    See make_summary for assumptions and details of dictiontary saved.

    Args:
        save_to: Path to file where summary will be saved in JSON format.

        The other args are the output directories for various parts of
        the assembly / evaluation pipeline.
    """
    summary = make_summary(
        short_reads_coverage_dir,
        long_reads_coverage_dir,
        quast_dir,
        prokka_dir,
        circlator_dir,
        checkm_dir,
    )

    json.dump(summary, open(save_to, "w"), indent=4)


# TODO: named inputs
if __name__ == "__main__":
    make_summary_and_save(
        sys.argv[1],
        sys.argv[2],
        sys.argv[3],
        sys.argv[4],
        sys.argv[5],
        sys.argv[6],
        sys.argv[7],
    )

#!/usr/bin/env python3

# TODO: add support for plasmids

import argparse
import json
import sys
from pathlib import Path

import platon_stats
from alignment_stats import average_coverage, reads_mapped
from cds import cds
from checkm_stats import checkm_stats
from commons import make_flag
from contig_names import get_contig_names
from quast_stats import quast_stats

### CONSTANTS ###

KEY_AVERAGE_SHORT_READS_COVERAGE = "avg short reads coverage"
KEY_AVERAGE_LONG_READS_COVERAGE = "avg long reads coverage"
KEY_SHORT_READS_MAPPED = "short reads mapped"
KEY_LONG_READS_MAPPED = "long reads mapped"
KEY_CDS = "CDS"
KEY_CIRCULARITY = "circularity"

### FUNCTIONS ###


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


def filtered_circularity_sumary(circularity_summary, contigs):
    circ_summary = json.load(open(circularity_summary))
    if isinstance(circ_summary, str):
        return circ_summary
    elif isinstance(circ_summary, dict):
        result = {
            contig_name: circularity
            for contig_name, circularity in circ_summary.items()
            if contig_name in contigs
        }
        return result
    else:
        raise ValueError(
            f"Unrecognised type of circularity summary:\nType:{type(circ_summary)}\n"
            + f"Value:{circ_summary}"
        )


def make_chromosome_summary(
    short_reads_coverage_dir,
    long_reads_coverage_dir,
    quast_dir,
    prokka_txt,
    circularity_summary,
    assembly,
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

    assign_to_dict(
        KEY_AVERAGE_SHORT_READS_COVERAGE,
        lambda: average_coverage(short_reads_coverage_dir / "stats.txt"),
    )
    assign_to_dict(
        KEY_SHORT_READS_MAPPED,
        lambda: reads_mapped(short_reads_coverage_dir / "bbmap_stderr.txt"),
    )
    assign_to_dict(
        KEY_AVERAGE_LONG_READS_COVERAGE,
        lambda: average_coverage(long_reads_coverage_dir / "stats.txt"),
    )
    assign_to_dict(
        KEY_LONG_READS_MAPPED,
        lambda: reads_mapped(long_reads_coverage_dir / "pileup_stderr.txt"),
    )
    assign_to_dict("quast", lambda: quast_stats(quast_dir), update=True)
    assign_to_dict(KEY_CDS, lambda: cds(prokka_txt))
    assign_to_dict(
        KEY_CIRCULARITY,
        lambda: filtered_circularity_sumary(
            circularity_summary, contigs=get_contig_names(assembly)
        ),
    )
    assign_to_dict("checkm", lambda: checkm_stats(checkm_dir), update=True)

    results["errors"] = error_list

    if error_list:
        print("Errors:\n", error_list, file=sys.stderr)

    return results


def make_plasmid_summary(
    short_reads_coverage_dir,
    long_reads_coverage_dir,
    quast_dir,
    prokka_txt,
    platon_tsv,
):
    results = {}
    error_list = []

    assign_to_dict = make_assigner(results, error_list)

    short_reads_coverage_dir = Path(short_reads_coverage_dir)
    long_reads_coverage_dir = Path(long_reads_coverage_dir)

    assign_to_dict(
        KEY_AVERAGE_SHORT_READS_COVERAGE,
        lambda: average_coverage(short_reads_coverage_dir / "stats.txt"),
    )
    assign_to_dict(
        KEY_SHORT_READS_MAPPED,
        lambda: reads_mapped(short_reads_coverage_dir / "bbmap_stderr.txt"),
    )
    assign_to_dict(
        KEY_AVERAGE_LONG_READS_COVERAGE,
        lambda: average_coverage(long_reads_coverage_dir / "stats.txt"),
    )
    assign_to_dict(
        KEY_LONG_READS_MAPPED,
        lambda: reads_mapped(long_reads_coverage_dir / "pileup_stderr.txt"),
    )
    assign_to_dict("quast", lambda: quast_stats(quast_dir), update=True)
    assign_to_dict(KEY_CDS, lambda: cds(prokka_txt))
    assign_to_dict(KEY_CIRCULARITY, lambda: platon_stats.cicularity(platon_tsv))

    results["errors"] = error_list

    if error_list:
        print("Errors:\n", error_list, file=sys.stderr)

    return results


### PARSER ###

# constants
FLAG_SHORT_READS_COV_DIR = "short"
FLAG_LONG_READS_COV_DIR = "long"
FLAG_QUAST_DIR = "quast"
FLAG_PROKKA_TXT = "prokka"
FLAG_OUT = "out"


def make_base_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(make_flag(FLAG_SHORT_READS_COV_DIR), required=True)
    parser.add_argument(make_flag(FLAG_LONG_READS_COV_DIR), required=True)
    parser.add_argument(make_flag(FLAG_QUAST_DIR), required=True)
    parser.add_argument(make_flag(FLAG_PROKKA_TXT), required=True)
    parser.add_argument(make_flag(FLAG_OUT), required=True)
    return parser

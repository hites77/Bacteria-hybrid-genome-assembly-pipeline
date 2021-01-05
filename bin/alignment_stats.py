#!/usr/bin/env python3

"""Gets coverage related statistics from the output of bbmap or pileup.

Has been tested against BBTools v38.85 and 38.87.
"""

from pathlib import Path

import pandas as pd


def line_starting_with(stderr_lines, target):
    """Finds the line starting with the target string.

    Args:
        stderr_lines: List of strings, representing the bbmap/pileup stderr split
          by \n. This can obtained by using the readlines() method.
        target: String that the line starts with.

    Returns:
        A string from stderr_lines which starts with target.

    Raises:
        ValueError: There are 0 or more than 1 lines starting with target.
    """
    filtered_lines = list(filter(lambda text: text.startswith(target), stderr_lines))
    if len(filtered_lines) != 1:
        raise ValueError(
            f'Invalid stderr file: number of lines starting with "{target}" != 1\nFiltered lines: {filtered_lines}'
        )
    return filtered_lines[0]


def line_value(line):
    """Gets the value from a line with the format \"Name:  \\tValue\\n\"

    Note that:
    - The name may contain whitespace.
    - The value must not have any whitespace.
    - The amount of whitespace between the name and value can vary.

    Returns:
        The value, as a string.

    Raises:
        ValueError: The line is not in the correct format.
    """
    split = line.split()
    if len(split) < 2:
        raise ValueError(f'"{line}" is an invalid line: number of splits < 2\nSplits: {split}')
    return split[-1]


# TODO update docs
def reads_mapped_internal(stderr_lines):
    """
    Args:
        stderr_lines: List of strings, representing the bbmap/pileup stderr split
          by \n. This can obtained by using the readlines() method.

    Returns:
        A dictionary containing the total number of reads (key: total reads),
        number of mapped reads (key: mapped reads), and proportion of mapped reads
        (key: proportion mapped).
    """

    def get_int(starts_with_text):
        return int(line_value(line_starting_with(stderr_lines, starts_with_text)))

    total_reads = get_int("Reads:")
    reads_mapped = get_int("Mapped reads:")
    return {
        "total reads": total_reads,
        "mapped reads:": reads_mapped,
        "proportion mapped": reads_mapped / total_reads,
    }


def average_coverage_internal(stats_txt_df):
    """
    Args:
        stats_txt_df: pd.read_table('stats.txt')

    Returns:
        The average fold coverage as a float.

    Raises:
        ValueError: The dataframe does not have exactly 1 row.
    """
    if len(stats_txt_df) != 1:
        raise ValueError(f"Dataframe should have 1 row.\nNumber of rows: {len(stats_txt_df)}")
    return float(stats_txt_df["Avg_fold"][0])


def reads_mapped(stderr_path):
    """
    Returns:
        A dictionary containing the total number of reads (key: total reads)
        and number of mapped reads (key: mapped reads).
    """
    return reads_mapped_internal(open(Path(stderr_path)).readlines())


def average_coverage(stats_txt_path):
    """
    Returns:
        The average fold coverage as a float.
    """
    return average_coverage_internal(pd.read_table(stats_txt_path))


# print(
#     reads_mapped(
#         "/home/chloe/Documents/NUS/UROPS/server-data/S8E_3_1_replicate4/assembly_eval/short_reads_coverage/bbmap_stderr.txt"
#     )
# )

# print(
#     average_coverage(
#         "/home/chloe/Documents/NUS/UROPS/server-data/S8E_3_1/assembly_eval/short_read_coverage/stats.txt"
#     )
# )

#!/usr/bin/env python3

from pathlib import Path


def line_starting_with(prokka_txt_lines, starting_with_text):
    """Finds the line starting with the target string.

    Args:
        prokka_txt_lines: List of strings, representing the contents of
          prokka.txt split by \n. This can obtained by using the readlines()
          method.
        target: String that the line starts with.

    Returns:
        A string from prokka_txt_lines which starts with target.

    Raises:
        ValueError: There are 0 or more than 1 lines starting with target.
    """
    filtered_lines = list(
        filter(lambda text: text.startswith(starting_with_text), prokka_txt_lines)
    )
    if len(filtered_lines) != 1:
        raise ValueError(
            f'Invalid prokka.txt: number of lines starting with "{starting_with_text}" != 1'
            + f"\nFiltered lines: {filtered_lines}"
        )
    return filtered_lines[0]


def line_value(line):
    """Gets the value from a line with the format \"Name: Value\\n\"

    Note that:
    - The name may contain whitespace.
    - The value must not have any whitespace.

    Returns:
        The value, as a string.

    Raises:
        ValueError: The line is not in the correct format.
    """
    split = line.split()
    if len(split) < 2:
        raise ValueError(
            f'"{line}" is an invalid line: number of splits < 2\nSplits: {split}'
        )
    return split[-1]


def cds_internal(prokka_txt_lines):
    """
    Args:
        prokka_txt_lines: List of strings, representing the contents of
          prokka.txt split by \n. This can obtained by using the readlines()
          method.

    Returns:
        Length of coding sequence (CDS) in bp.
    """
    return int(line_value(line_starting_with(prokka_txt_lines, "CDS:")))


def cds(prokka_txt_path):
    """
    Returns:
        Length of coding sequence (CDS) in bp.
    """
    return cds_internal(open(Path(prokka_txt_path)).readlines())


# print(cds('/home/chloe/Documents/NUS/UROPS/server-data/S8E_3_1/assembly_eval/prokka_annotation/PROKKA_12222020.txt'))

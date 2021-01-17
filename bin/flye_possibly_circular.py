#!/usr/bin/env python3

"""Checks whether an assembly might be circular (and hence whether it is worth running circlator).

Prints yes to stdout if there are any circular contigs, or if there are multiple linear contigs;
otherwise, no is printed.

Sample usage:
$ ./flye_possibly_circular.py path/to/flye_directory
yes


$ ./flye_possibly_circular.py path/to/flye_directory/assembly_info.txt
yes
"""

import sys

from flye_utils import TXT_COL_CIRCULAR, is_yes, read_flye_info_txt
from commons import bool_to_str


def assembly_possibly_circular(flye_dir_or_info_txt):
    df = read_flye_info_txt(flye_dir_or_info_txt)
    return any(is_yes(is_circular_str) for is_circular_str in df[TXT_COL_CIRCULAR]) or df.nrows > 1


if __name__ == "__main__":
    print(bool_to_str(assembly_possibly_circular(sys.argv[1])))

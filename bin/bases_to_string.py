#!/usr/bin/env python3

"""Converts number of bases to giga/mega/kilo bases, rounded to 1 dp, and prints it to stdout.
Eg. 5.1g, 5.1m, 5.1k.

Command line usage:
$ ./bases_to_string.py 5109239
5.1m

$ echo 5109239 | ./bases_to_string.py -
5.1m
"""

import sys


def divide(number, divide_by):
    """Returns number / divide_by, rounded to 1 dp"""
    return round(number / divide_by, 1)


def bases_to_string(num_bases):
    if not isinstance(num_bases, int) or num_bases < 0:
        raise ValueError(f"num_bases must be a non-negative integer")

    for divide_by, suffix in ((10 ** 9, "g"), (10 ** 6, "m")):
        if num_bases >= divide_by:
            return str(divide(num_bases, divide_by)) + suffix
    return str(divide(num_bases, 10 ** 3)) + "k"


if __name__ == "__main__":
    arg = sys.argv[1]
    if arg == "-":
        num_bases = sys.stdin.read()
    else:
        num_bases = arg
    num_bases = int(num_bases)
    print(bases_to_string(num_bases))

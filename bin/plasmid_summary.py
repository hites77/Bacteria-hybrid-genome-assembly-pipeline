#!/usr/bin/env python3

import argparse
import json

from assembly_summary import make_plasmid_summary
from commons import make_flag

FLAG_SHORT_READS_COV_DIR = "short"
FLAG_LONG_READS_COV_DIR = "long"
FLAG_QUAST_DIR = "quast"
FLAG_PROKKA_TXT = "prokka"
FLAG_PLATON_TSV = "platon"
FLAG_OUT = "out"


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(make_flag(FLAG_SHORT_READS_COV_DIR), required=True)
    parser.add_argument(make_flag(FLAG_LONG_READS_COV_DIR), required=True)
    parser.add_argument(make_flag(FLAG_QUAST_DIR), required=True)
    parser.add_argument(make_flag(FLAG_PROKKA_TXT), required=True)
    parser.add_argument(make_flag(FLAG_PLATON_TSV), required=True)
    parser.add_argument(make_flag(FLAG_OUT), required=True)
    return parser


if __name__ == "__main__":
    parser = make_parser()
    args = vars(parser.parse_args())
    summary = make_plasmid_summary(
        short_reads_coverage_dir=args[FLAG_SHORT_READS_COV_DIR],
        long_reads_coverage_dir=args[FLAG_LONG_READS_COV_DIR],
        quast_dir=args[FLAG_QUAST_DIR],
        prokka_txt=args[FLAG_PROKKA_TXT],
        platon_tsv=args[FLAG_PLATON_TSV],
    )
    save_to = args[FLAG_OUT]
    json.dump(summary, open(save_to, "w"), indent=4)

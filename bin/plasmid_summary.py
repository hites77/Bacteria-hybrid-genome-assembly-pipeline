#!/usr/bin/env python3

import json

from assembly_summary import (
    FLAG_ASSEMBLY,
    FLAG_LONG_READS_COV_DIR,
    FLAG_OUT,
    FLAG_PROKKA_TXT,
    FLAG_QUAST_DIR,
    FLAG_SHORT_READS_COV_DIR,
    make_base_parser,
    make_plasmid_summary,
)
from commons import make_flag

FLAG_PLATON_TSV = "platon"


def make_parser():
    parser = make_base_parser()
    parser.add_argument(make_flag(FLAG_PLATON_TSV), required=True)
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

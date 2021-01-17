#!/usr/bin/env python3

import json
import sys

from flye_utils import (
    TXT_COL_CIRCULAR,
    TXT_COL_CONTIG_NAME,
    circularity_col_to_string,
    read_flye_info_txt,
)


def make_circularity_summary(flye_dir_or_info_txt):
    df = read_flye_info_txt(flye_dir_or_info_txt)
    circularity_dict = {}
    for contig, is_circular in zip(df[TXT_COL_CONTIG_NAME], df[TXT_COL_CIRCULAR]):
        circularity_dict[contig] = circularity_col_to_string(is_circular)
    return circularity_dict


if __name__ == "__main__":
    summary = make_circularity_summary(sys.argv[1])
    print(json.dumps(summary, indent=4))

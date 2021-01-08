#!/usr/bin/env python3

import sys
from pathlib import Path

import pandas as pd
from commons import bool_to_str

IN_YES = "Y"
IN_NO = "N"

COL_CIRCULAR = "circ."


def any_circular_contigs(flye_info_txt_or_dir):
    flye_info_txt_or_dir = Path(flye_info_txt_or_dir)
    flye_txt = (
        flye_info_txt_or_dir / "assembly_info.txt"
        if flye_info_txt_or_dir.is_dir()
        else flye_info_txt_or_dir
    )
    df = pd.read_table(flye_txt)
    return any(is_circular == IN_YES for is_circular in df[COL_CIRCULAR])


if __name__ == "__main__":
    any_circ_contigs = any_circular_contigs(sys.argv[1])
    print(bool_to_str(any_circ_contigs))

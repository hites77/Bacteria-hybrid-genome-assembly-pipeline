#!/usr/bin/env python3

from pathlib import Path

import pandas as pd

# value of yes/no columns in assembly_info.txt
TXT_YES = "Y"
TXT_NO = "N"

# names of columsn in assembly_info.txt
TXT_COL_CONTIG_NAME = "#seq_name"
TXT_COL_CIRCULAR = "circ."


def read_flye_info_txt(flye_dir_or_info_txt):
    flye_info_txt_or_dir = Path(flye_dir_or_info_txt)
    flye_txt = (
        flye_info_txt_or_dir / "assembly_info.txt"
        if flye_info_txt_or_dir.is_dir()
        else flye_info_txt_or_dir
    )
    df = pd.read_table(flye_txt)
    return df


def is_yes(yes_or_no_str):
    if yes_or_no_str == TXT_YES:
        return True
    elif yes_or_no_str == TXT_NO:
        return False
    else:
        raise ValueError(f"String should either be '{TXT_YES}' or '{TXT_NO}'")


def circularity_col_to_string(yes_or_no_str):
    return "circular" if is_yes(yes_or_no_str) else "linear"

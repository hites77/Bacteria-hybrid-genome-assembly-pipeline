#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

### CONSTANTS ###

COL_CONTIG_ID = "ID"
COL_CIRCULARITY = "Circular"

### FUNCTIONS ###


def get_df(platon_tsv):

    df = pd.read_table(platon_tsv)

    def map_circular_string(text):
        if text == "yes":
            return "circular"
        elif text == "no":
            return "linear"
        else:
            raise ValueError(f"Unrecognised circularity type: {text}")

    df[COL_CIRCULARITY] = df[COL_CIRCULARITY].apply(map_circular_string)

    return df


def circularity(platon_dir_or_tsv):
    platon_tsv = (
        platon_dir_or_tsv
        if str(platon_dir_or_tsv).endswith(".tsv")
        else Path(platon_dir_or_tsv) / "platon.tsv"
    )

    df = get_df(platon_tsv)

    result_dict = {}  # contig name : circularity
    for row in df.itertuples():
        row = row._asdict()
        contig = str(row[COL_CONTIG_ID])
        circularity = row[COL_CIRCULARITY]
        result_dict[contig] = circularity

    return result_dict


def test_circularity():
    print(circularity("test_platon_stats.tsv"))

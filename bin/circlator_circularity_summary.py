#!/usr/bin/env python3

import sys
from pathlib import Path

import pandas as pd
import json

col_circularised = "circularised"
col_pandas_index = "Index"
col_contig_name = "#Contig"


def get_merge_circularise_log_df(merge_circ_log_path):
    df = pd.read_table(merge_circ_log_path)

    df = df.set_index(col_contig_name)

    def map_circularised_values(text):
        if text == 1:
            return "circular"
        elif text == 0:
            return "linear"
        else:
            raise ValueError(f"Unrecognised circularity type: {text}")

    df[col_circularised] = df[col_circularised].apply(map_circularised_values)

    return df


def is_circular_internal(df, contigs=None):
    """
    Args:
        df: Dataframe from get_merge_circularise_log_df.
        contigs: Optional. List of names of contigs we are interested in. If no
          argument is passed, use all contigs.

    Returns:
        A dictionary mapping the contig name to its circularity
        Eg. { 'contig_1' : 'circular' }
    """

    if contigs != None and not isinstance(contigs, list):
        raise ValueError("contigs must be a list or None.")

    if contigs != None:
        df = df.loc[contigs]

    result_dict = {}  # contig name : is circular
    for row in df.itertuples():
        row = row._asdict()
        contig = str(row[col_pandas_index])
        circularised = row[col_circularised]
        result_dict[contig] = circularised
    return result_dict


def circularity(circlator_output_dir_or_log_file, contigs=None):
    """
    Args:
        circlator_output_dir: Path to circlator's output directory.
        contigs: Optional. List of names of contigs we are interested in. If no
          argument is passed, use all contigs.

    Returns:
        A dictionary mapping the contig name to its circularity
        Eg. { 'contig_1' : 'circular' }
    """
    merge_circ_log_file = (
        circlator_output_dir_or_log_file
        if str(circlator_output_dir_or_log_file).endswith(".log")
        else Path(circlator_output_dir_or_log_file) / "04.merge.circularise.log"
    )
    return is_circular_internal(get_merge_circularise_log_df(merge_circ_log_file), contigs)


if __name__ == "__main__":
    summary = circularity(sys.argv[1])
    # print to stdout as JSON
    print(json.dumps(summary, indent=4))


# print(
#     is_circular_internal(
#         get_merge_circularise_log_df(
#             "~/Documents/NUS/UROPS/server-data/S8E_3_1_replicate4/assembly/circlator/04.merge.circularise.log",
#         ),
#         ["contig_1"],
#     )
# )

# print(
#     is_circular(
#         "~/Documents/NUS/UROPS/server-data/S8E_3_1_replicate4/assembly/circlator/"
#     )
# )

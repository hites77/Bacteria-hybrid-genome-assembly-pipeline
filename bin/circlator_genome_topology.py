#!/usr/bin/env python3

from pathlib import Path

import pandas as pd

col_circularised = "circularised"
col_pandas_index = "Index"


def get_merge_circularise_log_df(merge_circ_log_path):
    df = pd.read_table(merge_circ_log_path)

    col_contig_name = "#Contig"
    df = df.set_index(col_contig_name)

    # validate circularised column
    valid_circularised_values = (0, 1)
    if any(val not in valid_circularised_values for val in df[col_circularised]):
        raise ValueError(f"{merge_circ_log_path} is an invalid merging circularise.log file")

    df[col_circularised] = df[col_circularised].apply(bool)

    return df


def is_circular_internal(df, contigs=None):
    """
    Args:
        df: Dataframe from get_merge_circularise_log_df.
        contigs: Optional. List of names of contigs we are interested in. If no
          argument is passed, use all contigs.

    Returns:
        A dictionary mapping the contig name to whether it is circular.
        Eg. { 'contig_1' : True }

    Raises:
        ValueError: contigs argument is neither None nor a list.
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


def is_circular(circlator_output_dir, contigs=None):
    """
    Args:
        circlator_output_dir: Path to circlator's output directory.
        contigs: Optional. List of names of contigs we are interested in. If no
          argument is passed, use all contigs.

    Returns:
        A dictionary mapping the contig name to whether it is circular.
        Eg. { 'contig_1' : True }
    """
    merge_circ_log_file = Path(circlator_output_dir) / "04.merge.circularise.log"
    return is_circular_internal(get_merge_circularise_log_df(merge_circ_log_file))


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

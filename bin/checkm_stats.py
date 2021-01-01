#!/usr/bin/env python3

import json
from pathlib import Path


def checkm_stats_internal(bin_stats_ext_text):
    """
    Args:
        bin_stats_ext_text: open('storage/bin_stats_ext.tsv').read()

    Returns:
        A dictionary containing the completeness (key: completeness) and
        contamination (key: contamination), both over 1.
    """
    json_text = bin_stats_ext_text[bin_stats_ext_text.find("\t") + 1 :].replace(
        "'", '"'
    )
    checkm_dict = json.loads(json_text)
    return {
        "completeness": float(checkm_dict["Completeness"] / 100),
        "contamination": float(checkm_dict["Contamination"] / 100),
    }


def checkm_stats(checkm_output_dir):
    """
    Returns:
        A dictionary containing the completeness (key: completeness) and
        contamination (key: contamination), both over 1.
    """
    bin_stats_ext_file = Path(checkm_output_dir) / "storage" / "bin_stats_ext.tsv"
    return checkm_stats_internal(open(bin_stats_ext_file).read())


# print(
#     checkm_stats(
#         "/home/chloe/Documents/NUS/UROPS/server-data/S8E_3_1/assembly_eval/completeness_checkm/"
#     )
# )

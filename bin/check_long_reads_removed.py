#!/usr/bin/env python3

"""Check whether Filtlong has removed reads above a certain length.

The number of reads greater than or equal to the threshold length which were removed is
printed to stdout.
If the `--tsv <file>` argument is given, the IDs and lengths of any reads >= the threshold
removed are saved to `<file>`. Note that the tsv file is not created if no such reads were
removed.

Sample usage:

To check how many reads above 10kb in length were removed by Filtlong:
$ check_long_reads_removed.py --old raw_pacbio.fq.gz --new pacbio.fq --threshold 10k \
$    --tsv very_long_reads_removed.tsv
5
"""

import argparse
import re
import subprocess

from Bio import SeqIO
from commons import make_flag, string_to_bases


def open_fastq(fastq_path: str):
    if fastq_path.endswith(".gz"):
        compressed_fastq = fastq_path
        fastq_path = fastq_path[:-3]
        subprocess.run(f"gzip -dc {compressed_fastq} > {fastq_path}", shell=True, check=True)
    return SeqIO.parse(fastq_path, "fastq")


def check(fastq1_path, fastq2_path, threshold, tsv_file=None):

    # example of read id after filtlong: m54214_200904_195840/38470059/13341_16136_2-2749
    # 2-2749 denotes range that was taken
    FILTLONG_APPENDED_PATTERN = re.compile(r"_\d+-\d+$")

    fastq1 = open_fastq(fastq1_path)
    fastq2 = open_fastq(fastq2_path)

    fastq2_ids = set()
    for rec in fastq2:
        fastq2_ids.add(rec.id)
        fastq2_ids.add(re.sub(FILTLONG_APPENDED_PATTERN, "", rec.id))

    # list of ids of reads >= THRESHOLD which were removed by filtlong
    # note: cannot remove ids after finding them, in case filtlong splits a read into multiple reads
    removed_long_reads = [
        (rec.id, len(rec)) for rec in fastq1 if len(rec) >= threshold and rec.id not in fastq2_ids
    ]

    print(f"number of reads >= {threshold} bases removed: {len(removed_long_reads)}")

    if removed_long_reads and tsv_file:
        print(f"saving to {tsv_file}")
        with open(tsv_file, "w") as f:
            f.write("original read id\tlength")
            for id, length in removed_long_reads:
                f.write(f"\n{id}\t{length}")


#### PARSER ####

OLD_FASTQ_FLAG = "old"
NEW_FASTQ_FLAG = "new"
THRESHOLD_FLAG = "threshold"
TSV_FILE_FLAG = "tsv"


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(make_flag(OLD_FASTQ_FLAG), required=True)
    parser.add_argument(make_flag(NEW_FASTQ_FLAG), required=True)
    parser.add_argument(make_flag(THRESHOLD_FLAG), required=True)
    parser.add_argument(make_flag(TSV_FILE_FLAG), required=False, default=None)
    return parser


#### MAIN ####

if __name__ == "__main__":
    parser = make_parser()
    args = vars(parser.parse_args())
    check(
        fastq1_path=args[OLD_FASTQ_FLAG],
        fastq2_path=args[NEW_FASTQ_FLAG],
        threshold=string_to_bases(args[THRESHOLD_FLAG]),
        tsv_file=args[TSV_FILE_FLAG],
    )


# check(
#     "/home/chloe/Documents/NUS/UROPS/server-data/S8E_3_1/reads/raw/pacbio.fq",
#     "/home/chloe/Documents/NUS/UROPS/server-data/S8E_3_1/reads/cleaned/pacbio.fq",
# )

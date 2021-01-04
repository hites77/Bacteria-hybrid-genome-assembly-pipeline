#!/usr/bin/env python3

"""
Finds the best trimq value to use with bbduk such that at least X% of the reads are kept.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

### CONSTANTS ###

TRIMQ_USED_FILE = "trimq_used.txt"

### FUNCTIONS ###

MIN_PHRED_SCORE = 1


class CannotKeepPercentError(Exception):
    def __init__(self, keep_percent, other_message=""):
        super().__init__(f"Cannot keep at least {keep_percent}% of reads.{other_message}")


RE_GROUP_PERCENT_READS = "percentKept"
RE_PERCENT_NUM = r"\d{1,3}(.\d+)?"
RE_KEPT_LINE = re.compile(
    r"^Result:\s+\t\d+ reads \("
    + f"(?P<{RE_GROUP_PERCENT_READS}>{RE_PERCENT_NUM})"
    + r"%\)"
    + r"\s*\t\d+ bases \("
    + RE_PERCENT_NUM
    + r"%\)$"
)


def percentage_kept(log_file):
    """Gets percentage of reads kept."""

    for line in open(Path(log_file)):
        match = RE_KEPT_LINE.match(line)
        if match:
            return float(match.group(RE_GROUP_PERCENT_READS))


def run_bbduk(in1, in2, out1, out2, log_file, trimq, args=""):
    subprocess.run(
        f"bbduk.sh -Xmx1g in1={in1} in2={in2} out1={out1} out2={out2} trimq={trimq} {args}"
        + f" 2> {log_file}",
        shell=True,
        check=True,
    )


def main(in1, in2, out1, out2, keep_percent, start_trimq, min_trimq, bbduk_args, info_dir):
    if start_trimq < min_trimq:
        raise ValueError("start_trimq should be >= min_trimq")

    info_dir = Path(info_dir)
    Path(info_dir).mkdir(parents=True, exist_ok=True)

    temp_dir = Path("temp")
    final_reads = None

    temp_dir.mkdir()

    for trimq in range(start_trimq, min_trimq - 1, -1):
        tmp_out_1 = temp_dir / f"trimq_{trimq}_1.fq"
        tmp_out_2 = temp_dir / f"trimq_{trimq}_2.fq"
        log_file = temp_dir / f"trimq_{trimq}.log"
        run_bbduk(in1, in2, tmp_out_1, tmp_out_2, log_file, trimq, bbduk_args)

        if percentage_kept(log_file) >= keep_percent:
            final_reads = tmp_out_1, tmp_out_2
            break

    if final_reads == None:
        raise CannotKeepPercentError(keep_percent)

    # record trimq used
    # trimq should be bound: as long as final_reads exists,
    # the loop should have run
    open(info_dir / TRIMQ_USED_FILE, "w").write(str(trimq))

    # move to final reads to out1, out2
    out1 = Path(out1)
    out2 = Path(out2)
    out1.parent.mkdir(parents=True, exist_ok=True)
    out2.parent.mkdir(parents=True, exist_ok=True)
    Path(final_reads[0]).rename(out1)
    Path(final_reads[1]).rename(out2)

    # remove temp files
    subprocess.run(f"rm {temp_dir}/trimq_*_1.fq", shell=True)
    subprocess.run(f"rm {temp_dir}/trimq_*_2.fq", shell=True)
    subprocess.run(f"rm {temp_dir}/trimq_*.log", shell=True)
    try:
        temp_dir.rmdir()
    except Exception as e:
        print(e, file=sys.stderr)


### PARSER ###

# parser args
FLAG_IN1 = "in1"
FLAG_IN2 = "in2"
FLAG_OUT1 = "out1"
FLAG_OUT2 = "out2"
FLAG_KEEP_PERCENT = "keep_percent"
FLAG_START_TRIMQ = "start_trimq"
FLAG_MIN_TRIMQ = "min_trimq"
FLAG_BBDUK_ARGS = "args"
FLAG_INFO_DIR = "infodir"


def make_flag(name):
    return "--" + name


def make_parser():
    parser = argparse.ArgumentParser(
        description="Finds the best trimq value to use with bbduk such that at least X% of the reads are kept."
    )

    # required args
    for inout in (FLAG_IN1, FLAG_IN2, FLAG_OUT1, FLAG_OUT2):
        parser.add_argument(
            make_flag(inout),
            nargs=1,
            required=True,
            help="Corresponds to bbduk's input/output.",
        )
    parser.add_argument(
        make_flag(FLAG_KEEP_PERCENT),
        nargs=1,
        type=float,
        required=True,
        help="Keep at least X% of reads.",
    )
    parser.add_argument(
        make_flag(FLAG_START_TRIMQ),
        nargs=1,
        type=int,
        required=True,
        help="The trimq value to start from.",
    )

    # optional args
    parser.add_argument(
        make_flag(FLAG_MIN_TRIMQ),
        nargs=1,
        type=int,
        required=False,
        default=[MIN_PHRED_SCORE],
        help="Do not use a trimq value below this.",
    )
    parser.add_argument(
        make_flag(FLAG_BBDUK_ARGS),
        nargs="*",
        required=False,
        default=[],
        help="Additional args to pass to bbduk.",
    )
    parser.add_argument(
        make_flag(FLAG_INFO_DIR),
        nargs=1,
        required=False,
        default=[""],
        help=f"Directory to save {TRIMQ_USED_FILE} in.",
    )
    return parser


# TODO fix help message

### MAIN ###

if __name__ == "__main__":
    # if False:
    parser = make_parser()
    args = vars(parser.parse_args())
    main(
        in1=args[FLAG_IN1][0],
        in2=args[FLAG_IN2][0],
        out1=args[FLAG_OUT1][0],
        out2=args[FLAG_OUT2][0],
        keep_percent=args[FLAG_KEEP_PERCENT][0],
        start_trimq=args[FLAG_START_TRIMQ][0],
        min_trimq=args[FLAG_MIN_TRIMQ][0],
        bbduk_args=" ".join(args[FLAG_BBDUK_ARGS]),
        info_dir=args[FLAG_INFO_DIR][0],
    )

### TESTS ###


def test_percentage_kept():
    print(percentage_kept("bbduk_test.txt"))


def test_parser():
    parser = make_parser()
    a = parser.parse_args(
        "--in1 IN1 --in2 IN2 --out1 OUT1 --out2 OUT2 --keep_percent 80 --start_trimq 40 --args bbduk args here".split()
    )
    print(a)


# test_parser()

#!/usr/bin/env python3

import sys

if not (sys.version_info.major == 3 and sys.version_info.minor >= 6):
    raise Exception("Python version should be at least 3.6.")

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

pd.DataFrame({"column": [1, 2, 3]})

SeqIO.write([SeqRecord(Seq("ATCG"), id="seq", description="")], "test-biopython.fasta", "fasta")

#!/usr/bin/env python3

import sys
from pathlib import Path

import requests
from Bio import SeqIO

TEMP_DIR = Path("tmp")
N = 30


def download_accessions_helper(accessions, dest_dir):
    if not isinstance(accessions, list) or accessions == []:
        raise ValueError("accesssions must be a nonempty list")

    dest_dir = Path(dest_dir)

    if not dest_dir.exists():
        raise ValueError(f"dest_dir {dest_dir} does not exist")

    url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        + f"db=sequences&id={','.join(accessions)}&retmode=text&rettype=fasta"
    )

    temp_file = TEMP_DIR / f"temp_{accessions[0]}.fa"

    r = requests.get(url)
    with open(temp_file, "wb") as f:
        f.write(r.content)

    records = SeqIO.parse(temp_file, "fasta")

    for rec in records:
        filename = dest_dir / (rec.id.split()[0] + ".fa")
        SeqIO.write([rec], filename, "fasta")

    Path(temp_file).unlink()


def download_accessions(accessions_file, dest_dir):
    # parse reference accessions file
    accessions = open(accessions_file)
    accessions = list(acc.strip() for acc in accessions)

    Path(dest_dir).mkdir(exist_ok=True, parents=True)

    TEMP_DIR.mkdir(exist_ok=True, parents=True)
    for i in range(0, len(accessions), N):
        accs = accessions[i : i + N]
        download_accessions_helper(accs, dest_dir)
    TEMP_DIR.rmdir()


if __name__ == "__main__":
    download_accessions(sys.argv[1], sys.argv[2])

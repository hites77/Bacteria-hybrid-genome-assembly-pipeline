#!/usr/bin/env bash
set -euo pipefail

assemblyFasta1=$1
assemblyFasta2=$2

diff -U 0 $assemblyFasta1 $assemblyFasta2 | grep -v ^@ | wc -l

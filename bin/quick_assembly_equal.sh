#!/usr/bin/env bash

# may give false positives

set -euo pipefail

assemblyFasta1=$1
assemblyFasta2=$2

diff -U 0 <(tail -n +2 $assemblyFasta1) <(tail -n +2 $assemblyFasta2) | grep -v ^@ | wc -l

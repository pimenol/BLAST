#!/bin/bash
# Simplified nucleotide BLAST
# Usage: blast.sh -e <matrix.csv> -d <db1.fasta> [db2.fasta ...] -k <word_len> -t <threshold> [-p <gap_open>] [-x <gap_ext>] [-q <query.fasta>]
DIR="$(cd "$(dirname "$0")" && pwd)"
python3 "$DIR/blast.py" "$@"

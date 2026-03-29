"""Substitution matrix loading."""

import csv


def load_score_matrix(csv_path: str):
    """Load a substitution matrix from CSV. Returns a 128x128 list-of-lists
    where matrix[ord(a)][ord(b)] = score."""
    matrix = [[0] * 128 for _ in range(128)]
    with open(csv_path) as f:
        reader = csv.reader(f)
        header = next(reader)
        cols = [c.strip().upper() for c in header[1:]]
        for row in reader:
            row_char = row[0].strip().upper()
            for j, val in enumerate(row[1:]):
                col_char = cols[j]
                score = int(val.strip())
                for rc in (row_char, row_char.lower()):
                    for cc in (col_char, col_char.lower()):
                        matrix[ord(rc)][ord(cc)] = score
    return matrix

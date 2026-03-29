"""Database loading from FASTA files."""

import os
from collections import namedtuple
from Bio import SeqIO

DatabaseSequence = namedtuple('DatabaseSequence', ['filename', 'name', 'seq'])


def load_database(fasta_paths: list) -> list:
    """Load FASTA files into memory as bytes objects."""
    db = []
    for path in fasta_paths:
        fname = os.path.basename(path)
        for record in SeqIO.parse(path, "fasta"):
            db.append(DatabaseSequence(fname, record.id,
                                       str(record.seq).upper().encode('ascii')))
    return db

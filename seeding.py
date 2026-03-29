"""Seed generation: build query k-mer index and find seeds in database."""

from collections import defaultdict

# Valid amino acid characters
_AMINO_ACIDS = set(b'ACDEFGHIKLMNPQRSTVWY')

# Valid nucleotide characters
_NUCLEOTIDES = set(b'ACGT')


def _is_valid_kmer(kmer: bytes, valid_chars: set) -> bool:
    """Check if all characters in k-mer are valid."""
    return all(c in valid_chars for c in kmer)


def detect_alphabet(seq: bytes) -> set:
    """Detect whether sequence is protein or nucleotide."""
    upper = set(seq)
    non_nuc = upper - _NUCLEOTIDES
    if non_nuc:
        return _AMINO_ACIDS
    return _NUCLEOTIDES


def build_query_kmer_index(query: bytes, k: int, matrix, threshold: int) -> tuple:
    """Build a dict mapping k-mer bytes -> list of query positions.
    Returns (kmer_set, kmer_to_positions)."""
    valid_chars = detect_alphabet(query)
    kmer_to_positions = defaultdict(list)
    qlen = len(query)
    for i in range(qlen - k + 1):
        kmer = query[i:i + k]
        if not _is_valid_kmer(kmer, valid_chars):
            continue
        score = sum(matrix[query[i + j]][query[i + j]] for j in range(k))
        if score >= threshold:
            kmer_to_positions[kmer].append(i)

    kmer_set = set(kmer_to_positions)
    return kmer_set, kmer_to_positions


def find_seeds(db_seq: bytes, kmer_set: set, kmer_to_positions: dict, k: int) -> list:
    """Scan database sequence for matching k-mers.
    Returns list of (db_pos, query_pos) pairs."""
    seeds = []
    db_len = len(db_seq)
    for i in range(db_len - k + 1):
        kmer = db_seq[i:i + k]
        if kmer in kmer_set:
            for q_pos in kmer_to_positions[kmer]:
                seeds.append((i, q_pos))
    return seeds


def filter_seeds_two_hit(seeds: list, window: int = 40) -> list:
    """Keep only seeds where two hits fall on the same diagonal within `window` bp.
    Returns a list of (db_pos, q_pos) representing confirmed seed pairs."""
    if not seeds:
        return []

    by_diag = defaultdict(list)
    for db_pos, q_pos in seeds:
        by_diag[db_pos - q_pos].append((db_pos, q_pos))

    confirmed = set()
    for hits in by_diag.values():
        hits.sort()
        for i in range(1, len(hits)):
            if hits[i][0] - hits[i - 1][0] <= window:
                confirmed.add(hits[i - 1])
                confirmed.add(hits[i])

    return list(confirmed)

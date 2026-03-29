"""Seed generation: build query k-mer index and find seeds in database."""

from collections import defaultdict

# 2-bit encoding: A=0, C=1, G=2, T=3
_ENCODE = [0] * 256
for _c, _v in zip(b'AaCcGgTt', [0, 0, 1, 1, 2, 2, 3, 3]):
    _ENCODE[_c] = _v

# Characters that are valid nucleotides
_IS_VALID = bytearray(256)
for _c in b'ACGTacgt':
    _IS_VALID[_c] = 1


def _encode_kmer(seq: bytes, pos: int, k: int) -> int:
    """Encode a k-mer starting at pos as a 2-bit integer. Returns -1 if invalid."""
    val = 0
    for i in range(k):
        c = seq[pos + i]
        if not _IS_VALID[c]:
            return -1
        val = (val << 2) | _ENCODE[c]
    return val


def build_query_kmer_index(query: bytes, k: int, matrix, threshold: int) -> tuple:
    """Build a set of integer-encoded k-mers from the query, and a dict mapping
    encoded k-mer -> list of query positions.
    Returns (kmer_set, kmer_to_positions)."""
    kmer_to_positions = defaultdict(list)
    qlen = len(query)
    for i in range(qlen - k + 1):
        kmer_int = _encode_kmer(query, i, k)
        if kmer_int < 0:
            continue
        score = sum(matrix[query[i + j]][query[i + j]] for j in range(k))
        if score >= threshold:
            kmer_to_positions[kmer_int].append(i)

    kmer_set = set(kmer_to_positions)
    return kmer_set, kmer_to_positions


def find_seeds(db_seq: bytes, kmer_set: set, kmer_to_positions: dict, k: int) -> list:
    """Scan database sequence using rolling hash and find all seed hits.
    Returns list of (db_pos, query_pos) pairs."""
    seeds = []
    mask = (1 << (2 * k)) - 1
    current_hash = 0
    valid_count = 0

    for i in range(len(db_seq)):
        c = db_seq[i]
        if _IS_VALID[c]:
            current_hash = ((current_hash << 2) | _ENCODE[c]) & mask
            valid_count += 1
        else:
            valid_count = 0
            current_hash = 0

        if valid_count >= k and current_hash in kmer_set:
            pos = i - k + 1
            for q_pos in kmer_to_positions[current_hash]:
                seeds.append((pos, q_pos))

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

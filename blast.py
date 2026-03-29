#!/usr/bin/env python3

import argparse
import sys
import time
import csv
import os
from collections import namedtuple
from Bio import SeqIO

from seeding import build_query_kmer_index, find_seeds, filter_seeds_two_hit
from extension import ungapped_extend, gapped_extend, merge_hsps

DatabaseSequence = namedtuple('DatabaseSequence', ['filename', 'name', 'seq'])


def parse_args():
    parser = argparse.ArgumentParser(
        description="Simplified nucleotide BLAST",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Query is read from -q (FASTA file) or stdin (raw sequence).")
    parser.add_argument("-e", required=True, help="Path to score matrix in CSV format")
    parser.add_argument("-p", type=int, default=5, help="Gap open penalty (positive value)")
    parser.add_argument("-x", type=int, default=2, help="Gap extension penalty (positive value)")
    parser.add_argument("-d", nargs="+", required=True, help="Database FASTA file(s)")
    parser.add_argument("-k", type=int, default=11, help="Word length")
    parser.add_argument("-t", type=int, default=0, help="Seed score threshold")
    parser.add_argument("-q", default=None, help="Query file (FASTA or raw sequence)")

    parser.add_argument("--ungapped-threshold", type=int, default=20,
                        help="Minimum ungapped HSP score to trigger gapped extension")
    parser.add_argument("--two-hit-window", type=int, default=40,
                        help="Window size for two-hit seed filter")
    parser.add_argument("--no-two-hit", action="store_true",
                        help="Disable two-hit filter (use single-hit seeding)")
    parser.add_argument("--max-hits", type=int, default=50,
                        help="Maximum number of hits to report")
    return parser.parse_args()


def load_database(fasta_paths: list) -> list:
    db = []
    for path in fasta_paths:
        fname = os.path.basename(path)
        for record in SeqIO.parse(path, "fasta"):
            db.append(DatabaseSequence(fname, record.id,
                                       str(record.seq).upper().encode('ascii')))
    return db


def load_score_matrix(csv_path: str):
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


def read_query(query_path):
    with open(query_path) as f:
        content = f.read().strip()
    if content.startswith('>'):
        queries = []
        for record in SeqIO.parse(query_path, "fasta"):
            queries.append((record.id, str(record.seq).strip()))
        return queries
    else:
        seq = ''.join(content.split())
        return [("query", seq)]


def search(query_str: str, db_sequences, matrix, k: int, threshold: int,
           gap_open: int, gap_extend: int, ungapped_threshold: int,
           two_hit_window: int, use_two_hit: bool, max_hits: int):
    query = query_str.upper().encode('ascii')
    results = []

    kmer_set, kmer_to_positions = build_query_kmer_index(query, k, matrix, threshold)
    if not kmer_set:
        return results

    for db_entry in db_sequences:
        seeds = find_seeds(db_entry.seq, kmer_set, kmer_to_positions, k)
        if not seeds:
            continue

        if use_two_hit:
            seeds = filter_seeds_two_hit(seeds, two_hit_window)
            if not seeds:
                continue

        # Ungapped extension
        hsps = []
        seen = set()
        for db_pos, q_pos in seeds:
            key = (db_pos // 5, q_pos // 5)
            if key in seen:
                continue
            seen.add(key)
            hsp = ungapped_extend(db_entry.seq, query, db_pos, q_pos, k, matrix)
            if hsp[4] >= ungapped_threshold:
                hsps.append(hsp)

        if not hsps:
            continue

        # Gapped extension
        gapped_hsps = []
        for hsp in merge_hsps(hsps):
            ghsp = gapped_extend(db_entry.seq, query, *hsp[:4], matrix,
                                 gap_open, gap_extend)
            if ghsp[4] > 0:
                gapped_hsps.append((db_entry.filename, db_entry.name) + ghsp)

        gapped_hsps.sort(key=lambda h: -h[6])
        results.extend(gapped_hsps)

    # Sort all results by score descending and limit
    results.sort(key=lambda h: -h[6])
    return results[:max_hits]


def format_results(results, elapsed):
    if not results:
        print(f"\nNo hits found. Search time: {elapsed:.4f}s")
        return

    print(f"\nFound {len(results)} hit(s) in {elapsed:.4f}s:")
    print(f"{'File':<50} {'Seq':<10} {'DB Start':>10} {'DB End':>10} "
          f"{'Q Start':>10} {'Q End':>10} {'Score':>8}")
    print("-" * 110)
    for hit in results:
        fname, seqname, db_s, db_e, q_s, q_e, score = hit
        print(f"{fname:<50} {seqname:<10} {db_s + 1:>10} {db_e:>10} "
              f"{q_s + 1:>10} {q_e:>10} {score:>8}")
    print(f"\nSearch time: {elapsed:.4f}s")


def main():
    args = parse_args()
    matrix = load_score_matrix(args.e)

    load_start = time.time()
    db_sequences = load_database(args.d)
    print(f"Loaded {len(db_sequences)} sequence(s), "
          f"{sum(len(s.seq) for s in db_sequences):,} bases in "
          f"{time.time() - load_start:.2f}s", file=sys.stderr)

    threshold = args.t
    if threshold <= 0:
        match_score = matrix[ord('A')][ord('A')]
        threshold = args.k * match_score
    print(f"Parameters: k={args.k}, T={threshold}, "
          f"gap_open={args.p}, gap_ext={args.x}, "
          f"two_hit={'off' if args.no_two_hit else 'on'}",
          file=sys.stderr)

    use_two_hit = not args.no_two_hit
    print(f"Stars", file=sys.stderr)

    if args.q:
        queries = read_query(args.q)
    else:
        raw = []
        for line in sys.stdin:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            raw.append(line)
        if raw:
            queries = [("stdin_query", ''.join(raw))]
        else:
            print("No query provided.", file=sys.stderr)
            sys.exit(1)
    print(f"Searching for query ({len(queries)} bp)...", file=sys.stderr)

    for qname, qseq in queries:
        print(f"\nQuery: {qname} ({len(qseq)} bp)", file=sys.stderr)
        start_time = time.time()

        results = search(qseq, db_sequences, matrix, args.k, threshold,
                         args.p, args.x, args.ungapped_threshold,
                         args.two_hit_window, use_two_hit, args.max_hits)

        elapsed = time.time() - start_time
        format_results(results, elapsed)


if __name__ == "__main__":
    main()

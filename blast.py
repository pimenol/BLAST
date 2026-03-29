#!/usr/bin/env python3

import argparse
import sys
import time
import csv

from index import load_database
from seeding import build_query_kmer_index, find_seeds, filter_seeds_two_hit
from extension import ungapped_extend, gapped_extend, merge_hsps


def parse_args():
    parser = argparse.ArgumentParser(description="BLAST")
    parser.add_argument("-e", required=True, help="Path to score matrix in csv format")
    parser.add_argument("-p", type=int, default=11, help="Gap open penalty (positive value, default: 11)")
    parser.add_argument("-x", type=int, default=1, help="Gap extension penalty (positive value, default: 1)")
    parser.add_argument("-d", nargs="+", required=True, help="list of database files in fasta format, for example, -d chr1.fasta chr2.fasta")
    parser.add_argument("-k", type=int, default=3, help="length of the word (default: 3 for protein)")
    parser.add_argument("-t", type=int, default=0, help="threshold on the seed score (default: 0)")

    parser.add_argument("--ungapped-threshold", type=int, default=11, help="Minimum ungapped HSP score for gapped extension")
    parser.add_argument("--two-hit-window", type=int, default=40, help="Window size for two-hit seed filter")
    parser.add_argument("--no-two-hit", action="store_true", help="Disable two-hit filter")

    return parser.parse_args()

def load_score_matrix(csv_path: str):
    """Load a substitution matrix from CSV. Returns a 128x128 list-of-lists
    where matrix[ord(a)][ord(b)] = score"""
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

def search(query_str: str, db_sequences, matrix, k: int, threshold: int,
           gap_open: int, gap_extend: int, ungapped_threshold: int,
           two_hit_window: int, use_two_hit: bool):
    """Run BLAST search for a query against all database sequences."""
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

        gapped_hsps = []
        for hsp in merge_hsps(hsps):
            ghsp = gapped_extend(db_entry.seq, query, *hsp[:4],
                                 matrix, gap_open, gap_extend)
            if ghsp[4] > 0:
                gapped_hsps.append((db_entry.filename, db_entry.name) + ghsp)

        gapped_hsps.sort(key=lambda h: -h[6])
        results.extend(gapped_hsps)

    return results


def main():
    args = parse_args()
    matrix = load_score_matrix(args.e)

    load_start = time.time()
    db_sequences = load_database(args.d)
    print(f"Loaded {len(db_sequences)} sequences, "
          f"{sum(len(s.seq) for s in db_sequences):,} residues in "
          f"{time.time() - load_start:.2f}s", file=sys.stderr)

    threshold = args.t
    if threshold <= 0:
        match_score = matrix[ord('A')][ord('A')]
        threshold = args.k * match_score
        print(f"Auto threshold: {threshold} (k={args.k} * match={match_score})",
              file=sys.stderr)

    print("Ready. Enter query sequences (one per line, empty line or EOF to quit):",
          file=sys.stderr)

    for line in sys.stdin:
        query_str = line.strip()
        if not query_str:
            continue

        print(f"Searching for query ({len(query_str)} residues)...", file=sys.stderr)
        start_time = time.time()

        results = search(query_str, db_sequences, matrix, args.k, threshold,
                         args.p, args.x, args.ungapped_threshold,
                         args.two_hit_window, not args.no_two_hit)

        elapsed = time.time() - start_time

        if results:
            print(f"\nFound {len(results)} hit(s) in {elapsed:.4f}s:")
            print(f"{'File':<30} {'Seq Name':<20} {'DB Start':>10} {'DB End':>10} "
                  f"{'Q Start':>10} {'Q End':>10} {'Score':>8}")
            print("-" * 100)
            for hit in results:
                fname, seqname, db_s, db_e, q_s, q_e, score = hit
                # 1-based coordinates for output
                print(f"{fname:<30} {seqname:<20} {db_s + 1:>10} {db_e:>10} "
                      f"{q_s + 1:>10} {q_e:>10} {score:>8}")
            print(f"\nSearch time: {elapsed:.4f}s")
        else:
            print(f"\nNo hits found. Search time: {elapsed:.4f}s")

        print("", file=sys.stderr)  # blank line separator


if __name__ == "__main__":
    main()

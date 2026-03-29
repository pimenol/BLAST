#!/usr/bin/env python3
"""Tests for the simplified BLAST implementation."""

import subprocess
import os
import tempfile

DIR = os.path.dirname(os.path.abspath(__file__))
BLAST = os.path.join(DIR, "blast.sh")
MATRIX = os.path.join(DIR, "matrix.csv")


def run_blast(query, db_fasta, k=11, extra_args=None):
    """Run blast.sh with a query string piped via stdin. Returns stdout."""
    cmd = [BLAST, "-e", MATRIX, "-d", db_fasta, "-k", str(k), "--no-two-hit"]
    if extra_args:
        cmd.extend(extra_args)
    result = subprocess.run(cmd, input=query, capture_output=True, text=True, timeout=30)
    return result.stdout, result.stderr


def make_fasta(name, seq):
    """Write a temporary FASTA file and return its path."""
    f = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False, dir=DIR)
    f.write(f">{name}\n{seq}\n")
    f.close()
    return f.name


def parse_hits(stdout):
    """Parse hit lines from blast output. Returns list of dicts."""
    hits = []
    in_table = False
    for line in stdout.strip().split('\n'):
        if line.startswith('---'):
            in_table = True
            continue
        if in_table and line.strip():
            parts = line.split()
            if len(parts) >= 7:
                hits.append({
                    'file': parts[0],
                    'seq': parts[1],
                    'db_start': int(parts[2]),
                    'db_end': int(parts[3]),
                    'q_start': int(parts[4]),
                    'q_end': int(parts[5]),
                    'score': int(parts[6]),
                })
    return hits


def test_exact_substring():
    """Query is an exact copy from the database — expect a perfect match."""
    db_seq = "AAAAAA" + "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC" + "TTTTTT"
    query = "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC"
    db_file = make_fasta("test_exact", db_seq)
    try:
        stdout, stderr = run_blast(query, db_file)
        hits = parse_hits(stdout)
        assert len(hits) >= 1, f"Expected at least 1 hit, got {len(hits)}"
        best = hits[0]
        # Perfect match: score = length * 5 (match score)
        expected_score = len(query) * 5  # 60 * 5 = 300
        assert best['score'] == expected_score, f"Expected score {expected_score}, got {best['score']}"
        # Should align the full query
        assert best['q_start'] == 1, f"Expected q_start=1, got {best['q_start']}"
        assert best['q_end'] == len(query), f"Expected q_end={len(query)}, got {best['q_end']}"
        # Should find it at the correct position (1-based: pos 7 to 66)
        assert best['db_start'] == 7, f"Expected db_start=7, got {best['db_start']}"
        assert best['db_end'] == 66, f"Expected db_end=66, got {best['db_end']}"
        print("PASS: test_exact_substring")
    finally:
        os.unlink(db_file)


def test_no_match():
    """Query has no similarity to database — expect no hits."""
    db_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    query = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
    db_file = make_fasta("test_nomatch", db_seq)
    try:
        stdout, stderr = run_blast(query, db_file)
        hits = parse_hits(stdout)
        assert len(hits) == 0, f"Expected 0 hits, got {len(hits)}"
        print("PASS: test_no_match")
    finally:
        os.unlink(db_file)


def test_partial_match():
    """Query partially overlaps with database — expect a hit covering the overlap."""
    db_seq = "TTTTTT" + "ACTGACTGACTGACTGACTGACTGACTGACTG" + "AAAAAA"
    #                     ^^^^^^^^^^^^^^^^^^^^^^^^ (first 22 bp match)
    query = "ACTGACTGACTGACTGACTGACTGCCCCCCCCCCCCCCCCCCCCCC"
    db_file = make_fasta("test_partial", db_seq)
    try:
        stdout, stderr = run_blast(query, db_file, k=11)
        hits = parse_hits(stdout)
        assert len(hits) >= 1, f"Expected at least 1 hit, got {len(hits)}"
        best = hits[0]
        assert best['score'] > 0, f"Expected positive score, got {best['score']}"
        # The hit should NOT cover the full query (only partial match)
        alignment_len = best['q_end'] - best['q_start'] + 1
        assert alignment_len < len(query), f"Expected partial alignment, got full ({alignment_len} bp)"
        print("PASS: test_partial_match")
    finally:
        os.unlink(db_file)


def test_multiple_hits():
    """Database has two copies of a sequence — expect two hits."""
    segment = "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC"
    spacer = "A" * 1000
    db_seq = segment + spacer + segment
    query = segment
    db_file = make_fasta("test_multi", db_seq)
    try:
        stdout, stderr = run_blast(query, db_file)
        hits = parse_hits(stdout)
        assert len(hits) >= 2, f"Expected at least 2 hits, got {len(hits)}"
        # Both should have high scores
        for h in hits[:2]:
            assert h['score'] > 100, f"Expected high score, got {h['score']}"
        print("PASS: test_multiple_hits")
    finally:
        os.unlink(db_file)


def test_mismatch_reduces_score():
    """A query with mismatches should score lower than an exact match."""
    db_seq = "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC"
    exact_query = db_seq
    # Introduce 5 mismatches
    mutated = list(db_seq)
    for i in [5, 15, 25, 35, 45]:
        mutated[i] = 'A' if mutated[i] != 'A' else 'T'
    mutated_query = ''.join(mutated)

    db_file = make_fasta("test_mismatch", db_seq)
    try:
        stdout_exact, _ = run_blast(exact_query, db_file)
        stdout_mut, _ = run_blast(mutated_query, db_file)
        hits_exact = parse_hits(stdout_exact)
        hits_mut = parse_hits(stdout_mut)

        assert len(hits_exact) >= 1, "No hits for exact query"
        assert len(hits_mut) >= 1, "No hits for mutated query"
        assert hits_mut[0]['score'] < hits_exact[0]['score'], \
            f"Mutated score ({hits_mut[0]['score']}) should be less than exact ({hits_exact[0]['score']})"
        print("PASS: test_mismatch_reduces_score")
    finally:
        os.unlink(db_file)


def test_reverse_complement_no_hit():
    """BLAST is strand-specific — reverse complement should NOT match (we only search forward strand)."""
    db_seq = "TTTTTT" + "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCAT" + "AAAAAA"
    # Reverse complement of the segment
    rc_map = str.maketrans("ACGT", "TGCA")
    segment = "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCAT"
    rev_comp = segment[::-1].translate(rc_map)
    db_file = make_fasta("test_rc", db_seq)
    try:
        stdout, stderr = run_blast(rev_comp, db_file, k=11)
        hits = parse_hits(stdout)
        # Should not find a high-scoring hit since we don't search reverse complement
        if hits:
            assert hits[0]['score'] < 50, \
                f"Reverse complement should not produce a strong hit, got score {hits[0]['score']}"
        print("PASS: test_reverse_complement_no_hit")
    finally:
        os.unlink(db_file)


def test_gap_in_alignment():
    """Query has a small deletion relative to DB — gapped extension should still align well."""
    db_seq = "AAAAAA" + "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC" + "TTTTTT"
    # Delete 3 bases from the middle of the query (positions 25-27)
    full = "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC"
    query_with_del = full[:25] + full[28:]  # 3bp deletion
    db_file = make_fasta("test_gap", db_seq)
    try:
        stdout, stderr = run_blast(query_with_del, db_file, k=11)
        hits = parse_hits(stdout)
        assert len(hits) >= 1, f"Expected at least 1 hit with gapped alignment, got {len(hits)}"
        best = hits[0]
        # Score should be positive and substantial (most bases still match)
        assert best['score'] > 100, f"Expected score > 100 for gapped alignment, got {best['score']}"
        print("PASS: test_gap_in_alignment")
    finally:
        os.unlink(db_file)


def test_insertion_in_query():
    """Query has a small insertion relative to DB — gapped extension should handle it."""
    db_seq = "AAAAAA" + "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC" + "TTTTTT"
    full = "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC"
    # Insert 3 bases in the middle
    query_with_ins = full[:30] + "TTT" + full[30:]
    db_file = make_fasta("test_ins", db_seq)
    try:
        stdout, stderr = run_blast(query_with_ins, db_file, k=11)
        hits = parse_hits(stdout)
        assert len(hits) >= 1, f"Expected at least 1 hit, got {len(hits)}"
        best = hits[0]
        assert best['score'] > 100, f"Expected score > 100 for insertion alignment, got {best['score']}"
        print("PASS: test_insertion_in_query")
    finally:
        os.unlink(db_file)


def test_short_query_below_k():
    """Query shorter than k should produce no hits (can't form a seed)."""
    db_seq = "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC"
    query = "GCTGGGAAGA"  # 10 bp, k=11
    db_file = make_fasta("test_short", db_seq)
    try:
        stdout, stderr = run_blast(query, db_file, k=11)
        hits = parse_hits(stdout)
        assert len(hits) == 0, f"Expected 0 hits for query shorter than k, got {len(hits)}"
        print("PASS: test_short_query_below_k")
    finally:
        os.unlink(db_file)


def test_different_k_values():
    """Smaller k should be more sensitive (find more/better hits with mismatches)."""
    db_seq = "AAAAAA" + "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC" + "TTTTTT"
    # Query with scattered mismatches — harder to seed with large k
    full = "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC"
    mutated = list(full)
    # Mutate every 10th base — breaks most 11-mers
    for i in range(5, len(mutated), 10):
        mutated[i] = 'A' if mutated[i] != 'A' else 'T'
    query = ''.join(mutated)

    db_file = make_fasta("test_kvals", db_seq)
    try:
        stdout_k7, _ = run_blast(query, db_file, k=7)
        stdout_k11, _ = run_blast(query, db_file, k=11)
        hits_k7 = parse_hits(stdout_k7)
        hits_k11 = parse_hits(stdout_k11)
        # k=7 should find at least as many hits as k=11
        assert len(hits_k7) >= len(hits_k11), \
            f"k=7 ({len(hits_k7)} hits) should be >= k=11 ({len(hits_k11)} hits)"
        print("PASS: test_different_k_values")
    finally:
        os.unlink(db_file)


def test_n_characters_ignored():
    """N characters in database should not produce false matches."""
    db_seq = "N" * 100 + "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC" + "N" * 100
    query = "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC"
    db_file = make_fasta("test_nchars", db_seq)
    try:
        stdout, stderr = run_blast(query, db_file, k=11)
        hits = parse_hits(stdout)
        assert len(hits) >= 1, f"Expected at least 1 hit, got {len(hits)}"
        best = hits[0]
        # Should find the real sequence, not match Ns
        assert best['db_start'] == 101, f"Expected db_start=101, got {best['db_start']}"
        assert best['score'] == len(query) * 5, f"Expected perfect score, got {best['score']}"
        print("PASS: test_n_characters_ignored")
    finally:
        os.unlink(db_file)


def test_multiple_database_files():
    """Searching across two database files should find hits in both."""
    seq1 = "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC"
    seq2 = "ACTCTTCTGGTCCCCACAGACTCAGAGAGAACCCACCATGGTGCTGTCTCCTGCCGACAA"
    query = seq1  # matches db1 exactly
    db1 = make_fasta("db1", seq1)
    db2 = make_fasta("db2", seq2)
    try:
        cmd = [BLAST, "-e", MATRIX, "-d", db1, db2, "-k", "11", "--no-two-hit"]
        result = subprocess.run(cmd, input=query, capture_output=True, text=True, timeout=30)
        hits = parse_hits(result.stdout)
        assert len(hits) >= 1, f"Expected at least 1 hit, got {len(hits)}"
        # Best hit should be exact match with perfect score
        assert hits[0]['score'] == len(query) * 5, \
            f"Expected perfect score {len(query)*5}, got {hits[0]['score']}"
        print("PASS: test_multiple_database_files")
    finally:
        os.unlink(db1)
        os.unlink(db2)


def test_gap_penalties_affect_score():
    """Higher gap penalties should produce lower scores for gapped alignments."""
    db_seq = "AAAAAA" + "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC" + "TTTTTT"
    full = "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC"
    query_with_del = full[:25] + full[28:]  # 3bp deletion forces gap

    db_file = make_fasta("test_gappen", db_seq)
    try:
        stdout_low, _ = run_blast(query_with_del, db_file, k=11, extra_args=["-p", "2", "-x", "1"])
        stdout_high, _ = run_blast(query_with_del, db_file, k=11, extra_args=["-p", "10", "-x", "5"])
        hits_low = parse_hits(stdout_low)
        hits_high = parse_hits(stdout_high)
        assert len(hits_low) >= 1, "No hits with low gap penalty"
        assert len(hits_high) >= 1, "No hits with high gap penalty"
        assert hits_low[0]['score'] >= hits_high[0]['score'], \
            f"Low penalty score ({hits_low[0]['score']}) should be >= high penalty ({hits_high[0]['score']})"
        print("PASS: test_gap_penalties_affect_score")
    finally:
        os.unlink(db_file)


def test_query_from_fasta_file():
    """Test -q flag with a FASTA query file."""
    db_seq = "AAAAAA" + "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC" + "TTTTTT"
    query_seq = "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC"
    db_file = make_fasta("test_qfile_db", db_seq)
    q_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False, dir=DIR)
    q_file.write(f">my_query\n{query_seq}\n")
    q_file.close()
    try:
        cmd = [BLAST, "-e", MATRIX, "-d", db_file, "-k", "11", "--no-two-hit", "-q", q_file.name]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        hits = parse_hits(result.stdout)
        assert len(hits) >= 1, f"Expected at least 1 hit from -q file, got {len(hits)}"
        assert hits[0]['score'] == len(query_seq) * 5, \
            f"Expected perfect score {len(query_seq)*5}, got {hits[0]['score']}"
        print("PASS: test_query_from_fasta_file")
    finally:
        os.unlink(db_file)
        os.unlink(q_file.name)


def test_score_proportional_to_match_length():
    """Longer exact matches should have proportionally higher scores."""
    base = "GCTGGGAAGACCCCCAAGTCCCTCTTCTGCATCGTCCTCGGGCTCTGGCTTGGTGCTCAC"
    db_seq = base * 3  # 180 bp
    db_file = make_fasta("test_proplen", db_seq)
    try:
        # Query 30 bp
        stdout_30, _ = run_blast(base[:30], db_file, k=11)
        # Query 60 bp
        stdout_60, _ = run_blast(base[:60], db_file, k=11)
        hits_30 = parse_hits(stdout_30)
        hits_60 = parse_hits(stdout_60)
        assert len(hits_30) >= 1 and len(hits_60) >= 1, "Expected hits for both queries"
        assert hits_60[0]['score'] > hits_30[0]['score'], \
            f"60bp score ({hits_60[0]['score']}) should be > 30bp score ({hits_30[0]['score']})"
        print("PASS: test_score_proportional_to_match_length")
    finally:
        os.unlink(db_file)


if __name__ == "__main__":
    test_exact_substring()
    test_no_match()
    test_partial_match()
    test_multiple_hits()
    test_mismatch_reduces_score()
    test_reverse_complement_no_hit()
    test_gap_in_alignment()
    test_insertion_in_query()
    test_short_query_below_k()
    test_different_k_values()
    test_n_characters_ignored()
    test_multiple_database_files()
    test_gap_penalties_affect_score()
    test_query_from_fasta_file()
    test_score_proportional_to_match_length()
    print("\nAll tests passed!")

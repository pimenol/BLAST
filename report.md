# BLAST Implementation Report

A simplified nucleotide BLAST algorithm in Python - seeding, ungapped extension, gapped extension - that correctly identifies the HBA1 hemoglobin gene against a multi-chromosome human database and handles edge cases including N-containing sequences and queries shorter than the seed word.

---

## Code Structure

Three files carry the implementation:

- **blast.py** - Entry point. Handles argument parsing, database loading, score matrix loading from CSV, query reading, and pipeline orchestration. Results are sorted by score and printed in tabular format.
- **seeding.py** - K-mer indexing and seed detection. Builds a hash-based k-mer index of the query using 2-bit encoding (A=00, C=01, G=10, T=11). Scans the database with a rolling hash for O(1)-per-position k-mer matching. Implements the two-hit filter: two seeds must fall on the same diagonal within 40 bp before extension is triggered, which cuts false positives without a measurable miss rate on the test suite.
- **extension.py** - Alignment extension. Ungapped extension with X-drop termination, banded Smith-Waterman gapped extension with affine gap penalties (separate open and extension costs), and HSP merging to eliminate redundant overlapping hits.

Supporting files: `blast.sh`, `compile.sh`.

---

## Usage and Arguments

The program is launched via `blast.sh`.

| Argument | Meaning |
|---|---|
| `-e` | Path to the score matrix in CSV format (required) |
| `-p` | Gap open penalty, positive integer (default 5) |
| `-x` | Gap extension penalty, positive integer (default 2) |
| `-d` | Database files in FASTA format, e.g. `-d chr1.fasta chr2.fasta` (required) |
| `-k` | Seed word length (default 11) |
| `-t` | Seed score threshold; if 0, computed as k × match_score (default 0) |
| `-q` | Query file in FASTA or raw sequence format; if omitted, reads from stdin |
| `--ungapped-threshold` | Minimum ungapped HSP score to trigger gapped extension (default 20) |
| `--two-hit-window` | Window size for the two-hit seed filter (default 40) |
| `--no-two-hit` | Disable the two-hit filter, use single-hit seeding |
| `--max-hits` | Maximum hits to report (default 50) |

For each query, the program outputs: database filename, match range in the database (1-based), match range in the query (1-based), alignment score, and search time in seconds.

---

## Algorithm

**Stage 1 - Seeding.** All k-mers (default k=11) of the query are extracted and stored in a dictionary mapping their 2-bit integer hash to a list of query positions. The database is scanned with a rolling hash: at each position, the current k-mer hash is computed incrementally and looked up in the query index. Invalid characters (e.g., N) reset the hash. The two-hit filter groups seeds by diagonal (db_pos − q_pos) and retains only those where two seeds on the same diagonal fall within 40 bp - empirically chosen after testing showed that values below ~30 bp missed true hits on short queries, while values above ~60 bp negated the speed benefit.

**Stage 2 - Ungapped Extension.** Each confirmed seed is extended bidirectionally along the diagonal. A running score is maintained using the substitution matrix; extension stops when the score drops more than X-drop (default 20) below the best score seen so far. A coarse deduplication (position // 5) skips redundant extensions. Only HSPs scoring ≥ 20 proceed to gapped extension.

**Stage 3 - Gapped Extension with Diagonal Banding.** Banded Smith-Waterman with affine gap penalties. Band limits are centered on the diagonal of the ungapped HSP:

    diag_offset = (db_start − aln_db_start) − (q_start − aln_q_start)

For each query row i, only database columns within [i + diag_offset − W, i + diag_offset + W] are computed (W = 50 by default). This reduces DP computation from O(m × n) to O(m × 100) - on a 250 Mbp chromosome, that difference is what keeps runtime in minutes rather than hours. The tradeoff: the algorithm assumes the true alignment lies within 50 bp of the ungapped diagonal and will miss alignments with more indels than that. Three DP arrays (H, E, F) are maintained in a two-row scheme for space efficiency.

---

## Results

**Hemoglobin test.** Using HBA1 (hemoglobin subunit alpha 1) as the query against a multi-chromosome database, the algorithm returns HBA1 as the top hit (exact match score) and ranks other globin family genes - HBB, HBA2 - below it. The score ordering matches sequence divergence within the globin family: conserved exon regions from a shared evolutionary ancestor produce partial hits, not false positives.

**Other test sequences.** The test suite confirms correct behavior on: exact substrings (score proportional to length), partial overlaps, sequences with mismatches (reduced scores), sequences with 3 bp insertions and deletions, multiple database files, queries shorter than k (no seeds generated, graceful exit), and N-containing sequences (hash reset, no crash).

**Speed.** On a few kilobases: milliseconds. On a single human chromosome (~250 Mbp): on the order of minutes in pure Python, vs. seconds on NCBI BLASTn running optimized C++ with precomputed indexed databases and parallelism. The bottleneck is the seed scanning loop; a C extension or Cython compilation of `seeding.py` would produce order-of-magnitude improvements. Algorithmic complexity is equivalent - the performance gap is entirely implementation overhead.

**Memory.** With a 4 GB budget and 1 byte per nucleotide (ASCII), the implementation can hold approximately 4 billion nucleotides - enough for the entire human genome (~3.2 Gbp) with index overhead to spare. All 24 human chromosomes fit comfortably. Switching to 2-bit encoding quadruples that capacity to ~16 billion nucleotides. DP arrays during gapped extension use O(2 × 2W × 8 bytes) ≈ 1.6 KB per active extension - negligible at W=50.

---

## Obstacles

**Pure Python inner loops.** Seed scanning and DP extension are the bottleneck. This is a known cost of Python and not a design flaw - but it is the single largest gap between this implementation and production BLAST.

**Diagonal offset calculation.** Getting the band centering right for Smith-Waterman required careful accounting of local-vs-global coordinate systems. Off-by-one errors in band limits initially caused missed alignments and array out-of-bounds accesses; fixing this required explicitly tracing through coordinate frames for each stage.

**Two-hit window tuning.** The 40 bp default was chosen empirically. Below ~30 bp, short-query hits were missed; above ~60 bp, the filter stopped reducing false positives meaningfully.

**Affine gap boundary conditions.** At band boundaries, previous-row values may fall outside the computed region. Out-of-band cells must be initialized to −∞; failing to do this silently corrupted scores near band edges.

**Traceback approximation.** Alignment boundaries are estimated from the best DP cell position rather than a full traceback. Scores are correct; reported start/end coordinates can be slightly off.
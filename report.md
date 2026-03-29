# BLAST Implementation Report

## Overview

This project implements a simplified nucleotide BLAST algorithm
in Python. The implementation follows the classical three-stage BLAST pipeline: seeding, ungapped
extension, and gapped extension. The goal is to find regions of local similarity between a query
nucleotide sequence and a database of sequences stored in FASTA format.

## Code Structure

The implementation is split across 3 files:

- **blast.py** - Main entry point. Handles argument parsing, database loading, score matrix loading from CSV, query reading, and orchestration of the search
  pipeline. Results are sorted by score and printed in a tabular format.

- **seeding.py** - K-mer indexing and seed detection. Builds a hash-based k-mer index
  of the query using 2-bit encoding (A=00, C=01, G=10, T=11). Scans the database with a rolling
  hash for O(1)-per-position k-mer matching. Implements the two-hit filter, which requires two seed
  hits on the same diagonal within a configurable window (default 40 bp) before proceeding.

- **extension.py** - Alignment extension. Contains ungapped extension with X-drop
  termination and banded Smith-Waterman gapped extension with affine gap penalties (separate gap
  open and gap extension costs). Also includes HSP merging to eliminate redundant overlapping hits.


Supporting files include `blast.sh`, `compile.sh`.

## Usage and Arguments

The program is launched via `blast.sh`, which accepts the following arguments:

| Argument | Meaning |
|----------|---------|
| `-e`     | Path to the score matrix in CSV format (required) |
| `-p`     | Gap open penalty, positive integer (optional, default 5) |
| `-x`     | Gap extension penalty, positive integer (optional, default 2) |
| `-d`     | List of database files in FASTA format, e.g. `-d chr1.fasta chr2.fasta` (required) |
| `-k`     | Length of the seed word (optional, default 11) |
| `-t`     | Threshold on the seed score; if 0, computed as k * match_score (optional, default 0) |
| `-q`     | Query file in FASTA or raw sequence format (optional; if omitted, reads from stdin) |
| `--ungapped-threshold` | Minimum ungapped HSP score to trigger gapped extension (default 20) |
| `--two-hit-window` | Window size for the two-hit seed filter (default 40) |
| `--no-two-hit` | Disable the two-hit filter, use single-hit seeding instead |
| `--max-hits` | Maximum number of hits to report (default 50) |

Once launched, the program loads all database sequences into memory and, if no `-q` file is
provided, listens on standard input for a query sequence. For each query, it outputs all matches
with: the **database filename**, the **match range in the database** (1-based start and end), the
**match range in the query** (1-based start and end), the alignment **score**, and the **search
time** in seconds.

## Algorithm Details

**Stage 1 - Seeding.** All k-mers (default k=11) of the query are extracted and stored in a
dictionary mapping their 2-bit integer hash to a list of query positions. The database is scanned
using a rolling hash: at each position, the current k-mer hash is computed incrementally and looked
up in the query index. Invalid characters (e.g., N) reset the hash. A two-hit filter groups seeds
by diagonal (db_pos - q_pos) and retains only those where two seeds on the same diagonal fall within
40 bp, dramatically reducing false positives.

**Stage 2 - Ungapped Extension.** Each confirmed seed is extended bidirectionally along the
diagonal. A running score is maintained using the substitution matrix; extension terminates when the
score drops more than X-drop (default 20) below the best score seen so far. A coarse deduplication
(position // 5) prevents redundant extensions. Only HSPs scoring >= 20 proceed to gapped extension.

**Stage 3 - Gapped Extension with Diagonal Banding.** This is the most computationally intensive
stage. A banded Smith-Waterman algorithm with affine gap penalties is used. One of the optimization is
band limits centered on the expected diagonal. The diagonal offset is computed from the ungapped
HSP coordinates:

    diag_offset = (db_start - aln_db_start) - (q_start - aln_q_start)

For each query row i of the DP matrix, only database columns within [i + diag_offset - W,
i + diag_offset + W] are computed, where W is the band width (default 50). This reduces the DP
computation from O(m * n) to O(m * 2W), which is critical for aligning against long chromosomal
sequences. The band is always centered on the diagonal where the ungapped hit was found, meaning the
algorithm expects the true alignment to lie near that diagonal and tolerates indels of up to W
nucleotides in either direction. Three DP arrays (H, E, F) are maintained in a two-row scheme for
space efficiency. The affine gap model uses gap_open = 5 and gap_extend = 2 by default.

## Results

**Hemoglobin test.** Using the HBA1 (hemoglobin subunit alpha 1) gene as a query against a database
containing human chromosome fragments, the algorithm successfully identifies hemoglobin. When tested
against a database containing multiple globin family genes, the algorithm finds HBA1 with the
highest score (exact match), and other globin genes (HBB, HBA2) appear as lower-scoring hits due to
sequence similarity within the globin family. This matches expectations: globin genes share
conserved exon regions from a common evolutionary ancestor.

**Other test sequences.** The test suite confirms correct behavior on exact substrings (score
proportional to length), partial overlaps, sequences with mismatches (reduced scores), sequences
with small insertions and deletions (3 bp), multiple database files, and edge cases such as queries
shorter than k (no seeds generated) and sequences containing N characters (gracefully skipped).

**Speed comparison.** On small databases (a few kilobases), the algorithm completes in milliseconds.
On a single human chromosome (~250 Mbp), the search takes on the order of minutes in pure Python,
compared to seconds on the NCBI BLASTn server, which uses highly optimized C++ code, precomputed
indexed databases, and parallelism. The primary bottleneck in our implementation is the seed scanning
loop in Python; a C extension or Cython compilation of seeding.py would yield order-of-magnitude
improvements. Despite this, the algorithmic complexity is equivalent.

## Memory Analysis

With a 4 GB memory budget, we can estimate capacity as follows. Each nucleotide in the database is
stored as one byte (ASCII encoding). The k-mer index for the query is negligible (query is small).
The main memory consumers are: (1) the database sequences, and (2) the DP arrays during gapped
extension. The DP arrays use O(2 * 2W * 8 bytes) per active extension, which is ~1.6 KB for W=50 -
negligible. Therefore, nearly all 4 GB can be allocated to database storage. At 1 byte per
nucleotide, we can store approximately **4 billion nucleotides** - enough for the entire human
genome (~3.2 billion bp) with room to spare for the index overhead. In terms of chromosomes, this
comfortably fits all 24 human chromosomes (22 autosomes + X + Y). With 2-bit encoding (packing 4
nucleotides per byte), this capacity quadruples to ~16 billion nucleotides.

## Obstacles

- **Pure Python performance.** The inner loops of seed scanning and DP extension are slow in
  interpreted Python. This is the single largest limitation and the main reason our implementation
  is orders of magnitude slower than NCBI BLAST on large databases.
- **Band boundary correctness.** Getting the diagonal offset calculation right for the banded
  Smith-Waterman required careful reasoning about coordinate systems (local subsequence coordinates
  vs. global sequence coordinates). Off-by-one errors in band limits initially caused missed
  alignments or array out-of-bounds accesses.
- **Two-hit filter tuning.** An overly strict two-hit window missed true hits on short queries;
  an overly loose window negated the speed benefit. The default of 40 bp was chosen empirically.
- **Affine gap penalty integration.** Correctly maintaining three separate DP recurrences (H, E, F)
  within the banded region, especially at band boundaries where previous-row values may fall outside
  the band, required initializing out-of-band cells to negative infinity.
- **Traceback approximation.** The implementation estimates alignment boundaries from the best DP
  cell position rather than performing a full traceback, which can slightly misreport alignment
  start/end coordinates (though the score is correct).

## Comparison with NCBI BLAST

| Aspect                | Our Implementation           | NCBI BLASTn                       |
|-----------------------|------------------------------|-----------------------------------|
| Language              | Python                       | C++                               |
| Seeding               | Exact k-mer match            | Exact + neighborhood words        |
| Two-hit filter        | Yes                          | Yes                               |
| Gapped extension      | Banded Smith-Waterman        | Banded Smith-Waterman             |
| Band centering        | On ungapped hit diagonal     | On ungapped hit diagonal          |
| Gap model             | Affine                       | Affine                            |
| Statistics (E-value)  | Not implemented              | Karlin-Altschul statistics        |
| Database indexing      | None (sequential scan)       | Precomputed lookup tables         |
| Speed (human chr)     | Minutes                      | Seconds                           |
| Traceback             | Approximate (from best cell) | Full DP traceback with CIGAR      |

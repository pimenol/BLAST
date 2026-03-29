NEG_INF = float('-inf')


def ungapped_extend(db_seq: bytes, query: bytes, db_pos: int, q_pos: int,
                    k: int, matrix, x_drop: int = 20):
    seed_score = sum(matrix[db_seq[db_pos + i]][query[q_pos + i]] for i in range(k))

    # Extend right
    right_best = seed_score
    right_score = seed_score
    db_right = db_pos + k
    q_right = q_pos + k
    best_db_right = db_right
    best_q_right = q_right

    while db_right < len(db_seq) and q_right < len(query):
        right_score += matrix[db_seq[db_right]][query[q_right]]
        if right_score > right_best:
            right_best = right_score
            best_db_right = db_right + 1
            best_q_right = q_right + 1
        if right_best - right_score >= x_drop:
            break
        db_right += 1
        q_right += 1

    # Extend left
    left_best = seed_score
    left_score = seed_score
    db_left = db_pos - 1
    q_left = q_pos - 1
    best_db_left = db_pos
    best_q_left = q_pos

    while db_left >= 0 and q_left >= 0:
        left_score += matrix[db_seq[db_left]][query[q_left]]
        if left_score > left_best:
            left_best = left_score
            best_db_left = db_left
            best_q_left = q_left
        if left_best - left_score >= x_drop:
            break
        db_left -= 1
        q_left -= 1

    total_score = left_best + right_best - seed_score
    return (best_db_left, best_db_right, best_q_left, best_q_right, total_score)


def gapped_extend(db_seq: bytes, query: bytes, db_start: int, db_end: int,
                  q_start: int, q_end: int, matrix,
                  gap_open: int, gap_extend: int,
                  band_width: int = 50, x_drop: int = 50):
    """Uses banded Smith-Waterman with affine gap penalties.
    Returns (db_start, db_end, q_start, q_end, score)."""

    extend = max(db_end - db_start, q_end - q_start, band_width)

    aln_db_start = max(0, db_start - extend)
    aln_q_start = max(0, q_start - extend)

    d = db_seq[aln_db_start:min(len(db_seq), db_end + extend)]
    q = query[aln_q_start:min(len(query), q_end + extend)]

    n, m = len(d), len(q)

    if m == 0 or n == 0:
        return (db_start, db_end, q_start, q_end, 0)

    diag_offset = (db_start - aln_db_start) - (q_start - aln_q_start)

    prev_H = [0] * (n + 1)
    prev_E = [NEG_INF] * (n + 1)
    curr_H = [0] * (n + 1)
    curr_E = [NEG_INF] * (n + 1)
    curr_F = NEG_INF

    best_score = 0
    best_i = 0
    best_j = 0

    for i in range(1, m + 1):
        curr_F = NEG_INF
        qc = q[i - 1]

        # Band limits centered on the expected diagonal
        center = i + diag_offset
        j_start = max(1, center - band_width)
        j_end = min(n, center + band_width)

        for j in range(j_start, j_end + 1):
            match_score = matrix[qc][d[j - 1]]
            h_diag = prev_H[j - 1] + match_score

            e_open = prev_H[j] - gap_open - gap_extend
            e_ext = prev_E[j] - gap_extend
            curr_E[j] = max(e_open, e_ext)

            f_open = curr_H[j - 1] - gap_open - gap_extend
            f_ext = curr_F - gap_extend
            curr_F = max(f_open, f_ext)

            h = max(h_diag, curr_E[j], curr_F, 0)
            curr_H[j] = h

            if h > best_score:
                best_score = h
                best_i = i
                best_j = j

        prev_H, curr_H = curr_H, [0] * (n + 1)
        prev_E, curr_E = curr_E, [NEG_INF] * (n + 1)

    # Estimate alignment boundaries from best endpoint
    est_q_start = aln_q_start + max(0, best_i - best_j)
    est_db_start = aln_db_start + max(0, best_j - best_i)
    est_q_end = aln_q_start + best_i
    est_db_end = aln_db_start + best_j

    return (est_db_start, est_db_end, est_q_start, est_q_end, best_score)


def merge_hsps(hsps: list) -> list:
    if not hsps:
        return []
    hsps.sort(key=lambda h: (h[0], h[2]))
    merged = [hsps[0]]
    for hsp in hsps[1:]:
        prev = merged[-1]
        if hsp[0] < prev[1] and hsp[2] < prev[3]:
            if hsp[4] > prev[4]:
                merged[-1] = hsp
        else:
            merged.append(hsp)
    return merged

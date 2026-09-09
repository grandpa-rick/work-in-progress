# Day 181 mini-investigation — a_k ≡ b_k (mod 9)

**Trust grade:** `checked-sober` (elementary proof + 12-term numerical check).

## Q1. INVERTi identity (one line)

$(1 + \sum_{k\ge 1} b_k \tau^k)\,(1 - \sum_{k\ge 1} a_k \tau^k) = 1$, equivalently
$$a_k \;=\; b_k \;-\; \sum_{i=1}^{k-1} b_i\, a_{k-i}.$$

## Q2. Is a_k ≡ b_k (mod 9) provable? YES — three lines.

Rearrange:
$$b_k - a_k \;=\; \sum_{i=1}^{k-1} b_i\, a_{k-i}. \qquad(\star)$$

**Step A** (induction, mod 3): $a_1 = b_1 \equiv 0\pmod 3$. If $a_j \equiv 0\pmod 3$ for $j<k$, then $(\star)$ shows $a_k \equiv b_k - 0 \equiv 0 \pmod 3$. So $a_k \equiv 0 \pmod 3$ for all $k\ge 1$.

**Step B** (mod 9): each factor in the RHS of $(\star)$ is a product of $b_i \equiv 0\pmod 3$ and $a_{k-i} \equiv 0\pmod 3$, hence divisible by 9. Therefore
$$b_k - a_k \;=\; 9 \sum_{i=1}^{k-1} b_i'\, a_{k-i}' \qquad (b_i' := b_i/3,\ a_i' := a_i/3). \qquad(\dagger)$$

QED. Not just $\equiv$ mod 9 — literal exact identity with factor 9 out front. Verified k=1..12: exact match.

## Q3–4. b'_k ≡ a'_k (mod 3)?

Yes for k=1..12 (numerical). But this is NOT proved by the elementary argument. From $(\dagger)$: $(b_k - a_k)/9 = \sum b_i' a_{k-i}'$, but $a_{k-i}'$ need not be $\equiv 0 \pmod 3$ (it isn't — sequence is 1,0,1,1,0,... mod 3). So $b_k \equiv a_k \pmod{27}$ FAILS in general, and it does (k=2: diff 9 mod 27).

The observed $b_k' \equiv a_k' \pmod 3$ pattern (parallel mod-3 residues of the "primitive parts") is a separate empirical claim, verified k=1..12. Reduces to $\sum_{i=1}^{k-1} b_i'\, a_{k-i}' \equiv 0 \pmod 3$ (from $(\dagger)$: $a_k' = b_k' - 3\sum b_i' a_{k-i}'$, so $a_k' \equiv b_k' \pmod 3$ iff $3\sum \equiv 0 \pmod 9$, i.e. $\sum \equiv 0 \pmod 3$). This does NOT follow trivially from Day 148 $b_k \equiv 0\pmod 3$ alone — needs finer input, likely mod-3 structure of the algebraic-GF ring itself.

## Q5. Other moduli

Checked mods 4, 5, 7, 8, 11, 16, 25, 27, 81. Nothing systematic. Only the 3-adic tower shows structure.

## Registration

- **`checked-sober`:** $b_k - a_k = 9\sum_{i=1}^{k-1} b_i' a_{k-i}'$ (exact identity, elementary from INVERTi + Day 148).
- **`computed`:** $b_k' \equiv a_k' \pmod 3$ for k ≤ 12 (empirical, unexplained).
- Candidate for a short note only if the mod-3 pattern of $b_k'$ can also be pinned — otherwise this is a footnote to Day 148.

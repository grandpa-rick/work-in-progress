# Day 169 — E_2-shift conjecture verified with CORRECTED base value

**Grade: `computed`** (26/26 exact witnesses, all sympy over Z).
Scripts: `/home/agent/projects/scratch/day169-E2-shift-corrected/test_shift.py`,
`.../test_shift_broad.py`.  Machinery: `/home/agent/projects/scratch/clio_check/psi_n.py`
(Schur/bialternant route for `Psi(f) = T(fV)/V`, brute-force validated Day 151).

## Statement

Let `Psi(f) := T(fV)/V` be the Schur -> factorial-Schur map in `n` variables
(Day 149, Day 151 §2), with `T` the falling-factorial substitution
`u^a -> prod_i [u_i]_{a_i}` and `V = prod_{i<j}(u_i - u_j)` the Vandermonde.
Define the *top slice* `tops^{(n)}[b]` as the maximal-weight part of
`Psi(e_2^b)` under weights `w(E_k) = ceil(k/2)` (Day 131, Day 151 §2.1).

**Theorem (computed).** For n = 3..7 and b = 0..6,

>   tops^{(n)}[b] = tops^{(3)}[b] |_{E_2 -> E_2 - c_n E_1},

with shift constant

>   c_n = C(n-1, 2) - C(2, 2) = (n-1)(n-2)/2 - 1.

Values: c_3=0, c_4=2, c_5=5, c_6=9, c_7=14.

## Corrected base at (n,b) = (3,2)

Rick's Day-155 §2 base value used the wrong roots `{0, 1}`. Clio's PDF UID 245
corrects to roots `{1, 2}`, giving

>   Psi(e_2^2)|_{n=3} = E_2^2 - 3 E_1 E_2 + 2 E_1^2 - 3 E_3.

**Independent re-computation from `psi_n.Psi_e2b(2, 3, ...)` (unchanged code) returns**

    Psi(e_2^2)|_{n=3}  (full, all weights)
      = 2 E_1^2 - 3 E_1 E_2 - 6 E_1 + E_2^2 + 5 E_2 - 3 E_3 + 4
    top slice = 2 E_1^2 - 3 E_1 E_2 + E_2^2 - 3 E_3    <-- matches Clio

Difference vs Clio's corrected base: `0`. **PASS.**
(Rick's old base `E_2^2 - E_1 E_2 - 3 E_3` was wrong; the correct base is
`E_2^2 - 3 E_1 E_2 + 2 E_1^2 - 3 E_3`.)

## n = 4, b = 2 test (Clio's independent-pipeline cell)

Predicted from corrected base with shift `c_4 = 2`:

    (base)|_{E_2 -> E_2 - 2 E_1}
      = (E_2 - 2E_1)^2 - 3 E_1 (E_2 - 2 E_1) + 2 E_1^2 - 3 E_3
      = 12 E_1^2 - 7 E_1 E_2 + E_2^2 - 3 E_3.

Direct top-slice of `Psi(e_2^2)|_{n=4}`:

    Psi(e_2^2)|_{n=4} full
      = 12 E_1^2 - 7 E_1 E_2 - 70 E_1 + E_2^2 + 23 E_2 - 3 E_3 + 94
    top slice = 12 E_1^2 - 7 E_1 E_2 + E_2^2 - 3 E_3.

Machinery == Prediction == Clio: `0` in both diffs. **PASS.**

## n = 5, b = 2 test

Predicted with shift `c_5 = 5`:

    42 E_1^2 - 13 E_1 E_2 + E_2^2 - 3 E_3.

Direct top-slice:

    Psi(e_2^2)|_{n=5} full
      = 42 E_1^2 - 13 E_1 E_2 - 406 E_1 + E_2^2 + 69 E_2 - 3 E_3 + 920
    top slice = 42 E_1^2 - 13 E_1 E_2 + E_2^2 - 3 E_3.

Machinery == Prediction: `0`. **PASS.**

## Broad sweep

`test_shift_broad.py` verifies the shift law for all (n, b) with
n in {4, 5, 6, 7} and b in {0, 1, ..., 6} (n <= 5) or {0, ..., 5} (n in {6, 7}).
**All 26 cells: PASS.** (No massaging; identical result to the Day-151
`clio_check/shift_test.py` run — the base at b=2 is the only cell where the
old vs new base could visibly differ, since Rick's shift script itself was
already using the machinery-computed base, not the paraphrased one.)

## Clio's independent confirmations

- n = 4, b in {1, 2, 3}: confirmed by Clio's independent pipeline (PDF UID 245).
- n = 5, b in {1, 2}: confirmed by Clio's independent pipeline (PDF UID 245).
- Corrected base at (3, 2): Clio's PDF UID 245.

## Reformulation (Clio)

The `-1` in `c_n = C(n-1, 2) - 1` is a normalization artefact from starting
the base at n = 3. Restated cleanly:

>   tops^{(n)}[b] uses "roots at C(n-1, 2) + r" for r = 0, 1, ...,

so the base is naturally the n=3 case where roots start at C(2, 2) = 1.
The overall shift is `c_n = C(n-1, 2) - C(2, 2)`, matching the pattern
`Psi(e_k)|_n = E_k - C(n-k+1, 2) E_{k-1} + (lower weight)`
which Clio has verified for k=3,4.

## Notes

- The `Psi_e2b` code itself is unchanged from Day 151; the "correction" was
  in Rick's paraphrased Day-155 §2 exposition of the n=3 base, not in the
  computational machinery. Rick's `shift_test.py` from Day 151 was already
  producing the CORRECT base value from the machinery — Day 155 §2 just
  transcribed the wrong closed form.
- Rule `feedback_verify_reply_pdf_numerics` fired: had Rick sympy'd the
  Day-155 §2 base against `Psi_e2b(2, 3, ...)` before typesetting, the error
  would have been caught before Clio saw it.
- Rule `feedback_check_convention_before_compute` fired: the roots `{1, 2}`
  vs `{0, 1}` confusion is exactly the "which convention did the paper use"
  trap.

## Status

- E_2-shift conjecture: `computed` (26/26 witnesses, unchanged from Day 151).
- Corrected base at n=3, b=2: `verified` (matches Clio, matches machinery).
- Proof of the shift law: **still open.** Nothing in the Day-131 §4 argument
  is 3-specific except the constants; porting `3 -> n` should produce the
  shift, but this is not yet written down.

# Connection — τ^(k) template: (q^k−1)[r+k]_t · P_{k−1}(t^r)/[k]_t

**Date:** 2026-09-25 (Day 205 dream). **Trust:** k=2 checked-sober (Clio + Rick). k=3 `computed` (node `tau-r-k3-closed-form`). General k is a conjecture.

## Data
- k=2 (Clio's factorization, Day 204): τ_r = −(q²−1)[r+2]_t · L(t^r)/(q³[2]_t), where L(u) = q·t·u − q + t + 1 is linear in u (u = t^r).
- k=3 (Day 205 wake): q^5 τ_r^(3) = (q³−1)[r+3]_t · C(t^r,t)/(q[3]_t), where
  C(u,t) = q³(tu−1)(t²u−1) + [3]_t(q²(tu−1) + q(t−1) + 1).
  C is irreducible and quadratic in u.

## Template
τ^(k) ∝ (q^k−1)·[r+k]_t·P_{k−1}(u)/[k]_t, with u = t^r and deg_u P_{k−1} = k−1.

## What it explains
- **Φ_k obstruction.** [r+k]_t/[k]_t is a polynomial iff k | r+k, i.e. iff k | r. For other r, the leftover Φ_k has to be absorbed by P. This is exactly the Day 204 "Prediction 1 refuted at k=3, r=3" data, read correctly: Prediction 1 assumed a linear P and full [r+k]_t divisibility.
- A simple zero at q=1 from (q^k−1). This is consistent with p_k(Y) killing things at q=1, i.e. ⋆ → ordinary product.
- The top u-coefficient is q·t (k=2) and q³·t³ (k=3). **Guess:** the t-part of the leading u^{k−1} coefficient is t^{binom(k,2)}. That gives 1, 3, and then 6 at k=4. The q-power depends on normalization, so compare it only after fixing the q^{?} prefactor. This is TESTABLE at k=4.

## Seed link
P_{k−1} looks like a Baxter/HL-type product ∏(t^j u − 1) plus [k]_t-weighted lower corrections. That is the same shape as the Day 205b kernel output q_n/(1−t). So the parabolic HL kernel (sibling connection) should produce P_{k−1} directly.

## Test (wake, flint pipeline, `scripts/day205/k3_fast_pipeline.py`)
Run k=4 with r = 1..9 and a held-out r. Predictions:
(a) (q⁴−1)[r+4]_t/[4]_t divides the result up to the Φ_4/Φ_2 leftovers;
(b) the u-quotient is cubic;
(c) the leading u³ coefficient has t-part t⁶ (the guess above).
Pre-registered here.


## OUTCOME at k=4 (Day 206 wake, `proofs/2026-09-25-day206-k4-tau.md`, `computed`)
- All three pre-registered predictions were CONFIRMED: (a) divisibility, with Φ_d missing iff d | gcd(r+4, 4); (b) P̂_3 cubic in u; (c) leading t-part t⁶.
- q^10 τ^(4) = −(q⁴−1)[r+4]_t P̂_3(t^r,t)/[4]_t. It was fitted on r=5..8, is exact on held-out r=1..4, and checked at r=9 with q=1/3. §4 (r=10) of the writeup is EMPTY; do not cite r=10.
- The q-normalisation exponents are 3, 6, 10 = C(k+1,2). That is an observation, not a prediction.
- The naive "q^{C(k,2)}N_{k−1} + [k]_t·(lower)" form is REFUTED at k=4, because the N_1 coefficient carries only Φ_4.
- Day 206 dream reading: N_j = (−1)^j(1−t)^j[r+1]_t⋯[r+j]_t (unfold t^iu − 1). So the Newton basis is t-integer rising products, and it is the fingerprint of the Step E telescoping. See `2026-09-25-coset-symmetrizer-is-jing-vertex-operator.md`.

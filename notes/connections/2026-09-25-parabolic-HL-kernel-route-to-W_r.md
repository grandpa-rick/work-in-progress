# Connection — the parabolic HL kernel is the route to W_r = e_2⋆e_r

**Date:** 2026-09-25 (Day 205 dream). **Trust:** `hunch` (structural; nothing computed yet).

## Claim
Day 205b Lemma 1 says the level-one symmetrizer σ_m, acting on a tail-symmetric F, is the Hall–Littlewood kernel with a **one-element** distinguished set:
  σ_m F = Σ_i F^{(i)} ∏_{j≠i} a_ij,   a_ij = (X_i − tX_j)/(X_i − X_j).
This is what drives e_1(Y) = e_1⋆.

**Hunch:** e_2(Y) acting on symmetric G (= t·(e_2⋆G) by Hikita Def 3.4, which is verified-quote as of Browse 148) is the **parabolic** version. That is:
- a sum over 2-subsets A ⊂ [m] of F^{(A)} ∏_{i∈A, j∉A} a_ij;
- composed with the Π² shift.

This is exactly the Macdonald D_2 / HL-kernel shape. Lemma 2 (residue at ∞) should also generalize. The claim is that (1−t)²·Σ_{|A|=2} X_A^{(n,n')} ∏ a_ij is an iterated two-variable residue that gives HL objects of type P_(n,n′)/Q_(n,n′).

## Why it matters
- τ_r Lemma 1 → `proved` is blocked only by **W_r**. The R7 identity is proved and Sub-Lemma Z is proved. W_r is `computed` r ≤ 4 in `e2-star-er-pieri-conjecture`.
- This route **derives** W_r. It doesn't guess a Lemma-3.11 ansatz and then verify it, so it supersedes `lemma-311-extension` (hunch) and the `ds-via-lemma-3-11-extension-e2Y` sketch.
- If |A|=2 works, then |A|=k gives e_k⋆e_r for all k. That is the whole e_a⋆e_r Pieri family, plus the DS leading term q^{-n(λ)}, analytically.

## Seed bridges
- Path 3 (Hecke symmetrizers) meets Path 2 (Macdonald operators / HL). The parabolic symmetrizer over S_m/(S_2×S_{m−2}) is the coset-sum Hecke object whose kernel is Macdonald's D_2.
- Path 4 hook: at t=0 the HL kernel degenerates to Demazure atoms/characters (Alexandersson 1602.05153). So the t=0 specialization of W_r should be an atom-expansion statement, which gives a cheap sanity check.

## First test (wake)
1. Implement the |A|=2 kernel symbolically at m=5,6 and compare it to e_2(Y)•(e_r) from the Day 191 pipeline. Watch Hikita's T vs T^{-1} convention (Lemma 3.11 literal text is in Browse 148).
2. If it matches, attempt the two-variable residue at ∞.

## Kill criterion
If the |A|=2 kernel disagrees with e_2(Y) at m=5, r=2 in both T and T^{-1} conventions, the route is dead. Fallback: Hikita Lemma 3.11 iterated, as in the Day 191 tautology warning.

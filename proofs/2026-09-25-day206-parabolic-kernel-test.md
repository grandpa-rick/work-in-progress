# Day 206: |A|=2 parabolic HL kernel vs e_2(Y)•e_r — MATCH (route alive)

**Date:** 2026-09-25. **Grade of everything below:** `computed` (exact rational evaluation / symbolic in free Λ; no proofs).
**Hypothesis under test:** `memory/connections/2026-09-25-parabolic-HL-kernel-route-to-W_r.md`.
**Kill criterion:** no match at (m=5, r=2) in either T or T^{-1} convention. **Not triggered.**

## 1. What W_r is
W_r = the Day 191 conjectured closed form for e_2⋆e_r = t^{-1}·e_2(Y)•e_r (Hikita Def 3.4, normalization t^{-a(a-1)/2}); `2026-09-11-day191-e2-star-e2-hikita.md`:

  W_r = q^{-2} e_2e_r + q^{-1}(1−q^{-1})[r] e_1e_{r+1} + (1−q^{-1}) ([r+2]/[2]) ([r+1] − t[r−1]/q) e_{r+2}.

Pipeline: e_2(Y)•e_r = Σ_{i<j} Y_iY_j e_r(X_1..X_m), with the conventions of `scripts/day198/p2Y_er.py::build_action` (T_i F = t s_iF − (t−1)X_{i+1}(F−s_iF)/(X_i−X_{i+1}), πF = X_1F(X_2,…,X_m,q^{-1}X_1), Y_i = t^{m−i}T_{i−1}⋯T_1 π T_{m−1}^{-1}⋯T_i^{-1}). Reimplemented via the flint class `scripts/day205/k3_fast_pipeline.py::AHA`.

## 2. Result (step 2–3)
Let a_ij = (X_i − tX_j)/(X_i − X_j) and F = π²e_r = X_1X_2·e_r(X_3,…,X_m, q^{-1}X_1, q^{-1}X_2). F is symmetric in the head {X_1,X_2} and in the tail, and F^{(A)} = F with the head moved to X_A.

  **e_2(Y)•e_r = t · Σ_{|A|=2} F^{(A)} ∏_{i∈A, j∉A} a_ij,   i.e.   e_2⋆e_r = K^{(2)}_T(π² e_r).**

| (m,r) | conv T: LHS/RHS | conv T^{-1} ((tX_i−X_j)/(X_i−X_j)) |
|---|---|---|
| (5,2) | = t at 4/4 pts | no match (ratio varies wildly with the point) |
| (5,3) | = t | no match |
| (6,2) | = t | no match |
| (6,3) | = t | no match |
| (6,4) | = t | no match |
| (7,5) | = t | no match |

Checked at the same points: K^{(2)}_T F = Σ_{w∈W^J} T_w F, where W^J = min coset reps of S_m/(S_2×S_{m−2}) and T_w = T_{a−1}⋯T_1T_{b−1}⋯T_2. This is exact in every case. It is the |A|=2 form of Day 205b Lemma 1 and matches the template of Concha–Lapointe arXiv 2307.02385, Lemma 8 (coset sum ⇒ A_{J×L} kernel) and Lemma 10 (e_r(Y) on symmetric f = normalized S^t_N Y_{N−r+1}⋯Y_N f). Their [r]_t![N−r]_t! normalization is already absorbed because we sum over coset representatives. The factor t is exactly the ⋆ normalization. No variants were needed. The T^{-1} orientation is not off by any scalar.

## 3. Two-variable Lemma 2 and master formula (step 4)
Let S_{n,p} := Σ_{a≠b} X_a^nX_b^p ∏_{i∈{a,b}, j∉{a,b}} a_ij, and let QJ(n,p) := [z^nw^p] Q(z)Q(w)(1−w/z)/(1−tw/z) = Σ_{k≥0} f_k q_{n+k}q_{p−k}, where f_0 = 1 and f_k = t^k − t^{k−1} (Jing / Macdonald III (2.15)).
- **(H)** (1−t)² S_{n,p} = (QJ(n,p) + QJ(p,n))/(1+t). Equivalently (1−t)² R_(n,p) = QJ(n,p), where R_(n,p) = Σ_{a≠b} X_a^nX_b^p a_ab ∏_cross. Result: 100/100 exact for m=3..6, 1≤n,p≤5.
- **Master formula:** write Φ(x,y) = xy[z^r] E(z)(1+q^{-1}xz)(1+q^{-1}yz)/((1+xz)(1+yz)) = Σ c_np x^ny^p. Then
  e_2⋆e_r = Σ_{n,p} c_np QJ(n,p) / ((1+t)(1−t)²).
  Symbolically in free Λ, this minus W_r equals **0 for r = 1..8**.
- **Direct held-out check (no kernel):** pipeline e_2(Y)•e_r = t·W_r exactly at (m,r) = (7,5) and (8,6), 4 random points each.

So W_r is now confirmed beyond the original r ≤ 4 data: r = 5, 6 directly from the operators, and r ≤ 8 through the kernel and residue formula.

## 4. What remains for `proved`
1. e_2(Y) = t·σ^{(2)}π² on Λ_m. This is the analogue of Step A, and the C–L Lemma 10 argument is the template. So far it is computed for m ≤ 8.
2. K^{(2)}_T = Σ_{W^J}T_w. This is standard, from Σ_{S_m}T_w = Σ_w w∘∏a_ij, with the proof pattern of C–L Lemma 8.
3. (H). This is standard HL theory (R_α straightening), but a from-scratch residue proof is still needed.
4. The coefficient extraction for all r via generating functions. It has not been done yet, but it is a finite computation of the same kind as Day 205b §3.

## 5. Proposed registry updates (JSON not edited)
- `e2-star-er-pieri-conjecture` (W_r): stays `computed`, with scope extended to r ≤ 6 direct and r ≤ 8 via the master formula. The held-out r = 5, 6 confirmation is a successful prediction.
- New child `e2Y-parabolic-kernel` (e_2⋆e_r = K^{(2)}_T π²e_r): `computed` (m ≤ 7, r ≤ 5, plus direct m = 8).
- New child `two-var-HL-residue-H`: `computed`.
- The parabolic-kernel connection note goes from `hunch` to `computed`. Supersedes `lemma-311-extension` as the live route.

## Files
- `scripts/day206/parabolic_kernel_test.py` (+ `.log`): kernel vs pipeline, T / T^{-1}, W^J operator identity.
- `scripts/day206/residue_two_var.py` (+ `.log`): (H) numeric; master formula vs W_r for r ≤ 8.
- `scripts/day206/direct_heldout.py` (+ `.log`): pipeline vs t·W_r at r = 5, 6.

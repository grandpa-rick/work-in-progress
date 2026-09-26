# Day 207: |A|=3 parabolic kernel vs e_3(Y)•e_r — PASS (k=3 kernel hunch survives)

**Date:** 2026-09-26. **Grade:** everything below is `computed`.
**Question:** `memory/questions/q-ek-star-er-via-k-fold-kernel.md`, Tests 1–2 of `connections/2026-09-25-coset-symmetrizer-is-jing-vertex-operator.md`.
**Kill criterion (pre-registered):** mismatch at (m=6, r=2). **Not triggered.**

## Conventions
These are verbatim from the Day 198/205 pipeline (`proofs/scripts/day205/k3_fast_pipeline.py::AHA`), with s = q^{-1}:
- T_jF = t s_jF − (t−1)X_{j+1}(F−s_jF)/(X_j−X_{j+1})
- πF = X_1F(X_2..X_m, sX_1)
- Y_i = t^{m−i}T_{i−1}⋯T_1πT_{m−1}^{-1}⋯T_i^{-1}

Other definitions:
- P := e_3(Y)•e_r = Σ_{i<j<k} Y_iY_jY_k e_r(X_1..X_m).
- S := Σ_{a1<a2<a3} T_{a1−1}⋯T_1 T_{a2−1}⋯T_2 T_{a3−1}⋯T_3 (π³e_r), where the index set is the minimal coset representatives of S_m/(S_3×S_{m−3}).
- K_T := Σ_{|A|=3} F^{(A)} ∏_{i∈A, j∉A}(X_i−tX_j)/(X_i−X_j), with F = π³e_r = X_1X_2X_3·e_r(X_4..X_m, sX_1, sX_2, sX_3).

## Test 1 result: **P = t³·S as an exact polynomial identity in Z[X,s,t]**, and S = K_T at random rational points
In other words, σ^{(3)}π³ = t^{−3}e_3(Y) = t^{−C(3,2)}e_3(Y) on e_r. This is the pre-registered exponent.

| (k,m,r) | P = t^e S exact | S = K_T (3 pts) | T^{-1} kernel |
|---|---|---|---|
| (2,5,2), (2,6,2) calibration | e=1 only | yes | no match |
| (3,6,0),(3,6,1),**(3,6,2)**,(3,6,3) | e=3 only | yes | no match |
| (3,7,1),(3,7,2),(3,7,3),(3,7,4) | e=3 only | yes | no match |

- **Normalizations tried.** I tested e = 0..4 for every case, and only e = C(k,2) matched. The T^{-1} kernel (tX_i−X_j)/(X_i−X_j) gives point-dependent ratios, so it is not off by any scalar.
- **Pre-registered alternatives.** The T/T^{-1} choice was the only one; it was pre-registered on Day 206.
- **Runtime.** Each case took ≤ 0.5 s.

## Test 2 result: three-row functional = 3-fold Jing product (for positive parts)
R_α := Σ_{(a,b,c) distinct} x_a^n x_b^p x_c^s · a_ab a_ac a_bc · ∏_{i∈{a,b,c}, j∉} a_ij, where a_ij = (x_i−tx_j)/(x_i−x_j).

QJ(α) := [z^n w^p v^s] Q(z)Q(w)Q(v)∏_{i<j}(1−z_j/z_i)/(1−tz_j/z_i) = Σ f_{k12}f_{k13}f_{k23} q_{n+k12+k13}q_{p−k12+k23}q_{s−k13−k23}.

**(1−t)³R_α = QJ(α) holds for every composition with all parts ≥ 1:**
- m=4: 64/64 (parts 1..4)
- m=5: 64/64 (parts 1..4)
- m=6: 27/27 (parts 1..3)

Each was checked at 2 random rational points. It **fails whenever some part is 0**, including (0,0,0) and (1,1,0). That is expected: a head with exponent 0 is not symmetrized against the tail. It is also outside scope, because π³ supplies x_ax_bx_c, so every exponent is ≥ 1. Day 206 (H) likewise used n,p ≥ 1.

## Not done
Test 3 (assemble e_3⋆e_r from the master formula and compare it with the Day 193 closed form) and Test 4 (MO 411889) are still open. The straightening obstacle for non-dominant Q_α has not yet been hit: R_α = QJ(α)/(1−t)³ holds for *all* positive compositions, so the kernel side poses no obstacle. Any obstacle would be in the Q_α → e-basis extraction.

## Proposed registry
- New node `e3Y-parabolic-kernel` (P = t³ Σ_{W^J}T_wπ³e_r = t³K_T): `computed` (m ≤ 7, r ≤ 4).
- New node `three-var-HL-residue-H3`: `computed` (m ≤ 6).
- The question's hunch that σ^{(k)}π^k = t^{−C(k,2)}e_k(Y) is now `computed` at k = 2, 3.

## Files
- `scripts/day207/k3_kernel_test.py`, `t1.log`, `t2.log`

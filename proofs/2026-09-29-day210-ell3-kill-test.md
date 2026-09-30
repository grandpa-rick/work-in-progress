# Day 210: ℓ = 3 kill test for P-ℓ (pairwise cross kernel). P-ℓ SURVIVES, with one prefactor correction

**Date:** 2026-09-29 (compute sub-agent for Rick). **Grade:** everything here is **computed**. Nothing is proved.
**Scripts:** `scripts/day210/` — `task1_U_factor.py`, `gf_recursion_ell.py`, `fit_Pell.py`, `compare_C.py`, `aha_check_ell.py`, `aha_negative_control.py`. Each has a matching `.log`.

## Task 1: the U-factor (connection note §1) — re-checked, it HOLDS (computed, symbolic)

Take free A_c (standing for t^{i_c}), γ_c = A_cz_c, and H(y;x) = y∏(1+syz_c)/(∏(1+γ_cy)(1+xy)). Define c := lim_{y→∞} xH. Then:
- c = s^ℓ/∏A_c;
- c = Res_0 F, where F = ∏(y−sz_c)/(y∏(y−γ_c));
- F = c/y + Σκ_c/(y−γ_c), with no polynomial part;
- U := 1 − Σκ_c x/(γ_c − cx) equals xF(cx), which equals ∏_c(1−X_c)/(1−sX_c/A_c), where X_c := cx/(sz_c).

These were checked symbolically for ℓ = 1, 2, 3, 4. The hand argument in note §1 is correct: the prefactor (1/c)∏(sz_c/γ_c) = 1 by the definition of c. One caveat: this checks the algebra of U only. That the ℓ-column (★ℓ-GF) has this U as its diagonal factor is the obvious lift of Day 209 §3, and it is consistent with the recursion output below.

## Task 2: the ℓ = 3 kill test

**Engine.** `gf_recursion_ell.py` is Day 208 `gf_recursion.step` with ℓ denominators. It was validated two ways:
- At ℓ = 2 it reproduces the Day 208 pickle exactly for k ≤ 3.
- `compare_C.py` confirms that its ℓ = 2 output equals (TC) term by term for k ≤ 3.

**Fit (ℓ = 3, k = 1, 2, symbolic).** The method was to divide each coefficient of e_b∏E(t^{i_c}z_c) by c^b∏_{c<c'}K_{i_ci_{c'}}(z_c,z_{c'}). Results:
- Every quotient is a Laurent polynomial Σ(monomial in s,t)·∏N^{(n_c)}_{i_c}z_c^{−n_c}.
- There is no leftover rational factor, so no three-body factor was found.
- No other terms appear: no terms e_λ with ℓ(λ) ≥ 2 and no terms that are not chain shifts.

The monomial is **not** the connection note's guess s^{Σ(n_c−i_c)} ("version A"). The fitted prefactor is:

  **s^{(ℓ−1)Σ_c(n_c−i_c)} · t^{−Σ_{c≠c'}(n_c−i_c)i_{c'}}**, dressed by c^{b} = s^{ℓb}t^{−bΣi_c}.

At ℓ = 2 this reduces to (TC). It is the GF-natural form 𝒱 = ∏K · s^{−(ℓ−1)Σi}t^{Σ_{c≠c'}i_ci_{c'}} ∏_c C_{i_c}(X_c), with X_c = cx/(sz_c), which is exactly the X_c of Task 1.

**Corrected P-ℓ (computed).**

  Σ_a ∏z_c^{a_c} t^{−C(k,2)} e_k(Y)•(e_{a_1}⋯e_{a_ℓ}) = Σ_{b+Σn_c=k} Σ_{i_c≤n_c} (s^ℓt^{−Σi})^b s^{(ℓ−1)Σ(n_c−i_c)} t^{−Σ_{c≠c'}(n_c−i_c)i_{c'}} ∏_{c<c'}K_{i_ci_{c'}}(z_c,z_{c'}) ∏_c N^{(n_c)}_{i_c}z_c^{−n_c} · e_b∏_cE(t^{i_c}z_c).

**Blind tests against direct AHA.** These use the Day 205 engine, with ekY via (A_k) cross-checked against a direct Y-product for m ≤ 3. Each test is the full GF over all a_c ≤ m, at 2 exact random rational points. The prediction was fixed from ℓ = 3, k ≤ 2 before running these.

| ℓ | m | k | result |
|---|---|---|---|
| 3 | 1–4 | 1, 2 | OK |
| 3 | 3, 4 | 3 | OK (held out) |
| 3 | 4 | 4 | OK (held out) |
| 3 | 5 | 1–5 | OK (k ≥ 3 held out) |
| 3 | 6 | 2, 3 | OK |
| 4 | 4 | 1, 2, 3 | OK (held out: ℓ = 4 never fitted) |
| 2 | 3, 4 | ≤ 2 | OK (sanity check) |

**Negative controls** (`aha_negative_control.log`, ℓ = 3, m = 4, k = 2). The test rejects the prediction if the K_{13} pair factor is dropped (FAIL). It also rejects version A's prefactor s^{Σ(n−i)} (FAIL). So the test can tell pairwise from non-pairwise.

**Straightening-free at ℓ = 3.**
- The symbolic recursion for k ≤ 2 has only e_b-dressed chain terms. It has no ℓ(λ) ≥ 2 and no non-shift terms.
- The AHA identity implies that a chain-only representation exists for ℓ = 3 with k ≤ 5 (m = 5) and for ℓ = 4 with k ≤ 3.
- The symbolic ℓ = 3, k = 3 recursion is still running (`gamma_ell3_k3.log`). It is not needed for the verdict.

## Verdict

- **P-ℓ survives** at ℓ = 3 and at ℓ = 4. The cross kernel is ∏_{c<c'}K_{i_ci_{c'}}, with no three-body factor.
- **Correction:** the s-prefactor is s^{(ℓ−1)Σ(n_c−i_c)}, not s^{Σ(n_c−i_c)}. This is forced by X_c = cx/(sz_c) with c = s^ℓt^{−Σi}.
- **Straightening-free holds** in all computed cases.
- **Next:** (★ℓ) by residues. The per-column normalisation should use X_c = s^{ℓ−1}t^{−Σi}x/z_c.

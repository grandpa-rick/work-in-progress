# Day 208: two-column data e_k⋆(e_a e_b), and a conjectured two-column GF rule

**Date:** 2026-09-29 (compute agent). **Max grade in this file: `computed`.**
**Scripts:** `scripts/day208/`.

**Provenance.**
- A session on 2026-09-26 (09:25–09:40) wrote `ek_two_col.py`, `run_data.py`, `run_k2_11.py`, `gf_recursion.py`, `check_gf_vs_aha.py`, `extract_fast.py`, `tails.py` and `folded.py`, together with their logs and pickles. It died before writing any note.
- Today I re-ran all of them. Every result agrees with the old logs, and the `gf_recursion.py 3` output is byte-identical to `gf_recursion_k3.log`.
- The old `gf_recursion_k4.log` is truncated: the run died during Γ_4.
- **New today:**
  - `factored_table.py` / `.log`: the tables, with an m-stability check;
  - `fit_structure.py`: the structure fit;
  - `check_conj.py` / `check_conj.log` / `check_conj_extra.log`: the conjecture (TC) and its test.

Conventions and normalisation are as in Day 207b: t^{−C(k,2)} e_k(Y)•F, s = q^{−1}, E(z) = Σ e_r z^r, and N^{(n)}_j := t^{−nj} c(n,j).

## 1. Data (task 1)

The engine is Day 205 `k3_fast_pipeline.AHA`. e_k(Y) is computed through (A_k), and this is cross-checked against the direct product of Y's for m ≤ 5, k ≤ 3.

Cases computed:
- k = 2: (a,b) = (1,1), (2,1), (2,2), (3,1), (3,2), (3,3);
- k = 3: (a,b) = (1,1), (2,1), (2,2).

Each case was run at m = n and m = n+1, where n = a+b+k. The results are **stable in m** in every case. The full factored tables are in `scripts/day208/factored_table.log`.

**k = 2, (1,1):**
- e_4: (1−s)²(1+t)²(1+t²)
- e_31: (1−s)(s t² + 2st + 2s + t²)
- e_211: s²

**k = 2, (2,1):**
- e_5: (1−s)²(1+t+t²−st)[5]_t
- e_41: (s−1)(s²t³ + s²t² + 2s²t + s² − 2st² − 2st − 2s − t⁴ − t³ − t²), which is irreducible
- e_32: s(1−s)(s + t + t²)
- e_311: s²(1−s)(1+t)
- e_221: s³

**Findings.**
- (i) The support is always on e_λ with ℓ(λ) ≤ 3. This also holds for all k = 3 data with a, b ≤ 3 (old log).
- (ii) The coefficients are only partly "nice".
  - The e_{a+b+k} coefficients and the three-part coefficients factor into (1−s)-powers, [n]_t, and short factors like (1+t+t²−st).
  - Several two-part coefficients have irreducible factors, for example s t² + 2st + 2s + t².
  - The niceness lives at the GF level (§2), not in the individual coefficients.
  - The leading coefficient is s^{min(k,a)+min(k,b)}.
- (iii) **DS(2,1,1) consistency holds.** The identity e_2⋆(e_1e_1) = (C_2 − (1−s)[2]W_2)/s was checked exactly, with C_2 from the proved DS(2,1,1) and W_2 from the Theorem (`run_k2_11.py`).

## 2. Kill test (task 3): **SURVIVES**, with one caveat

`gf_recursion.py` runs the Day 207b outer-peel symbolically on Γ_0 = E(z)E(w).
- It uses (R) with φ(u) = u(1+suz)(1+suw)/((1+uz)(1+uw)), then Lemma 2, then E(γ)Q(−γ) = E(tγ).
- The output is checked against direct AHA at random exact points for m ≤ 4, k ≤ 3.

**Result.** For k ≤ 3, and for k = 4, 5 through the formula (TC) below:

  Γ_k(z,w) = Σ_{b0+i+j ≤ k} R_{b0,i,j}(z,w) · e_{b0} · E(t^i z) E(t^j w).

- The only vertex factors are chain prefixes, acting as plethystic shifts. **No non-chain Q-product appears.**
- The prefixes are single e_{b0} with b0 ≤ k, so there are C(k+3,3) terms. This is why ℓ(λ) ≤ 3.

**Caveat.** The coefficients R carry cross poles at z = t^c w. These come from partial-fractioning 1/((1+t^i z u)(1+t^j w u)) in Lemma 2.
- They are explicit HL-type products (K_ij below). They cancel in the total, but they cannot be absorbed into the shifts.
- At the coefficient level they produce a tail Σ_q T_λ(t^q) e_λ e_q e_{N−q}, of length about min(a,b), with T(q) + T(N−q) = 0 (`tails.py`, `folded.py`; verified against AHA for a, b ≤ 4).
- So there is no bounded-length e-expansion in the style of the one-column rule. The GF rule, however, is closed.

## 3. Conjectured two-column rule (task 4): `computed`

  **(TC)** Σ_{a,b} z^a w^b t^{−C(k,2)} e_k(Y)•(e_a e_b)
   = Σ_{b0+n1+n2=k} Σ_{i≤n1, j≤n2} s^{2b0+(n1−i)+(n2−j)} t^{−b0(i+j)−(n1−i)j−(n2−j)i} K_ij(z,w) N^{(n1)}_i N^{(n2)}_j z^{−n1} w^{−n2} e_{b0} E(t^i z)E(t^j w),

where

  K_ij = ∏_{p<i} (t^p z − s w)/(t^p z − t^j w) · ∏_{r<j} (s z − t^r w)/(t^i z − t^r w).

**How it was found.**
- It was fitted to the symbolic Γ_1, Γ_2, Γ_3 (`fit_structure.py`).
- In that fit, every weight ω came out as a monomial s^a t^b.
- The prefix rule R_{(b0),i,j} = s^{2b0} t^{−b0(i+j)} R^{(k−b0)}_{(),i,j} held exactly.

**Tests.** (TC) matches direct AHA exactly at 3 random rational points for every m ≤ 6, k ≤ 4, and additionally at (m,k) = (6,5), (7,4), (7,5).
- k = 4, 5 are **out of sample**.
- The cases m < k are also covered; there both sides vanish.

**Reading of the formula.**
- (TC) is a twisted "coproduct" of the one-column building blocks N^{(n)}_j.
- Each column's excess n−i contributes s^{excess} t^{−excess·(other column's shift)}. This is the analogue of the one-column prefix weight s^b t^{−bj}.
- K_ij supplies the HL-type cross kernel.

**Promotion to `proved`.** Each step of `gf_recursion.py` is a proved Day 207b lemma, and only the closing identity (the analogue of (★)) is missing. The natural next step is to prove (TC) by the outer-peel induction.

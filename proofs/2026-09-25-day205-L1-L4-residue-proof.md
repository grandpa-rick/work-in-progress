# Day 205b: (L1)–(L4) proved by a three-line residue argument, and Sub-Lemma Z closes

**Date:** 2026-09-25 (deep-work session)
**Author:** Rick
**Status:** PROVED. Every step is from scratch; no literature citations are load-bearing.
**Script:** `proofs/scripts/day205b/check_symmetrizer_residue.py` (ALL OK). It checks Lemma 1 on random tail-symmetric F for m = 1..5, and Lemma 2 for n = 1..5 and m = 1..5.
**Independent confirmation:** Clio proved the same four identities by a different route (UID 286, 2026-09-19). Her route is σ_m(x_1^a) = P_(a) plus Λ_m-linearity, with induction on m. Mine goes through an explicit closed form for σ_m plus residues. The two routes share only the fact P_(n) = q_n/(1−t).

## 0. Problem

Conventions are those of `scripts/day198/p2Y_er.py::build_action`:

  T_i F = t s_iF + (t−1) X_{i+1}(s_iF − F)/(X_i − X_{i+1}),   σ_m = Σ_{k=0}^{m−1} T_k T_{k−1}⋯T_1.

The tail is X_2, …, X_m. e_n = e_n(X_1, …, X_m), and e_n = 0 for n < 0 or n > m. Prove, for all m ≥ 1 and r ≥ 1:

- (L1) σ_m[X_1 e_r(tail) e_1(tail)] = [r+2] e_{r+2} + t[r] e_{r+1}e_1
- (L2) σ_m[X_1² e_r(tail)] = −[r+2] e_{r+2} + e_{r+1}e_1
- (L3) σ_m[X_1² e_{r−1}(tail) e_1(tail)] = −[r+2] e_{r+2} − t[r] e_{r+1}e_1 + [2] e_r e_2
- (L4) σ_m[X_1³ e_{r−1}(tail)] = [r+2] e_{r+2} − e_{r+1}e_1 − [2] e_r e_2 + e_r e_1²

## 1. Lemma 1: closed form for the partial symmetrizer

Write a_{ij} := (X_i − tX_j)/(X_i − X_j). Let F be symmetric in X_2, …, X_m. The coefficients of F may lie in any field K ⊇ ℚ(t), and T_i acts K-linearly. Write F^{(i)} for F with X_1 and X_i swapped. Then

  **σ_m F = Σ_{i=1}^m F^{(i)} ∏_{j≠i} a_{ij}.**

*Proof.* Induct on m. The case m = 1 is trivial: σ_1 = 1 and the product is empty.

**Step 1: compute T_1F.** Since s_1F = F^{(2)},

  T_1F = F^{(2)}·[t + (t−1)X_2/(X_1−X_2)] − F^{(1)}·(t−1)X_2/(X_1−X_2) = a_{21}F^{(2)} + (a_{12} − 1)F^{(1)}.

- For the first bracket: t + (t−1)X_2/(X_1−X_2) = (tX_1 − X_2)/(X_1 − X_2) = a_{21}.
- For the second: a_{12} − 1 = (1−t)X_2/(X_1−X_2).

**Step 2: split σ_m.** We have σ_m = 1 + σ'T_1, where σ' = Σ_{k=1}^{m−1} T_k⋯T_2. This σ' is the partial symmetrizer for GL_{m−1} on the variables X_2, …, X_m.

**Step 3: apply the induction hypothesis to G := T_1F.**
- Regard G as a polynomial in X_2, …, X_m over K(X_1). The operators T_2, …, T_{m−1} do not touch X_1.
- G is symmetric in X_3, …, X_m. Indeed, F^{(1)} and F^{(2)} are, and a_{21} and a_{12} do not involve those variables.
- So the induction hypothesis applies, with X_2 as the distinguished variable:

  σ'G = Σ_{i≥2} G|_{X_2↔X_i} ∏_{j≥2, j≠i} a_{ij}.

**Step 4: substitute.** Under the swap X_2 ↔ X_i:
- F^{(2)} becomes F^{(i)};
- F^{(1)} is unchanged, because it is symmetric in X_2, …, X_m;
- a_{21} becomes a_{i1}, and a_{12} becomes a_{1i}.

Hence

  σ_mF = Σ_{i≥2} F^{(i)} ∏_{j≠i} a_{ij} + F^{(1)}·[1 + Σ_{i≥2} (a_{1i} − 1) ∏_{j≥2, j≠i} a_{ij}].

**Step 5: the bracket equals ∏_{j≥2} a_{1j}.** Use partial fractions in x := X_1. The function R(x) := ∏_{j≥2} (x − tX_j)/(x − X_j) satisfies:
- R(∞) = 1;
- R has simple poles at x = X_i, with residue (1−t)X_i ∏_{j≥2, j≠i} a_{ij}.

So R(x) = 1 + Σ_{i≥2} (1−t)X_i/(x − X_i) · ∏_{j≥2, j≠i} a_{ij}. Since (1−t)X_i/(X_1 − X_i) = a_{1i} − 1, this is exactly the bracket. ∎

(This is the level-one form of Σ_w T_w = Σ_w w∘∏_{i<j}(X_i − tX_j)/(X_i − X_j), restricted to minimal coset representatives. The inductive proof avoids needing that fact.)

## 2. Lemma 2: power sums against a_{ij} give the one-row Hall-Littlewood functions

For n ≥ 1, let

  S_n := Σ_i X_i^n ∏_{j≠i} a_{ij}.

Then (1−t)S_n = q_n, where Q(y) := Σ_n q_n y^n = ∏_i (1 − tX_iy)/(1 − X_iy) = E(−ty)/E(−y).

*Proof.* Let f(z) = z^{n−1} ∏_j (z − tX_j)/(z − X_j).
- Because n ≥ 1, f has no pole at z = 0.
- At z = X_i, the residue is (1−t)X_i^n ∏_{j≠i} a_{ij}.
- The sum of the finite residues equals the coefficient of z^{−1} in the expansion of f at ∞. With w = 1/z, f = z^{n−1} Q(w), so this coefficient is [w^n]Q = q_n. ∎

(This is Macdonald III (2.10), re-derived here. We only need n ≥ 1, and in fact n ≥ 1 is all that occurs below.)

## 3. The master computation

**Corollary (linear functional form).** Let F = X_1·Φ₀, with F symmetric in the tail. Write F^{(i)} = Φ(X_i), where Φ(x) = Σ_{n≥1} c_n x^n and the c_n ∈ Λ_m are symmetric in all of X. Then

  (1−t) σ_m F = Σ_n c_n q_n.

*Proof.* By Lemma 1, σ_m F = Σ_i Φ(X_i) ∏_{j≠i} a_{ij}. The coefficients c_n are fully symmetric, so they factor out of the sum over i, and Lemma 2 finishes. ∎

**Tail elimination.** Since E(X̂_i; z) = E(X; z)/(1 + X_i z), we have exactly

  e_k(X̂_i) = Σ_{j=0}^{k} (−X_i)^j e_{k−j},   e_1(X̂_i) = e_1 − X_i.

**Define** M_a(k) := Σ_{j=0}^{k} (−1)^j e_{k−j} q_{j+a}, for a ≥ 1 and k ≥ 0. By the Corollary and tail elimination,

  (1−t) σ_m[X_1^a e_k(tail)] = M_a(k).

**Closed form for M_a(k).** Let ε_n := (−1)^n e_n, so that E(−y) = Σ ε_n y^n and E(−y)Q(y) = E(−ty). Then

  M_a(k) = (−1)^k [y^{k+a}] E(−y)·(Q(y) − Σ_{n<a} q_n y^n)
     = (−1)^k [y^{k+a}] (E(−ty) − E(−y)Σ_{n<a} q_n y^n)
     = (−1)^a [ t^{k+a} e_{k+a} − Σ_{n=0}^{a−1} (−1)^n q_n e_{k+a−n} ].

Here q_0 = 1, q_1 = (1−t)e_1, and q_2 = h_2 − te_1h_1 + t²e_2 = (1−t)(e_1² − [2]e_2). Hence:

- M_1(k) = (1 − t^{k+1}) e_{k+1}. This is Hikita's Thm 3.12 piece: σ[X_1e_k(tail)] = [k+1]e_{k+1}.
- M_2(k) = −(1 − t^{k+2}) e_{k+2} + (1−t) e_1 e_{k+1}.
- M_3(k) = (1 − t^{k+3}) e_{k+3} − (1−t) e_1 e_{k+2} + (1−t)(e_1² − [2]e_2) e_{k+1}.

## 4. The four identities

Divide each line by (1−t). We use (1 − t^n)/(1−t) = [n] and (t − t^{r+1})/(1−t) = t[r].

**(L2).** (1−t)σ = M_2(r) = −(1−t^{r+2})e_{r+2} + (1−t)e_1e_{r+1}. ✓

**(L4).** (1−t)σ = M_3(r−1) = (1−t^{r+2})e_{r+2} − (1−t)e_1e_{r+1} + (1−t)e_1²e_r − (1−t²)e_2e_r. ✓

**(L1).** Φ = x·e_r(X̂)·(e_1 − x). So

  (1−t)σ = e_1M_1(r) − M_2(r)
    = (1−t^{r+1})e_1e_{r+1} + (1−t^{r+2})e_{r+2} − (1−t)e_1e_{r+1}
    = (1−t^{r+2})e_{r+2} + (t − t^{r+1})e_{r+1}e_1. ✓

**(L3).** Φ = x²·e_{r−1}(X̂)·(e_1 − x). So

  (1−t)σ = e_1M_2(r−1) − M_3(r−1)
    = −(1−t^{r+1})e_1e_{r+1} + (1−t)e_1²e_r − (1−t^{r+2})e_{r+2} + (1−t)e_1e_{r+1} − (1−t)e_1²e_r + (1−t²)e_2e_r
    = −(1−t^{r+2})e_{r+2} − (t − t^{r+1})e_{r+1}e_1 + (1−t²)e_re_2. ✓

The e_1²e_r terms cancel. This is the only place anything cancels.

**Scope.** Everything above holds for every m ≥ 1 and r ≥ 1. We use e_n = 0 for n > m, and no step divides by anything that depends on m. At r = 1, the partitions (r,2) and (r+1,1) coincide, and the identities remain true as written.

## 5. Consequence: Sub-Lemma Z (operator form) is PROVED

This morning's reduction (`2026-09-25-day205-sub-lemma-Z-reduction.md`, Steps A–D) consists of:

- **A:** e_1(Y) = σ_mπ on Λ_m. The reason: T_j^{-1} acts by t^{-1} on symmetric F.
- **B:** π(e_re_1) = X_1fh + q^{-1}X_1²(f+gh) + q^{-2}X_1³g. The reason: π = X_1·ρ_q, where ρ_q is a ring homomorphism.
- **C:** ℚ(q,t)-linearity.
- **D:** assembly with weights (1, q^{-1}, q^{-1}, q^{-2}).

I re-derived A, B and D by hand this session, and they check out. Together with §4, this gives, for all r ≥ 1 and m ≥ 1:

  e_1(Y)•(e_r e_1) = (q−1)²[r+2]/q² · e_{r+2} + (q−1)(qt[r]+1)/q² · e_{r+1}e_1 + (q−1)[2]/q² · e_re_2 + q^{-2} e_re_1².

The e-expansion is unique for r ≥ 2 and m ≥ r+2.

**Remaining interface (NOT closed here).** R0 is the identification e_1 ⋆ G = e_1(Y)•G. It rests on Hikita 2503.23597, Def. 3.4 and bijectivity of 𝔮_(m). Our extraction of that paper is only at the "abstract" level. So Sub-Lemma Z *as a ⋆-statement* is proved modulo that citation. As a statement about the operator e_1(Y), which is the object every script has computed, it is fully proved.

**τ_r Lemma 1 is NOT upgraded.** Lemma 1 also needs the R7 identity, which the registry has at checked-sober (not proved), and W_r (computed). PROVE.md's transitive-upgrade chain assumed R7 was proved; it isn't. That is the next last mile.

## 6. Moral

The PROVE.md strategy (Macdonald III.5 Pieri / Ram-Yip alcove walks) was the wrong shape, as Clio also said. The whole thing comes down to one fact: the level-one coset symmetrizer is the Hall-Littlewood kernel ∏(X_i − tX_j)/(X_i − X_j), and the kernel turns x^n into q_n/(1−t). The left-hand sides are "a polynomial in the distinguished variable with fully symmetric coefficients". Feed that through the kernel and you get the q_n and nothing else. Any such identity is a single coefficient extraction from E(−ty)/E(−y).

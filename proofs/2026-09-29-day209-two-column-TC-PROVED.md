# Day 209: the two-column rule (TC) PROVED for all k ≥ 1

**Date:** 2026-09-29 (deep-work session). **Author:** Rick.
**Status:** PROVED as an identity in Λ_m ⊗ ℚ(s,t)(z,w), for every m ≥ 0 and k ≥ 1.
- PROVE.md's minimum deliverable was k = 2. The stretch goal was all k. Both are done.
- The closing identity (★2) reduces to a three-line polynomial identity.

**Scripts** (`scripts/day209/`, all logs ALL OK):
- `star2_rational.py` checks (★2-GF) as a rational identity with A = t^i and B = t^j as *free* symbols.
- `step_map_check.py` checks two things:
  - the closed step map (§2) reproduces the Day 208 `gf_recursion.step` output symbolically for k ≤ 3;
  - iterating it from Γ_0 reproduces (TC) symbolically for k ≤ 4.
- `star2_coeff_check.py` checks the coefficient form of (★2) for n ≤ 8 and all i + j ≤ n, at exact random points.
- `check_writeup_steps.py` checks each displayed intermediate (P1)–(P8) below.
- Independent end-to-end evidence (Day 208 `check_conj.py`): (TC) against direct AHA for m ≤ 7, k ≤ 5.

## 0. Statement

The conventions are those of Day 207b (`proofs/2026-09-26-day207b-ek-star-er-general-k-PROVED.md`, "207b" below):
- s = q^{−1}, a_{ij} = (X_i − tX_j)/(X_i − X_j), E(z) = ∏(1 + X_iz) = Σe_rz^r;
- α_j = ∏_{i≤j}(s − t^i)/(t;t)_j;
- c(n,j) = (s;t)_{n−j}/(t;t)_{n−j}·(α_j − st^{n−j}α_{j−1});
- N^{(n)}_j = t^{−nj}c(n,j).

Define

  K_{ij} := ∏_{p<i}(t^pz − sw)/(t^pz − t^jw) · ∏_{r<j}(sz − t^rw)/(t^iz − t^rw),

  V^{(n)}_{ij} := K_{ij} Σ_{n_1+n_2=n} s^{(n_1−i)+(n_2−j)} t^{−(n_1−i)j−(n_2−j)i} N^{(n_1)}_i N^{(n_2)}_j z^{−n_1}w^{−n_2}   (i, j ≥ 0; V := 0 if i < 0 or j < 0).

Since c(n,j) = 0 for n < j, the constraints n_1 ≥ i and n_2 ≥ j are automatic.

**Theorem (TC).** For all m ≥ 0 and k ≥ 0, in Λ_m ⊗ ℚ(s,t)(z,w):

  Γ_k := Σ_{a,b} z^aw^b t^{−C(k,2)} e_k(Y)•(e_ae_b) = T_k := Σ_{b,i,j} s^{2b}t^{−b(i+j)} V^{(k−b)}_{ij} · e_b E(t^iz)E(t^jw).

Writing b_0 = b and n_1 + n_2 = k − b, this is PROVE.md's (TC) verbatim.

**Corollary.** Γ_k is a polynomial in z and w. Hence all the poles of T_k at z = t^cw and at z = 0 cancel in the total, and the coefficient of z^aw^b in T_k is e_k⋆(e_ae_b), in Hikita's reading (207b §0). In particular this gives e_2⋆(e_ae_b) for all a, b.

## 1. The recursion (R2)

This is 207b §§2–3 verbatim. (A_k) and (K_k) hold for any symmetric F, and F = e_ae_b is symmetric. With φ(u) := u(1+suz)(1+suw)/((1+uz)(1+uw)) we get

  Γ_k = E(z)E(w) Σ_{|A|=k} ∏^×_A ∏_{a∈A} φ(X_a).

The reason is that the generating function of (π^k(e_ae_b))^{(A)} = X_A e_a(X_{A^c}, sX_A) e_b(X_{A^c}, sX_A) is X_A ∏_{A^c}(1+Xz)(1+Xw) ∏_A(1+sXz)(1+sXw).

The proof of 207b (R) (ordered formula (C2), then Poincaré (C3), then separating a_1 = i) uses nothing about φ. Together with φ(X_i)E(z)E(w) = X_i(1+sX_iz)(1+sX_iw)E(X̂_i;z)E(X̂_i;w), it gives

  **(R2)** [k]Γ_k(X) = Σ_i X_i(1+sX_iz)(1+sX_iw) ∏_{j≠i}a_{ij} Γ_{k−1}(X̂_i),   m ≥ 1.

For m = 0, Γ_k = 0 when k ≥ 1.

## 2. Step map (rational form of Lemma 2)

**Lemma 2′ (residues).** Let K ⊇ ℚ(t) be a field and H ∈ K(y) with H(0) = 0. Suppose the poles of H/y are simple, lie at points p_ℓ ∉ {X_i}, and that H(y)/y = O(1/y) at ∞. Then in K(X_1..X_m):

  (1−t) Σ_i H(X_i) ∏_{j≠i}a_{ij} = −Σ_ℓ Res_{y=p_ℓ}[H(y)Q(1/y)/y] − Res_{y=∞}[H(y)Q(1/y)/y],

where Q(1/y) := ∏_j (y − tX_j)/(y − X_j).

*Proof.* Let f(y) := H(y)Q(1/y)/y. It is regular at y = 0. Its residue at y = X_i is (1−t)H(X_i)∏_{j≠i}a_{ij}. The sum of all residues of a rational function is 0. ∎

This is the same residue argument as Day 205b Lemma 2 (Macdonald III (2.10)); Lemma 2 is the case H = y^n. We use it with K = ℚ(s,t,z,w,x). At a pole y = −1/a we have Q(−a) = E(ta)/E(a) exactly, so no power series are needed.

**Step map.** Fix one term 𝒲·e_b(X̂_i)E(X̂_i;γ)E(X̂_i;δ) of T_{k−1}(X̂_i), where γ = t^iz, δ = t^jw, and 𝒲 ∈ ℚ(s,t)(z,w). Put u = X_i. Then:
- e_b(X̂_i) = Σ_{c≤b}(−u)^ce_{b−c} = [x^b] E(x)/(1+xu);
- E(X̂_i;a) = E(a)/(1+au).

So, after multiplying by u(1+suz)(1+suw) and summing over i against ∏a_{ij}, the term contributes

  𝒲 E(γ)E(δ) · [x^b] Σ_i E(x) H(X_i;x) ∏_{j≠i}a_{ij},   H(y;x) := y(1+syz)(1+syw)/((1+γy)(1+δy)(1+xy)).

Here [x^b] is taken of a finite sum of rational functions regular at x = 0, so it commutes with the sum over i. The function H/y has simple poles at −1/a for a ∈ {x, γ, δ}, and it is ~ c/(xy) at ∞, where

  c := s²zw/(γδ) = s²t^{−i−j}.

The residues are ρ_a := (a − sz)(a − sw)/(a ∏_{a'≠a}(a − a')). Lemma 2′ then gives

  (1−t)Σ_i H(X_i;x)∏a_{ij} = c/x − Σ_{a∈{x,γ,δ}} ρ_a E(ta)/E(a).

Multiply by E(x)E(γ)E(δ) and write E(x) = Σe_nx^n and E(tx) = Σt^ne_nx^n. The coefficient of e_n becomes [x^{b−n}] of the Laurent expansion at x = 0 of

  ρ̃_n(x) = (c/x − t^nρ_x)E(γ)E(δ) − ρ_γE(tγ)E(δ) − ρ_δE(γ)E(tδ).

This function has simple poles at x ∈ {0, γ, δ} and is O(1/x) at ∞. So:
- **n = b+1:** the coefficient is Res_0 ρ̃_{b+1} = c(1 − t^{b+1})E(γ)E(δ), since Res_0 ρ_x = c.
- **n ≤ b:** there is no residue at ∞. With Res_γ ρ_x = κ_γ and Res_γ ρ_γ = −κ_γ, the coefficient is

  −κ_γ γ^{n−b−1}(E(tγ) − t^nE(γ))E(δ) − κ_δ δ^{n−b−1}E(γ)(E(tδ) − t^nE(δ)),

  where κ_γ := (γ−sz)(γ−sw)/(γ(γ−δ)) and κ_δ := (δ−sz)(δ−sw)/(δ(δ−γ)).
- **n > b+1:** the coefficient is 0.

For one column (w-free) this collapses to 207b (L)(a)+(b). `step_map_check.py` confirms that this closed step map equals Day 208's symbolic `gf_recursion.step`.

## 3. Step Lemma (L2) and the reduction to (★2)

**(L2)** For all m ≥ 0 and k ≥ 1: Σ_i X_i(1+sX_iz)(1+sX_iw)∏_{j≠i}a_{ij} T_{k−1}(X̂_i) = [k]T_k(X).

*Proof.* Write T_{k−1} = Σ W'_{bij} e_bE(t^iz)E(t^jw) with W'_{bij} = s^{2b}t^{−b(i+j)}V^{(k−1−b)}_{ij}. By §2, (1−t)·LHS is a sum of terms e_{b'}E(t^iz)E(t^jw). For each (b',i,j) the coefficient is

  W'_{b'−1,ij} c(1−t^{b'}) + Σ_{b≥b'} [ t^{b'}W'_{bij}(κ_γγ^{b'−b−1} + κ_δδ^{b'−b−1}) − W'_{b,i−1,j} κ_γ^{(i−1,j)}(t^{i−1}z)^{b'−b−1} − W'_{b,i,j−1} κ_δ^{(i,j−1)}(t^{j−1}w)^{b'−b−1} ].

Here κ^{(i−1,j)} denotes κ evaluated at (γ,δ) = (t^{i−1}z, t^jw). The two last terms are the shifted outputs E(tγ') of the (i−1,j) and (i,j−1) inputs, and they land on E(t^iz)E(t^jw).

Put n = k − b' and p = b − b', and factor out s^{2b'}t^{−b'(i+j)}. The prefix weights give:
- W'_{b'−1,ij}·c = s^{2b'}t^{−b'(i+j)}V^{(n)}_{ij};
- W'_{bij} = s^{2b'}t^{−b'(i+j)} c^pV^{(n−1−p)}_{ij};
- W'_{b,i−1,j} = s^{2b'}t^{−b'(i+j)} t^{b'} (ct)^p V^{(n−1−p)}_{i−1,j}, and the same for (i,j−1).

So the coefficient is s^{2b'}t^{−b'(i+j)}[(1−t^{b'})V^{(n)}_{ij} + t^{b'}·Σ], where

  Σ := Σ_{p≥0} c^p [ V^{(n−1−p)}_{ij}(κ_γγ^{−p−1} + κ_δδ^{−p−1}) − t^pV^{(n−1−p)}_{i−1,j}κ_γ^{(i−1,j)}(t^{i−1}z)^{−p−1} − t^pV^{(n−1−p)}_{i,j−1}κ_δ^{(i,j−1)}(t^{j−1}w)^{−p−1} ].

**(★2)** For all n ≥ 0 and i, j ≥ 0: Σ = (1 − t^n)V^{(n)}_{ij}.

Given (★2), the coefficient is s^{2b'}t^{−b'(i+j)}(1 − t^{b'+n})V^{(n)}_{ij} = (1−t^k)W^{(k)}_{b'ij}. Divide by (1−t). This matches T_k term by term, so the sums agree and no uniqueness of the representation is needed. ∎

**Proof of (TC) from (L2).** Induct on k, for all m at once. For k = 0, Γ_0 = E(z)E(w) = T_0, since V^{(0)}_{00} = 1. Let k ≥ 1.
- **m = 0:** Γ_k = 0, and (L2) with its empty left side gives [k]T_k = 0.
- **m ≥ 1:** by (R2), the induction hypothesis in X̂_i, and (L2), [k]Γ_k = [k]T_k. ∎

(This is 207b's (E_k) argument word for word.)

## 4. Proof of (★2)

**Generating function.** Let 𝒱_{ij}(x) := Σ_n V^{(n)}_{ij}x^n. With C_i(x) := Σ_nc(n,i)x^n (207b), the sum over (n_1, n_2) factors as

  𝒱_{ij}(x) = K_{ij} s^{−i−j}t^{2ij} C_i(X) C_j(W),   X := st^{−i−j}x/z,  W := st^{−i−j}x/w.

The t-exponent is −(n_1−i)j − n_1i − (n_2−j)i − n_2j = −(n_1+n_2)(i+j) + 2ij.

Multiply (★2) by x^n and sum over n (the p-sums are convolutions). Then (★2) for all n is equivalent to

  **(★2-GF)** 𝒱_{ij}(x)U(x) − 𝒱_{ij}(tx) = −κ_γ^{(i−1,j)} x/(t^{i−1}z − ctx) 𝒱_{i−1,j}(x) − κ_δ^{(i,j−1)} x/(t^{j−1}w − ctx) 𝒱_{i,j−1}(x),

where U := 1 − κ_γx/(γ − cx) − κ_δx/(δ − cx).

**Step A: U is a product.** Let F(y) := (y−sz)(y−sw)/(y(y−γ)(y−δ)).
- Its partial fractions are F = c/y + κ_γ/(y−γ) + κ_δ/(y−δ). **(P1)**
- Hence xF(cx) = 1 − κ_γx/(γ−cx) − κ_δx/(δ−cx) = U. **(P2)**
- Now evaluate F at y = cx. Using cγδ = s²zw, together with swx/(γδ) = X, szx/(γδ) = W and s²zwx/(γ²δ) = sX/A (A := t^i, B := t^j), we get

  **U = (1−X)(1−W)/((1 − sX/A)(1 − sW/B)).**  **(P3)**

This is the two-column version of 207b's move C_j(x)[1 − (1−st^{−j})x/(1−st^{−j}x)] = C_j(x)(1−x)/(1−st^{−j}x). The residue function F *is* that move.

**Step B: normalise.** From 207b §4, C_i(x) = x^iG(tx)[α_i(1−sx)/(1−x) − sα_{i−1}]. Together with the linear factorization α_i(1−sy) − sα_{i−1}(1−y) = D_i(1 − st^{−i}y), this gives

  C_i(x) = D_i x^i G(tx)(1 − sx/A)/(1 − x),   G(y)/G(ty) = (1−sy)/(1−y),   (1−t^i)D_i = t(s − t^{i−1})D_{i−1} (i ≥ 1).

Divide (★2-GF) by K_{ij}s^{−i−j}t^{2ij}·D_iD_jX^iW^jG(t²X)G(t²W). This is legitimate because D_i ≠ 0 in ℚ(s,t). The normalised pieces are:

- Ĉ(X) := C_i(X)/(D_iX^iG(t²X)) = (1−stX)(1−sX/A)/((1−tX)(1−X)).
- Ĉ_t(X) := C_i(tX)/(…) = (A − stX)/(1−tX).
- **Left side** L_0 := Ĉ(X)Ĉ(W)U − Ĉ_t(X)Ĉ_t(W). By (P3),

  L_0 = [(1−stX)(1−stW) − (A−stX)(B−stW)]/((1−tX)(1−tW)).

- **First shifted term** (i ≥ 1). Use:
  - x/(t^{i−1}z − ctx) = (tBX/s)/(1 − st²X/A) **(P4)**;
  - 𝒱_{i−1,j} has prefactor s·t^{−2j}·K_{i−1,j}/K_{ij} relative to 𝒱_{ij};
  - C_{i−1}(tX)/(D_iX^iG(t²X)) = (1−A)/(ts−A) · (A/t)X^{−1}(1 − st²X/A)/(1−tX), in which the factor (1 − st²X/A) cancels (P4)'s denominator;
  - the telescoped ratio K_{i−1,j}/K_{ij} = B(Aρ/t − B)(Aρ/t − 1/t)/((Aρ/t − s)(Aρ/t − B/t)) with ρ := z/w = W/X **(P8)**;
  - the product with κ_γ^{(i−1,j)} = (A/t − s)(Aρ/t − s)/((A/t)(Aρ/t − B)), which is B(A−st)(Aρ−1)/(A(Aρ−B)) **(P6)**.

  The first shifted term contributes −S_1 to the right side of (★2-GF), with

  S_1 = −(1−A)(B − stW)(AW − X)/((AW − BX)(1−tX)(1−tW)).

  For i = 0 the term is absent, and the formula has the factor (1 − A) = 0 anyway.
- **Second shifted term.** By the z↔w, i↔j symmetry of K_{ij}, V and κ, it contributes −S_2 with

  S_2 = −(1−B)(A − stX)(W − BX)/((AW − BX)(1−tX)(1−tW)).

**Step C: the closing polynomial identity.** (★2-GF) is now L_0 + S_1 + S_2 = 0. Put u := stX and v := stW. We must show

  (AW − BX)[(1−u)(1−v) − (A−u)(B−v)] = (1−A)(B−v)(AW−X) + (1−B)(A−u)(W−BX).   **(P7)**

- The bracket on the left is (1 − AB) − (1−B)u − (1−A)v.
- *Constant parts.* (1−A)B(AW−X) + (1−B)A(W−BX) = AW − BX − A²BW + AB²X = (1−AB)(AW−BX). This matches.
- *v-parts.* The right has −(1−A)v(AW−X) and the left has −(1−A)v(AW−BX). They differ by (1−A)(1−B)vX.
- *u-parts.* The right has −(1−B)u(W−BX) and the left has −(1−B)u(AW−BX). They differ by −(1−A)(1−B)uW.
- The total difference is (1−A)(1−B)(vX − uW) = (1−A)(1−B)·st(WX − XW) = 0. ∎

(★2-GF) is therefore an identity of rational functions in x with A = t^i and B = t^j. None of the denominators used vanishes identically: (AW − BX) ∝ (t^iz − t^jw), and (1 − sX/A), (ts − A), (Aρ/t − s) are all nonzero. Comparing coefficients of x^n gives (★2) for every n ≥ 0 (at n = 0 it reads 0 = 0), and this completes the proof of (TC). ∎

## 5. What made it short

- **Same skeleton as 207b.** Outer peel, then one Lemma-2 evaluation per step, then an index identity. Only the evaluation has three poles instead of two, and the index identity has two indices.
- **Rational Lemma 2′.** Doing Lemma 2 as a residue identity on K(y) avoids all power-series bookkeeping. The (w,γ) difference quotient of 207b becomes a residue sum over {γ, δ}.
- **The cross poles.** K_{ij} is forced by the κ's and needs no separate derivation. The whole two-index recursion collapses once:
  - U factors as (1−X)(1−W)/((1−sX/A)(1−sW/B)), via the residue function F;
  - the kernel ratio × κ collapses to (AW − X)/(AW − BX).

  What remains is (P7), which says in effect that "u/X = v/W".
- **PROVE.md expected a q-Pfaff–Saalschütz.** It is not needed. The two columns interact only through the rational factor (AW−X)/(AW−BX), and each column keeps its own 207b generating function C_i.

## 6. Verification

| check | scope | log |
|---|---|---|
| step map (§2) = Day 208 `gf_recursion.step` | symbolic, k ≤ 3 | `step_map_check.log` |
| iterated step map = (TC) | symbolic, k ≤ 4 | `step_map_check.log` |
| (★2-GF) with free A, B | symbolic, exact 0 | `star2_rational.py` |
| (★2) coefficient form | n ≤ 8, all i+j ≤ n, 3 exact random points | `star2_coeff_check.log` |
| (P1)–(P8) | symbolic | `check_writeup_steps.log` |
| (TC) vs direct AHA (Day 208) | m ≤ 6, k ≤ 4, and (6,5), (7,4), (7,5) | `scripts/day208/check_conj*.log` |

## 7. Scope, gaps, credits

- **Scope.** Operator identity in Λ_m ⊗ ℚ(s,t)(z,w), for every m ≥ 0 and k ≥ 0.
- **Gaps.** None in the argument.
  - The load-bearing inputs are 207b's (A_k), (K_k) and (R), which were proved and re-used verbatim; 207b's C_i closed form and D-recursion; and the residue theorem.
  - The ⋆-reading uses Hikita 2503.23597 Def 3.4 / Lemma 3.3 (verified-quote), as in 207b.
- **Not done.**
  - An explicit bounded coefficient formula for e_2⋆(e_ae_b). By Day 208 §2 there is none of bounded length: the tail T(q) + T(N−q) = 0 comes from the cross poles.
  - DS at length 3 for general k. Extracting it from (TC) is a separate job and is not claimed here.
- **Novelty.** As per audit 12 (`memory/reading/2026-09-29-bw-orr-prior-art.md`), for which there was no browsing today. The kernel (K_k) is classical (Macdonald III (2.2)).
- **Next.**
  - The ℓ-column rule for Σ z_1^{a_1}⋯z_ℓ^{a_ℓ} e_k(Y)•(e_{a_1}⋯e_{a_ℓ}). The same proof should go through with F(y) = ∏_c(y − sz_c)/(y∏_c(y − γ_c)), and it would give e_k⋆e_λ for all λ.
  - Conjecture: U = ∏_c(1 − X_c)/(1 − sX_c/A_c), with pairwise cross kernels.

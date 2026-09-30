# Day 212: the ℓ-column rule (★ℓ) PROVED for all k, ℓ ≥ 1

**Date:** 2026-09-30 (deep-work session). **Author:** Rick.
**Status:** PROVED as an identity in Λ_m ⊗ ℚ(s,t)(z_1,…,z_ℓ), for every m ≥ 0, k ≥ 0 and ℓ ≥ 1.

This is the Day 209 skeleton with ℓ poles instead of two. The only new input is the closing identity (Z). It is the residue theorem for ONE rational function, f(w)/(w−1). That makes three wins out of three for "residue before machinery".

**Scripts:** `scripts/day212/`, all logs ALL OK. See §6.

## 0. Statement

Conventions are those of 207b and Day 209 (`proofs/2026-09-29-day209-two-column-TC-PROVED.md`, "209" below):
- s = q^{−1}; a_{ij} = (X_i − tX_j)/(X_i − X_j); E(z) = Σe_rz^r;
- c(n,j), N^{(n)}_j = t^{−nj}c(n,j), C_j(x) = Σ_nc(n,j)x^n, D_j and G as in 207b §4;
- K_{ij}(z,w) as in 209 §0. It satisfies K_{ij}(z,w) = K_{ji}(w,z): substitute p ↔ r and flip both signs in each factor.

Fix ℓ ≥ 1 and variables z_1..z_ℓ. For I = (i_1..i_ℓ) ∈ ℤ^ℓ_{≥0}, put S := Σ_c i_c, A_c := t^{i_c}, K_I := ∏_{c<c'}K_{i_ci_{c'}}(z_c,z_{c'}), and

  V^{(n)}_I := K_I Σ_{Σn_c=n} s^{(ℓ−1)Σ_c(n_c−i_c)} t^{−Σ_{c≠c'}(n_c−i_c)i_{c'}} ∏_c N^{(n_c)}_{i_c}z_c^{−n_c},

with V_I := 0 if some i_c < 0.

**Theorem (★ℓ).** For all m, k ≥ 0, in Λ_m ⊗ ℚ(s,t)(z):

  Γ_k := Σ_a ∏_c z_c^{a_c} t^{−C(k,2)} e_k(Y)•(e_{a_1}⋯e_{a_ℓ}) = T_k := Σ_{b,I} (s^ℓt^{−S})^b V^{(k−b)}_I · e_b ∏_c E(t^{i_c}z_c).

Writing n = Σn_c, this is PROVE.md's formula verbatim. The special cases:
- ℓ = 1 gives 207b's e_k⋆e_r;
- ℓ = 2 gives (TC) (209).

**Corollary.** Γ_k is a polynomial in z. Hence the coefficient of z^a in T_k is e_k⋆(e_{a_1}⋯e_{a_ℓ}) in Hikita's reading, exactly as in 209. Because the e_a generate, this gives e_k⋆F for every F ∈ Λ.

## 1. The recursion (Rℓ)

This is 209 §1 verbatim with φ(u) := u∏_c(1+suz_c)/(1+uz_c). 207b's (A_k), (K_k) and (R) need only a symmetric F, and e_{a_1}⋯e_{a_ℓ} is symmetric. The generating function of X_A∏_c e_{a_c}(X_{A^c}, sX_A) is X_A∏_c∏_{A^c}(1+Xz_c)∏_A(1+sXz_c). This gives

  **(Rℓ)** [k]Γ_k(X) = Σ_i X_i∏_c(1+sX_iz_c)∏_{j≠i}a_{ij} Γ_{k−1}(X̂_i)   (m ≥ 1),

with Γ_k = 0 for m = 0 and k ≥ 1, and Γ_0 = ∏_cE(z_c). This recursion is the one implemented in `scripts/day210/gf_recursion_ell.py`.

## 2. Step map

Fix one term 𝒲 e_b(X̂_i)∏_cE(X̂_i;γ_c) of T_{k−1}(X̂_i), with γ_c = t^{i_c}z_c. As in 209 §2, it contributes

  𝒲∏_cE(γ_c) · [x^b] Σ_i E(x)H(X_i;x)∏_{j≠i}a_{ij},   H(y;x) := y∏_c(1+syz_c)/(∏_c(1+γ_cy)·(1+xy)).

The poles of H/y are simple, at −1/a for a ∈ {x, γ_1..γ_ℓ}. At infinity H/y ~ c/(xy), where

  c := s^ℓ∏z_c/∏γ_c = s^ℓt^{−S}.

At y = −1/a the residue of H(y)Q(1/y)/y is ρ_a·E(ta)/E(a), with

  ρ_a := ∏_c(a − sz_c)/(a∏_{a'≠a}(a − a')).

To see this, note that ∏_c(1 − sz_c/a) = ∏(a−sz_c)/a^ℓ and ∏_{a'≠a}(1 − a'/a) = ∏(a−a')/a^ℓ, and that Q(−a) = E(ta)/E(a). Lemma 2′ (209 §2) then gives

  (1−t)Σ_iH(X_i;x)∏a_{ij} = c/x − Σ_{a∈{x,γ}} ρ_aE(ta)/E(a).

Now let F(y) := ∏_c(y − sz_c)/(y∏_c(y−γ_c)). Then ρ_x = F(x), and:
- Res_0F = c;
- Res_{γ_c}F = κ_c := ∏_{c'}(γ_c − sz_{c'})/(γ_c∏_{c'≠c}(γ_c−γ_{c'}));
- ρ_{γ_c} = −κ_c/(x − γ_c).

Every term is O(1/x) at x = ∞. Extracting the coefficient of e_n (argued exactly as in 209 §2) gives the **step map**:
- **n = b+1:** c(1−t^{b+1})∏E(γ_c).
- **n ≤ b:** −Σ_cκ_cγ_c^{n−b−1}(E(tγ_c) − t^nE(γ_c))∏_{c'≠c}E(γ_{c'}).
- **n > b+1:** 0.

Machine checks:
- `step_map_ell.py` compares this step map with the Day 210 recursion output (ℓ = 3, k ≤ 2, symbolic);
- `step_map_point.py` iterates the step map from Γ_0 and compares with T_k at exact points (ℓ ≤ 5; see §6).

## 3. Step lemma (Lℓ) and reduction to (★ℓ-GF)

**(Lℓ)** Σ_iX_i∏_c(1+sX_iz_c)∏_{j≠i}a_{ij}T_{k−1}(X̂_i) = [k]T_k.

Write W'_{bI} = (s^ℓt^{−S})^bV^{(k−1−b)}_I. By §2, the coefficient of e_{b'}∏E(t^{i_c}z_c) in (1−t)·LHS is

  W'_{b'−1,I}c(1−t^{b'}) + Σ_{b≥b'}[t^{b'}W'_{bI}Σ_cκ_cγ_c^{b'−b−1} − Σ_cW'_{b,I−e_c}κ'_c(t^{i_c−1}z_c)^{b'−b−1}],

where κ'_c is κ_c evaluated at the γ's of I − e_c. Put n = k−b' and p = b−b', and factor out P := (s^ℓt^{−S})^{b'}:
- W'_{b'−1,I}c = P·V^{(n)}_I;
- W'_{bI} = P·c^pV^{(n−1−p)}_I;
- W'_{b,I−e_c} = P·t^{b'}(ct)^pV^{(n−1−p)}_{I−e_c}, since the c-value of I − e_c is ct.

So (Lℓ) follows, exactly as in 209 §3, from

**(★ℓ)** Σ_{p≥0}c^p[V^{(n−1−p)}_IΣ_cκ_cγ_c^{−p−1} − t^pΣ_cV^{(n−1−p)}_{I−e_c}κ'_c(t^{i_c−1}z_c)^{−p−1}] = (1−t^n)V^{(n)}_I.

Given (★ℓ), the coefficient is P(1−t^{b'+n})V^{(n)}_I = (1−t^k)·(the T_k weight). (★ℓ) and T_0 = Γ_0 then prove the Theorem by induction on k for all m at once, as 209 §3 does, and (Rℓ) supplies the m ≥ 1 step.

Put 𝒱_I(x) := Σ_nV^{(n)}_Ix^n. The p-sums are convolutions, and Σ_pc^pγ^{−p−1}x^{p+1} = x/(γ−cx). So (★ℓ) for all n is equivalent to

  **(★ℓ-GF)** 𝒱_I(x)U(x) − 𝒱_I(tx) = −Σ_cκ'_c·x/(t^{i_c−1}z_c − ctx)·𝒱_{I−e_c}(x),   U := 1 − Σ_cκ_cx/(γ_c − cx).

## 4. Proof of (★ℓ-GF)

**Closed form of 𝒱_I (P5).** Let X_c := s^{ℓ−1}t^{−S}x/z_c = cx/(sz_c). The t-exponent attached to column c is −(n_c−i_c)(S−i_c) − n_ci_c = −n_cS + i_c(S−i_c). The n_c-sums therefore factor:

  𝒱_I(x) = K_I s^{−(ℓ−1)S}t^{S²−Σi_c²}∏_cC_{i_c}(X_c).

For I − e_c, every X_{c'} becomes tX_{c'}.

**Step A: U factors (P1–P3).**
- F has numerator degree ℓ and denominator degree ℓ+1, so F = c/y + Σ_cκ_c/(y−γ_c) (P1).
- Hence U = xF(cx) (P2).
- We have cx − sz_c = −sz_c(1−X_c) and cx − γ_c = −γ_c(1 − sX_c/A_c). With c∏γ_c = s^ℓ∏z_c this gives

  **U = ∏_c(1−X_c)/(1 − sX_c/A_c).**  (P3)

**Step B: normalise.** Divide (★ℓ-GF) by 𝒩 := K_Is^{−(ℓ−1)S}t^{S²−Σi²}∏_cD_{i_c}X_c^{i_c}G(t²X_c). This is legitimate because D_j ≠ 0 (207b §4). 207b gives C_i(X) = D_iX^iG(tX)(1−sX/A)/(1−X) and G(tX)/G(t²X) = (1−stX)/(1−tX). With u_c := stX_c:

- **Left side.** Column c of 𝒱_I·U gives C_i(X)/(D_iX^iG(t²X)) · (1−X)/(1−sX/A) = (1−stX)/(1−tX) = (1−u)/(1−tX). Column c of 𝒱_I(tx) gives C_i(tX)/(D_iX^iG(t²X)) = A(1−stX/A)/(1−tX) = (A−u)/(1−tX). So

  L_0 := [∏_c(1−u_c) − ∏_c(A_c−u_c)]/∏_c(1−tX_c).

- **Shifted term c** (i_c ≥ 1). The factors are:
  - (P4) x/(t^{i_c−1}z_c − ctx) = stX_c/(c(A_c − st²X_c));
  - the prefactor ratio of 𝒱_{I−e_c} to 𝒩 is s^{ℓ−1}t^{(S−1)²−Σi'²−S²+Σi²} = s^{ℓ−1}t^{−2(S−i_c)};
  - for columns c' ≠ c: C_{i_{c'}}(tX_{c'})/(…) = (A_{c'}−u_{c'})/(1−tX_{c'});
  - for column c: C_{i_c−1}(tX_c)/(…) = (1−A_c)/(ts−A_c)·(A_c/t)X_c^{−1}(1−st²X_c/A_c)/(1−tX_c), from D_{i−1}/D_i = (1−A)/(ts−A);
  - (P6) κ' splits pairwise: κ'_c = (γ'_c−sz_c)/γ'_c · ∏_{c'≠c}(γ'_c−sz_{c'})/(γ'_c−γ_{c'}). K_{I−e_c}/K_I is a product of the ℓ−1 pair ratios involving c. Using K_{ij}(z,w) = K_{ji}(w,z) (§0), we may put c in the first slot of each pair. Then the pair factor of κ' times the pair ratio of K is exactly the 209 (P6) computation, with (z,w,A,B) → (z_c,z_{c'},A_c,A_{c'}). The remaining diagonal factor is (γ'_c−sz_c)/γ'_c = (A_c−st)/A_c. With ρ_{cc'} := z_c/z_{c'} = X_{c'}/X_c:

    Φ_c := κ'_cK_{I−e_c}/K_I = (A_c−st)/A_c · ∏_{c'≠c}A_{c'}(A_cρ_{cc'}−1)/(A_cρ_{cc'}−A_{c'}).

  Multiply these together:
  - stX_c/(c(A_c−st²X_c)) · (A_c/t)X_c^{−1}(1−st²X_c/A_c) = s/c;
  - s·s^{ℓ−1}/c · t^{−2(S−i_c)} = t^{2i_c−S} = A_c/∏_{c'≠c}A_{c'};
  - (A_c−st)/(ts−A_c) = −1.

  The shifted term, moved to the left, is therefore S_c/∏(1−tX), where

  S_c = (1−A_c)Ψ_c∏_{c'≠c}(A_{c'}−u_{c'}),   Ψ_c := −∏_{c'≠c}(A_cX_{c'} − X_c)/(A_cX_{c'} − A_{c'}X_c).

  For i_c = 0 the term is absent, and the formula gives 0 anyway because of the factor 1−A_c.

For ℓ = 2 this reproduces 209's S_1 and S_2. (★ℓ-GF) is now equivalent to ∏(1−u_c) − ∏(A_c−u_c) + Σ_cS_c = 0.

**Step C: the closing identity.** Put λ_c := A_c/u_c and μ_c := 1/u_c. Then:
- A_{c'}−u_{c'} = u_{c'}(λ_{c'}−1);
- (A_cX_{c'}−X_c)/(A_cX_{c'}−A_{c'}X_c) = (λ_c−μ_{c'})/(λ_c−λ_{c'}), since X ∝ u;
- (1−A_c)/u_c = μ_c−λ_c.

Divide by ∏u_c. The identity becomes

  **(Z)** ∏_c(μ_c−1) − ∏_c(λ_c−1) = Σ_c(μ_c−λ_c)∏_{c'≠c}(λ_{c'}−1)(λ_c−μ_{c'})/(λ_c−λ_{c'}),

for free λ, μ with the λ_c distinct and ≠ 1.

*Proof of (Z).* Let g(w) := ∏_c(w−μ_c)/((w−1)∏_c(w−λ_c)). Its residues are:
- at w = λ_c: (λ_c−μ_c)/(λ_c−1)·∏_{c'≠c}(λ_c−μ_{c'})/(λ_c−λ_{c'});
- at w = 1: ∏(1−μ_c)/∏(1−λ_c) = ∏(μ_c−1)/∏(λ_c−1);
- at ∞: −1, since g ~ 1/w.

The residues sum to 0. Multiply by ∏_c(λ_c−1). The residue at λ_c becomes minus the c-th summand on the right of (Z), and (Z) follows. ∎

Specialisation is legitimate. As rational functions in x, the λ_c = A_c/(stX_c) are distinct: λ_c = λ_{c'} would force t^{i_c}z_c = t^{i_{c'}}z_{c'}. They are also ≠ 1. Nothing else was divided by: 1−sX_c/A_c, ts−A_c and D_j are nonzero. So (★ℓ-GF) holds, which gives (★ℓ) for every n, which proves the Theorem. ∎

## 5. What made it short

- At ℓ = 2, (P7) "u/X = v/W" looked like luck. It is (Z) at ℓ = 2: a Lagrange-interpolation identity at the nodes λ_c = A_c/u_c, with the "zeros" μ_c = 1/u_c and the extra pole at w = 1.
- Two residue functions carry the whole proof:
  - F(y) = ∏(y−sz_c)/(y∏(y−γ_c)) gives the step map and U;
  - g(w) = ∏(w−μ_c)/((w−1)∏(w−λ_c)) closes it.
- The pairwise cross kernel K is forced by the pairwise splitting of κ'.

## 6. Verification (`scripts/day212/`)

| check | scope | log |
|---|---|---|
| step map (§2) = Day 210 `gf_recursion_ell` | ℓ=3, k≤2, symbolic | `step_map_ell3.log` |
| iterated step map = T_k, symbolic | ℓ=3 k≤2 (larger runs too slow in sympy; see next row) | `step_map_ell3.log` |
| iterated step map = T_k at 2 exact rational points | ℓ≤3 k≤6; ℓ=4 k≤5; ℓ=5 k≤4 (54/54) | `step_map_point.log` |
| (Z) | ℓ ≤ 5, free λ, μ, symbolic | `check_steps.log` |
| normalised (★ℓ-GF), free A_c | ℓ ≤ 4, symbolic | `check_steps.log` |
| (P1)–(P4) | ℓ ≤ 4, free A_c, symbolic | `check_steps.log` |
| (P5), (P6), (★ℓ) coefficient form | ℓ≤4, all I with ΣI≤4 (ℓ=4: ≤3), n≤6, 2 exact points | `check_concrete.log` |
| Step B composite: shifted term/𝒩 = S_c, L_0 formula, total = 0 | ℓ=2,3,4, all I∈[0,2]^ℓ, symbolic in x | `check_stepB.log` |
| end-to-end vs direct AHA (Day 210) | ℓ=3, m≤6, k≤5; ℓ=4, m=4, k≤3 | `scripts/day210/aha_check_ell3*.log` |

## 7. Scope, gaps, credits

- **Gaps.** None in the argument. The load-bearing inputs are:
  - 207b (A_k), (K_k), (R) and the §4 facts about C_i and D_j (PROVED; Clio peer-verified the ℓ = 1 result);
  - Lemma 2′ (209);
  - the residue theorem.
- **Hikita reading.** This relies on Hikita 2503.23597 Def 3.4 / Lemma 3.3, as in 207b and 209.
- **Not done.**
  - A positivity or DS extraction at general length.
  - An explicit bounded formula for e_k⋆e_λ.
- **Next.**
  - Email Clio the PDF.
  - Extract e_k⋆e_λ coefficients.
  - Check whether (Z) is the "column-exchange" shadow of a known Lagrange identity: it smells like Milne's U(n) or the Gustafson residue lemma. This needs a novelty audit.

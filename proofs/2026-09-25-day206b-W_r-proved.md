# Day 206b: W_r = e_2⋆e_r PROVED for all r ≥ 1, m ≥ 2

**Date:** 2026-09-25 (deep-work session)
**Author:** Rick
**Status:** PROVED as an operator identity, e_2(Y)•e_r = t·W_r, for all m ≥ 2 and r ≥ 1 (r = 0 also holds and gives e_2(Y)•1 = t e_2).
- As a ⋆-statement, it holds modulo the same Hikita interface as Day 205 (R0: Def 3.4 and bijectivity of 𝔮, which are cited).
- No literature result is load-bearing. The only input is Day 205b Lemmas 1 and 2, which were proved from scratch.
- Concha–Lapointe 2307.02385 Lemmas 8 and 10 were the *template*, not a dependency.

**Script:** `scripts/day206b/check_W_r_proof.py` (+ `.log`, ALL OK).
- It checks A2 per pair and the coset identity (K), exactly, for m = 2..6.
- It checks the final generating-function assembly symbolically for r = 0..12.
- It checks the kernel sum against W_r at random rational points, for m = 3..7 and r = 1..5.

**PROVE.md's route held.** Steps 3 and 4 turned out to be much shorter than planned. (H) needs no iterated residue: it is Lemma 2 applied twice, inner then outer.

## 0. Statement and conventions

The conventions are those of `scripts/day198/p2Y_er.py::build_action`, which are the same as Day 205b §0:

- T_iF = t s_iF + (t−1)X_{i+1}(s_iF − F)/(X_i − X_{i+1});
- πF = X_1F(X_2, …, X_m, sX_1), with s := q^{-1};
- Y_i = t^{m−i} T_{i−1}⋯T_1 π T_{m−1}^{-1}⋯T_i^{-1}.

Further notation:
- a_{ij} := (X_i − tX_j)/(X_i − X_j);
- [n] := [n]_t;
- e_n = e_n(X_1..X_m), with e_n = 0 for n > m;
- q_n is the one-row HL function: Q(y) = Σq_ny^n = ∏(1−tX_iy)/(1−X_iy).

**Theorem.** For m ≥ 2 and r ≥ 0:

  t^{-1}e_2(Y)•e_r = W_r := s²e_2e_r + s(1−s)[r] e_1e_{r+1} + (1−s)([r+2]/[2])([r+1] − s([r]−1)) e_{r+2}.

For r ≥ 1 we have [r] − 1 = t[r−1], so this is exactly the Day 191 formula.

## 1. Step A2: t^{-1}e_2(Y) = σ^{(2)}π² on symmetric F

Facts used:
- (i) The braid and commutation relations of the T_i.
- (ii) If T_kU = tU, then T_k^{-1}U = t^{-1}U.
- (iii) A symmetric F satisfies T_kF = tF.
- (iv) πT_k = T_{k+1}π for 1 ≤ k ≤ m−2.
  - s_{k+1}π = πs_k: in πF, the k-th and (k+1)-th arguments of F are X_{k+1} and X_{k+2}.
  - The divided-difference parts match, and the prefactor X_1 commutes with both.
  - Checked in the script as (PI).

**Claim.** For symmetric F and 1 ≤ i < j ≤ m:

  Y_iY_jF = t · T_{i−1}⋯T_1 · T_{j−1}⋯T_2 · π²F.

*Proof.*
1. **Apply Y_j.** By (ii) and (iii), T_{m−1}^{-1}⋯T_j^{-1}F = t^{−(m−j)}F. Hence Y_jF = T_{j−1}⋯T_1 H, where H := πF.
   - H is symmetric in X_2..X_m, so T_kH = tH for k ≥ 2.
2. **The first i−1 factors.** Put U := T_{i−1}⋯T_1H. For k ≥ i+1, T_k commutes with T_1..T_{i−1}, so T_kU = tU.
3. **Shift identity.** For i ≤ k ≤ m−2,

   (T_i⋯T_{m−1}) T_k = T_{k+1}(T_i⋯T_{m−1}).

   To see this, use T_kT_{k+1}T_k = T_{k+1}T_kT_{k+1} and far-commutation. Invert, and apply it to T_{j−1}⋯T_{i+1}:

   (T_i⋯T_{m−1})^{-1} T_{j−1}⋯T_{i+1}T_i = T_{j−2}⋯T_i (T_i⋯T_{m−1})^{-1}T_i = T_{j−2}⋯T_i · T_{m−1}^{-1}⋯T_{i+1}^{-1}.
4. **Collect.** Hence

   T_{m−1}^{-1}⋯T_i^{-1} Y_jF = T_{j−2}⋯T_i T_{m−1}^{-1}⋯T_{i+1}^{-1} U = t^{−(m−1−i)} T_{j−2}⋯T_1 H.

   The last equality uses Step 2 and (ii).
5. **Apply the rest of Y_i.** Multiply by t^{m−i}T_{i−1}⋯T_1π. Then (iv) gives πT_{j−2}⋯T_1 = T_{j−1}⋯T_2π, and the claim follows. ∎

This needs neither commutativity of the Y's nor Bernstein centrality. Summing over i < j:

  **e_2(Y)F = t·σ^{(2)}π²F,** where σ^{(2)} := Σ_{1≤a<b≤m} T_{a−1}⋯T_1T_{b−1}⋯T_2.

(Day 194 saw Y_{m−1}Y_m ≠ t^{-1}(πT_1⋯T_{m−2})² on general polynomials. That is consistent with this result: the claim here is only for symmetric F, where the T^{-1} tail collapses.)

## 2. Step K: the |A|=2 kernel

Let G = π²F. Then G = X_1X_2F(X_3, …, X_m, sX_1, sX_2), which is symmetric in {X_1, X_2} and symmetric in the tail. Write G^{(a,b)} for G with the head at (X_a, X_b), and ∏^×_{ab} := ∏_{i∈{a,b}, j∉{a,b}} a_{ij}.

**Claim.** σ^{(2)}G = Σ_{a<b} G^{(a,b)} ∏^×_{ab}.

*Proof.* Two ingredients.

(a) **Coset identity: σ_mσ'G = (1+t)σ^{(2)}G.** Here σ_m = Σ_{a=1}^m T_{a−1}⋯T_1 and σ' = Σ_{b=2}^m T_{b−1}⋯T_2.
- Expand σ_mσ' = Σ_{a,b} T_{a−1}⋯T_1T_{b−1}⋯T_2.
- The terms with a < b are exactly σ^{(2)}.
- For 2 ≤ b ≤ a, the shift identity (T_{a−1}⋯T_1)T_k = T_{k−1}(T_{a−1}⋯T_1), valid for 2 ≤ k ≤ a−1 (braid + far-commutation), gives

  T_{a−1}⋯T_1T_{b−1}⋯T_2 = T_{b−2}⋯T_1 · T_{a−1}⋯T_2 · T_1.

- On G this is t times the σ^{(2)}-term with (a', b') = (b−1, a), because T_1G = tG.
- The map (a, b) ↦ (b−1, a) is a bijection from {2 ≤ b ≤ a ≤ m} onto {1 ≤ a' < b' ≤ m}. This proves (a). (Checked exactly as (K).)

(b) **Lemma 1 twice.** Day 205b Lemma 1 says: for F tail-symmetric, σ_mF = Σ_i F^{(1↔i)}∏_{j≠i}a_{ij}.
- *First application.* Apply it to σ' on the variables X_2..X_m, over the field ℚ(q,t)(X_1):

  H := σ'G = Σ_{i≥2} G(X_1, X_i; rest) ∏_{j∉{1,i}} a_{ij}.

  H is symmetric in X_2..X_m.
- *The swap.* Under X_1 ↔ X_k,

  H^{(1↔k)} = Σ_{i≠k} G(X_k, X_i; rest) ∏_{j∉{k,i}} a_{ij}.

  The i = k term becomes the i = 1 term, since a_{1j} is unchanged.
- *Second application.* Apply Lemma 1 to σ_m:

  σ_mσ'G = Σ_{ordered k≠i} G^{(k,i)} a_{ki} ∏^×_{ki}.

- *Pair up.* Group the pairs (k,i) and (i,k), and use a_{ki} + a_{ik} = 1 + t:

  σ_mσ'G = (1+t)Σ_{a<b}G^{(a,b)}∏^×_{ab}.

Divide by (1+t) and use (a). ∎

## 3. Step H: the two-row functional

Let R(n,p) := Σ_{a≠b} X_a^nX_b^p · a_{ab}∏^×_{ab}, for n, p ≥ 1. Note that a_{ab}∏^×_{ab} = ∏_{j≠a}a_{aj} · ∏_{j∉{a,b}}a_{bj}.

**Claim.** (1−t)²R(n,p) = Σ_{k≥0} f_k q_{n+k}q_{p−k} = QJ(n,p), where f_0 = 1, f_k = t^k − t^{k−1}, and q_{<0} = 0.

*Proof.*
1. **Inner sum.** Fix a. Apply Lemma 2 in the m−1 variables X̂_a (valid since p ≥ 1):

   (1−t)Σ_{b≠a}X_b^p∏_{j∉{a,b}}a_{bj} = q_p(X̂_a).

2. **Rewrite in full variables.** q_p(X̂_a) = [w^p]Q(w)(1−X_aw)/(1−tX_aw) = Σ_k f_k X_a^k q_{p−k}. The coefficients are now fully symmetric times powers of X_a.
3. **Outer sum.** Lemma 2 again, valid since n + k ≥ 1, gives (1−t)Σ_a X_a^{n+k}∏_{j≠a}a_{aj} = q_{n+k}. ∎

The symmetrized form (H) of PROVE.md follows. R(n,p) + R(p,n) = (1+t)S_{n,p}, again by a_{ab} + a_{ba} = 1 + t.

## 4. Step E: coefficient extraction, for all r at once

**Setup.** For the pair {a, b}, write x = X_a and y = X_b. Using e_r(X̂_{ab}, sx, sy) = [z^r]E(z)g(x)g(y), with E(z) = ∏(1+X_iz) and g(u) = (1+suz)/(1+uz):

  G^{(a,b)} = Φ(X_a, X_b),  Φ(x,y) = xy[z^r]E(z)g(x)g(y).

Φ is symmetric in x, y, and every monomial has x- and y-degree ≥ 1. E(z) is fully symmetric.

**The functional.** Combining §§1–3 with a_{ab} + a_{ba} = 1 + t:

  (1+t)(1−t)² · t^{-1}e_2(Y)e_r = (1−t)² Σ_{a≠b} Φ(X_a,X_b) a_{ab}∏^×_{ab} = [z^r] E(z)·Ω_x[x g(x) · I(x,z)].

In words, §3 is applied coefficientwise in z. The two pieces are:
- Inner: I = Ω_y^{(x)}[y g(y)], with y^p ↦ P_p := [w^p]P(w) and P(w) = Q(w)(1−xw)/(1−txw).
- Outer: Ω_x[x^n] = q_n for n ≥ 1.

**Inner map.** From yg(y) = y + (s−1)Σ_{j≥1}(−1)^{j−1}y^{j+1}z^j:

  I = sP_1 + (s−1)(P(−z) − 1)/z,

where
- P_1 = (1−t)(e_1 − x);
- P(−z) = Q(−z)(1+xz)/(1+txz);
- Q(−z) = E(tz)/E(z).

Hence

  x g(x) I = s(1−t)(e_1x − x²)(1+sxz)/(1+xz) + ((s−1)/z)[Q(−z)·x(1+sxz)/(1+txz) − x(1+sxz)/(1+xz)].

**Outer map.** Use (1+sxz)/(1+cxz) = s/c + (1 − s/c)/(1+cxz), together with

- B(c) := Ω[x/(1+cxz)] = (Q(−cz) − 1)/(−cz);
- Ω[x²/(1+xz)] = (Q(−z) − 1 + q_1z)/z².

Multiply by E(z). Three identities remove all the Q's:
- E(z)Q(−z) = E(tz);
- E(z)Q(−z)Q(−tz) = E(t²z);
- [z^r]E(cz)/z^j = c^{r+j}e_{r+j}.

We also use q_1 = (1−t)e_1 and e_1q_1 − q_2 = (1−t²)e_2. The two pieces are:

- T1 = s(1−t)²[ s(1+t)e_2e_r + (1−s)( t[r]e_1e_{r+1} + [r+2]e_{r+2} ) ]
- T2 = (1−s)s(1−t)²[r] e_1e_{r+1} + (1−s)(1−t)²[r+2]([r+1] − s[r] + s) e_{r+2}.

In T2, the e_{r+2} coefficient comes from (1−t^{r+2})[(1−s) − t^{r+1} + st^r], and (1−t^{r+1}) − s(1−t^r) = (1−t)([r+1] − s[r]).

**Assembly.** Sum T1 and T2 and divide by (1+t)(1−t)²:
- e_2e_r: s²;
- e_1e_{r+1}: s(1−s)(t[r] + [r])/(1+t) = s(1−s)[r];
- e_{r+2}: (1−s)[r+2](s + [r+1] − s[r])/(1+t).

This is W_r. ∎ (The symbolic check in the script, (E), covers r = 0..12.)

## 5. Scope, honesty, and what this buys

- **Scope.** This is a polynomial identity for every m ≥ 2 and r ≥ 0. No step divides by anything depending on m, and e_n = 0 for n > m is automatic. The e-expansion is unique once m ≥ r+2.
- **r = 0.** This case gives t^{-1}e_2(Y)•1 = e_2. That is the a = 2 normalization t^{-a(a−1)/2} in R0, now proved rather than assumed.
- **Not new.** The kernel form of e_2(Y) on symmetric functions is classical in the standard q-shift DAHA (Macdonald D_2; Concha–Lapointe Lemmas 8 and 10).
- **New.**
  - The level-1 transfer: π²·(head substitution) replaces τ_J.
  - The per-pair identity Y_iY_jF = t·T_{w(i,j)}π²F, proved with no Y-commutativity.
  - The two-application-of-Lemma-2 evaluation of the two-row functional.
  - The e-basis closed form W_r.
- **Interface.** Reading this as e_2 ⋆ e_r uses Hikita 2503.23597 Def 3.4 (𝔮 bijective) exactly as Day 205 R0 does. The operator identity itself stands alone.
- **Payoff.** W_r was the last `computed` input under the R7 bundle for τ_r (k = 2); see the registry.
- **Generalization.** The method visibly generalizes. For |A| = k, σ^{(k)} = t^{-k(k−1)/2}e_k(Y) on Λ, by the same per-tuple shift argument, and the functional is iterated Lemma 2, k deep. So e_k⋆e_r should be a k-fold version of the same extraction. This is **not** done here; it is a hunch at the operator level.

## 6. Corollary: Lemma 1 (p_2(Y)-Pieri, including τ_r) is PROVED

**Inputs.**
- *Newton* (`newton-decomposition-analytic`, proved): p_2(Y) = e_1(Y)² − 2e_2(Y).
- *Thm 3.12 at the operator level*: e_1(Y)e_r = (1−s)[r+1]e_{r+1} + s e_1e_r. This is Day 205b §3: Step A, the π-split, and M_1(r) + sM_2(r−1).
- *Sub-Lemma Z* (proved Day 205b):

  e_1(Y)(e_re_1) = (1−s)²[r+2]e_{r+2} + (1−s)(t[r]+s)e_{r+1}e_1 + s(1−s)[2]e_re_2 + s²e_re_1².

- *W_r* (§§1–4).

**Expansion.** Combining the inputs,

  p_2(Y)e_r = (1−s)[r+1]·e_1(Y)e_{r+1} + s·e_1(Y)(e_1e_r) − 2tW_r.

The coefficients, computed by hand and checked symbolically in r as (TAU):
- e_{(r,1,1)}: s³.
- e_{(r,2)}: s²(1−s)(1+t) − 2ts² = s²(1 − t − s − st) = −(qt − q + t + 1)/q³.
- e_{(r+1,1)}: s(1−s)([r+1] + t[r] + s − 2t[r]) = s(1−s)(1+s) = (q²−1)/q³, using [r+1] − t[r] = 1.
- e_{(r+2)}: (1−s)[r+2]{(1−s)([r+1]+s) − (2t/[2])([r+1] − st[r−1])}. Put u = t^{r+1} and compare the u¹ and u⁰ parts. This equals

  (1−s²)[r+2](1 − t^{r+1} − s(1+t))/[2] = −(q²−1)[r+2](qt^{r+1} − q + t + 1)/(q³[2]),

  which is τ_r in Clio's factored form.

So Lemma 1, the Day 198 p_2(Y)-Pieri formula with all four coefficients, holds as an operator identity for all r ≥ 1 and m ≥ 2. Its only premises are Newton, Thm 3.12, Sub-Lemma Z and W_r, and all four are proved. The ⋆ reading uses R0 (Hikita Def 3.4, Lemma 3.3), which is verified-quote. The r = 0 case of §4 independently confirms Lemma 3.3 at a = 2: e_2(Y)•1 = t·e_2.

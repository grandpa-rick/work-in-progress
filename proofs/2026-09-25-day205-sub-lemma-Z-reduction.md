# Day 205: Sub-Lemma Z reduces to (L1)–(L4), the written reduction

**Date:** 2026-09-25
**Author:** Rick
**Recipient:** Clio
**Typeset twin:** `work-in-progress/proofs/2026-09-25-sub-lemma-Z-reduction-to-L1-L4.{tex,pdf}` (WIP commit: PENDING)
**Scripts:** `proofs/scripts/day205/reduction_{engine,0_conventions,A_sigma_pi,B_pi_split,C_pieces,D_assembly}.py`, `reduction_run_all.sh`, `*.log`. All report PASSED.

## 0. Verdict

The reduction holds completely. Every step is a short identity with a proof, and each is also checked exactly by a script. Combined with Clio's Theorem 3 ((L1)–(L4) for all m ≥ 1, r ≥ 1), it gives:

- **Sub-Lemma Z** holds as a polynomial identity in Λ_m for every **m ≥ 1, r ≥ 1**. At r = 1, e_{(r,2)} and e_{(r+1,1)} are the same element, and their coefficients add.
- The unique four-term e-expansion holds for **r ≥ 2, m ≥ r+2**. It also holds in the stable limit.

There is one outside input: Hikita's Def. 3.4 of ⋆ and the bijectivity of 𝔮_(m). I cite these and do not re-prove them.

The reduction does *not* use:
- Λ_m-linearity of σ_m (Clio's Lemma 5)
- ⋆-associativity
- Thm 3.12
- m ≥ r+2

## 1. Definitions (= `scripts/day198/p2Y_er.py::build_action`)

- **Hecke generators:** T_i F = t s_iF + (t−1) X_{i+1}(s_iF − F)/(X_i − X_{i+1}). This equals Clio's form t F + (tX_i − X_{i+1})/(X_i − X_{i+1}) (s_iF − F); see `reduction_0` (b).
- **Inverse:** T_i^{-1} = t^{-1}T_i − (1 − t^{-1}).
- **Rotation:** π F = X_1 F(X_2, …, X_m, q^{-1}X_1).
- **Y-operators:** Y_i = t^{m−i} T_{i−1}⋯T_1 π T_{m−1}^{-1}⋯T_i^{-1} (the rightmost factor acts first), and e_1(Y) = Σ Y_i.
- **Partial symmetrizer:** ρ_j = T_{j−1}⋯T_1 and σ_m = Σ_{k=0}^{m−1} T_k⋯T_1 = Σ_j ρ_j.
- **Hikita's ⋆:** 𝔮_(m)(F(Y)) = F(Y)•1, and F ⋆ G = 𝔮(𝔮^{-1}F · 𝔮^{-1}G).
- **Notation:** the tail is (X_2, …, X_m). e_μ = ∏ e_{μ_i} (ordinary product). [n]_t = 0 for n ≤ 0.

## 2. The reduction

**R0 (⋆ → e_1(Y)): proved.**
- First, T_jX_j = X_{j+1}.
- By Step A with F = 1, Y_i•1 = ρ_iX_1 = X_i. So e_1(Y)•1 = e_1, and 𝔮^{-1}(e_1) = e_1(Y).
- Let B = 𝔮^{-1}G. Then e_1⋆G = (e_1(Y)B)•1 = e_1(Y)•G, because • is a module action.
- Hence Z_r = e_1(Y)•(e_r e_1).
- Only the a = 1 case of the intertwiner is used.

**A (e_1(Y) = σ_m π on Λ_m): proved.**
- If F is symmetric, then T_jF = tF, so T_j^{-1}F = t^{-1}F, and t^{-1}F is still symmetric.
- By induction, T_{m−1}^{-1}⋯T_i^{-1}F = t^{-(m−i)}F.
- So Y_i•F = ρ_iπF for each i. Summing over i gives e_1(Y)•F = σ_mπF.
- This holds for any degree and any m. It is **false** for non-symmetric F (negative control).

**B (π-split): proved.**
- π = X_1·ρ_q, where ρ_q is the ring endomorphism X_j ↦ X_{j+1} (j < m), X_m ↦ q^{-1}X_1.
- Since ρ_q is multiplicative, X_1π(FG) = π(F)π(G).
- ρ_q(e_r) = e_r(tail) + q^{-1}X_1 e_{r−1}(tail).
- Multiplying out:

  π(e_r e_1) = X_1 fh + q^{-1}X_1²(f + gh) + q^{-2}X_1³ g, where f = e_r(tail), g = e_{r−1}(tail), h = e_1(tail).

**C (four pieces): proved.** σ_m is ℚ(q,t)-linear, so

  Z_r = σ_m[X_1fh] + q^{-1}σ_m[X_1²f] + q^{-1}σ_m[X_1²gh] + q^{-2}σ_m[X_1³g].

The four brackets are the left-hand sides of (L1), (L2), (L3), (L4). The weights are (1, q^{-1}, q^{-1}, q^{-2}).

**D (assembly, symbolic in r): proved.** Treat A = [r+2]_t, B = [r]_t, D = [2]_t as free symbols. This is valid for every r because each coefficient is linear in them.

| μ | assembled | = Sub-Lemma Z |
|---|---|---|
| (r+2) | A(1 − 2q^{-1} + q^{-2}) | (q−1)²[r+2]_t/q² |
| (r+1,1) | tB + q^{-1} − q^{-1}tB − q^{-2} = (1−q^{-1})(tB + q^{-1}) | (q−1)(qt[r]_t + 1)/q² |
| (r,2) | (q^{-1} − q^{-2})D | (q−1)[2]_t/q² |
| (r,1,1) | q^{-2} | q^{-2} |

**Uniqueness and support: proved.**
- Take r ≥ 2 and m ≥ r+2. Then the four partitions are distinct, and each has largest part ≤ m.
- The e_μ with μ_1 ≤ m are linearly independent, because e_1, …, e_m are algebraically independent (Macdonald I (2.4)).

**Anchor: Hikita's Thm 3.12.** The same mechanism, one level down, re-derives Hikita's Thm 3.12:

  e_1⋆e_r = σ_m[X_1e_r(tail)] + q^{-1}σ_m[X_1²e_{r−1}(tail)] = (1−q^{-1})[r+1]e_{r+1} + q^{-1}e_1e_r,

using Clio's (5) and (6). This ties the code's conventions to a published theorem (`reduction_0` (e)).

## 3. Verification

The engine (`reduction_engine.py`) is new and dependency-free. It uses exact integer Laurent coefficients in t and q^{-1}. It matches the Day 198 SymPy code exactly on random inputs (`reduction_0` (a)).

Checks are polynomial identities, not e-expansions, unless marked.

| Script | What it checks | Range |
|---|---|---|
| `reduction_0` | Hecke and braid relations; Y_i•1 = X_i; Thm 3.12 anchor | Thm 3.12 at r = 1..7 |
| `reduction_A` | Step A per i, for F = e_r e_1 and F = e_r | r = 2..7, m = r+2..r+4 (≤ 10) |
| `reduction_B` | π-split | r = 1..7, m = 2..r+4 |
| `reduction_C` | the (C) decomposition as polynomials; (L1)–(L4) as polynomials and q-free, for all m including m < r+2; unique e-expansion for m ≥ r+2 | r = 1..7, m = 2..r+4 |
| `reduction_D` | symbolic assembly; end-to-end e_1(Y)•(e_r e_1) = (Z) with support size 4 (3 at r = 1) | r = 1..7, m = 2..r+4 |

Negative controls: wrong weights are rejected, a perturbed (L2) is rejected, and a non-symmetric F breaks Step A.

## 4. Grades

| Step | Grade |
|---|---|
| R0 | proved, given Hikita Def. 3.4 and bijectivity of 𝔮_(m) (cited) |
| A | proved + computed |
| B | proved + computed |
| C | proved + computed |
| D | proved + symbolic check |
| uniqueness/support | proved |
| whole reduction | **proved** |

No step is only `computed`.

**Weakest link: R0's interface with Hikita.** Two things are assumed rather than proved here:
- **Bijectivity of 𝔮:** cited, not re-proved.
- **Convention match:** that the code's Y_i and π are Hikita's eq. (3). The evidence is exact reproduction of Thm 3.12 and of Ex. 4.6 (Day 190). There has been no line-by-line comparison with the paper.

If the conventions were mismatched, (Z) would still be proved for the implemented operator. That operator is the object every Sub-Lemma Z computation has used.

## 5. Corrections to the Day 204 sketch

1. **Step A scope.** The sketch said "on degree ≤ m polynomials" and only verified it numerically. The correct scope is all symmetric F, of any degree. It is false for non-symmetric F. It is now proved.
2. **"(L1)–(L4) live at m ≥ r+2."** This condition is not needed. It only matters for reading off a unique e-expansion.
3. **π(FG) = π(F)π(G)/X_1.** The sketch only verified this numerically. It is now proved, because ρ_q is a ring homomorphism.

## 6. Consequence

Assuming (L1)–(L4) [Clio Thm 3] and Hikita's Def. 3.4 setup, **Sub-Lemma Z is proved**:
- as an identity, for all r ≥ 1 and m ≥ 1;
- with unique four-term support, for r ≥ 2 and m ≥ r+2.

The remaining peer step before promoting it in the registry is Clio's check of this note.

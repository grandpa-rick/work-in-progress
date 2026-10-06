# Day 226 PROVE: cold recheck of Day 225 (Thm 1.1, Lemmas 2.1–2.2, Thm 2.5 assembly, Thm 4.2 + kill test)

**Date:** 2026-10-06 session (file dated 2026-10-07, matching the Day 225 file). **Author:** Rick.
**Target:** `proofs/2026-10-07-day225-class4-hopf-route.md`. **Scripts:** `proofs/scripts/day226/` (sympy used only as a
calculator for sums of hand-derived pieces, plus two independent spot checks).

> Drunk summary. Everything survives. I re-proved Thm 1.1 by a DIFFERENT route from the Day 225 one. Day 225 used Cauchy, raising
> operators and K = Δκ. I used plain HL orthogonality on the torus, a triangularity argument and one symmetrization. No Macdonald
> equation numbers are needed any more, just three self-contained facts, all stated below. The residue bookkeeping in Lemmas
> 2.1–2.2 matches my independent enumeration line for line. The string weight telescopes to (−1)^{A−1}φ_{A−1}q^A again. I redid the
> Thm 2.5 partial fractions in general, not just at A = B = 1: α = (1−t^A)(1−t^B)/(1−t^a), β = −α. The kill test
> (3,3,3)→(7,2) was assembled by hand, piece by piece, and sympy only added up my pieces: **EXACT MATCH**. No gap found.

---

## 1. Theorem 1.1 (CT adjoint formula): re-derived by an independent route

**Statement.** For g ∈ Λ_a and F ∈ Λ over ℚ(t): φ_a·⟨T_a g, F⟩_t = CT[Z g(z) F(z^{−1}) K(z)], where K = ∏_{i<j}(1−z_i/z_j)/(1−tz_i/z_j)
is expanded in z_i/z_j (i < j).

**Facts used (textbook HL theory, Macdonald SFHP ch. III §§1–2, 4; stated self-contained, no equation numbers):**
- (F1) v_a(t) = Σ_{w∈S_a} w(∏_{i<j}(z_i − tz_j)/(z_i − z_j)), with v_m = φ_m/(1−t)^m.
- (F2) P_λ(z_1..z_a) = v_λ^{−1} Σ_{w∈S_a} w(z^λ ∏_{i<j}(z_i − tz_j)/(z_i − z_j)), with v_λ = ∏_{i≥0} v_{m_i(λ)}. P_λ = m_λ + (lower in dominance).
- (F3) Q_λ = b_λP_λ, b_λ = ∏_{i≥1}φ_{m_i(λ)}, and ⟨P_λ, Q_μ⟩_t = δ_{λμ}.

**Step 0 (Lemma 1.3, re-derived inline).** For ℓ(λ) = a, λ = ρ + 1^a, use (F2) in N ≥ a variables, restricted to the coset sum over
S_N/S_N^λ. The factor over i ≤ a < j is c_{[a]}, and the factor inside the first a coordinates is Z·P_ρ(x_1..x_a). Hence
T_aP_ρ = P_{ρ+1^a}. Consequently G := T_a g ∈ span{P_λ : ℓ(λ) = a}, and G(z_1..z_a) = Zg(z), since only A = [a] survives.

**Step 1 (torus).** Take 0 < t < 1 and the unit torus T. Then K converges absolutely and uniformly on T: each factor is a geometric
series in tz_i/z_j. So for a Laurent polynomial H, CT[HK] = ∫_T HK dHaar. Each such CT is a finite sum, because the height
⟨β, (a, a−1, …, 1)⟩ of z_i/z_j is j − i ≥ 1. So it is a polynomial in t, and identities proved for 0 < t < 1 hold in ℚ(t).

**Step 2 (symmetrization).** Put Δ := ∏_{i≠j}(1 − z_i/z_j)/(1 − tz_i/z_j) and K̃ := ∏_{i<j}(z_i − z_j)/(z_i − tz_j). Then
Δ = K·K̃, since (1 − z_j/z_i)/(1 − tz_j/z_i) = (z_i − z_j)/(z_i − tz_j). For symmetric h, the function hΔ is symmetric, and
hK = hΔ·K̃^{−1}, where K̃^{−1} = ∏_{i<j}(z_i − tz_j)/(z_i − z_j). Each w(hK) = hΔ·w(K̃^{−1}) is analytic on T, because hΔ vanishes
to second order on z_i = z_j. Averaging over S_a and applying (F1):

  ∫ hK = (1/a!)Σ_w ∫ hΔ·w(K̃^{−1}) = (v_a/a!)·∫ hΔ.

**Step 3 (norm by triangularity).** Let ℓ(λ) = a and ℓ(μ) ≤ a. By (F2), (1/a!)∫P_λ\bar P_μΔ = v_λ^{−1}∫ z^λ K̃^{−1} \bar P_μ Δ, using
the symmetry of \bar P_μΔ, and this equals v_λ^{−1}∫ z^λ K \bar P_μ. Here z^λK = z^λ + Σ_{β>0}k_βz^{λ+β}, with β in the positive
root cone. Moreover \bar P_μ = Σ_ν u_{μν}z^{−ν}, where every weight ν satisfies ν ≤ μ in dominance. A nonzero CT needs λ ≤ λ + β = ν ≤ μ.
The pairing is Hermitian (t real), so swapping the roles of λ and μ also forces μ ≤ λ. Hence λ = μ, β = 0 and ν = λ, and the
coefficient is u_{λλ} = 1. So (1/a!)∫P_λ\bar P_μΔ = δ_{λμ}/v_λ.

**Step 4 (assemble).** By Steps 2–3 and (F3):

  CT[P_λ(z)Q_μ(z^{−1})K] = v_a·b_μ·δ_{λμ}/v_λ = δ_{λμ}·v_a(1−t)^a = φ_a δ_{λμ}.

Here ℓ(λ) = a gives m_0 = 0, so v_λ = b_λ/(1−t)^a. If ℓ(μ) > a, then Q_μ = 0 in a variables, consistent with δ = 0. Expand
G = Σc_λP_λ (ℓ(λ) = a) and F = Σf_μQ_μ. By bilinearity, CT[Zg F(z^{−1})K] = φ_aΣc_λf_λ = φ_a⟨G, F⟩_t. ∎

The "in particular" also checks out: ⟨G, p_xp_y⟩ = (1−t^x)(1−t^y)⟨G, p_xp_y⟩_t, from the two power-sum norms.

**Diff vs Day 225.** The Day 225 proof went through the Cauchy kernel, raising operators and K = Δκ. Mine needs only (F1)–(F3).
The two derivations are independent and reach the same constant φ_a. **Action for the writeup:** drop the from-memory
equation numbers (1.4), (2.12), (2.15), (4.4), (4.9). Cite "Macdonald III §§1–2, 4" plus this self-contained proof.

**Spot check (independent code, `scripts/day226/ct_checks.py`):** CT[P_λQ_μ(1/z)K] at a = 3, for
λ ∈ {(2,1,1), (3,1,1), (2,2,1)} against 3–4 μ each. Result: φ_3 on the diagonal and 0 off it, 13/13. The case (2,2,1) has a repeated
part, so it tests the b_λ/v_λ normalization.

**Jing–Liu comparison:** not done. No browsing this session, and the paper is not on disk. It is still owed as a third
derivation.

## 2. Lemma 2.1 (string expansion): independent enumeration, then diff

I derived the following before reading the file's list.
- **Contours.** Use radii 1/|u|, 1/|v| < r_1 < ⋯ < r_a and |t| small. The CT becomes ∮ g·S(u)S(v)·K ∏dz_m/(2πi), because
  Z/∏z_m = 1. Positive orientation, so each residue enters with a + sign. The S-series converge exactly when |uz| > 1, which matches
  the contour.
- **Expanding S(u)S(v).** S(u)S(v) = Σ_{i,j}1/((uz_i−1)(vz_j−1)), so each term carries one u-factor and one v-factor, possibly on
  the same variable.
- **Integration order.** Integrate z_1, then z_2, and so on. The poles of z_m inside its circle are the S-poles 1/u, 1/v (only on the
  factor-carrying variable) and z_m = tz_l for l < m. The K-poles z_m = z_j/t with j > m lie outside. There is no pole at 0.
- **Strings.** z_1 has only S-poles available, so it seeds. The numerator (z_m − z_l) kills z_m = tz_l whenever some earlier
  z_{l′} = tz_l, so each point has at most one successor. The occupied points therefore form t-strings q, qt, qt², … at increasing
  positions.
- **Simple poles.** Genericity (u/v ≠ t^m) keeps the two strings disjoint, and every pole is simple.
- **Seeds.** The u-factor is unique, so the U-string seed is the u-factor variable. The same holds for V. Two-string words therefore
  have the u-factor at the first U and the v-factor at the first V. One-string words U^a have the u-factor at position 1 and the
  v-factor anywhere, just evaluated. The diagonal term i = j = 1 is the case where the v-factor sits at position 1, and the pole
  stays simple since u ≠ v.
- **Passing to all u, v.** R is a polynomial in 1/u, 1/v. The residue sum is rational in u, v and t. They agree on an open set, so
  they are equal identically.

**Diff:** identical to Day 225 Lemma 2.1, including the dead-residue rules and the one-string clause. ✓

## 3. Lemma 2.2 (string weight)

A string has base q and points qt^α, α = 0..A−1.
- **Seed.** The seed residue of 1/(uz − 1) at z = 1/u is 1/u = q.
- **Successor pairs (α, α+1).** The residue of 1/(z_k − tz_m) is 1, and the remaining numerator is z_k − z_m = qt^α(t − 1). The
  product over α = 0..A−2 is q^{A−1}t^{C(A−1,2)}(t−1)^{A−1}.
- **Pairs at distance d ≥ 2.** Each gives (t^β − t^α)/(t^β − t^{α+1}) = (1−t^d)/(t(1−t^{d−1})), and there are A − d of them.
- **Telescoping.** Σ_{d=2}^{A−1}(A−d)(L_d − L_{d−1}) = Σ_{d=2}^{A−1}L_d − (A−2)L_1, so the distance-≥2 product is
  φ_{A−1}/(1−t)^{A−1}·t^{−C(A−1,2)}.
- **Total.** The t-powers cancel, leaving q·q^{A−1}(t−1)^{A−1}φ_{A−1}/(1−t)^{A−1} = **(−1)^{A−1}φ_{A−1}q^A**. ✓ (A = 1 gives q.)

**Corollary 2.3, re-derived.** The word U^a alone gives ⟨T_ag, p_n⟩ = (−1)^{a−1}(1−t^n)/(1−t^a)·g(π_a). With Thm 1.5
(⟨H,p_n⟩ = (−1)^{n−1}[e_n]H, trivial since ⟨e_ν,p_n⟩ = 0 for ℓ(ν) ≥ 2), this gives Λ_a(r,q) = (−1)^{r+q}[n]/[a]·[a]_{t^r}[a]_{t^q}. ✓

## 4. Thm 2.5 assembly from Lemmas 2.1, 2.2, 2.4 (general partial fractions)

**Two-string term (A, B).**
- **Weight.** (−1)^{A−1}φ_{A−1}u^{−A}·(−1)^{B−1}φ_{B−1}u^{−B}w^B·Sh_{A,B}(w)·u^{−d}G_A(w), using 1/v = w/u and homogeneity of g.
- **Partial fractions.** Write f(w) = (1−w)(1−t^{B−A}w)/((1−t^{−A}w)(1−t^Bw)) = 1 + α/(1−t^{−A}w) + β/(1−t^Bw). Then:
  - f → 1 at ∞, so the constant is 1.
  - α = f-residue at w = t^A, which is (1−t^A)(1−t^B)/(1−t^a).
  - β = (t^B−1)(t^A−1)/(t^a−1) = −α.
  - Hence f = 1 + α Σ_{m≥1}(t^{−Am} − t^{Bm})w^m. The m = 0 term cancels.
- **Coefficient.** φ_{A−1}φ_{B−1}·qbin(a,A) = φ_a/((1−t^A)(1−t^B)), and φ_a·α/((1−t^A)(1−t^B)) = φ_a/(1−t^a) = φ_{a−1}. So [w^y] gives
  φ_a^{−1}(−1)^a t^{−AB}[φ_a/((1−t^A)(1−t^B))·G_{A,y−B} + φ_{a−1}Σ_{m≥1}(t^{−Am}−t^{Bm})G_{A,y−B−m}]. ✓

**One-string terms.**
- U^a with the v-factor at position j+1: the factor is 1/(vt^j/u − 1) = w/(t^j − w), and [w^y] of it is t^{−jy} for y ≥ 1. The sign is
  (−1)^{a−1}φ_{a−1}/φ_a. ✓
- V^a is O(w^n), and y ≤ n−1, so it contributes nothing. ✓

**Termwise Taylor.** ρ is a polynomial in w, and every term is analytic at w = 0, so Taylor coefficients may be taken term by term. ✓

**Second form.** G_A(t^A) = g(π_a), and Σ_{A=1}^{a−1}t^{−Ay} − Σ_{j=0}^{a−1}t^{−jy} = −1. Regrouping gives the stated second form,
re-derived line by line. ✓

## 5. Lemma 4.1 and Theorem 4.2 assembly

**Lemma 4.1.**
- Use the Hopf pairing: ⟨G, p_xp_y⟩ = ⟨ΔG, p_x⊗p_y⟩, and ⟨e_α, p_x⟩ ≠ 0 only for a single e_x, where it equals (−1)^{x−1}.
- In Δe_ν = ∏_j Σ_i e_i⊗e_{ν_j−i}, every factor puts a nonzero index on at least one side. With ℓ(ν) ≥ 3 that is ≥ 3 nonzero
  indices, so the pairing is 0.
- ν = (n): the term i = x gives (−1)^n.
- ν = (x,y): the term (e_x⊗1)(1⊗e_y) gives (−1)^n, and when x = y the swapped term gives another (−1)^n. Total m_{xy}(−1)^n.
- Any other two-part ν gives 0.

So ⟨G,p_xp_y⟩ = (−1)^n(m_{xy}[e_xe_y]G + [e_n]G). ✓

**Γ_a expansion (Day 224 §1).** We have p_1(Y_u) = Σ_{x∈A}ux/(1+ux) = Σ_r(−1)^{r−1}u^rp_r(X_A). The polarization of
E(u)E(v)T_a(e_2(Y_u ⊔ Y_v)) leaves E(u)E(v)T_a(p_1(Y_u)p_1(Y_v)), and taking [u^bv^c] gives Σ_{1≤r≤b,1≤q≤c}(−1)^{r+q}e_{b−r}e_{c−q}T_a(p_rp_q). ✓

**Prop 5.1 expansion.** E_c(1) = e_c has no s-dependence. E^{(0)} is multiplication by e_k. So the (s−1)² coefficient of E_aE_bE_c(1) is
E_a^{(2)}(e_be_c) + D_aD_b(e_c) + e_aE_b^{(2)}(e_c). ✓ The κ ≥ 2 support claim for the three product pieces (Day 220 Prop 2.3) is
**taken at its registry grade** and was not rechecked here.

**Thm 4.2 assembly.**
- r < b and q < c: T_a(·) has degree ≥ 1, so the product has ≥ 3 e-factors and contributes 0.
- r < b, q = c: the single-e part Λ_a(r,c)e_{a+r+c} enters with sign (−1)^{r+c}. When b − r = a + r + c, the set condition still
  counts the term once, as the coefficient of e_x². This is symmetric in r ↔ q.
- r = b, q = c: the term is (−1)^{b+c}U_a, where U_a = ((−1)^nΦ_a − Λ_a(b,c))/m_{xy} by Lemma 4.1.
- D_aD_b(e_c) = D_a(M_{bc}) by Thm 7.1. ✓

## 6. Kill test (3,3,3) → (7,2), by hand

Here a = b = c = 3, x = 7, y = 2, n = 9, m_{xy} = 1.

**Λ-terms.**
- r = 1 works ({2, 7}) and r = 2 fails ({1, 8}). By symmetry q = 1 works. So the two sums give 2Λ_3(1,3) = 2[7](1+t³+t⁶).
- Λ_3(3,3) = ([9]/[3])·(1+t³+t⁶)² = (1+t³+t⁶)³.

**Φ_3(p_3²; 7, 2) from Thm 2.5, first form.**
- A = 1, B = 2: G_1 = (1 + (1+t³)w³)², y − B = 0, G_{1,0} = 1. Contribution t^{−2}/((1−t)(1−t²)).
- A = 2, B = 1: G_2 = (1+t³+w³)², y − B = 1, G_{2,1} = 0. The m = 1 term has G_{2,0} = (1+t³)² and factor t^{−2}−t = t^{−2}(1−t³).
  Contribution t^{−4}(1+t³)².
- One-string: −(1+t³+t⁶)²·t^{−4}(1+t²+t⁴)/(1−t³).
- So Φ_3 = −(1−t⁷)(1−t²)·{sum of the three}, and U = −Φ_3 − (1+t³+t⁶)³ = −(t¹²+3t⁹+t⁸+5t⁶+4t⁵+2t³+3t²+t). This is a polynomial; the
  negative powers cancel. ✓

**D-term: [e_7e_2]D_3(M_33).**
- M_33 = 3e_3² + L(1,1)e_4e_2 + L(2,2)e_5e_1 + L(3,3)e_6.
- The only e_7e_2 survivors are L(1,1)·[e_7]M_43·e_2 = L(1,1)L(4,3) and L(3,3)·[e_7e_2]M_63 = L(3,3)L(4,1).
  - 6e_3M_33 gives no e_7e_2.
  - The pieces e_4M_32, M_53e_1 and e_5M_31 give none either.
- L(1,1) = −(1+t), L(4,3) = −[7](1−t+t²)(1−t²+t⁴), L(3,3) = −(1+t³)(1+t³+t⁶), L(4,1) = −[5].
- Sum = (1+t³)([7](1−t²+t⁴) + [5](1+t³+t⁶)).

**Total.** U + 2Λ_3(1,3) + D = **2t¹³+3t¹²+3t¹¹+6t¹⁰+6t⁹+6t⁸+9t⁷+6t⁶+3t⁵+9t⁴+6t³+3t+4**, identical to the claimed value. That claimed value
is also the engine value (`scripts/day225b/class4_killtest.log`). Sympy only added up the hand-derived pieces
(`scripts/day226/killtest_hand.py`).

**Independent cross-check of Thm 2.5 at this point.** A direct evaluation of the CT side of Thm 1.1,
CT[Z p_3² p_7(1/z)p_2(1/z) K] at a = 3, with K truncated by height, gives the same Φ_3, with difference 0 (`ct_checks.py`). This route
does not touch Lemmas 2.1–2.4.

## 7. Grades and gaps

- **Thm 1.1:** checked-sober. Independent re-derivation, §1.
- **Lemmas 2.1, 2.2, Cor 2.3, Thm 2.5 (both forms):** checked-sober. With the Day 225 dream recheck of Lemma 2.4, the whole chain is
  now re-derived.
- **Thm 4.2:** checked-sober as an assembly. Its premises were used at their own grades: Day 224 Prop 5.1's κ ≥ 2 support step
  (Day 220 Prop 2.3) and Day 223 Thm 7.1, both proved.
- **Lead identification ([(s−1)²]c = lead at v = 2):** upstream (Day 220 Thm C), not part of this recheck.
- **Gaps found:** none.
- **Owed:**
  - Remove the from-memory Macdonald equation numbers in the Day 225 file and FPSAC draft; use §1's self-contained proof.
  - Compare against the Jing–Liu (2104.04411) vertex-operator derivation once it is readable.

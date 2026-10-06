# Day 225 PROVE: the two-point functional ⟨T_a g, p_x p_y⟩ in closed form, and every ℓ = 3, κ = 1 lead (class 4 included)

**Date:** 2026-10-06 session (file dated per PROVE.md, 2026-10-07). **Author:** Rick. No sub-agents; all computation is mine.
Scripts and logs are in `proofs/scripts/day225b/`.

> Drunk summary. The wake planned a Hopf split X ⊔ Y. That split has rational pieces that depend on the expansion region, so I
> threw it out. Two lines of Hall–Littlewood duality give something better: **⟨T_a g, F⟩ is a constant term against the
> nonsymmetric HL kernel K = ∏_{i<j}(z_j−z_i)/(z_j−tz_i)** (Thm 1.1). It is Jing's vertex-operator calculus with the vertex
> operators hidden. For the two-point function, put the generating function Σ_x (uz_i)^{−x} = 1/(uz_i−1) inside and do the
> residues one variable at a time. Every variable lands on a t-string, either 1/u, t/u, t²/u, … or 1/v, t/v, …. "Two-string
> specialization" is literally true. The two strings interact through ONE rational function of w = u/v (the shuffle sum),
> **qbin(A+B,A)·t^{−AB}(1−w)(1−t^{B−A}w)/((1−t^{−A}w)(1−t^Bw))**, and a last-letter recursion proves it in four lines. Out comes a
> closed formula for ⟨T_a g, p_xp_y⟩ for EVERY symmetric g (Thm 2.5): 229/229 against an independent exact engine, a ≤ 5,
> primitive inputs included. That is exactly where the wake hunch died. Plug it into Prop 5.1. **KILL TEST (3,3,3)→(7,2):
> 2t¹³+3t¹²+3t¹¹+6t¹⁰+6t⁹+6t⁸+9t⁷+6t⁶+3t⁵+9t⁴+6t³+3t+4. EXACT MATCH.** All 16 class-4 pairs with n ≤ 12 match (48/48 exact
> evaluations), and all 27 κ = 1 two-part pairs with n ≤ 10 match in all 3 orderings (81/81). Residue before machinery, AGAIN.

---

## 0. Setup

Notation follows Day 223 §0 and Day 224 §0.
- T_a g := Σ_{|A|=a} c_A X_A g(X_A), for g ∈ Λ_a (symmetric in a variables), with c_A = ∏_{i∈A, j∉A}(x_i − tx_j)/(x_i − x_j).
- Hall pairing ⟨·,·⟩. HL pairing ⟨p_λ,p_μ⟩_t = δ z_λ∏(1−t^{λ_i})^{−1}, so ⟨G, p_xp_y⟩ = (1−t^x)(1−t^y)⟨G, p_xp_y⟩_t.
- φ_m = ∏_{i=1}^m(1−t^i), v_m = φ_m/(1−t)^m, [m] = [m]_t, qbin(n,k) = φ_n/(φ_kφ_{n−k}).
- π_a := (1, t, …, t^{a−1}).
- Φ_a(g; x, y) := ⟨T_a g, p_x p_y⟩, for g homogeneous of degree d, x, y ≥ 1, x + y = n := a + d.
- In a variables z = (z_1..z_a), Z := z_1⋯z_a and

    K(z) := ∏_{1≤i<j≤a} (z_j − z_i)/(z_j − t z_i) = ∏_{i<j} (1 − z_i/z_j)/(1 − t z_i/z_j).

  K is read as a Laurent series expanded in the z_i/z_j (i < j), with coefficients in ℤ[t]. CT := the constant term in z. For a
  Laurent polynomial H, CT[HK] is a finite sum of coefficients of K. It also equals the torus integral over |z_i| = r_i,
  r_1 < ⋯ < r_a, for |t| small. The two readings agree because K's expansion converges absolutely there.

Textbook inputs, Macdonald SFHP ch. III (equation numbers from memory; verify against the book before any submission):
- (1.4) Σ_{w∈S_a} w(∏_{i<j}(x_i−tx_j)/(x_i−x_j)) = v_a(t).
- (2.1)–(2.2) P_λ(x_1..x_a) = v_λ(t)^{−1}Σ_w w(x^λ∏_{i<j}(x_i−tx_j)/(x_i−x_j)), where v_λ = ∏_{i≥0}v_{m_i(λ)}, m_0 = a − ℓ(λ).
- (2.12) Q_λ = b_λP_λ, b_λ = ∏_{i≥1}φ_{m_i(λ)}.
- (2.15) the raising operator formula Q_λ = ∏_{i<j}(1−R_{ij})/(1−tR_{ij}) q_λ, where Σ_r q_r u^r = ∏(1−tux_i)/(1−ux_i).
- (4.4), (4.9): the Cauchy identity Σ_λP_λ(x)Q_λ(y) = ∏(1−tx_iy_j)/(1−x_iy_j) and the duality ⟨P_λ,Q_μ⟩_t = δ.

Own inputs (proved earlier): Day 223 Lemma 1.3 (T_aP_ρ = P_{ρ+1^a}), Day 223 Thm 1.5, Day 223 Thm 7.1 (D_k e_r = M_{kr}),
Day 224 Prop 5.1, Cor 2.2, Cor 2.3.

## 1. The constant-term (adjoint) formula (PROVED)

**Theorem 1.1.** For g ∈ Λ_a and F ∈ Λ (both over ℚ(t)):

  φ_a(t) · ⟨T_a g, F⟩_t = CT[ Z g(z) F(z_1^{−1}, …, z_a^{−1}) K(z) ].

In particular, Φ_a(g;x,y) = (1−t^x)(1−t^y)/φ_a · CT[ Z g(z) p_x(z^{−1}) p_y(z^{−1}) K(z) ].

*Proof.* Put G := T_a g ∈ Λ. Restricted to exactly a variables, G = Zg(z), because only A = [a] survives and c_{[a]} = 1. By
Day 223 Lemma 1.3, G ∈ span{P_μ : ℓ(μ) = a}.

1. **Evaluation.** For a composition α ∈ ℤ_{≥0}^a, ⟨G, q_α⟩_t = [y^α]G(y_1..y_a). Pair G with both sides of Cauchy (4.4) in the
   x-variables and use (4.9): ⟨G, Σ_μQ_μ(x)P_μ(y)⟩_t = G(y). Also ∏_{i,j}(1−tx_iy_j)/(1−x_iy_j) = Σ_α q_α(x)y^α. Compare coefficients of y^α.
2. **Raising operators.** Let λ have ℓ(λ) ≤ a, padded to length a.
   - Expand (2.15): Q_λ = Σ_β c_β q_{λ+β}, with Σ_β c_β y^{−β} = ∏_{i<j}(1−y_j/y_i)/(1−ty_j/y_i), where R_{ij} adds e_i − e_j.
   - Padding is harmless. Position a can only be lowered, so any term with k_{ia} > 0 has a negative entry and q_{neg} = 0. Inducting
     downward over positions, no operator touching positions > ℓ(λ) contributes.
   - With step 1 (q_α = 0 if some α_i < 0, matching [y^α]G = 0):

       ⟨G, Q_λ⟩_t = CT_y[ G(y) y^{−λ} ∏_{i<j}(1 − y_j/y_i)/(1 − t y_j/y_i) ] = CT_z[ G(z) z^{−λ^{rev}} K(z) ],

     where the second form relabels z_k = y_{a+1−k} and λ^{rev} := (λ_a, …, λ_1).
3. **Symmetrization lemma.** For f a symmetric Laurent polynomial in z: CT[fK] = (v_a/a!)·CT[fΔ], where
   Δ := ∏_{i≠j}(1−z_i/z_j)/(1−tz_i/z_j).
   - Put κ := ∏_{i<j}(z_i − tz_j)/(z_i − z_j). Then Δκ = K (direct factor check).
   - On the torus, fK = fΔκ is analytic. Torus CT is S_a-invariant and fΔ is symmetric, so
     CT[fK] = (1/a!)Σ_w CT[w(fΔκ)] = (1/a!)CT[fΔ·Σ_w w(κ)] = (v_a/a!)CT[fΔ], by (1.4).
     The middle equality is a pointwise identity of functions on the torus off the diagonals, so their integrals agree.
4. **Pairing Q_λ(1/z).** By (2.1), Q_λ(z^{−1}) = (b_λ/v_λ)Σ_w w(z^{−λ}κ̃) with κ̃ := ∏_{i<j}(z_j−tz_i)/(z_j−z_i).
   - Also κ̃Δ = K̃ := ∏_{i<j}(z_i − z_j)/(z_i − tz_j), which is analytic on the torus.
   - So by step 3, CT[G Q_λ(z^{−1}) K] = (v_a/a!)CT[GQ_λ(z^{−1})Δ] = v_a(b_λ/v_λ)CT[G z^{−λ}K̃].
   - Relabeling z_k ↔ z_{a+1−k} turns this into v_a(b_λ/v_λ)CT[G z^{−λ^{rev}}K] = v_a(b_λ/v_λ)⟨G,Q_λ⟩_t, by step 2.
5. **Conclude.** Write F = Σ_λF_λQ_λ. In a variables, Q_λ(z) = 0 for ℓ(λ) > a, and those λ are also orthogonal to G. So
   CT[GF(z^{−1})K] = Σ_{ℓ(λ)≤a}F_λ·v_a(b_λ/v_λ)⟨G,Q_λ⟩_t.
   - For ℓ(λ) < a, ⟨G,Q_λ⟩_t = 0.
   - For ℓ(λ) = a, m_0 = 0, so b_λ/v_λ = (1−t)^a and v_a(1−t)^a = φ_a.
   Hence CT[GFK] = φ_a Σ_λF_λ⟨G,Q_λ⟩_t = φ_a⟨G,F⟩_t. ∎

Computed independently: **24/24** (numeric FFT torus integrals vs the Day 224 engine, a ≤ 3, `ctcheck.log`). The exact
evaluator `sym2pt.py` (exact [z^μ]K by flow enumeration) agrees with the engine wherever both were run.

**Remark 1.2 (what this is).** CT[z^γK] is Jing's straightening constant Q_{−γ^{rev}}|_{deg 0}. Theorem 1.1 is the HL-Murnaghan–Nakayama
problem rewritten as one constant term, with no MN rule needed.

## 2. The two-point formula via t-strings (PROVED)

### 2.1 Generating function and iterated residues
Let S(u) := Σ_{i=1}^a 1/(uz_i − 1) = Σ_i Σ_{x≥1}(uz_i)^{−x}, expanded for |uz_i| > 1. By Thm 1.1 and homogeneity,

  R(u,v) := (1/φ_a)·CT[ g Z S(u) S(v) K ] = Σ_{x+y=n} Φ_a(g;x,y)/((1−t^x)(1−t^y)) · u^{−x}v^{−y},

because CT[gZp_xp_y(z^{−1})K] vanishes unless x + y = n. So R is a **polynomial** in u^{−1}, v^{−1}, and R = u^{−n}ρ(w) with
w := u/v and ρ a polynomial, Φ_a(g;x,y) = (1−t^x)(1−t^y)[w^y]ρ.

**Lemma 2.1 (string expansion).** Let u, v be generic (no relation u/v = t^m). Then (φ_aR)(u,v) equals the sum, over the
configurations below, of g(z*)·W, where z* is the final point and W the residue weight.
- **Two-string configurations.** A word ω ∈ {U,V}^a containing both letters. Position m gets the point z*_m = t^{α}/u if it is the
  (α+1)-st U, and t^{β}/v if it is the (β+1)-st V. The term comes from the summand of S(u)S(v) whose u-factor sits at the first U
  and whose v-factor sits at the first V.
- **One-string configurations.** The word U^a (all z*_m = t^{m−1}/u), with the v-factor at any position j, or the mirror V^a.

*Proof.* Choose |t| small and radii max(1/|u|,1/|v|) < r_1 < ⋯ < r_a. Then the CT is the torus integral of
g·S(u)S(v)·K ∏dz_m/(2πi), because Z/∏z_m = 1. Integrate z_1, z_2, …, z_a in turn. When z_m is integrated, z_1..z_{m−1} have
already been set to points of modulus ≤ max(1/|u|,1/|v|), and z_{m+1}.. still run over their circles.
- **Poles inside |z_m| = r_m.**
  1. z_m = 1/u or 1/v, if z_m carries the u- or v-factor (the residue of 1/(uz−1) at 1/u is 1/u);
  2. z_m = t z_l for l < m, from (z_m − z_l)/(z_m − tz_l).
  The poles z_m = z_j/t (j > m) have modulus r_j/|t| > r_m, so they lie outside. g is a polynomial and Z cancels ∏dz/z, so there
  is no pole at 0.
- **Dead residues.**
  - The numerator ∏_{l<m}(z_m − z_l) kills a residue at z_m = tz_l whenever some earlier z_{l′} = tz_l. So each point has at most one
    successor, and the points form t-strings.
  - It also kills z_m = 1/u if the U-string is already seeded. But the seed can only come from the unique u-factor anyway.
  - For m = 1 only type 1 is available, so position 1 seeds a string.
- **Genericity** makes all poles simple and distinct.

So each surviving iterated residue is one of the configurations. Two strings need two seeds, hence the u-factor at the first U
and the v-factor at the first V. With one string the unused factor is just evaluated. ∎

**Lemma 2.2 (string weight).** A string of length A with base q (points q, qt, …, qt^{A−1}, at increasing positions) contributes
the factor (seed residue)·(within-string factors) = **(−1)^{A−1}φ_{A−1}(t)·q^A**.

*Proof.*
- The seed residue is q.
- A successor pair (l → m) contributes (z_m − z_l) = qt^α(t−1).
- A non-successor pair at string distance δ ≥ 2 contributes (t^δ−1)/(t^δ−t) = (1−t^δ)/(t(1−t^{δ−1})).
- So the weight is q·∏_{α=0}^{A−2}qt^α(t−1)·∏_{δ=2}^{A−1}[(1−t^δ)/(t(1−t^{δ−1}))]^{A−δ}.
- The last product telescopes to t^{−C(A−1,2)}φ_{A−1}/(1−t)^{A−1}: the factor (1−t^δ) has net exponent 1 for 2 ≤ δ ≤ A−1, and
  (1−t) has net exponent −(A−2). With the successor factors q^{A−1}t^{C(A−1,2)}(t−1)^{A−1}, everything multiplies out to the claim. ∎

**Corollary 2.3 (Thm 1.5, third proof).** For the one-point function, replace S(u)S(v) by S(u). Only the configuration U^a
survives, and Σ_n⟨T_ag,p_n⟩u^{−n}/(1−t^n) = φ_a^{−1}(−1)^{a−1}φ_{a−1}u^{−a}g(π_a/u). Hence
⟨T_ag,p_n⟩ = (−1)^{a−1}[n]/[a]·g(π_a), which is Day 223 Thm 1.5 (⟨H,p_n⟩ = (−1)^{n−1}lin_eH). ✓

### 2.2 The interaction: one rational function
Fix a U-string of length A and a V-string of length B, with bases 1/u and 1/v, and positions interleaved according to ω. For
positions l < m in different strings, the cross factor is (z_m − z_l)/(z_m − tz_l). In terms of w = u/v:
- U_α before V_β gives (wt^{β−α}−1)/(wt^{β−α}−t);
- V_β before U_α gives (1−wt^{β−α})/(1−wt^{β−α+1}).

Let Sh_{A,B}(w) := Σ_{shuffles ω} ∏_{cross pairs}(cross factor).

**Lemma 2.4 (shuffle identity).** Sh_{A,B}(w) = qbin(A+B, A)·t^{−AB}·(1−w)(1−t^{B−A}w) / ((1−t^{−A}w)(1−t^Bw)).

*Proof.* Induct on A + B. For A = 0 or B = 0 both sides equal 1. Otherwise, condition on the last letter.
- **Last letter U_{A−1}.** It comes after every V_β, with factors ∏_{β=0}^{B−1}(1−wt^{β−A+1})/(1−wt^{β−A+2}) = (1−wt^{1−A})/(1−wt^{B+1−A})
  (telescoping). The rest is a shuffle of strings of lengths A−1 and B, with the same bases.
- **Last letter V_{B−1}.** The factors are ∏_{α=0}^{A−1}(wt^{B−1−α}−1)/(wt^{B−1−α}−t) = t^{−A}(1−wt^{B−1})/(1−wt^{B−1−A}).
- So Sh_{A,B} = (1−wt^{1−A})/(1−wt^{B+1−A})·Sh_{A−1,B} + t^{−A}(1−wt^{B−1})/(1−wt^{B−1−A})·Sh_{A,B−1}.
- **Check the closed form.** Put X = t^A, Y = t^B, s(X,Y) := Sh/(qbin·t^{−AB}). Use
  qbin(A+B−1,A−1)/qbin(A+B,A) = (1−X)/(1−XY) and qbin(A+B−1,A)/qbin(A+B,A) = (1−Y)/(1−XY).
  - The first term becomes Y(1−X)(1−w)/((1−XY)(1−Yw)), because the factor (1−wtY/X)/(1−wt/X) cancels the prefactor.
  - The second term becomes (1−Y)(1−w)/((1−XY)(1−w/X)), by the same cancellation.
  - So the claim reduces to Y(1−X)(1−w/X) + (1−Y)(1−Yw) = (1−XY)(1−wY/X). Both sides expand to 1 − XY − wY/X + wY². ∎

Computed: the closed form equals the brute shuffle sum for A, B ≤ 3 (`shuffle.log`), and sympy confirms the recursion identity in
ℚ(t,w,X,Y) (`shuffle_rec.log`).

**Partial fractions.** (1−w)(1−t^{B−A}w)/((1−t^{−A}w)(1−t^Bw)) = 1 + κ_{AB}Σ_{m≥1}(t^{−Am} − t^{Bm})w^m, with
κ_{AB} = (1−t^A)(1−t^B)/(1−t^{A+B}). The residues at w = t^A and w = t^{−B} are ±κ_{AB}, and the value at w = 0 is 1. Note that
qbin(a,A)φ_{A−1}φ_{B−1} = φ_a/((1−t^A)(1−t^B)), and multiplying that by κ_{AB} gives φ_{a−1}.

### 2.3 The theorem
**Theorem 2.5 (two-point formula).** Let g ∈ Λ_a be homogeneous of degree d, n = a + d, x, y ≥ 1, x + y = n. For 1 ≤ A ≤ a−1 put
B = a − A and

  G_A(w) = Σ_k G_{A,k}w^k := g(1, t, …, t^{A−1}, w, wt, …, wt^{B−1})   (two-string principal specialization).

Then

  Φ_a(g;x,y) = (−1)^a(1−t^x)(1−t^y)·{ Σ_{A=1}^{a−1} t^{−AB}[ G_{A,y−B}/((1−t^A)(1−t^B)) + (1/(1−t^a))·Σ_{m≥1}(t^{−Am} − t^{Bm})·G_{A,y−B−m} ]
                                     − g(π_a)/(1−t^a) · Σ_{j=0}^{a−1} t^{−jy} }.

Equivalently, regrouping the t^{−Am} sum against the last term so that no negative powers of t remain inside the brackets:

  (−1)^{a−1}Φ_a = (1−t^x)(1−t^y)/(1−t^a)·[ g(π_a) + Σ_A t^{−AB}( Σ_{k≥y−B}G_{A,k}t^{A(k−y+B)} + Σ_{k<y−B}G_{A,k}t^{B(y−B−k)} ) ]
                  − (1−t^x)(1−t^y)Σ_A t^{−AB}G_{A,y−B}/((1−t^A)(1−t^B)).

*Proof.* By Lemmas 2.1–2.4, ρ(w) = u^nR is the following sum of rational functions of w:
- **Two strings (A, B).** U-weight (−1)^{A−1}φ_{A−1}u^{−A}, V-weight (−1)^{B−1}φ_{B−1}v^{−B} = (−1)^{B−1}φ_{B−1}u^{−B}w^B, the
  interaction Sh_{A,B}(w), and g(z*) = u^{−d}G_A(w). Total φ_a^{−1}(−1)^aφ_{A−1}φ_{B−1}w^BSh_{A,B}(w)G_A(w).
- **U^a with the v-factor at position j+1.** φ_a^{−1}(−1)^{a−1}φ_{a−1}g(π_a)·1/(vt^j/u − 1) = …·w/(t^j − w).
- **V^a.** w^n·(…)·Σ_j 1/(t^jw − 1), which is O(w^n).

Lemma 2.1 holds for an open set of (u,v) with arbitrary generic ratio, so this is an identity of rational functions of w. Each term
is analytic at w = 0 (its poles are at t^{±k}), so [w^y]ρ can be taken termwise from Taylor expansions at 0. For 1 ≤ y ≤ n−1:
- V^a contributes nothing;
- U^a contributes φ_a^{−1}(−1)^{a−1}φ_{a−1}g(π_a)Σ_j t^{−jy};
- (A, B) contributes, by Lemma 2.4 and the partial fractions,
  φ_a^{−1}(−1)^at^{−AB}[φ_a/((1−t^A)(1−t^B))·G_{A,y−B} + φ_{a−1}Σ_{m≥1}(t^{−Am}−t^{Bm})G_{A,y−B−m}].

Multiplying by (1−t^x)(1−t^y) gives the first form, since φ_{a−1}/φ_a = 1/(1−t^a). For the second form:
- Σ_k G_{A,k}t^{Ak} = G_A(t^A) = g(π_a). So Σ_{m≥1}t^{−Am}G_{A,y−B−m} = t^{−A(y−B)}g(π_a) − Σ_{k≥y−B}G_{A,k}t^{A(k−y+B)}.
- t^{−AB}t^{−A(y−B)} = t^{−Ay}, and Σ_{j=0}^{a−1}t^{−jy} − Σ_{A=1}^{a−1}t^{−Ay} = 1.
- The t^{Bm} sum is already the truncated sum over k < y−B. ∎

**Gate.** At a = 1 the A-sum is empty and Φ_1 = (1−t^x)(1−t^y)/(1−t)·g(1) = (1−t)[x][y]g(1), the known `twopoint.log` value. ✓

**Key example (the wake's failing primitive input).** a = 2, g = p_2, x = y = 2. Then A = B = 1, G_1(w) = 1 + w², y − B = 1,
G_{1,1} = 0 and G_{1,0} = 1. Theorem 2.5 gives
- (1−t²)²{t^{−1}(t^{−1}−t)/(1−t²) − (1+t²)(1+t^{−2})/(1−t²)} = (1−t²)·(t^{−2} − 1 − 2 − t² − t^{−2}) = −(1−t²)(3+t²)
- = (t−1)(t+1)(t²+3),

which is the engine value. The wake hunch had only the leading two-string term. The interaction term (the w-Taylor series of
Sh) and the one-string term carry the rest.

**Verification:**
- Thm 2.5, first form: **229/229** exact agreement with the exact CT evaluator, for every p_ν with |ν| ≤ 5, 4, 4, 3, 2 at a = 1, 2, 3, 4, 5
  respectively, and every split x + y = n (`formula_raw.log`).
- Second form: 173/173, a ≤ 4 (`formula2pt.log`).
- The string expansion itself matches 7/7 complete (a, g) cases (`strings.log`).
- The CT evaluator agrees with the Day 224 engine (`ctcheck.log`, 24/24).

## 3. Bonus: nonsymmetric one-point lemma (PROVED)

**Lemma 3.1.** For any polynomial h(z_1..z_a), homogeneous of degree n − a:

  CT[ h(z) Z z_1^{−n} K ] = (−1)^{a−1}φ_{a−1}·h(1, t, …, t^{a−1}).

*Proof.* Induct on a. The case a = 1 is the CT of a monomial of degree 0. By linearity take h = z^β.
- **Eliminate the top variable.** z_a sits on the largest circle and appears as z_a^{1+β_a}∏_{i<a}(1−z_i/z_a)/(1−tz_i/z_a). Its CT
  is Q_{1+β_a}(z_{<a}), where Σ_eQ_eu^e := ∏_{i<a}(1−uz_i)/(1−tuz_i).
- **Apply induction** to h′ = z_1^{β_1}⋯z_{a−1}^{β_{a−1}}Q_{1+β_a} in a−1 variables. The degrees match.
- **Telescope.** Q_e(1,…,t^{a−2}) = [u^e](1−u)/(1−t^{a−1}u) = −(1−t^{a−1})t^{(a−1)(e−1)} for e ≥ 1. ∎

Computed: exact for a ≤ 4, d ≤ 3, all compositions β (`test1pt.log`). It shows that the principal specialization in Thm 1.5 already
holds before symmetrization: in the subset formula only the smallest variable can carry the p_n(z^{−1}) "source".

## 4. Every ℓ = 3, κ = 1 lead in closed form (class 4 CLOSED)

Let λ = (a,b,c) be any ordering of a length-3 partition, μ = (x,y) two-part, κ(λ,μ) = 1, n = x + y, and m_{xy} = 2 if x = y, else 1.
Put
- Λ_a(r,q) := lin_eT_a(p_rp_q) = (−1)^{r+q}[a+r+q]/[a]·[a]_{t^r}[a]_{t^q} (Thm 1.5);
- U_a(b,c;x,y) := [e_xe_y]T_a(p_bp_c) = ((−1)^nΦ_a(p_bp_c;x,y) − Λ_a(b,c))/m_{xy}, with Φ_a from Theorem 2.5.

**Lemma 4.1.** For G homogeneous of degree n, ⟨G,p_xp_y⟩ = (−1)^n(m_{xy}[e_xe_y]G + lin_eG).
*Proof.*
- ⟨e_ν, p_xp_y⟩ = ⟨Δe_ν, p_x⊗p_y⟩ vanishes for ℓ(ν) ≥ 3: Δe_ν = ∏(Σe_i⊗e_{ν_j−i}), and each tensor side must be a single e to pair
  with p_x or p_y.
- ⟨e_n,p_xp_y⟩ = (−1)^{n−2}.
- ⟨e_xe_y,p_xp_y⟩ = (−1)^{x−1}(−1)^{y−1}m_{xy}.
- Other two-part ν give 0. ∎

**Theorem 4.2 (closed ℓ = 3, κ = 1 lead).**

  [(s−1)²]c_{λμ} = (−1)^{b+c}U_a(b,c;x,y) + Σ_{1≤r<b, {b−r, a+r+c}={x,y}}(−1)^{r+c}Λ_a(r,c) + Σ_{1≤q<c, {c−q, a+b+q}={x,y}}(−1)^{b+q}Λ_a(b,q)
                   + [e_xe_y] D_a(M_{bc}),

where D_a is the derivation with D_a(e_j) = M_{aj} (Day 223 Thm 7.1). So [e_xe_y]D_a(M_{bc}) is the finite sum of products of
the L-weights read off from M_{bc} and M_{a·}.

*Proof.* Use Day 224 Prop 5.1: [(s−1)²]c = [e_μ](Γ_a(e_b,e_c) + D_aD_b(e_c)), and D_b(e_c) = M_{bc}. Expand
Γ_a(e_b,e_c) = Σ_{r≤b,q≤c}(−1)^{r+q}e_{b−r}e_{c−q}T_a(p_rp_q) (Day 224 §1).
- T_a(·) is homogeneous of degree ≥ a ≥ 1, so e_{b−r}e_{c−q}T_a(…) has ≥ 3 e-factors when r < b and q < c, and contributes nothing.
- For r < b, q = c it contributes e_{b−r}·(single-e part of T_a(p_rp_c)) = Λ_a(r,c)e_{b−r}e_{a+r+c}. The case r = b, q < c is symmetric.
- For r = b, q = c it contributes U_a.
- Lemma 4.1 converts Φ_a into U_a. ∎

**Corollary 4.3.** Every connected v = 2 lead with ℓ(λ) = 3 is now closed.
- μ with three padded parts reduces by the Column Lemma (Day 224 Cor 2.2) to μ_3 = 0.
- One-part μ is Thm G.
- Two-part μ is Thm 4.2.

With block multiplicativity (Day 224 Thm 6.1), **every v = 2 lead (κ = ℓ − 2, any ℓ) has a closed formula.** Classes 1–3 of
Day 224 §8 become special cases, and class 4 is no longer open.

**Verification (exact, all computed from the closed formulas only: Thm 2.5 + Thm 1.5 + Thm 7.1):**
- **KILL TEST:** (3,3,3)→(7,2) = 2t¹³+3t¹²+3t¹¹+6t¹⁰+6t⁹+6t⁸+9t⁷+6t⁶+3t⁵+9t⁴+6t³+3t+4, identical to the engine (`class4_killtest.log`).
- All **27** κ = 1 two-part pairs with n ≤ 10, in **each of the 3 orderings** (a = λ_1, λ_2, λ_3): **81/81** symbolic matches with the
  engine leads of `day224/newcases_n10.log` (`class4_all_n10.log`). Order independence is one more consistency check against Hikita
  commutativity.
- All **16** open class-4 pairs with n ≤ 12 (list of `day224/classify.py`): **48/48** exact matches against the engine `second()` at
  t ∈ {2, 3, −1/2} (`class4_n12.log`).

## 5. What died, and why the wake hunch failed
- **The X ⊔ Y Hopf split** is not the right decomposition. Its pieces are rational in u and depend on the expansion region, which is
  exactly what the wake saw. The t-string decomposition is the right one: the pieces are again rational in w, but the TOTAL is the
  sum of their Taylor coefficients at w = 0, with no region ambiguity.
- **The wake hunch** "out_i = t^{(a−i)(x−i)}Σ⟨T_if₁,p_x⟩⟨T_{a−i}f₂,p_y⟩" is roughly the G_{A,y−B}/((1−t^A)(1−t^B)) term. It misses the
  interaction series κ_{AB}Σ_m(t^{−Am}−t^{Bm})w^m and the one-string term. The wake saw the hunch fit on its e-type inputs. On
  primitive inputs the missing terms are visibly nonzero (the §2.3 example).

## 6. Gaps and grades
- **Thm 1.1:** proved, modulo textbook Macdonald III facts. The equation numbers are quoted from memory: (1.4), (2.1)–(2.2), (2.12),
  (2.15), (4.4), (4.9). They should be checked against the book before external use. Not deep-read this session.
- **Lemmas 2.1, 2.2, 2.4, Thm 2.5, Cor 2.3, Lemma 3.1:** proved, self-contained residue calculus plus algebra. Every identity was also
  computer-checked (above).
- **Thm 4.2, Cor 4.3:** proved, from Prop 5.1 (Day 224), Thm 7.1, Thm 1.5 (Day 223), Thm 2.5 and Lemma 4.1. Lead identification at
  v = 2 uses Thm C as before (Day 220; (N)-dependent via t = 0). This is unchanged from Day 224.
- **No cold recheck yet.** These are proved, not checked-sober. The dream cycle should re-derive Lemma 2.4 and the coefficient
  extraction cold.
- **Novelty flags (NOT checked; no browsing):**
  - Φ_a(P_ρ;x,y) = (1−t^x)(1−t^y)X^{ρ+1^a}_{(x,y)}(t)/b_{ρ+1^a}-type data, i.e. Thm 2.5 is a closed formula for **two-part Green
    polynomials** X^λ_{(x,y)}(t) at ℓ(λ) = a, for arbitrary λ. Green polynomials for two-row/two-part μ have a long history: Green,
    Morris 1963, Kirillov's "Ramifications of Hall–Littlewood", and Garsia–Procesi's charge formulas. This needs a first-hand
    search before ANY claim.
  - The t-string residue evaluation of HL constant terms is classical in spirit, e.g. Macdonald's proof of the norm formula and
    Cherednik-style constant terms.
  - The shuffle identity Lemma 2.4 is probably a known q-identity (a two-block parabolic symmetrizer at string points).
  - What is ours is the application to ⋆ leads (Thm 4.2 / Cor 4.3).

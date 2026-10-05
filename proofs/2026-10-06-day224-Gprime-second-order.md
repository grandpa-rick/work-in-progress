# Day 224 PROVE: G′ at second order — Box Complement, block multiplicativity, and closed forms for (n−1,1) and λ ∋ 1

**Date:** 2026-10-05 session (file dated per PROVE.md, 2026-10-06). **Author:** Rick. No sub-agents. All scripts are in
`proofs/scripts/day224/` (engine `hl.py`, `ops.py`, `lead2.py`, `second_gen.py`; closed formulas `closed.py`, `one_formula.py`).

> Drunk summary. PROVE.md asked for an order-2 Taylor formula and then a closed form for three n = 6 leads. The order-2
> formula is a one-liner (§1). Then the data showed something better. (2,2,2)→(4,1,1) has the SAME lead as
> (1,1,1)→(3), and (3,3,2)→(4,4) the same as (2,1,1)→(4). Stripping a column does nothing, and neither does complementing
> in a box. **c_{N^ℓ−λ, N^ℓ−μ} = s^{N·C(ℓ,2) − (ℓ−1)|λ|} c_{λμ}, EXACTLY, for all s and t.** It's x ↦ 1/x in N variables:
> E_{N−k} is E_k turned inside out. That kills two of the three targets on sight, since (2,2,2)→(3,3) and →(4,1,1) are the
> coarsening (1,1,1)→(3) in a mirror, so Thm G gives (t+2)[3]. For the third, (2,2,2)→(5,1), peel x_1 one order below the
> top (§4). The target becomes a *linear* e-coefficient, and Thm 1.5 already knows every linear coefficient in closed
> form for all s: **lin_e(e_k⋆G) = (−1)^d[k+d]/[k]·G[(s−1)[k]_t]**. Out comes a closed formula for c_{λ,(n−1,1)}(s,t),
> and as a bonus a ten-line re-proof of the 207b Pieri rule. Then the polarization bookkeeping of Day 220 Prop 2.3 shows
> that EVERY lead is a sum over block decompositions of products of connected leads (§6, generalizing Thm F). So G′
> reduces to connected leads, and at v = 2 to ℓ = 3. What is left open is an explicit family starting at
> (3,3,3)→(7,2). One more thing died: leads are NOT positive ((4,4,2)→(7,3) has a −t), so "G′ = a t-count of flows" is dead.

---

## 0. Setup
Notation as in Day 220 §0 and Day 223 §0:
- E_kF = Σ_{|A|=k} c_A X_A F(X_{A^c}, sX_A), c_A = ∏_{i∈A,j∉A}(x_i−tx_j)/(x_i−x_j).
- e^⋆_λ = E_{λ_1}⋯E_{λ_ℓ}(1) = Σ_μ c_{λμ}(s,t)e_μ, with E_0 = id. Zero parts are allowed throughout ("padded" λ).
- I = (e_1, e_2, …). lin_e G = [e_n]G for G homogeneous of degree n.
- κ(λ,μ) and v = (s−1)-valuation as in Day 220 §3.
- E_k = Σ_p (s−1)^pE_k^{(p)}, D_k = E_k^{(1)}, and B(e_k,e_r) = M_{kr} (Day 223 Thm 7.1).
- T_k g = Σ_A c_AX_A g(X_A) (Day 223 §0).
- L(a,b) = (1−t^{a+b})(t^{ab}−1)/((1−t^a)(1−t^b)), and [m] = [m]_t.

## 1. (A) The all-order Taylor lemma, and the bilinear piece (PROVED)
**Lemma 1.1.** Put y_{ia} := u_ix_a/(1+u_ix_a). Then for every p,

  E_k^{(p)}(∏_iE(u_i)) = ∏_iE(u_i) · T_k( e_p({y_{ia}}_{i, a∈A}) ).

*Proof.* (1+su x)/(1+ux) = 1 + (s−1)y. Hence ∏_{i,a∈A}R = ∏(1+(s−1)y_{ia}) = Σ_p(s−1)^pe_p(y). Insert this in Day 223 Lemma
1.1 and compare (s−1)-coefficients. T_k contains no s. ∎

**Order 2.** e_2(Y_u ⊔ Y_v) = e_2(Y_u) + e_2(Y_v) + p_1(Y_u)p_1(Y_v). So the polarization (Day 220 Lemma 2.1) of the order-2
operator on two factors is

  Γ_a(e_b,e_c) := E_a^{(2)}(e_be_c) − e_bE_a^{(2)}(e_c) − e_cE_a^{(2)}(e_b) = [u^bv^c] E(u)E(v)·T_a(p_1(Y_u)p_1(Y_v))
             = Σ_{r≤b, q≤c} (−1)^{r+q} e_{b−r}e_{c−q} T_a(p_rp_q),

using p_1(Y_u) = Σ_r(−1)^{r−1}u^rp_r(X_A). The diagonal x = x′ part of p_1(Y_u)p_1(Y_v) is
Σ_x uvx²/((1+ux)(1+vx)) = (u·p_1(Y_v) − v·p_1(Y_u))/(u−v). By D_a(E(u)) = E(u)T_a(p_1(Y_u)) (Lemma 1.1, p = 1),

  Γ_a^{diag}(E(u),E(v)) = (u E(u)·D_aE(v) − v E(v)·D_aE(u))/(u − v),

which is explicit by Thm 7.1. **For a = 1 there is no off-diagonal part**, since A is a single variable. Hence Γ_1 = Γ_1^{diag}, and
extracting [u^bv^c] (expand (u^{i+1}v^j − u^jv^{i+1})/(u−v)) gives

  (1.2) Γ_1(e_b,e_c) = Σ_{j=1}^{min(b,c)} e_{b+c−j}D_1(e_j) − Σ_{j=max(b,c)+1}^{b+c} e_{b+c−j}D_1(e_j),   D_1(e_j) = e_1e_j − [j+1]e_{j+1}.

Check: Γ_1(e_1,e_1) = e_1³ − (t+2)e_2e_1 + [3]e_3, which agrees with the engine (`gamfull_n6.log`).
**Gap (A):** the off-diagonal part of Γ_a for a ≥ 2 has no closed form here. It is the only missing piece of an explicit
E_a^{(2)}(e_be_c). E_a^{(2)}(e_b) itself is explicit as (s−1)²-coefficient of the Pieri rule (Cor. 3.2).

## 2. Box Complement Theorem (PROVED; 243/243 full (s,t) checks, n ≤ 6)

**Theorem 2.1.** Let ℓ ≥ 1, N ≥ 0, and let λ, μ be partitions with at most ℓ parts (padded with zeros to length ℓ) and all
parts ≤ N. Write N^ℓ − λ := (N−λ_ℓ, …, N−λ_1). Then

  c_{N^ℓ−λ, N^ℓ−μ}(s,t) = s^{N·C(ℓ,2) − (ℓ−1)|λ|} · c_{λμ}(s,t).

*Proof.*
1. **Stability for every N (incl. N < degree).** For F ∈ Λ_{N+1}: (E_kF)|_{x_{N+1}=0} = E_k(F|_{x_{N+1}=0}). Terms with N+1 ∈ A
   carry the factor x_{N+1} in X_A and are regular at x_{N+1} = 0, so they vanish. For N+1 ∉ A, c_A = c_A^{(N)}·∏_{i∈A}(x_i−tx_{N+1})/(x_i−x_{N+1})
   → c_A^{(N)}. Iterating, the N-variable subset-formula value of e^⋆_λ is the restriction of the stable one, namely
   Σ_μ c_{λμ}e_μ(x_1..x_N) (e_μ = 0 if μ_1 > N).
2. **Inversion lemma.** In N variables, let G be homogeneous of degree g and F(x) := e_N(x)^m G(1/x) (a rational function). Then

     E_{N−k}F = s^{m(N−k)−g} · e_N^{m+1} · (E_kG)(1/x).

   Put y = 1/x, sum over B = A^c with |A| = k, and check four things:
   - c_B(x) = ∏_{i∈B,j∈A}(y_j − t y_i)/(y_j − y_i) = c_A(y);
   - X_B = e_N(x)·Y_A;
   - F(x_A, s x_B) = s^{m(N−k)}e_N(x)^m G(y_A, s^{−1}y_B);
   - G(y_A, s^{−1}y_B) = s^{−g}G(y_B, s y_A) = s^{−g}G(Y_{A^c}, sY_A).
3. **Iterate.** Let G_j := E_{λ_{ℓ−j+1}}⋯E_{λ_ℓ}(1) (degree D_j) and F_j := E_{N−λ_{ℓ−j+1}}⋯E_{N−λ_ℓ}(1). By induction with step 2,
   F_j = s^{σ_j}e_N^jG_j(1/x) with σ_{j+1} = σ_j + j(N−λ_{ℓ−j}) − D_j. Summing the increments gives
   σ_ℓ = N·C(ℓ,2) − Σ_i(i−1)p_i − Σ_i(ℓ−i)p_i = N·C(ℓ,2) − (ℓ−1)|λ| (p_i := λ_{ℓ+1−i}), independent of the order.
4. **Compare coefficients.** e_{N−r}(x) = e_N(x)e_r(1/x) in N variables (also for r = 0). So e_N^ℓe_μ(1/x) = e_{N^ℓ−μ}(x) for ℓ(μ) ≤ ℓ.
   Every μ in the support has ℓ(μ) ≤ ℓ(λ) ≤ ℓ (DS, Day 214). Hence
   e^⋆_{N^ℓ−λ} = s^{σ}Σ_{μ_1≤N} c_{λμ}e_{N^ℓ−μ}(x). The map μ ↦ N^ℓ−μ is a bijection onto partitions with ≤ ℓ parts, all ≤ N, and
   {e_ν : ν_1 ≤ N} is a basis of Λ_N. ∎

**Corollary 2.2 (Column Lemma).** c_{λ+1^ℓ, μ+1^ℓ} = s^{C(ℓ,2)}c_{λμ} for λ, μ with ≤ ℓ parts (padded).
*Proof.* Complement with N, then with N+1: σ_N(λ) + σ_{N+1}(N^ℓ−λ) = (2N+1)C(ℓ,2) − (ℓ−1)ℓN = C(ℓ,2). Equivalently
(direct): [x_1^{m+1}]E_{k+1}F = s^mE_k([x_1^m]F) (Thm 4.1(a)), iterated ℓ times from F = 1. ∎

**Corollary 2.3 (lead invariance).** v and Lead are invariant under (λ,μ) ↦ (λ+1^ℓ, μ+1^ℓ) and (λ,μ) ↦ (N^ℓ−λ, N^ℓ−μ),
since s^σ is a unit of ℚ[[s−1]] with constant term 1. In particular, if μ padded to ℓ = ℓ(λ) has μ_2 = … = μ_ℓ, strip μ_ℓ
columns (μ ⊵ λ forces λ_ℓ ≥ μ_ℓ). If it has μ_1 = … = μ_{ℓ−1}, complement with N = μ_1. **Either way (λ,μ) is a coarsening in disguise,
and Thm G gives its lead in closed form.**

## 3. The plethystic linear coefficient, and a new proof of the ℓ = 2 Pieri rule (PROVED)

**Theorem 3.1.** For k ≥ 1 and G ∈ Λ ⊗ ℚ(t)[s] homogeneous of degree d:

  lin_e(e_k⋆G) = lin_e E_k(G) = (−1)^d ([k+d]_t/[k]_t) · G[(s−1)[k]_t].

Here G[(s−1)[k]_t] is the ring map e_j ↦ π_k(j) := [u^j]∏_{m=0}^{k−1}(1+sut^m)/(1+ut^m), equivalently p_r ↦ (s^r−1)[k]_{t^r}.
*Proof.* Take G = e_J by linearity. Day 223 Lemma 1.1 gives E_k(e_J) = Σ_{w≤J} e_{J−w}·T_k([u^w]∏_iR(u_i)).
- For w ≠ J, e_{J−w} ∈ I, and T_k(·) is homogeneous of degree ≥ k ≥ 1, so it is in I. The product is in I² and lin_e kills it.
- For w = J, Day 223 Thm 1.5 gives (−1)^d[k+d]/[k]·∏_i[u_i^{j_i}]R(u_i)|_{X_A=(1,…,t^{k−1})} = (−1)^d[k+d]/[k]∏_iπ_k(j_i). ∎

**Corollary 3.2 (Pieri, all s).** For k, r ≥ 0:

  e_k⋆e_r = Σ_{y=0}^{min(k,r)} s^y C_{k−y,r−y} e_{k+r−y}e_y,   C_{a,0} = C_{0,b} = 1,  C_{a,b} = (−1)^b[a+b]/[a]·π_a(b).

*Proof.*
- **Support.** By DS the support is μ ⊵ (k,r), i.e. μ = (k+r−y, y) with y ≤ min(k,r).
- **Coefficient.** Apply Cor 2.2 with ℓ = 2, y times. This works for compositions too: the x_1-top proof never uses sorting.
  It reduces the coefficient to lin_e E_{k−y}(e_{r−y}), and Theorem 3.1 evaluates that. ∎

This is an independent proof of the Day 207b Pieri rule. It uses neither the A_k key lemma nor the outer-peel induction; its inputs
are Thm 1.5 and the Column Lemma. **62/62** coefficients agree with the full (s,t) engine (`closed_check.py`).
Comparing with 207b yields the identity F_n(t^m) = (−1)^m[n+m]/[n]·π_n(m) for n, m ≥ 1, since both are proved formulas for the same coefficient.

## 4. The subleading x_1-coefficient, and a closed formula for c_{λ,(n−1,1)} (PROVED)

**Theorem 4.1.** Let F ∈ Λ_N (coefficients in ℚ(t)[s]) with deg_{x_1}F ≤ m, F = Σ_{j≤m}x_1^jF_j(x′), x′ = (x_2..x_N). Then
deg_{x_1}E_{k+1}F ≤ m+1, and

  (a) [x_1^{m+1}]E_{k+1}F = s^m E_k(F_m),
  (b) [x_1^{m}]E_{k+1}F = s^{m−1}E_k(F_{m−1}) + t^{k+1}E_{k+1}(F_m) + s^m(1−t)K_k(F_m),

with all operators acting in x′ and

  K_k(G) := Σ_{|A′|=k} c_{A′}X_{A′}p_1(X′_{A′^c})G(X′_{A′^c}, sX′_{A′}) = (s·e_1E_k(G) − E_k(e_1G))/(s−1)   (K_0 = e_1·).

*Proof.* E_{k+1}F is a polynomial, so its x_1-coefficients are those of its Laurent expansion at x_1 = ∞. Expand termwise.
- **1 ∈ A**, A = {1} ⊔ A′. Then c_A = R_1·c_{A′}(x′) with R_1 = ∏_{j∈x′∖A′}(x_1−tx_j)/(x_1−x_j) = 1 + (1−t)p_1(x′∖A′)/x_1 + O(x_1^{−2}).
  Also X_A = x_1X_{A′} and F(X_{A^c}, sX_A) = Σ_js^jx_1^jF_j(X′_{A′^c}, sX′_{A′}).
  - The coefficient of x_1^{m+1} is s^m F_m(…).
  - The coefficient of x_1^m is s^{m−1}F_{m−1}(…) + s^m(1−t)p_1(x′∖A′)F_m(…).
- **1 ∉ A.** The factors ∏_{i∈A}(x_i−tx_1)/(x_i−x_1) = t^{k+1} + O(1/x_1) have degree 0 in x_1, and F has x_1-degree ≤ m. So these terms
  contribute only t^{k+1}E_{k+1}(F_m) to the coefficient of x_1^m, and nothing to x_1^{m+1}.
- **Summing over A′.** This gives (a) and (b). For the K_k identity: E_k(e_1G) = Σ c X(p_1(A′^c) + s p_1(A′))G(…) and
  e_1E_k(G) = Σ c X(p_1(A′^c) + p_1(A′))G(…). ∎

**Theorem 4.2 (closed formula for μ = (n−1,1)).** Let λ = (λ_1,λ_2,λ_3) with all λ_i ≥ 1 (any order), n = |λ| ≥ 3,
ν = λ − 1³ and k = ν_1. Put
- Λ_k[G] := lin_eE_k(G) (Thm 3.1 for k ≥ 1; Λ_0[G] = lin_eG);
- F := e_{λ_2}⋆e_{λ_3} (Cor 3.2), F_1 := [x_1]F (via e_j = e′_j + x_1e′_{j−1});
- F_2 := s·e_{ν_2}⋆e_{ν_3}.
Then

  c_{λ,(n−1,1)} = s·Λ_k[F_1] + t^{λ_1}Λ_{λ_1}[F_2] + s²(1−t)·Λ^+_k[F_2] − ε_n s³·Λ_k[e_{ν_2}⋆e_{ν_3}],

where
- Λ^+_k[e_J] := lin_eK_k(e_J) = (−1)^{|J|}[k+|J|+1]·∏π_k(j_i) for k ≥ 1 (it is −Λ_k[e_1e_J]/(s−1), and Thm 3.1 gives a factor (s−1)[k] from π_k(1));
- Λ^+_0[G] = (constant term of G);
- ε_n = 1 for n > 3, and ε_3 = 3.

*Proof.*
1. **Apply Thm 4.1(b)** to e^⋆_λ = E_{λ_1}F with m = 2. By 4.1(a), F_2 = [x_1²]F = s·e^⋆_{(ν_2,ν_3)}. Then apply lin_e in x′.
2. **Left side.** [x_1²]e^⋆_λ = Σ_μ c_{λμ}[x_1²]e_μ. Here [x_1²]e_μ is
   - e_{μ_1−1}e_{μ_2−1} if ℓ(μ) = 2;
   - Σ_ie_{μ_i}∏_{j≠i}e_{μ_j−1} if ℓ(μ) = 3;
   - 0 if ℓ(μ) = 1.
3. **Which terms survive lin_e.** It keeps exactly μ = (n−1,1), coefficient 1, and μ = (n−2,1,1), coefficient ε_n.
4. **The (n−2,1,1) term.** c_{λ,(n−2,1,1)} = s³c_{ν,(n−3)} (Cor 2.2), and c_{ν,(n−3)} = Λ_k[e_{ν_2}⋆e_{ν_3}]. ∎

**Verification:**
- **20/20** full (s,t) agreement with the engine data, for every ordering of every ℓ = 3 λ with n ≤ 6 (`closed_check.py`). The order
  independence is a consistency check on Hikita commutativity.
- At n ≤ 10, for all 11 pairs with μ = (n−1,1) and κ = 1: the (s−1)^0 and (s−1)^1 coefficients vanish symbolically, and the (s−1)²
  coefficient equals the engine's lead (**11/11**, `closed_leads.py`).

## 5. ℓ = 3, κ = 1: reduction, and the closed form when λ ∋ 1 (PROVED)

**Proposition 5.1.** Let (a,b,c) be any ordering of λ and κ(λ,μ) = 1. Then

  [(s−1)²]c_{λμ} = [e_μ](Γ_a(e_b,e_c) + D_aD_b(e_c)).

*Proof.* [(s−1)²]E_aE_bE_c(1) = E_a^{(2)}(e_be_c) + D_aD_b(e_c) + e_aE_b^{(2)}(e_c), because E^{(p)}(1) = 0 for p ≥ 1. Expand
E_a^{(2)}(e_be_c) = Γ_a + e_bE_a^{(2)}e_c + e_cE_a^{(2)}e_b. Each of e_bE_a^{(2)}e_c, e_cE_a^{(2)}e_b, e_aE_b^{(2)}e_c is a product of
some e_{λ_i} with a function of type (other two parts) (Day 220 Prop 2.3). So it lies in span{e_{(λ_i)∪ν}: ν ⊵ …}. That gives κ ≥ 2 for
every μ in its support. ∎ (Engine: these three pieces vanish at every κ = 1 μ tested, `split3_n7.log`.)

**Theorem 5.2 (λ ∋ 1).** Let λ = (b,c,1), n′ = b+c, m = min(b,c), M = max(b,c), μ = (x,y) with x > y, x+y = n′+1, κ = 1. Then

  [(s−1)²]c_{λμ} = [x]·(𝟙_{x≥M+2} − 𝟙_{x≤m+1} − L(b−y,c−y)𝟙_{y<m}) − [y]·(𝟙_{2≤y≤m+1} + L(b−y+1,c−y+1)𝟙_{2≤y≤m}).

*Proof.* Use Prop 5.1 with a = 1. Then:
- D_1D_b(e_c) = D_1(M_{bc}) with M_{bc} = m·e_be_c + Σ_{i=0}^{m−1}L(b−i,c−i)e_ie_{n′−i} (Thm 7.1, reindexed).
- D_1 is a derivation with D_1(e_j) = e_1e_j − [j+1]e_{j+1}.
- Collect the two-factor terms:
  - those of D_1(e_ie_{n′−i}), i ≥ 1, are −[n′−i+1]e_ie_{n′−i+1} − [i+1]e_{n′−i}e_{i+1};
  - those of D_1(e_{n′}) and D_1(e_be_c) are e_1e_{n′}, e_be_{c+1} and e_ce_{b+1}, all of κ ≥ 2, so they drop.
- From (1.2), the two-factor terms of Γ_1 are
  −Σ_{j=1}^{m}[j+1]e_{n′−j}e_{j+1} + Σ_{j=M+1}^{n′−1}[j+1]e_{n′−j}e_{j+1}, plus −e_1e_{n′} (κ ≥ 2).
- Since x + y = n′+1: {i, n′+1−i} = {x,y} ⇔ i ∈ {x,y}, and {n′−j, j+1} = {x,y} ⇔ j+1 ∈ {x,y}.
- Bounding which of x, y lies in each range gives the displayed form. The sum form and the displayed form agree on 55/55 cases,
  n ≤ 15 (`one_formula.py`). ∎

**Verification:** it agrees with the engine leads on **10/10** pairs (8 distinct), n ≤ 10.

## 6. Block multiplicativity: every lead is a sum of products of connected leads (PROVED; generalizes Thm F)

**Theorem 6.1.** Let κ = κ(λ,μ) ≥ 1, ℓ = ℓ(λ), m = ℓ−κ. Then

  [(s−1)^m]c_{λμ} = Σ_{π} Σ_{(ν^C)_{C∈π}} ∏_{C∈π} [(s−1)^{|C|−1}]c_{λ_C, ν^C},

where
- π runs over set partitions of the positions [ℓ] into exactly κ blocks;
- (ν^C) runs over tuples of partitions with ⊔_Cν^C = μ (as multisets) and ν^C ⊵ λ_C;
- every factor then has κ(λ_C,ν^C) = 1.

All lower (s−1)-coefficients vanish (Thm A). With Thm C (v = ℓ − κ; (N)-dependent via t = 0) this is the Lead.

*Proof.*
1. **Histories.** Expand [(s−1)^m]e^⋆_λ = Σ_{Σp_i=m}E^{(p_1)}_{λ_1}⋯E^{(p_ℓ)}_{λ_ℓ}(1) by Day 220 Lemma 2.1 at each step, as in Prop 2.3.
   The result is a sum over histories H, where step i chooses S_i among the current blocks with |S_i| ≤ p_i. Each history gives
   ∏_{C∈π(H)}g^H_C with g^H_C of type λ_C.
2. **Tightness.** |π(H)| = ℓ − Σ|S_i| ≥ ℓ − m. Equality holds iff |S_i| = p_i for all i (tight).
3. **Non-tight terms drop.** A non-tight term has > κ blocks of types λ_C. Its e_μ-coefficient is nonzero only if μ splits
   compatibly into > κ blocks, contradicting maximality of κ.
4. **Tight histories factor.** In a tight history, step i merges only blocks inside the final block containing i, and
   δ_{S_i} depends only on E^{(p_i)}_{λ_i} and the g's in S_i. By induction over the steps, g^H_C equals the function produced by the
   restricted history on C. That restricted history is a single-block tight history of λ_C, since Σ_{i∈C}|S_i| = |C|−1.
   Conversely, single-block tight histories on the blocks of π interleave uniquely (by global position order) into a global tight
   history. Hence the tight part is Σ_{|π|=κ}∏_C Conn(λ_C), where Conn(λ_C) is the sum over single-block tight histories of λ_C.
5. **Identify the factors.** Apply steps 1–3 to λ_C itself with m = |C|−1. There, tight means single block. So
   [e_ν]Conn(λ_C) = [e_ν][(s−1)^{|C|−1}]e^⋆_{λ_C} whenever κ(λ_C,ν) = 1. Expanding [e_μ]∏_CConn(λ_C) = Σ_{(ν^C)}∏[e_{ν^C}]Conn(λ_C):
   - factors with ν^C ⋭ λ_C vanish (type);
   - if κ(λ_C,ν^C) ≥ 2, refining C gives > κ blocks, which contradicts maximality.

   **Convention:** c_{λ_C,ν} here means the coefficient in e^⋆ of the *composition* λ_C (positions of C in their original
   order), so the proof uses no commutativity. By Hikita's commutativity of ⋆ it equals the coefficient for the sorted partition.
   Hikita 2503.23597 is at verified-quote extraction only (memory/reading/sources.json); the order-independence is also confirmed by the engine
   (all orderings agree in `run_n6.py` and `closed_check.py`). ∎

**Verification:** **99/99** (t = 3, n ≤ 8) and **205/205** each at t = 2 and t = 5 (n ≤ 9), over all ℓ ≥ 4 pairs with κ = ℓ−2
(`decomp_test.py`). The general-ℓ (s−1)² engine was cross-checked **52/52** against the full (s,t) data for ℓ ≥ 4, n ≤ 6 (`gen_vs_cst.py`).

**Consequences.**
- Theorem F (coarsenings) is the special case in which every ν^C is a single part.
- **G′ reduces to connected (κ = 1) leads.** At v = 2 (κ = ℓ−2):
  Lead = Σ_{two disjoint pairs}L·L-type products (pair leads are the v = 1 weights of Thm 7.1/Cor 7.2) + Σ_{triples}(ℓ = 3, κ = 1 lead).

## 7. The three PROVE.md targets — all PROVED
- **(2,2,2)→(4,1,1):** Cor 2.2 strips 1³ to give (1,1,1)→(3). Thm G: pref·K = [3](t+2). ✓
- **(2,2,2)→(3,3):** Thm 2.1 with N = 3 gives (1,1,1)→(0,0,3), and σ = 9 − 12 = −3 (a unit). Lead = [3](t+2). ✓
  The PROVE.md hunch "(2,2,2)→(3,3) *is* a full merge of three unit moves" is literally true: it is a mirror of one.
- **(2,2,2)→(5,1):** Thm 4.2, evaluated exactly in ℚ(s,t) (`closed_leads.py`). The (s−1)^0 and (s−1)^1 coefficients vanish, and the
  (s−1)² coefficient is 2t⁷+3t⁶+6t⁵+6t⁴+9t³+6t²+3t+4, the engine value. ✓ (Its mirror (3,3,3)→(5,4) has the same lead,
  `newcases_n10.log`.)

## 8. What is covered at v = 2, ℓ = 3, and what is open
Pairs with ℓ(λ) = 3 and κ = 1, with μ padded to (μ_1,μ_2,μ_3) and canonical μ_3 = 0 (Cor 2.2), fall into these classes:
1. μ has a repeated padded part: coarsening in disguise (Cor 2.3), so Thm G applies.
2. μ = (x,y,0) with y = 1 or x−y = 1: the (n−1,1) family up to mirror. Thm 4.2 applies.
3. λ_3 = μ_3+1 or λ_1 = μ_1−1: the orbit contains a λ with a part 1. Thm 5.2 applies.
4. **OPEN:** everything else. Among two-part μ, n ≤ 12: 25 pairs fall in class 2, 30 in class 3 and **16 are open** (`classify.py`).
   The smallest open pairs are (3,3,3)→(7,2), then (4,4,2)→(7,3) and (4,3,3)→(8,2).
   For these, Thm 4.1(b) still gives an exact *recursion* in n: it reduces to ℓ ≤ 4 coefficients at two-part μ of degree n−2,
   but not a closed form. The missing ingredient is the off-diagonal part of Γ_a, a ≥ 2. Equivalently, it is the 2-point
   functional f ↦ ⟨T_af, p_xp_y⟩. That functional is (1−t)[x][y]f(1) at a = 1 (computed, `twopoint.log`) and has no product
   form at a ≥ 2.

**Escalation (three strikes):** see `memory/for-collaborator/2026-10-05-day224-escalation-class4.md`. The three attempts were:
HL Green polynomials (the two-part Green polynomials are sums, not products); commutativity relations (these fix Γ only up to its
symmetric part, which is circular); and iterated x_1-extraction (the operator family does not close). Lead-level mirror consistency
holds on 10/10 data pairs, n ≤ 10.

**Positivity is dead.** (4,4,2)→(7,3) = (t+1)(t²+1)(2t⁸+t⁷+t⁶+3t⁴+t³+t²−t+2) and (4,3,3)→(8,2) also have a negative coefficient. So
no "t-count of minimal flow forests" (the questions/ file's G′ guess) can hold in general. Lead(0) ∈ {2, 4} on all 27 two-part data points.
Lead(1) = 2n²−6n+3 for every μ = (n−1,1) (n = 6..10) is observed only, not proved.

## 9. Grades
- Lemma 1.1, (1.2) Γ_1: **proved**.
- Theorem 2.1 (Box Complement), Cor 2.2, Cor 2.3: **proved** ((N)-free; inputs are the subset formula and DS).
- Theorem 3.1, Cor 3.2 (Pieri re-proof): **proved** (inputs: Day 223 Thm 1.5 + Lemma 1.1, Cor 2.2, DS).
- Theorem 4.1, Theorem 4.2: **proved**; 4.2 is also computed 20/20 full and 11/11 leads.
- Prop 5.1, Theorem 5.2: **proved** (inputs: Thm 7.1, Prop 2.3).
- Theorem 6.1: **proved** (inputs: Day 220 Lemma 2.1, Prop 2.3, Thm A, Hikita commutativity). Lead identification needs Thm C.
- G′ at v = 2: **reduced** to connected ℓ = 3 leads. Those are closed in classes 1–3, and **open** in class 4.
- **Novelty flags (not checked; no browsing this session).**
  - Thm 2.1 is presumably the shadow, under (N), of Macdonald's complementation P_ν(x^{−1})e_N^ℓ = P_{ν^c}(x) (Macdonald VI). The
    direct (N)-free proof is new to me. Before claiming anything, check Hikita 2503.23597 and DFK for an inversion symmetry of E_k.
  - Theorem 6.1 is a cluster/multiplicativity statement of the kind familiar from Mayer expansions; the ⋆-specific content is
    the connected-lead input.
  - Cor 3.2 re-derives 207b; the identity F_n(t^m) = (−1)^m[n+m]/[n]π_n(m) is a ₂φ₁ evaluation, probably classical.

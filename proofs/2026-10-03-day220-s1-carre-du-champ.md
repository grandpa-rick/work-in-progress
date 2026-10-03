# Day 220 PROVE: the s = 1 edge of ⋆ — biderivation, the block valuation law, and exact valuations on coarsenings

**Date:** 2026-10-03. **Author:** Rick (the reasoning thread is mine; sub-agents were not used; scripts listed in §6).
**Scope.** e^⋆_λ := E_{λ_1}E_{λ_2}⋯E_{λ_ℓ}(1) = Σ_μ c_{λμ}(s,t) e_μ, with Hikita's e_k⋆ = E_k given by the subset formula
(Day 207b (A_k)+(K_k), Day 214 (0.1); registry `subset-formula-Ak-Kk`, proved):

  E_k F = Σ_{|A|=k} c_A · X_A · F(X_{A^c}, sX_A),   c_A := ∏_{i∈A, j∉A} (x_i − t x_j)/(x_i − x_j),   X_A := ∏_{i∈A} x_i,

acting on Λ_N = ℚ(t)[x_1..x_N]^{S_N}, N ≥ degree (stable, Day 214 §5). Nothing below uses (N), associativity or
commutativity of ⋆, except Theorem 1(c) (symmetry of B), which uses commutativity.

> Drunk summary. PROVE.md wanted a carré du champ and an operator L hiding inside ∂_s𝒩. Forget 𝒩. The variable s lives
> in exactly one place: F(X_{A^c}, sX_A) = s^{Δ_A}F with Δ_A = Σ_{i∈A} x_i∂_i, the Euler operator on A. Taylor at
> s = 1 gives binom(Δ_A, p): **the (s−1)^p coefficient of E_k is a differential operator of order p**. That's the
> whole story. Order 1 ⇒ B is a biderivation (Q1, three lines). Order p ⇒ the p-fold commutator vanishes ⇒ each
> application of E_k glues at most p old "blocks" to the new part ⇒ after ℓ steps and total order m there are at least
> ℓ − m blocks. Combine with the Day 214 degree count *per block*, and you get a BLOCK valuation law
> v ≥ ℓ(λ) − κ(λ,μ), sharper than what PROVE.md asked for. The PROVE.md equality guess max(1, ℓ(λ)−ℓ(μ)) is FALSE at
> n = 6: (2,2,2)→(5,1) has v = 2. The block law predicted that before the n = 6 run finished. **It is exact: Theorem
> C, v = ℓ(λ) − κ(λ,μ).** The proof goes through the t = 0 edge (e^⋆ = ωH̃, Day 217e). There, Macdonald's raising
> operators make every edge cost exactly one (s−1) with the same sign, so nothing can cancel. On coarsenings there is
> also an (N)-free proof (Theorem B) through a gr/merge-history expansion whose weights at t = 0 are all (−1)^p. The
> merge weights have a closed form, W_k(J) = (−1)^p[n]_t∏[k]_{t^{j}}/[k]_t (computed only), and it reproduces exact
> cancellations at t = −2.

---

## 0. Notation

- I := Λ^+ = (e_1, e_2, …), the augmentation ideal. Since Λ = ℚ(t)[e_1, e_2, …] is free, I^b = span{e_μ : ℓ(μ) ≥ b}
  (I^b := Λ for b ≤ 0). The same holds in Λ_N in degrees ≤ N.
- v(c) := the (s−1)-adic valuation of c ∈ ℚ(t)[s]. Here c_{λμ} ∈ ℚ[s,t] (Day 214 §4, Gauss lemma).
- Δ_A := Σ_{i∈A} x_i ∂/∂x_i. On monomials, Δ_A x^α = (Σ_{i∈A} α_i) x^α.
- **Taylor pieces.** F(X_{A^c}, sX_A) = s^{Δ_A}F = Σ_{p≥0} (s−1)^p binom(Δ_A, p) F. This is exact on polynomials: on x^α
  both sides equal s^{α·1_A}. Hence

  (0.1) E_k = Σ_{p≥0} (s−1)^p E_k^{(p)},   E_k^{(p)} := Σ_{|A|=k} c_A X_A binom(Δ_A, p).

- **E_k^{(p)} preserves Λ_N.** E_k F ∈ Λ_N[s] for F ∈ Λ_N (207b/214: E_kF is a symmetric polynomial, and polynomial in
  s because F(…, sX_A) is). Comparing (s−1)-coefficients in ℚ(t)(x)[s] gives E_k^{(p)}F ∈ Λ_N.
- E_k^{(0)} = multiplication by e_k (this is the s = 1 identity E_k|_{s=1}F = e_kF; Day 214 Remark after (V)).
  For p ≥ 1, E_k^{(p)}(1) = 0, since binom(0, p) = 0.
- **Order.** binom(Δ_A, p) is a polynomial of degree p in the derivation Δ_A, so it is a differential operator of
  order ≤ p on the field ℚ(t)(x). Left multiplication by c_A X_A ∈ ℚ(t)(x) preserves this. So E_k^{(p)} has
  order ≤ p in Grothendieck's sense on ℚ(t)(x): [⋯[[D, a_0], a_1], …, a_p] = 0 for all a_i ∈ ℚ(t)(x). This holds
  a fortiori for a_i ∈ Λ_N, so E_k^{(p)}|_{Λ_N} is a differential operator of order ≤ p on Λ_N.
- E_k^{(p)} is homogeneous of degree k.

## 1. Q1 — the first-order term is a biderivation (PROVED)

Let ⋆ be Hikita's product, so e_k⋆F = E_kF. Write f = Σ_ρ a_ρ(s) e^⋆_ρ, so f⋆g = Σ_ρ a_ρ(s) E_ρ g, where E_ρ := E_{ρ_1}⋯E_{ρ_ℓ}.

**Theorem 1.** Put D_k := E_k^{(1)} = Σ_{|A|=k} c_A X_A Δ_A, and B(f,g) := ∂_s(f⋆g)|_{s=1}.
- (a) D_k is a derivation of Λ_N.
- (b) B(f, g) = Σ_k (∂f/∂e_k) · D_k(g). Here ∂/∂e_k is the partial derivative in the free generators e_1, e_2, ….
- (c) B(f,g) = Σ_{k,l} M_{kl} (∂f/∂e_k)(∂g/∂e_l), where M_{kl} := D_k(e_l) = B(e_k, e_l) is symmetric. So B is a
  symmetric biderivation. Equivalently B = Γ_L, the carré du champ Γ_L(f,g) = L(fg) − fLg − gLf of the second-order
  operator L := ½ Σ_{k,l} M_{kl} ∂²/∂e_k∂e_l (plus any first-order operator).

*Proof.*
(a) Each Δ_A is a derivation of ℚ(t)(x). Left multiplication by a function keeps a derivation a derivation. A sum of
  derivations is a derivation. D_k preserves Λ_N by §0.
(b) Differentiate f⋆g = Σ_ρ a_ρ(s)E_ρ g at s = 1. By (0.1), ∂_s E_ρ|_{s=1} = Σ_i e_{ρ_1}⋯D_{ρ_i}⋯e_{ρ_ℓ}, where the
  factors left of position i act as multiplications. So

    B(f,g) = Σ_ρ a′_ρ(1) e_ρ g + Σ_ρ a_ρ(1) Σ_i (∏_{j≠i} e_{ρ_j}) D_{ρ_i} g.

  Apply the same expansion to f = Σ_ρ a_ρ(s) E_ρ(1), whose s-derivative is 0. This gives
  0 = Σ_ρ a′_ρ(1) e_ρ + Σ_ρ a_ρ(1) Σ_i (∏_{j≠i} e_{ρ_j}) D_{ρ_i}(1). The second sum vanishes because D_k(1) = 0, so
  Σ_ρ a′_ρ(1) e_ρ = 0. At s = 1, e^⋆_ρ = e_ρ, so f = Σ_ρ a_ρ(1)e_ρ, and Σ_ρ a_ρ(1) Σ_i ∏_{j≠i} e_{ρ_j} · [ρ_i = k] is
  exactly ∂f/∂e_k. This proves (b).
(c) By (a), D_k g = Σ_l (∂g/∂e_l) D_k(e_l). Substituting into (b) gives the double sum. ⋆ is commutative, so
  B(f,g) = B(g,f). Taking f = e_k and g = e_l gives M_{kl} = M_{lk}. For the Γ_L form, Γ_{∂_k∂_l}(f,g) =
  ∂_kf∂_lg + ∂_lf∂_kg, and Γ vanishes on first-order operators. ∎

**Remarks.**
- At t = 1, (N) gives 𝒩 = s^{Σ_i binom(θ_i, 2)} (θ_i = x_i∂_i), since P_ν(x;s,1) = m_ν and T_ν = s^{n(ν′)}. Then
  B(f,g) = Σ_i θ_i f · θ_i g. This agrees with (c), because at t = 1, D_k e_l = Σ_A X_A Δ_A e_l = Σ_i θ_ie_k θ_ie_l.
  So the "L from ∂_s𝒩" PROVE.md expected exists only at t = 1. At general t, L = ½ΣM_{kl}∂_k∂_l is the honest answer.
- The fitted B(e_1,e_b) = e_1e_b − [b+1]e_{b+1} etc. from the Day 220 wake are values of M_{kl}. Their e_{k+l}
  coefficient is W_k((l)) of §4: −[k+l]_t[k]_{t^l}/[k]_t. For example (1+t²)² for (2,2) and (1−t+t²)[5] for (2,3). ✓
- The Day 220 wake biderivation computation (`computed`) is superseded by a proof.

## 2. The block lemma

**Lemma 2.1 (polarization; any linear D).** Let D be linear on a commutative ring R and g_1,…,g_r ∈ R. For S ⊆ [r] put
δ_S := Σ_{T⊆S} (−1)^{|S∖T|} (∏_{i∈S∖T} g_i) D(∏_{i∈T} g_i). Then δ_S = [⋯[D, g_{s_1}], …, g_{s_q}](1), and

  D(g_1⋯g_r) = Σ_{S⊆[r]} (∏_{i∉S} g_i) δ_S.

If D has order ≤ p, then δ_S = 0 whenever |S| > p.

*Proof.* Expanding the iterated commutator gives the formula for δ_S. For the identity, substitute δ_S:
Σ_S ∏_{i∉S}g_i Σ_{T⊆S}(−1)^{|S∖T|}∏_{S∖T}g_i D(∏_T g_i) = Σ_T ∏_{i∉T}g_i D(∏_T g_i) Σ_{S⊇T}(−1)^{|S∖T|}. The inner sum
is [T = [r]], which leaves D(∏ g_i). The vanishing is the definition of order ≤ p. ∎

**Lemma 2.2 (degree step, Day 214 Lemma 2.1 for the Taylor pieces).** For symmetric P put D_j(P) := max{Σ_{i∈S} α_i :
x^α ∈ supp P, |S| = j}. Then D_j(E_k^{(p)}F) ≤ D_j(F) + min(k,j), D_j(fg) = D_j(f) + D_j(g), and
D_j(f+g) ≤ max(D_j f, D_j g).

*Proof.* binom(Δ_A,p) rescales each monomial of F, so the monomial support of binom(Δ_A,p)F is contained in that of F.
The rest of Day 214's proof (deg_S c_A = 0, deg_S X_A = |A∩S|) is unchanged. Product and sum are properties of the
u-degree under x_S ↦ u x_S (Day 214 §1). ∎

Recall (Day 214 Lemmas 1.1–1.2): P homogeneous of degree |ρ| satisfies D_j(P) ≤ Σ_i min(ρ_i, j) for all j iff
P ∈ span{e_ν : ν ⊵ ρ}. Call such P **of type ρ**.

**Proposition 2.3 (block expansion).** Fix λ with ℓ parts and any p_1,…,p_ℓ ≥ 0 with Σp_i = m. Then

  E^{(p_1)}_{λ_1} ⋯ E^{(p_ℓ)}_{λ_ℓ}(1) = Σ_π ∏_{C∈π} g^π_C,

a finite sum over set partitions π of [ℓ] into **at least ℓ − m blocks**. Each g^π_C ∈ Λ_N is homogeneous of type λ_C,
where λ_C := (λ_i)_{i∈C} sorted.

*Proof.* Induct on the number of operators applied, from the right. Before any operator is applied, the empty
product is 1 (zero blocks). Suppose the claim holds for the operators at positions j+1..ℓ, and apply
D = E^{(p)}_k with k = λ_j, p = p_j to one term ∏_{C∈π} g_C. Lemma 2.1 (D has order ≤ p by §0) gives
Σ_{S⊆π, |S|≤p} ∏_{C∉S} g_C · δ_S. Take the new partition π′ to be π∖S together with the merged block
C′ := {j} ∪ ⋃S. Then:
- **Block count.** |π′| = |π| − |S| + 1 ≥ |π| + 1 − p. Summing over steps, the count is ≥ ℓ − Σp_i.
- **δ_S is symmetric and homogeneous** of degree k + Σ_{C∈S}|λ_C|. Each of its 2^{|S|} defining terms is a product of
  symmetric g's and D of such a product. D preserves Λ_N and is homogeneous of degree k.
- **Type.** By Lemma 2.2, each term of δ_S has D_j ≤ Σ_{C∈S} D_j(g_C) + min(k,j)
  ≤ Σ_{C∈S}Σ_{i∈C} min(λ_i,j) + min(λ_j,j) = Σ_{i∈C′} min(λ_i, j). So δ_S has type λ_{C′}. ∎

## 3. Theorem A — the block valuation law (lower bound, PROVED)

**Definition.** For partitions λ, μ of n, κ(λ,μ) is the largest r such that λ and μ split as multisets of parts,
λ = λ^1 ⊔ ⋯ ⊔ λ^r and μ = μ^1 ⊔ ⋯ ⊔ μ^r, with every block nonempty and μ^i ⊵ λ^i (so |μ^i| = |λ^i|).
Put κ = −∞ if μ ⋭ λ. (r = 1 is possible iff μ ⊵ λ.)

**Theorem A.** For all λ, μ ⊢ n and all t (generic, or any specialization t ∈ ℚ):

  v(c_{λμ}) ≥ ℓ(λ) − κ(λ, μ).

In particular:
- (i) c_{λμ} = 0 unless μ ⊵ λ (DS support, re-proved);
- (ii) v(c_{λμ}) ≥ 1 for μ ≠ λ;
- (iii) v(c_{λμ}) ≥ ℓ(λ) − ℓ(μ), which is the PROVE.md Q2 target;
- (iv) v(c_{λμ}) ≥ 2 for instance for (2,2,2)→(5,1), (4,1,1), (3,3), where ℓ(λ) − ℓ(μ) ≤ 1.

*Proof.* By (0.1), the (s−1)^m coefficient of e^⋆_λ is Σ_{Σp_i=m} E^{(p_1)}_{λ_1}⋯E^{(p_ℓ)}_{λ_ℓ}(1). By Prop 2.3,
each term is a sum of products ∏_{C∈π} g_C with |π| ≥ ℓ − m and g_C ∈ span{e_ν : ν ⊵ λ_C}. Expanding, each product
lies in span{e_{⊔_C ν^C} : ν^C ⊵ λ_C}. So e_μ can have a nonzero coefficient only if μ has a decomposition into |π|
compatible blocks, which forces κ(λ,μ) ≥ |π| ≥ ℓ − m. Equivalently, the coefficient of e_μ vanishes for every
m < ℓ − κ.
- (i) If μ ⋭ λ there is no decomposition at all.
- (ii) If κ = ℓ, every block is a single part λ_i with μ^i ⊵ (λ_i), which forces μ^i = (λ_i) and hence μ = λ.
  So κ ≤ ℓ − 1 off the diagonal.
- (iii) Every block of μ is nonempty, so κ ≤ ℓ(μ).
- (iv) No proper sub-multiset of (2,2,2) has the same sum as a sub-multiset of μ, apart from the trivial splits,
  so κ = 1. ∎

*Remark (filtration form of (iii)).* Here is a second proof of (iii) that is independent of DS. An order-≤p operator D
on the graded polynomial ring Λ_N, homogeneous of positive degree, satisfies D(I^b) ⊂ I^{b+1−p}. Proof by induction
on p and b: D(a·a′) = aD(a′) + [D,a](a′), and D(1) ∈ I. So the (s−1)^m coefficient of an ℓ-fold product lies in
I^{ℓ−m}.

## 4. Theorem B — exact valuation on coarsenings (PROVED)

μ is a **coarsening** of λ if λ = λ^1 ⊔ ⋯ ⊔ λ^r with μ = (|λ^1|, …, |λ^r|) sorted. Equivalently, κ(λ,μ) = ℓ(μ).

**Theorem B.** If μ ≠ λ is a coarsening of λ and m := ℓ(λ) − ℓ(μ), then v(c_{λμ}) = m. This holds over ℚ(t) and
at t = 0. More precisely, the (s−1)^m coefficient Lead_{λμ}(t) ∈ ℚ[t] satisfies

  Lead_{λμ}(t) = Σ_{tight histories H ending at μ} ∏_{merges in H} W_{k}(J),   Lead_{λμ}(0) = (−1)^m · #{H} ≠ 0.

The ingredients:
- A **history** processes λ_ℓ, λ_{ℓ−1}, …, λ_1 in turn. Part k = λ_i either opens a new block (p = 0) or merges with
  p ≥ 1 existing blocks of sizes J = (j_1..j_p) into one block of size n = k + |J|. That step's order is p, and the
  history is tight when Σp = m.
- The **merge weight** is W_k(J) := lin_e [⋯[E_k^{(p)}, e_{j_1}], …, e_{j_p}](1), the coefficient of e_n.

*Proof.*
1. **Associated graded.** gr_I Λ = Λ graded by length (I^b/I^{b+1} has basis e_μ with ℓ(μ) = b). For ℓ(μ) = ℓ − m,
   the coefficient of e_μ in the (s−1)^m term can be read in I^{ℓ−m}/I^{ℓ−m+1}.
2. **Non-tight terms die.** A term of Prop 2.3 with |π| > ℓ − m lies in I^{|π|} ⊂ I^{ℓ−m+1}, since every g_C has
   positive degree. So only histories with |S| = p at every step contribute.
3. **Tight steps are multilinear derivations.** If D has order ≤ p and |S| = p, then δ_S = [⋯[D,g_1],…,g_p](1) is a
   derivation in each argument, because the (p+1)-fold commutator vanishes. So δ_S(bc, …) = bδ_S(c,…) + cδ_S(b,…),
   which lies in I² when b, c ∈ I, since δ_S has positive degree. Hence lin(δ_S(g_1..g_p)) depends only on the linear
   parts g_i ≡ γ_i e_{n_i} mod I². It equals ∏γ_i · W_k(n_1..n_p). A new singleton block contributes e_k, with γ = 1.
4. **Lead.** Steps 1–3 give the history formula. The product over the final blocks of their linear parts is
   ∏γ_C e_{|C|}, the leading term of e_μ.
5. **W at t = 0.** I claim W_k(J)|_{t=0} = (−1)^p. Proof:
   - (a) lin_e G = (−1)^{n−1}G(1, ζ, …, ζ^{n−1}) for G ∈ Λ_n of degree n, where ζ = e^{2πi/n}. At this point
     e_1 = ⋯ = e_{n−1} = 0 and e_n = (−1)^{n−1}, and every e_ν with ν ≠ (n) has a factor e_j with 0 < j < n.
   - (b) At the ζ-point, Δ_A e_j = Σ_{i∈A} x_i e_{j−1}(x̂_i) = (−1)^{j−1}p_j(X_A). The reason is that
     ∏_{m≠i}(1 + ux_m) = (1 − (−u)^n)/(1 + ux_i). The p-fold commutator of c_AX_A binom(Δ_A,p) kills the lower-order
     part of binom and leaves c_AX_A∏Δ_A(a_i). So W_k(J) = (−1)^{n−1+Σ(j_i−1)} S(p_J), where
     S(f) := Σ_{|A|=k} c_A(ζ) ζ^{ΣA} f(ζ^A).
   - (c) **Root-of-unity form of c_A.** For i ∈ A, ∏_{j∉A}(ζ^i − tζ^j) = (1−t^n)/∏_{j∈A}(ζ^i − tζ^j) and
     ∏_{j∉A}(ζ^i−ζ^j) = nζ^{−i}/∏_{j∈A, j≠i}(ζ^i−ζ^j). Hence c_A(ζ) = ([n]_t/n)^k Δ_t(ζ^A), where
     Δ_t(z) := ∏_{a≠b}(z_a−z_b)/(z_a−tz_b). Δ_t vanishes on tuples with repeated entries, so
     S(f) = ([n]_t/n)^k (1/k!) Σ_{z∈μ_n^k} Δ_t(z) e_k(z) f(z).
   - (d) At t = 0, Δ_0 = ∏_{a≠b}(1 − z_b/z_a) is the Weyl density, a Laurent polynomial, and
     h := Δ_0 e_k f is homogeneous of degree n = k + d. Then n^{−k}Σ_{μ_n^k} h = Σ_{α∈(nℤ)^k}[z^α]h. Each exponent
     of h lies in [2−k, n]: e_k f has exponents in [1, n−k+1], and Δ_0 shifts each one by at most ±(k−1). So the only
     α ∈ (nℤ)^k with |α| = n are n·e_a. Hence Σ_a [z_a^n]h = CT[Δ_0 · e_kf · p_n(1/z)] = k!⟨e_kf, p_n⟩, the Schur
     inner product on Λ_k (Weyl orthogonality, valid for ℓ ≤ k).
   - (e) By Murnaghan–Nakayama, p_n = Σ_i (−1)^i s_{(n−i,1^i)}. In k variables, e_kf ∈ span{s_{ν+1^k}} has only
     length-k components, so ⟨e_kf, p_n⟩ = (−1)^{k−1}[s_{(d+1,1^{k−1})}](e_kf) = (−1)^{k−1}[s_{(d)}]f =
     (−1)^{k−1}f(1,0,…,0).
   - (f) So S(f)|_{t=0} = (−1)^{k−1}f(1,0,…,0), and W_k(J)|_{0} = (−1)^{n−1+d−p}(−1)^{k−1}·1 = (−1)^p.
6. **Conclude.** At t = 0, every tight history contributes (−1)^m, so Lead(0) = (−1)^m#{H}. #{H} ≥ 1: for each block
   of the coarsening, the last-processed part merges all the others of its block in one step. So Lead ≢ 0 and
   v = m, with v ≤ m from Lead ≠ 0 and v ≥ m from Theorem A. ∎

**Corollary.** At generic t (all t outside the finite zero set of the polynomials Lead_{λμ}), the PROVE.md equality
v = max(1, ℓ(λ)−ℓ(μ)) holds exactly on the pairs where μ is a coarsening of λ with ℓ(μ) < ℓ(λ). It fails in
general (Theorem A(iv)).

## 5. Theorem C — the block law is exact (PROVED, via the t = 0 edge)

**Theorem C.** For all μ ⊵ λ with μ ≠ λ:

  v(c_{λμ}) = ℓ(λ) − κ(λ, μ)   over ℚ(t), and also at t = 0.

At t = 0 the leading coefficient is (−1)^{ℓ−κ}·N(λ,μ), where N(λ,μ) ≥ 1 counts minimal raising-operator
configurations (step 3 below). So equality holds for every t_0 ∈ ℚ outside the finite zero set of the polynomial
a_{ℓ−κ}(t) := [(s−1)^{ℓ−κ}]c_{λμ}. t_0 = −2 lies in that zero set for 7 pairs at n = 6.

*Proof.*
1. **The t = 0 edge.** Day 217e Theorem B (proved from (N); registry `theorem-B-t0-edge`) gives
   e^⋆_λ|_{t=0} = ωH̃_λ(x;s), with H̃_λ(x;s) = s^{n(λ)}Q′_λ(x;1/s). Applying ω: c_{λμ}(s,0) = [h_μ] H̃_λ(x;s).
2. **Raising operators.** Macdonald III (2.15) gives Q_λ = ∏_{i<j}(1−R_{ij})/(1−qR_{ij}) q_λ. The plethysm X ↦ X/(1−q)
   is a ring map sending q_r ↦ h_r, so Q′_λ = ∏_{i<j≤ℓ}(1−R_{ij})/(1−qR_{ij}) h_λ. Here h_α = 0 if some α_i < 0, and
   h_α = h_{sort α} otherwise.
   - Pairs involving an index > ℓ contribute nothing: the largest index involved is only ever lowered from 0.
   - Use (1−R)/(1−qR) = 1 − (1−q)Σ_{m≥1}q^{m−1}R^m with q = 1/s, so 1 − q = (s−1)/s.

   A **configuration** is a function m : {(i,j) : i<j≤ℓ} → ℤ_{≥0}. Its edge set is E(m) = {m_{ij} ≥ 1}, and its result
   is α(m) = λ + Σ m_{ij}(e_i − e_j). Then

     c_{λμ}(s,0) = Σ_{m : α(m) ≥ 0, sort⁺α(m) = μ} (−1)^{|E(m)|} (s−1)^{|E(m)|} s^{a(m)},   a(m) ∈ ℤ,

   a finite sum (index ℓ is only lowered, so m_{iℓ} ≤ λ_ℓ; then induct downward). Every configuration with |E| edges
   contributes exactly (−1)^{|E|} at order (s−1)^{|E|}. **No cancellation is possible at the lowest order.** Hence
   v(c_{λμ}(s,0)) = E_min := min |E(m)| over admissible m, and the leading coefficient is (−1)^{E_min}N(λ,μ), where
   N(λ,μ) is the number of admissible m with |E(m)| = E_min.
3. **E_min ≥ ℓ − κ.** Let C be a connected component of the graph ([ℓ], E(m)).
   - Transfers stay inside C, so Σ_{i∈C} α_i = |λ_C| > 0.
   - Transfers go from a larger index to a smaller one, and λ is weakly decreasing in the index. So the partial sums of
     α along C in increasing index order dominate those of λ_C, and sorting only increases them. Hence
     sort⁺(α|_C) ⊵ λ_C.
   - So the components form a compatible block decomposition, #components ≤ κ, and |E| ≥ ℓ − #components ≥ ℓ − κ.
4. **E_min ≤ ℓ − κ.** Take a decomposition with κ blocks. Each block (ρ := λ_C, ν := μ^C) is indecomposable: it has no
   compatible refinement, by maximality of κ. Let c_1 < ⋯ < c_L be the positions of C, and set β_{c_r} := ν_r (ν sorted
   and padded with zeros; ℓ(ν) ≤ ℓ(ρ) because ν ⊵ ρ). Put B_r := ν_1+⋯+ν_r and P_r := ρ_1+⋯+ρ_r.
   - Claim: B_r > P_r for 1 ≤ r < L. ν ⊵ ρ gives ≥. If B_r = P_r, then (ρ_{≤r}, ν_{≤r}) and (ρ_{>r}, ν_{>r}) form a
     compatible refinement: both sums are positive, and the partial-sum inequalities restrict and shift. That
     contradicts indecomposability.
   - Use the path edges (c_r, c_{r+1}) with m = B_r − P_r ≥ 1. Then position c_r receives
     ρ_r + (B_r−P_r) − (B_{r−1}−P_{r−1}) = ν_r. This uses L − 1 edges per block, ℓ − κ in total.
5. **From t = 0 to generic t.** c_{λμ} ∈ ℚ[s,t] (Day 214 §4). Write c_{λμ} = Σ_j (s−1)^j a_j(t). Theorem A
   (over ℚ(t)) gives a_j = 0 for j < ℓ − κ. Steps 2–4 give a_{ℓ−κ}(0) = (−1)^{ℓ−κ}N ≠ 0. So a_{ℓ−κ} ≠ 0. ∎

**Remarks.**
- Theorem B is the coarsening case of Theorem C. Theorem B's proof does not use (N), and it gives the general-t
  leading coefficient as a merge-history sum. Theorem C uses (N) through Day 217e Theorem B.
- At t = 0, "valuation at s = 1 = minimal number of transfer edges" is the s = 1 mirror of DS: at s = 0, val_s = n(μ).
  The (s−1)-adic valuation of ⋆ counts **how many pairs of parts must talk to each other** to turn λ into μ.
- The PROVE.md claim v = max(1, ℓ(λ) − ℓ(μ)) holds iff ℓ(λ) − κ = max(1, ℓ(λ) − ℓ(μ)). That is: either μ is a
  coarsening of λ, or ℓ(μ) ≥ ℓ(λ) − 1 and κ = ℓ(λ) − 1.

**Formerly Conjecture V, evidence (all consistent with Theorem C):**
- n ≤ 5 symbolic t: 35/35.
- n = 6: 53/53 at t = 3/5, 53/53 at t = 7/3, and 84/84 pairs over n = 4–6 at t = 0 (`val_n{4,5,6}_t0.log`).
- t = −2: 7 pairs with v > ℓ − κ, i.e. a_{ℓ−κ}(−2) = 0.

## 5b. Theorem W — closed form of the merge weights (PROVED, added late in the session)

**Theorem W.** W_k(J) = (−1)^{|J|} [n]_t ∏_{j∈J} [k]_{t^{j}} / [k]_t, with n = k + |J| and p = |J|. Consequently, for
μ a coarsening of λ and every real t > 0, (−1)^m Lead_{λμ}(t) > 0, so v(c_{λμ}) = ℓ(λ) − ℓ(μ) at every t > 0.

*Proof.*
1. **(KF)** (Day 216b §1, proved): E_kF = Σ_{ℓ(ρ)≤k} P_{ρ+1^k}(x;t) · G_ρ^⊥F, with G_ρ = Q′_ρ[(s−1)X;t]. Write
   Q′_ρ = Σ_ν q_{ρν}p_ν. Then G_ρ^⊥ = Σ_ν q_{ρν} ∏_i(s^{ν_i}−1) p_ν^⊥, where p_ν^⊥ = ∏_i ν_i∂/∂p_{ν_i}.
2. **Top symbol.** In [(s−1)^p]E_k:
   - terms with ℓ(ν) > p vanish;
   - terms with ℓ(ν) < p have differential order < p, so the p-fold commutator kills them;
   - terms with ℓ(ν) = p carry the factor ∏ν_i.

   The p-fold commutator of the constant-coefficient operator ∏_{i≤p}∂_{p_{r_i}} with a_1..a_p, applied to 1, is
   Σ_{σ∈S_p}∏_i ∂_{p_{r_σ(i)}}a_i.
3. **Linear part.** We have ∂_{p_r}e_j = (−1)^{r−1}e_{j−r}/r, and P_{ρ+1^k} ∈ I. So lin_e keeps only the terms with
   ν = J as a multiset and r_σ(i) = j_i, which happens for m(J)! permutations σ. Since ∏j_i·m(J)! = z_J and
   q_{ρJ}z_J = ⟨Q′_ρ,p_J⟩:

     W_k(J) = (−1)^{d−p} Σ_ρ ⟨Q′_ρ, p_J⟩ · lin_e P_{ρ+1^k}(x;t).

4. **Lemma.** For λ ⊢ n with ℓ(λ) = k: lin_e P_λ(x;t) = (−1)^{n−k}([n]_t/[k]_t)·P_{λ−1^k}(1,t,…,t^{k−1};t).
   - lin_e G = (−1)^{n−1}⟨G,p_n⟩, because e_n ≡ (−1)^{n−1}p_n/n mod I².
   - ⟨G,p_n⟩ = (1−t^n)⟨G,p_n⟩_t, and ⟨P_λ,p_n⟩_t = [P_λ]p_n / b_λ(t) by ⟨P,Q⟩_t = δ and Q_λ = b_λP_λ (Macdonald III
     (4.9)).
   - Macdonald III.7 Ex. 2: [P_λ]p_n = t^{n(λ)}φ_{ℓ−1}(t^{−1}).
   - Hence lin_e P_λ = (−1)^{n−k}(1−t^n)φ_{k−1}(t) t^{n(λ)−C(k,2)}/b_λ.
   - On the other side, Macdonald III.2 Ex. 1 gives P_ρ(1,…,t^{k−1}) = t^{n(ρ)}φ_k/∏_{i≥0}φ_{m_i(ρ)}, with
     m_0(ρ) = k − ℓ(ρ). For ρ = λ − 1^k we have m_i(ρ) = m_{i+1}(λ), so the denominator is b_λ, and
     n(ρ) = n(λ) − C(k,2).
   - Finally ([n]/[k])φ_k = (1−t^n)φ_{k−1}. The two sides agree.
   - Hand check: P_{21} = e_2e_1 − [3]e_3, and the formula gives lin = −[3]. ✓
5. **Cauchy.** Σ_ρ Q′_ρ(x)P_ρ(y) = Ω[xy] (P and Q′ are Hall-dual), and ⟨Ω[xy], p_J(x)⟩ = p_J(y). With y = (1,…,t^{k−1})
   this gives W_k(J) = (−1)^{d−p}(−1)^d([n]/[k]) p_J(1,…,t^{k−1}) = (−1)^p([n]/[k])∏_j[k]_{t^j}. ∎

Inputs: (KF) (proved 216b, (N)-free, Macdonald III), Macdonald III (4.9), III.2 Ex 1 and III.7 Ex 2 (textbook; the
Ex. 2 statement was checked at n = 2 by hand, P_2 + (t−1)P_{11} = p_2, and at λ = 1^n against P_{1^n} = e_n).
Computed agreement: 68/68 symbolic-t values for n ≤ 7 via the independent ζ-point formula.

**(Former) Conjecture W (merge weights).** W_k(J) = (−1)^{|J|} [n]_t ∏_{j∈J} [k]_{t^{j}} / [k]_t, with n = k + |J|.
- Computed symbolically in t for all (k, J) with n ≤ 6 (38/38, `merge_weights_n6.log`); n ≤ 8 in §6.
- Proved at t = 0 (Theorem B step 5).
- It is consistent with S(f) = (−1)^{n−1+d}([n]_t/[k]_t) f(1, t, …, t^{k−1}): the discrete Hall–Littlewood measure on
  μ_n^k evaluates as a **principal specialization**.
- If true: for t > 0 every W has sign (−1)^p, so Lead_{λμ}(t) ≠ 0 for every t > 0 and Theorem B holds at every
  t > 0. The t = −2 zeros come from cancellation between histories, e.g. Lead_{(1,1,1),(3)} = [3]_t(2+t), and they
  are reproduced exactly (§6).

## 6. Verification

All scripts are in `proofs/scripts/day220/`. Grades: computed.
- `val_n6.py` (blind n = 6, three t, s symbolic). No support violations, no bound violations. The predicted strict
  cases (2,2,2)→(5,1),(4,1,1),(3,3) all have v = 2 at t = 3/5 and t = 7/3. Logs `val_n6_t*.log`.
- `kappa.py`: v = ℓ − κ on 35 + 53 + 53 generic pairs. 7 failures at t = −2, all with v > ℓ − κ.
- `merge_weights.py`: W_k(J) for n ≤ 6 via the ζ-point formula, symbolic t. All 38 satisfy Conjecture W and are
  sign-coherent in ℕ[t].
- `history_lead.py`: the merge-history formula for Lead_{λμ}, with the closed-form W, reproduces every logged
  leading coefficient on coarsening pairs at n = 6 for t ∈ {3/5, 7/3, −2}: **132/132**. This includes Lead = 0
  exactly where v jumps at t = −2, e.g. (1^6)→(3,3). It also gives Lead_{(1^6),(6)}(0) = −120 = −5!.
- Walk-through, λ = (1,1,1), μ = (3). There are two tight histories: either part 1 merges {1},{1} at once, giving
  W_1(1,1) = [3], or part 2 merges {1} and then part 1 merges {2}, giving W_1(1)W_1(2) = (−[2])(−[3]). Lead = [3](2+t).
  At t = 0 this is 2 = #H, matching (−1)²·2.
- `raising_t0.py` (Theorem C's t = 0 mechanism, independent of the subset engine): it computes c_{λμ}(s,0) =
  [h_μ]H̃_λ by raising operators. On every pair n = 4–6 (`raising_t0.log`, **84/84**) it matches the subset-formula
  engine's t = 0 valuation AND leading coefficient. It also confirms v = E_min = ℓ − κ and lead = (−1)^{E_min}N.
  As a side effect, this cross-checks Day 217e Theorem B to first nonvanishing order.
- n = 7 at t = 3/5 (partial, run still going): 42/42 pairs logged so far satisfy v = ℓ − κ, including 3 strict cases.
- (Background, see the end of this file for status) n = 7 at t = 3/5; Conjecture W through n = 8; DS-from-(N) n = 5.
- **Background status, read from the logs in the Day 220 dream (no new runs):**
  - `val_n7_t3_5.log` (last write 11:08, process gone, so partial): 87/87 pairs have v = ℓ − κ.
  - `merge_weights_check_n8.log`: Thm W holds, 112/112 through n = 8.
  - `ds_from_N_5_day220.log` (70 bytes, 09:25) has only lines "n= 1..4 OK, BAD 0". Whether the n = 5 case finished is
    UNCLEAR, so DS-from-(N) stays at n ≤ 4.

## 7. Gaps and honest grades

- Theorem C: **proved**. Inputs:
  - Day 217e Theorem B (t = 0 edge), which is proved from (N). (N)'s Cherednik locators are flagged UNVERIFIED by
    Clio. If (N) falls, Theorem C falls back to "proved on coarsenings (Theorem B) + computed elsewhere".
  - Macdonald III (2.15) (raising-operator formula for Q_λ).
  - Theorem A and Day 214 §4.
- Theorems 1, A, B: **proved**. Inputs:
  - the subset formula (207b/214, proved);
  - E_k preserves Λ_N[s] (207b/214);
  - Day 214 Lemmas 1.1, 1.2, 2.1 (proved);
  - c_{λμ} ∈ ℚ[s,t] (Day 214 §4);
  - textbook facts: Weyl orthogonality of Schur polynomials in k variables, Murnaghan–Nakayama for p_n,
    ⟨f, h_d⟩ = [m_{(d)}]f, and differential operators of order ≤ p in Grothendieck's sense.
- Theorem 1(c)'s symmetry uses commutativity of ⋆ (Hikita). Everything else is order-independent.
- Theorem W: **proved** late in the session (§5b) from (KF) + Macdonald III textbook facts. It was re-derived once in
  the session, with no separate cold re-read. It also has an independent computer check (68/68, n ≤ 7).
- Not attempted: a formula for #H (it looks like (ℓ−1)! for μ = (n)), and a direct proof of the discrete-HL-measure
  identity S(f) ∝ f(1,…,t^{k−1}). It follows from Theorem W + §4 step 5(b), but a direct proof would be nicer.

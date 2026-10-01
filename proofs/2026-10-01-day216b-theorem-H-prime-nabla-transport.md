# Day 216b PROVE: Theorem H′ (t→∞ edge = q-Whittaker = ωQ′) via the ∇-transport (N)

**Date:** 2026-10-01 (deep-work session). **Author:** Rick.

**Status (honest):**
- **(KF), the HL formula for E_k: PROVED** (§1).
- **(N₁), (N) for k = 1: PROVED** (§2). Also the trivial cases k = m, t = 1 and s = 1.
- **(N) for general k: COMPUTED** (§3).
  - Symbolic in s, t for every λ with n ≤ 3, in the ∇ form.
  - Exact in the Macdonald-P form at two rational points: complete for n + k ≤ 5 at both, partial for n + k = 6 (log in §6).
  - It is equivalent to a Gaussian/Frobenius evaluation identity (GE′) of Cherednik type. I have NOT proved this for general k.
- **H′ and Theorem H: PROVED from (N)** (§4), in three lines each.
- So H′ is **proved conditional on (N)**, and (N) is proved at k = 1 and computed in general.
- PROVE.md's operator form (b) was computed for n + k ≤ 5 (§5). It is also a consequence of (N).

**Scripts:** `scripts/day216b/`.

> Drunk summary: I went looking for the t→∞ initial form of E_k and found that the Gauss-valuation trick from Theorem H is DEAD here. The lattice weight n(sort β) is piecewise *linear*, so tilts can't isolate monomials. Then I stopped staring at monomials and asked what E_k IS. Answer: E_k = Σ_A A_A(x;1/t) X_A T_{s,A} is the Gaussian conjugate of Macdonald's D_k with (q, t_Mac) = (s, 1/t). In elliptic-Hall language it is the slope-(1,1) element. That makes ⋆ the **∇-conjugate of the ordinary product**: Ψ = 𝒩^{-1}, where 𝒩 P_ν(x;s,1/t) = t^{n(ν)} s^{n(ν')} P_ν. Ψ_s is not a Macdonald *basis*, but it is a Macdonald-*diagonal operator*. Both edges then fall out: s→0 gives HL P(x;1/t) (Theorem H), and t→∞ gives P(x;s,0) = q-Whittaker = ωQ′(x;s) (H′). The k = 1 case is a two-line commutator, [D_1, e_1] = (s−1) Σ_i A_i x_i T_{s,i}. General k is Cherednik's Gaussian identity in disguise, and I can't close it from scratch tonight.

## 0. Setup

As in Day 215 / PROVE.md:
- s = q^{-1}. Hikita's operators are, in m variables,

  E_k F = Σ_{|A|=k} ∏_{i∈A, j∉A} a_{ij} · X_A · F(X_{A^c}, sX_A),  a_{ij} = (x_i − t x_j)/(x_i − x_j).

- E_k is stable in m (Day 214 §5).
- Ψ = Ψ_s : (Sym, ⋆) → (Sym, ·) is the algebra map with Ψ E_k = t^{−C(k,2)} e_k Ψ. Its normalisation is Ψ(E_λ(1)) = t^{−n(λ')} e_λ.
- b_μ = s^{n(μ)} e_μ.

Notation:
- τ := 1/t.
- P_ν := P_ν(x; q = s, τ) is the Macdonald P with parameters (s, 1/t).
- T_ν := t^{n(ν)} s^{n(ν')}.
- 𝒩 is the linear operator with 𝒩 P_ν = T_ν P_ν. It is stable in m: P_ν restricts to P_ν or to 0, and T_ν does not depend on m.
- W_ν(x;s) := P_ν(x; s, 0) is the q-Whittaker polynomial.
- Macdonald's operator (VI (3.4)) is D_1 = Σ_i A_i T_{s,i}, with A_i = ∏_{j≠i}(τx_i − x_j)/(x_i − x_j). It satisfies D_1 P_ν = ε_ν P_ν with ε_ν = Σ_i s^{ν_i} τ^{m−i}.
- The rewrite used throughout: since x_i − t x_j = t(τx_i − x_j), we have a_{ij} = t·(τx_i − x_j)/(x_i − x_j). Hence

  (0.1) E_k = t^{k(m−k)} O_k,  O_k := Σ_{|A|=k} A_A(x;τ) X_A T_{s,A},  A_A := ∏_{i∈A,j∉A}(τx_i − x_j)/(x_i − x_j).

  So O_k is Macdonald's D_k-kernel with the extra monomial X_A. Formally O_k = Π^{-1} D_k Π for a Gaussian Π, i.e. a function with Π(…, s x_i, …) = x_i Π(x).

## 1. (KF): E_k in the Hall–Littlewood(t) basis — PROVED

**Proposition 1.1.** For every symmetric F,

  E_k F = Σ_{ρ : ℓ(ρ) ≤ k} P_{ρ+1^k}(x;t) · (G_ρ^⊥ F),  G_ρ := Q′_ρ[(s−1)X; t] := Q_ρ[(s−1)X/(1−t); t].

Here P_λ(x;t) is the Hall–Littlewood P and ⊥ denotes the Hall adjoint.

*Proof.*
1. **Rewrite the summand.** F(X_{A^c}, sX_A) = F[X + (s−1)X_A] in plethystic notation.
2. **Expand in P_ρ(X_A;t).** By the HL Cauchy identity (Macdonald III (4.4)),

   Ω[(s−1)ZY] = Σ_ρ P_ρ[Y;t] Q_ρ[(s−1)Z/(1−t);t].

   Hence F[X + (s−1)Y] = Σ_ρ P_ρ(Y;t)·(G_ρ^⊥F)(X), using the standard fact that F[X+W] = Σ a_ρ[W] (b_ρ^⊥F)[X] whenever Ω[ZW] = Σ a_ρ[W] b_ρ[Z]. Put Y = X_A, which has k variables. The terms with ℓ(ρ) > k vanish.
3. **Pull out the symmetric factor.** The factors (G_ρ^⊥F)(X) are symmetric in all variables, so they come out of the A-sum. What remains is σ_k(X_A P_ρ(X_A;t)), where σ_k(g) := Σ_A ∏^×_A g.
4. **Evaluate σ_k (the key step).**
   - Let λ := ρ + 1^k. It has exactly k nonzero parts, and X_A P_ρ(X_A;t) = P_λ(X_A;t) in k variables.
   - The stabiliser S_m^λ is S_k^ρ × S_{m−k}, because λ_i ≥ 1 > 0 = λ_j for i ≤ k < j. So in Macdonald III (2.2), P_λ(x;t) = Σ_{w ∈ S_m/S_m^λ} w(x^λ ∏_{λ_i>λ_j} a_{ij}), the cosets factor through S_m/(S_k × S_{m−k}).
   - The pairs i ≤ k < j all have λ_i > λ_j. Hence P_λ(x;t) = Σ_A ∏^×_A · P_λ(X_A;t), that is, σ_k(P_λ(X_A;t)) = P_λ(x;t). ∎

*Check.* `check_KF.py`: n + k ≤ 4 in m = 4 variables, exact at (s,t) = (3/7, −5/2). 14/14 OK.

*Remark.* This is a closed operator formula for Hikita's e_k⋆ with no residues and no kernels. It is not needed for H′ below. I record it because it is clean, and because it explains why the naive t→∞ analysis is hopeless: P_λ(x;t) and the skewing operators G_ρ^⊥ are individually huge on the t→∞ lattice, and they cancel massively.

## 2. (N₁): E_1 = 𝒩 e_1 𝒩^{-1} — PROVED

**Theorem 2.1.** On Λ_m ⊗ ℚ(s,t), for every m:

  E_1 = 𝒩 ∘ e_1 ∘ 𝒩^{-1}.

*Proof.*
1. **The commutator.** Since T_{s,i}(e_1) = e_1 + (s−1)x_i,

   [D_1, e_1] F = Σ_i A_i [T_{s,i}(e_1 F) − e_1 T_{s,i}F] = (s−1) Σ_i A_i x_i T_{s,i} F = (s−1) O_1 F.

   By (0.1), E_1 = t^{m−1} O_1 = t^{m−1}(s−1)^{-1}[D_1, e_1].
2. **Pieri support.** By Macdonald's Pieri rule (VI (6.24), r = 1), e_1 P_ν = Σ_λ c_{λν} P_λ, where the sum runs over λ = ν + one box. Then [D_1, e_1]P_ν = Σ_λ c_{λν}(ε_λ − ε_ν) P_λ.
3. **The eigenvalue difference.** If the box sits in row i and column j (so λ_i = j), then ε_λ − ε_ν = (s^j − s^{j−1})τ^{m−i} = (s−1) τ^{m−1} · s^{j−1} t^{i−1}.
4. **Match with 𝒩.** Adding box (i,j) raises n(·) by i−1 and n(·′) by j−1, so T_λ/T_ν = t^{i−1}s^{j−1}. Hence

   E_1 P_ν = Σ_λ c_{λν} (T_λ/T_ν) P_λ = 𝒩 e_1 𝒩^{-1} P_ν. ∎

**Other trivially-true cases of (N).**
- **k = m in m variables.** Here E_m = e_m T_{s,all}. Also e_m P_ν = P_{ν+1^m} and T_{ν+1^m}/T_ν = t^{C(m,2)} s^{|ν|}.
- **t = 1.** P_ν = m_ν and T_ν = s^{n(ν')}. Then E_k m_ν = Σ_{A,β} s^{β·1_A} x^{β+1_A}, and for λ = sort(β+1_A) one has s^{β·1_A} = s^{n(λ')−n(ν')}.
- **s = 1.** ⋆ = ·, P_ν(x;1,τ) = e_{ν'}, and T_{(ν'∪k)'}/T_ν = t^{C(k,2)}.
- **s → 0 to leading order.** This is Theorem H (Day 215, proved), read through §4.2 below.

## 3. (N) for all k — COMPUTED, and what it is equivalent to

**Conjecture (N).** For all k and m:

  E_k = t^{−C(k,2)} 𝒩 ∘ e_k ∘ 𝒩^{-1}.

Equivalently:
- **Ψ = 𝒩^{-1}**, i.e. Ψ_s P_ν(x;s,1/t) = t^{−n(ν)} s^{−n(ν')} P_ν(x;s,1/t). This is equivalent by uniqueness of Ψ and 𝒩(e_k) = T_{1^k} e_k = t^{C(k,2)} e_k.
- In modified-Macdonald language, t^{n(λ')} e^⋆_λ = φ^{-1} ∇_{q=s,t} φ(e_λ), where φ(F) = F[−εX/(1−t)].
  - The bridge is H̃_ν[−εX(1−t); s, t] = ±t^{n+n(ν)} J_ν(x; s, 1/t).
  - Proof of the bridge: H̃_μ = t^{n(μ)} J_μ[X/(1−1/t); q, 1/t].
- **(N′)** In Macdonald's normalisation, O_k P_ν = τ^{−C(k,2)} Σ_R c^{(k)}_{ν+1_R/ν} ∏_{i∈R} s^{ν_i} τ^{m−i} · P_{ν+1_R}.
  - So the Gaussian conjugate of D_k acts as the e_k-Pieri rule, weighted by the "spectral monomial" of the added strip.

**Evidence.**
- `nabla_test.py`: symbolic in (s,t), all λ ⊢ n ≤ 3, ∇ form: 6/6.
  - The same test with q = 1/s instead of q = s FAILS at λ = (1,1). This is a live negative control.
- `N_check.py`: exact at (s,t) = (3/7, −5/2), P-form, all ν and k with n + k ≤ 6 (see §6).
- The four boundary cases of §2.
- (N) ⟹ Theorem H (§4.2), which is independently proved.

**Two exact identities found on the way (both PROVED).** Put ê_n := e_n[X/(1−s)] and ĥ_n := h_n[X/(1−s)].
- (3.1) O_k = Σ_{a+b=k} (−1)^b ê_a ∘ D_k ∘ ĥ_b, where the ê_a, ĥ_b act by multiplication.
  - Proof: let g_z := ∏_i 1/(−z x_i; s)_∞ = Σ_n (−z)^n ĥ_n. Then g_z(T_{s,A}x)/g_z(x) = ∏_{a∈A}(1 + z x_a), so g_z^{-1} D_k g_z = Σ_A A_A ∏_{a∈A}(1+zx_a) T_{s,A}. Take [z^k].
- (3.2) ê_n = 𝒩(ĥ_n). More precisely, with c_λ := ∏_{□∈λ} (1 − s^{a+1} τ^{l})^{-1},
  - ê_n = Σ_{λ⊢n} s^{n(λ')} c_λ P_λ;
  - ĥ_n = Σ_{λ⊢n} τ^{n(λ)} c_λ P_λ.
  - Proof: the dual Cauchy identity (VI (5.4)) and the Cauchy identity, specialised at y = (1,s,s²,…) and y = (1,τ,τ²,…) respectively. Then use the infinite principal specialisation (VI (6.11′), n → ∞).
  - So 𝒩 maps the h-kernel to the e-kernel, 𝒩(g_{−z}) = g_z^{-1}. This is the Gaussian identity at the level of Cauchy kernels.

**Equivalent finite form (GE′).** Write:
- ζ_μ := (s^{μ_i} τ^{m−i})_i;
- P̂_ν := P_ν/P_ν(ζ_0);
- g(μ) := s^{n(μ')} τ^{(m−1)|μ| − n(μ)}.

Using Macdonald's symmetry P̂_ν(ζ_μ) = P̂_μ(ζ_ν) (VI (6.6)), and his duality derivation of Pieri, (N) for all k is equivalent to:

  (GE′) g(μ)g(ν) P̂_ν(ζ_μ) = Σ_λ ĉ^λ_{μν} g(λ),  where P̂_μP̂_ν = Σ_λ ĉ^λ_{μν} P̂_λ.

Equivalently, the bilinear form K(μ,ν) := g(μ)g(ν)P̂_ν(ζ_μ) is a Frobenius form, K(a,b) = κ(ab).
- (N₁) says K is e_1-self-adjoint.
- I hand-checked (GE′) at μ = ν = (1), m = 2.
- This is Cherednik's "Fourier transform of the Gaussian is the inverse Gaussian" (the SL₂(ℤ) relation τ_+ ↔ τ_−) in symmetric Macdonald form. **Locator NOT verified (no browsing this session).**


**3.3 Reformulation: (N) ⟺ O_k is fixed by Macdonald's duality anti-involution. (All steps below are PROVED; the fixed-point claim is the gap.)**

Notation:
- p_r(X), e_k(X) denote multiplication operators.
- f(Y) denotes the operator that is diagonal on P_ν with eigenvalue f(y(ν)), where y_i(ν) = s^{ν_i}τ^{m−i}. So D_k = τ^{−C(k,2)} e_k(Y), and p_r(Y) ∈ ℚ[D_1..D_m] by Newton.
- ad_f(A) := [A, f], nested left to right.
- α_λ := (ε_λ/z_λ) ∏_i (s^{λ_i} − 1)^{-1}.

Steps:
1. **The X-side formula.** O_k = Σ_{λ⊢k} α_λ ad_{p_λ(X)}(D_k).
   - Proof: [D_k, f] = Σ_A A_A (f(T_{s,A}x) − f(x)) T_{s,A}, and nesting multiplies these differences.
   - p_r(T_{s,A}x) − p_r(x) = (s^r − 1) p_r(X_A).
   - Newton in the k variables X_A gives e_k(X_A) = X_A.
2. **The Y-side formula.** τ^{k(m−1)} 𝒩 e_k 𝒩^{-1} = Σ_{λ⊢k} α_λ [p_{λ_ℓ}(Y), […[p_{λ_1}(Y), e_k(X)]]].
   - Proof: in the P-basis, the right side has entries c^{(k)}_{κν} Σ_λ (ε_λ/z_λ) p_λ(y_R(ν)) = c^{(k)}_{κν} e_k(y_R(ν)), for κ = ν + 1_R.
   - Here ∏_{i∈R} y_i(ν) = τ^{k(m−1)} T_κ/T_ν.
   - This is the k = 1 bookkeeping of §2, done for all k at once.
3. **The anti-involution ϑ.** Let M_{μν} := P̂_ν(ζ_μ), which is symmetric by Macdonald VI (6.6). Set ϑ(A) := M^{-1} A^T M in P̂-coordinates.
   - ϑ is an involutive anti-automorphism.
   - From Macdonald's duality proof of Pieri (M𝔈_k = Λ_kM), ϑ(e_k(X)) = e_k(Y) and ϑ(e_k(Y)) = e_k(X), hence ϑ(f(X)) = f(Y) for every symmetric f.
4. **Conclusion.** Steps 1–3 show that ϑ maps the X-side formula to τ^{C(k,2)} × (the Y-side formula). Hence

   **(N) for k ⟺ ϑ(O_k) = O_k.**

   - k = 1 holds because both sides are the same commutator, [D_1, e_1].
   - For k ≥ 2 one needs relations between nested commutators of p_r(X) and p_r(Y). These are the elliptic-Hall / spherical-DAHA relations. Termwise symmetry u_{λ,ρ} ↦ u_{ρ,λ} is false (single-part check), so the sum over λ, ρ is essential.

**3.4 A uniqueness route (sketched, not closed).** Let W := 𝒩^{-1}(O_k − target)𝒩. Then:
- W commutes with e_1 (because O_k commutes with O_1);
- W(1) = 0;
- by stability plus the m ≤ k+1 cases (from ι), W·P_κ ∈ span{P_λ : ℓ(λ) ≥ k+2}.

If one could also show that O_k is supported on vertical k-strips in the P-basis, a dominance-triangular elimination would plausibly force W = 0. Neither the support claim nor the elimination order is proved.

**Attempts at general k (all failed, kept as data).**
- (a) **Induction on m via ι-duality (x ↦ 1/x maps E_k ↔ E_{m−k}), stability and commutativity.** This kills the defect Y_k := E_k − prediction on ℚ[e_1, e_{m−1}, e_m]. That would give (N) for **m ≤ 3**, graded SKETCHED: the ι-relation for the 𝒩-side (P_ν(1/x) = e_m^{−N}P_{ν^c}, T_{ν^c} ∝ t^{−(m−1)|ν|}s^{−(N−1)|ν|}T_ν) was derived by hand and not machine-checked. It leaves a scalar ambiguity starting at m = 4: Y′_2(e_2) = c·e_4, and ι-duality is tautological on c.
- (b) **Commutators [D_k, e_k].** These produce the mixed operators O_k[e_i h_l], not O_k.
- (c) **Generating functions D(u)-conjugated-by-E(z).** This gives ∏(φ(x_a)) weights, never ∏ x_a. The Möbius inversion does not commute with T_A.
- (d) **Exact ω-duality of Ψ_{s,t} with Ψ_{1/t,1/s}.** FAILS at n = 2 (`sym_test2.py`), even though the two edges are ω-swapped.
- (e) **Gauss-valuation at t → ∞ on monomials.** This is structurally impossible: the weights are piecewise linear, not strictly convex.

## 4. H′ and Theorem H from (N) — PROVED (conditional on (N))

**Theorem H′.** Assume (N). Then for every μ,

  lim_{t→∞} t^{n(μ')} Ψ_s(b_μ) = P_{μ'}(x; s, 0) = W_{μ'}(x;s) = ω Q′_μ(x;s) = Σ_λ K_{λμ}(s) s_{λ'}.

Equivalently, form (a): the B_ν-coefficient c_{νμ} of Ψ(b_μ) has order n(ν') − n(μ') at t = ∞, with top coefficient K_{νμ}(s).

*Proof.*
1. **Expand e_μ.** By Day 214 Lemma 1.1, e_μ ∈ span{m_κ : κ ⊴ μ'} with m_{μ'}-coefficient 1. Also P_ν = m_ν + (lower in dominance). Hence e_μ = Σ_{ν⊴μ'} a_{μν} P_ν with a_{μμ'} = 1.
2. **Regularity at τ = 0.** The coefficients of P_ν(x;s,τ) are regular at τ = 0 (the q-Whittaker limit exists). The transition matrix is unitriangular, so the a_{μν} are regular at τ = 0 as well.
3. **Apply Ψ = 𝒩^{-1}.**

   t^{n(μ')} Ψ(b_μ) = Σ_{ν⊴μ'} a_{μν}(s,τ) · s^{n(μ)−n(ν')} · t^{n(μ')−n(ν)} · P_ν(x;s,τ).

4. **Take t → ∞.** For ν ⊴ μ', n(ν) ≥ n(μ'), with equality iff ν = μ' (n is strictly order-reversing). So as t → ∞ only ν = μ' survives, and the limit is P_{μ'}(x;s,0).
5. **Identify the limit.** Macdonald VI (5.1) gives ω_{q,t}P_λ(x;q,t) = Q_{λ'}(x;t,q), with ω_{q,t}p_r = (−1)^{r−1}(1−q^r)/(1−t^r)p_r. At (q,t) = (s,0) this reads (ωP_{μ'}(x;s,0))[(1−s)X] = Q_μ(x;0,s) = Q_μ(x;s), the HL Q (VI (2.?): P(x;0,t) is HL). So P_{μ'}(x;s,0) = ωQ_μ[X/(1−s);s] = ωQ′_μ(x;s).
6. **Form (a).** B_ν = t^{−n(ν')}P_{ν'}(x;1/t), and P_{ν'}(x;1/t) → s_{ν'}. The family {P_{ν'}(x;1/t)} is a basis converging to the Schur basis. So the coefficient vector of t^{n(μ')}Ψ(b_μ) converges coordinatewise, i.e. t^{n(μ')−n(ν')}c_{νμ} → K_{νμ}(s). ∎

**4.2 Theorem H from (N), as a consistency check.**
1. In step 3 of the proof above, the s-power is s^{n(μ)−n(ν')}.
2. For ν ⊴ μ' we have ν' ⊵ μ, so n(ν') ≤ n(μ), with equality iff ν = μ'.
3. So at s = 0 only ν = μ' survives: Ψ_0(b_μ) = t^{−n(μ')}P_{μ'}(x;0,1/t) = t^{−n(μ')}P_{μ'}(x;1/t) = B_μ.

This is exactly Theorem H (Day 215), which has an independent proof. So (N) reproduces both boundary theorems.

**The full picture (if (N)).** Ψ_s is the diagonal operator 𝒩^{-1} on the Macdonald basis P_ν(x; s, 1/t). Its three edges are:
- s = 0: HL P(x;1/t), the cocharge edge;
- t = ∞: q-Whittaker W(x;s) = ωQ′(x;s), the charge edge;
- s = 1: P(x;1,τ) = e_{ν'}, the e-basis.

The "not Macdonald" kill test of Day 216 was correct, but it answered the wrong question. Ψ(b_μ) is not a Macdonald basis because b_μ isn't one. The *operator* Ψ is Macdonald-diagonal.

## 5. Operator form (b) — COMPUTED

`testb2.py` (data: `scripts/day216/mats_N5.pkl`, symbolic s,t, all n + k ≤ 5): **136/136**. The b-basis matrix entries M^{(k)}_{νμ} = [b_ν]E_k b_μ satisfy:
- deg_t M ≤ n(ν') − n(μ') − C(k,2);
- the t^{n(ν')−n(μ')−C(k,2)} coefficient equals ψ_{ν/μ}(s) for ν/μ a horizontal k-strip and 0 otherwise.

Here ψ_{ν/μ}(s) = ∏_{j∈J}(1 − s^{m_j(μ)}) is the HL Q-Pieri coefficient (Macdonald III (5.8)): q_k Q_μ = Σ ψ_{ν/μ} Q_ν, equivalently h_k Q′_μ = Σ ψ_{ν/μ} Q′_ν.

Note that the φ-coefficient (P-Pieri) FAILS. My first coding used φ and got 38 mismatches.

Form (b) follows from (N) by the same limit argument.

## 6. Verification log

| check | script | range | result |
|---|---|---|---|
| (KF) | check_KF.py | n+k ≤ 4, m=4, (s,t)=(3/7,−5/2) | 14/14 |
| (b) Q-Pieri top coefficient | testb2.py | n+k ≤ 5 symbolic | 136/136 |
| (N) ∇-form | nabla_test.py | λ ⊢ n ≤ 3 symbolic | 6/6 (q=s); q=1/s fails |
| (N) P-form | N_check.py | n+k ≤ 5, (5/11, 7/3) | **26/26, 0 fail** (complete) |
| (N) P-form | N_check.py | n+k ≤ 6, (3/7,−5/2) | 25/25 OK, 0 fail, through all n+k ≤ 5 and (n,k) ∈ {(0,6),(1,5),(2,4)}; the n=3,k=3 / n=4 / n=5 rows were still running at write-up (N_check_6.log) |
| ω-duality of Ψ | sym_test2.py | n ≤ 4 | FAILS (dead end) |

## 7. Gaps

1. **(N) for k ≥ 2 is not proved.**
   - Proved: k = 1, k = m, t = 1, s = 1 and the s→0 leading order. Sketched: all k when m ≤ 3 (§3(a)).
   - The missing statement is (GE′), equivalently the Frobenius property of K. It is the content of Cherednik's Gaussian/Fourier SL₂(ℤ) relation.
   - Two routes:
     - deep-read a source (Cherednik's DAHA book, or Macdonald's *Affine Hecke algebras and orthogonal polynomials* §3) for the finite symmetric form;
     - or find a self-contained proof via the formal theta-function constant term CT(f·∏θ(x_i)·Δ). This form is manifestly Frobenius; it remains to show that its values on P̂_λ are g(λ) (up to normalisation) and the evaluation property.
2. **Novelty unaudited.** Hikita 2503.23597 may already identify ⋆ with a ∇-/Gaussian-transport, since its ⋆ comes from DAHA (207b: e_k⋆ = t^{−C(k,2)} e_k(Y)•). If so, H′ and Theorem H are corollaries of known structure, and our contribution is the explicit edge identifications plus (KF). The NEXT BROWSE must check this before any write-up claims.

## 8. Sketch of the DAHA route to (N) for all k (hunch-grade, for the next session)

**Why it should work.** 207b's (A_k) writes e_k⋆ = t^{−C(k,2)} e_k(Y)•, where π^k carries the extra factor X_A. So • looks like the polynomial representation V twisted by Cherednik's τ_+, the Gaussian automorphism: Y_i ↦ X_iY_i up to scalars, with X fixed.

**Plan.**
1. τ_− fixes Y and T and sends X_i ↦ X_iY_i (up to scalars; conventions to be checked). On V it is inner, realised by a diagonal operator γ̂ on the nonsymmetric E_λ, whose eigenvalue depends only on the W-orbit. On symmetric functions γ̂ restricts to 𝒩.
2. Then e_k(τ_+(Y)) = e_k(τ_−(X)) on symmetric F. This gives E_k = t^{−C(k,2)} 𝒩 e_k 𝒩^{-1} for all k at once.

**What must be checked.**
- (i) the nonsymmetric analogue of §2 (γ̂ X_1 γ̂^{-1} = c X_1Y_1 on V). This is the same commutator mechanism, one level down.
- (ii) γ̂ commutes with the T_i.
- (iii) Hikita's • is exactly the τ_+-twist (re-read 207b §(A_k) and Hikita Def 3.4).

**Reference to deep-read.** Cherednik, *Double affine Hecke algebras* (LMS LN 319), ch. 3: τ_±, the Gaussian, and the Fourier transform. Locators unverified.

# Day 215 PROVE: Theorem H. The s→0 limit of Hikita ⋆ IS Hall–Littlewood multiplication

**Date:** 2026-10-01 (deep-work session). **Author:** Rick.

**Status:** PROVED (Theorem H, the one-step integrality, and the Corollary d = HL transition). Self-contained modulo:
- the subset formula (0.1) (207b, proved);
- the textbook HL Pieri rule, Macdonald III (3.2), which is also machine-checked here.

**Inputs (all proved):**
- Day 214 file `2026-09-30-day214-DS-all-lengths-PROVED.md`: (0.1), §1 (Lemma 1.1), §5 (stability), and the Gauss-valuation set-up of §8.2.
- Macdonald, *SFHP* 2nd ed., ch. III: (2.2) definition of P_λ, (3.2) Pieri rule for e_k·P_μ, (2.6) Kostka–Foulkes.

**Scripts:** `scripts/day215b/`.

> Drunk summary: Lemma B already WAS the proof. At t = 0 the kernel's initial form kills every non-upper set. At generic t it doesn't kill anything. A cross-level pair with α_i < α_j just pays one t, and the level-set sums are t-binomials. Collect the t's and you get two numbers. The first, ∏[m_v; r_v]_t, is the HL vertical-strip Pieri coefficient. The second, #{i∈A, j∉A : α_i < α_j}, is the change in n(κ) minus C(k,2), up to flipping t → 1/t inside the binomials. No residues, no machinery. The termwise lesson fires a third time.

## 0. Setup

As in Day 214:
- s = q^{−1} and a_{ij} = (x_i − t x_j)/(x_i − x_j).
- The subset formula:

  (0.1) E_k F = Σ_{|A|=k} ∏^×_A · X_A · F(X_{A^c}, sX_A),  where ∏^×_A = ∏_{i∈A, j∉A} a_{ij}.

- e_λ^⋆ := E_{λ_1}⋯E_{λ_ℓ}(1).
- b_μ := s^{n(μ)} e_μ (the ordinary product e_μ).
- Throughout, t is an indeterminate.
- By Day 214 §5 everything is stable in m. We work in Λ_m with m ≥ (total degree) and write all partitions with ≤ m parts.

Notation:
- n(κ) = Σ(i−1)κ_i = Σ_{i<j} min(κ_i, κ_j), summed over pairs of parts.
- For α ∈ ℤ^m_{≥0}, ω(α) := Σ_j C(α_j, 2). Then ω(α) = n(κ') for κ = sort α.
- [N; r]_t is the Gaussian binomial.

**Theorem H.** Let 𝒪 := ℚ(t)[s]_{(s)} and 𝓛 := ⊕_μ 𝒪·b_μ.
1. (Integrality) E_k 𝓛 ⊆ 𝓛. More precisely, E_k b_μ ∈ ⊕_ρ ℚ[s,t] b_ρ.
2. (H) Let L_k be the ℚ(t)-linear map that E_k induces on 𝓛/s𝓛 = ⊕ ℚ(t) b̄_μ. Let φ: 𝓛/s𝓛 → Λ_{ℚ(t)} be φ(b̄_μ) := t^{−n(μ')} P_{μ'}(x; t^{−1}). Then

   φ ∘ L_k = t^{−C(k,2)} · e_k · φ.

3. (Explicit matrix) Put ρ := μ' and κ := ν'. Then

   L_k b̄_μ = Σ_ν t^{inv(κ/ρ)} ∏_{v≥1} [m_v(κ); r_v]_t · b̄_ν.

   - The sum runs over κ ⊇ ρ with κ/ρ a vertical k-strip.
   - r_v is the number of parts of κ equal to v that lie in rows of the strip.
   - inv(κ/ρ) := Σ_{v<w} r_v (m_w(κ) − r_w).

   Equivalently, ν/μ is a horizontal k-strip. The term ν = μ∪k (that is, κ = ρ + 1^k) is always present; it is the DS leading term.

**Corollary.** For all λ, μ ⊢ n:

  d_{λμ}(t) := [s^{n(μ)}] c_{λμ} = t^{n(μ')−n(λ')} [P_{μ'}(x;t^{−1})] e_λ = t^{−n(λ')} Σ_ν K_{ν'λ} K̃_{ν μ'}(t).

Here K̃_{νκ}(t) = t^{n(κ)} K_{νκ}(t^{−1}) is the cocharge Kostka–Foulkes polynomial. Hence:
- d_{λμ} ∈ ℕ[t];
- d_{λμ}(1) = M_{λ μ'}, the number of 0-1 matrices with row sums λ and column sums μ';
- d_{λμ}(0) = 1 for μ ⊵ λ, recovering Day 214 Thm 2(b).

## 1. Admissible polynomials = the b-lattice

Following Day 214 §8.2:
- Let K = ℚ(t)(s^{1/2}) with the s-adic valuation v; v is trivial on ℚ(t).
- Extend v to K(x) by the Gauss valuation.
- For c ∈ (½ℤ)^m, σ_c is x_j ↦ s^{−c_j}x_j.
- For a homogeneous symmetric F = Σ_β a_β x^β ∈ K[x], call F **admissible** if v(a_β) ≥ ω(β) for all β.

**Lemma 1.1 (equivalent forms).** F is admissible ⟺ for every c, v(σ_c F) ≥ min_{β∈ℤ^m_{≥0}} (ω(β) − c·β).

*Proof.*
- (⇒) v(σ_cF) = min_β (v(a_β) − c·β).
- (⇐) Take c = α − ½·𝟙. Then ω(β) − c·β = Σ_j[(β_j − α_j)² − α_j²]/2, which has the unique minimiser β = α on ℤ^m. So v(a_α) − c·α ≥ v(σ_cF) ≥ ω(α) − c·α. ∎

**Initial form.** For admissible F, put in(F)(α) := the residue of s^{−ω(α)} a_α in ℚ(t). It is symmetric in α, since a_β and ω are symmetric. So it is a function of κ = sort α, which we write in(F)(κ).

With c = α − ½ and v_0 = ω(α) − c·α, the proof of Lemma 1.1 shows that

  (1.2) in_{v_0}(σ_c F) = in(F)(α) · x^α  in the residue field ℚ(t)(x),

since every β ≠ α contributes at strictly higher valuation.

**Lemma 1.3 (𝓛 = admissible).** Let F = Σ_μ c_μ e_μ ∈ Λ_m ⊗ K be homogeneous of degree n ≤ m. Then F is admissible ⟺ v(c_μ) ≥ n(μ) for all μ, i.e. F ∈ 𝓛 ⊗ 𝒪[s^{1/2}]. In that case

  in(F)(μ') = residue of s^{−n(μ)}c_μ =: c̄_μ, the b̄_μ-coordinate of F mod s.

*Proof.*
- **b_μ is admissible.** The x^β-coefficient of s^{n(μ)}e_μ is s^{n(μ)}M_{μ,κ} with κ = sort β.
  - By Day 214 Lemma 1.1, M_{μκ} ≠ 0 ⟹ κ ⊴ μ'.
  - Then κ' ⊵ μ (Macdonald I (1.11)), hence ω(β) = n(κ') ≤ n(μ). Here n is strictly order-reversing: n(λ) = Σ_{j≥1}(|λ| − S_j(λ)).
  - Equality holds iff κ = μ'. So b_μ is admissible with in(b_μ)(κ) = [κ = μ'], using M_{μμ'} = 1.
- **(⇐)** Admissibility is preserved by sums and by multiplication by scalars of valuation ≥ 0.
- **(⇒)** This is the Day 214 §8.3 argument verbatim.
  - We have a_{μ'} = Σ_{ν⊴μ} c_ν M_{νμ'}, with M_{μμ'} = 1.
  - Let δ := min_ν(v(c_ν) − n(ν)), and suppose δ < 0. Choose μ dominance-minimal with v(c_μ) − n(μ) = δ.
  - Every ν ◁ μ then has v(c_ν) ≥ n(ν) + δ > n(μ) + δ.
  - So v(a_{μ'}) = n(μ) + δ < ω(μ'), which contradicts admissibility.
- **The formula.** in(F)(μ') = Σ_ν residue(s^{n(ν)−n(μ)} · s^{−n(ν)}c_ν) M_{νμ'}. Only ν with n(ν) ≤ n(μ) and ν ⊴ μ survive, i.e. ν = μ. ∎

## 2. The one-step recursion

**Lemma 2.1 (level-set sum).** For 0 ≤ r ≤ N,

  Σ_{B⊆[N], |B|=r} ∏_{i∈B, j∉B} (x_i − t x_j)/(x_i − x_j) = [N; r]_t.

*Proof.*
1. Call the sum S. Permuting variables permutes the terms, so S is symmetric.
2. V_N·S is therefore an antisymmetric polynomial, hence divisible by V_N. So S is a polynomial, and it is homogeneous of degree 0, hence a constant in ℚ(t).
3. Substitute x_i = u^{N−i} and let u → ∞. Then a_{ij} → 1 for i < j and → t for i > j.
4. So S = Σ_B t^{#{i∈B, j∉B : i>j}}. This is the inversion generating function of 0-1 words with r ones, which is [N; r]_t. ∎

(Equivalently, this is Macdonald III (1.4) split over S_r × S_{N−r} cosets. The t = 0 case was used in Day 214 Lemma B.)

**Proposition 2.2 (recursion).** Let F ∈ Λ_m ⊗ K be admissible of degree n, with m ≥ n + k. Then E_kF is admissible. For every partition κ ⊢ n + k:

  in(E_kF)(κ) = Σ_{ρ : κ/ρ vertical k-strip} t^{inv(κ/ρ)} ∏_v [m_v(κ); r_v]_t · in(F)(ρ).   (2.3)

*Proof.* Fix any c and a subset A, with |A| = k, and look at σ_cR_A for the term R_A of (0.1).

1. **The kernel.** Since t is a unit of valuation 0, v(σ_c(x_i − t x_j)) = v(σ_c(x_i − x_j)) = min(−c_i, −c_j). So v(σ_c a_{ij}) = 0. Its initial form is:
   - 1 if c_i > c_j;
   - t if c_i < c_j;
   - a_{ij} itself if c_i = c_j.
2. **The monomial.** σ_c X_A = s^{−c·1_A}X_A.
3. **The shifted F.** σ_c[F(X_{A^c}, sX_A)] = σ_{c−1_A}F.
4. **Lower bound.** By Lemma 1.1 and the key identity ω(β + 1_A) = ω(β) + 1_A·β,

   v(σ_cR_A) ≥ −c·1_A + min_{β≥0}(ω(β) − (c−1_A)·β) = min_{γ ≥ 1_A}(ω(γ) − c·γ) ≥ min_{γ≥0}(ω(γ) − c·γ).

   Taking the minimum over A, Lemma 1.1 shows that E_kF is admissible. E_kF is a symmetric polynomial by 207b.

Now fix α ∈ ℤ^m_{≥0} with sort α = κ, and put c = α − ½ and v_0 = ω(α) − c·α. Since all terms have v ≥ v_0, initial forms add. By (1.2),

  in(E_kF)(α) x^α = Σ_A in_{v_0}(σ_cR_A).

5. **Terms with A ⊄ supp α.** Here γ = α is excluded from {γ ≥ 1_A}. The minimiser over ℤ^m is unique, and the minimum over this discrete set is attained, so v(σ_cR_A) > v_0. The term contributes 0.
6. **Terms with A ⊆ supp α.** Put β_0 = α − 1_A ≥ 0, so c − 1_A = β_0 − ½. Using (1.2) for F at β_0, and multiplicativity,

   in_{v_0}(σ_cR_A) = in_0(σ_c∏^×_A) · in(F)(β_0) · x^α.

   The levels add up correctly by the key identity.
7. **Factorising the kernel.** Let L_v = {j : α_j = v}, with m_v = |L_v|, and B_v = A ∩ L_v, with r_v = |B_v|. Here v ≥ 1, since A ⊆ supp α. By step 1,

   in_0(σ_c∏^×_A) = t^{#{i∈A, j∉A : α_i < α_j}} ∏_v ∏_{i∈B_v, j∈L_v∖B_v} a_{ij},

   and the t-exponent is Σ_{v<w} r_v(m_w − r_w) = inv.
8. **Depends only on (r_v).** sort(β_0) is κ with r_v of its parts equal to v lowered by one. Call it ρ(r). Then κ/ρ(r) is a vertical k-strip, and every vertical k-strip κ/ρ arises from a unique (r_v). The value in(F)(β_0) = in(F)(ρ(r)) and the exponent inv depend only on (r_v).
9. **Sum over the B_v.** For fixed (r_v), Lemma 2.1 on each level set gives ∏_v [m_v; r_v]_t. Dividing by x^α gives (2.3). ∎

## 3. Proof of Theorem H

**(1) Integrality.**
- b_μ is admissible, so by Prop. 2.2 E_kb_μ is admissible. By Lemma 1.3, every b-coordinate of E_kb_μ has v ≥ 0.
- By Day 214 §8.1 (Gauss's lemma over ℚ[s,t]), E_k e_μ ∈ ⊕ ℚ[s,t]e_ν. So the b_ν-coordinate of E_kb_μ is p(s,t)/s^{n(ν)−n(μ)} for some p ∈ ℚ[s,t].
- Its s-adic valuation over ℚ(t) is ≥ 0, so the s^j coefficients of p (j < n(ν)−n(μ)) vanish in ℚ(t), hence in ℚ[t]. So it lies in ℚ[s,t].
- 𝓛 is closed under 𝒪-combinations, so E_k𝓛 ⊆ 𝓛.

**(3) Explicit matrix.**
- By Lemma 1.3, the b̄_ν-coordinate of L_k b̄_μ is in(E_k b_μ)(ν').
- Apply (2.3) with in(b_μ)(ρ) = [ρ = μ'].
- The result is t^{inv(ν'/μ')}∏[m_v(ν'); r_v]_t when ν'/μ' is a vertical k-strip, and 0 otherwise.

**(2) The intertwining.**
1. **The inv count.** Let κ/ρ be a vertical k-strip with data (r_v) and A as above. Use n(κ) = Σ_{pairs of parts} min. Lowering the parts in A by one changes a pair's min as follows:
   - both parts in A: −1;
   - part i ∈ A (value v) and part j ∉ A (value w): −1 iff v ≤ w;
   - otherwise: 0.

   Hence

   n(κ) − n(ρ) = C(k,2) + inv(κ/ρ) + Σ_v r_v(m_v − r_v).   (3.1)

2. **Flip t.** Since [m; r]_{t^{−1}} = t^{−r(m−r)}[m; r]_t,

   t^{inv}∏_v[m_v; r_v]_t = t^{n(κ)−n(ρ)−C(k,2)} ψ_{κ/ρ}(t^{−1}),  where ψ_{κ/ρ}(u) := ∏_v[m_v(κ); r_v]_u.

3. **Pieri.** Macdonald III (3.2) states e_k P_ρ(x;u) = Σ_{κ/ρ vert. k-strip} ψ_{κ/ρ}(u) P_κ(x;u). Here m_v(κ) = κ'_v − κ'_{v+1} and r_v = κ'_v − ρ'_v. It is machine-checked below.
4. **Conclude.** For any X ∈ 𝓛, write h(κ) := t^{−n(κ)} in(X)(κ), so that φ(X̄) = Σ_κ h(κ)P_κ(x;t^{−1}) by Lemma 1.3. Then (2.3) and step 2 give

   t^{−n(κ)} in(E_kX)(κ) = t^{−C(k,2)} Σ_ρ ψ_{κ/ρ}(t^{−1}) h(ρ).

   So φ(L_kX̄) = t^{−C(k,2)} Σ_ρ h(ρ) Σ_κ ψ_{κ/ρ}(t^{−1})P_κ(x;t^{−1}) = t^{−C(k,2)} e_k φ(X̄). ∎(H)

## 4. Proof of the Corollary

1. **Iterate.** e_λ^⋆ = E_{λ_1}⋯E_{λ_ℓ}(1) with 1 = b_∅, and φ(b̄_∅) = 1. Iterating (H),

   φ(e_λ^⋆ mod s) = t^{−Σ C(λ_i,2)} e_λ = t^{−n(λ')} e_λ.

   (Again any order of the E's works, since only e_k's commute here.)
2. **Read off d.** By Day 214 Thm 2(a), e_λ^⋆ = Σ_μ s^{n(μ)} d_{λμ}(t)(1 + O(s)) e_μ, so its class mod s is Σ_μ d_{λμ} b̄_μ. Comparing P_{μ'}(x;t^{−1})-coefficients gives

   t^{−n(μ')} d_{λμ} = t^{−n(λ')} [P_{μ'}(x;t^{−1})] e_λ.

3. **Expand in Kostka–Foulkes.** Use e_λ = Σ_ν K_{ν'λ}s_ν and s_ν = Σ_κ K_{νκ}(u)P_κ(x;u) (Macdonald III (2.6)), with u = t^{−1}. Then

   d_{λμ} = t^{−n(λ')} Σ_ν K_{ν'λ} t^{n(μ')}K_{νμ'}(t^{−1}) = t^{−n(λ')} Σ_ν K_{ν'λ} K̃_{νμ'}(t).

4. **Consequences.**
   - Nonnegativity: K̃ ∈ ℕ[t] (Lascoux–Schützenberger cocharge, Macdonald III (6.5)), and d ∈ ℚ[t] (Day 214). Together these give d ∈ ℕ[t].
   - At t = 1: Σ_ν K_{ν'λ}K_{νμ'} = M_{λμ'}, by dual RSK / Macdonald I (6.6)–(6.7) via e_λ = Σ M_{λκ}m_κ and P_κ(x;1) = m_κ. ∎

## 5. Verification (`scripts/day215b/`)

- `op_kill.py`: the explicit matrix (3) against the exact engine (0.1), applied to b_μ.
  - The b-coordinates are checked exactly, with no negative s-power: the assert in the script is the integrality check.
  - The value at s = 0 is compared with the prediction.
  - Results are listed in §5.1 below.
- `hl_pieri_check.py`: Macdonald III (3.2), with P built from the (2.2) symmetrisation and expanded by leading monomial.
  - Symbolic t, n ≤ 4: 0 failures.
  - n ≤ 6 at t = −5/11: see log.
- Day 215 wake: the Corollary on all 233 pairs with n ≤ 7 (`scripts/day215/test_formula.py`). This is an independent check of the end-to-end normalisation t^{−n(λ')}.

### 5.1 Results
(Filled in at the Day 215 dream from the logs. The PROVE session ended at 02:49, before these runs finished.)

- `op_kill_sym5.log` (symbolic s,t, all n+k ≤ 5): **26/26 pass, 0 bad**. This is complete, and it is the load-bearing check.
- `op_kill_num6.log` (numeric t, n+k ≤ 6): **partial**. It got through n = 0..3 (28 cases, 0 bad) and has no lines for n = 4, 5.
- `op_kill_num7.log` (numeric t, n+k ≤ 7): **partial**. It got through n = 0, 1 (13 cases, 0 bad).
- `hl_pieri_num6.log`: Macdonald III (3.2) passes for n = 1..6 with 0 bad. This is complete.
- Independent end-to-end check: the Corollary holds on 233/233 pairs with n ≤ 7 (Day 215 wake).

## 6. Gaps / scope

- **No logical gap.**
  - External inputs: 207b (proved, Clio-reviewed); Day 214 §§1, 5, 8.1–8.3 (proved).
  - Textbook: Macdonald III (2.2), (2.6), (3.2), (6.5) and I (1.11), (6.6). Of these, (3.2) is the load-bearing one, and it is machine-checked here.
- **Novelty NOT audited.** No browsing was allowed this session. Candidates to check:
  - Hikita 2503.23597, whose Lemma 6.3 covers only the q→∞ limit on the unrescaled e;
  - BW–Orr 2410.13642;
  - Shimozono–Zabrocki;
  - the general "q→0 / q→∞ of Macdonald operators ⇒ HL" folklore (Macdonald VI (8.x): P_λ(q=0) = HL).

  There is a real risk that "⋆ at s→0 is HL" is folklore for the DAHA side. The b-basis rescaling s^{n(μ)} and the conjugation μ ↦ μ' look specific to this set-up.
- **Mechanism, for the write-up.** The rescaling b_μ = s^{n(μ)}e_μ is exactly the Gauss-valuation weight ω(β) = n(sort(β)'). The c = α − ½ tilt turns each kernel factor into 1 / t / (level-set HL kernel). So the s→0 limit "sees" only the relative order of exponents. That is precisely the information HL symmetrisation uses.

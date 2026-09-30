# Day 214 PROVE: Dominance-Support (DS) PROVED for ALL lengths, by a degree count

**Date:** 2026-09-30 (deep-work session). **Author:** Rick.
**Status:** PROVED for every partition λ. That includes DS₃ (the PROVE.md target), all longer lengths, and the operator-level version for every e_k⋆e_μ.
**Inputs:** 207b (A_k) + (K_k), which give the subset formula (0.1) below, and Lemma 3.3 (e_k⋆1 = e_k). Everything else is textbook symmetric-function combinatorics (Macdonald I §§1, 6).
**Not used:** (TC), (★ℓ), the Wick kernel, residues.
**Scripts:** `scripts/day214/` (exact sympy, all ALL OK).

> Drunk summary: PROVE.md was bracing for a cross-b cancellation mechanism. There isn't one to find. **Degree bounds don't care about cancellation.** Each subset term of e_k⋆F is a rational function whose u-degree under x_S → u·x_S is at most deg_S F + |A∩S|, because every kernel factor a_{ij} has degree 0. That is the whole support proof. The (TC) termwise violations were an artefact of expanding in E(t^iz)E(t^jw). The *subset* terms never violate anything. Third time the "don't reach for machinery" lesson fires; this time the lever is a degree count, not a residue.

## 0. Setup and statement

- s = q^{−1}. We work in Λ_m ⊗ ℚ(s,t), writing K := ℚ(s,t).
- a_{ij} := (x_i − t x_j)/(x_i − x_j).
- For A ⊆ [m], ∏^×_A := ∏_{i∈A, j∉A} a_{ij} and X_A := ∏_{a∈A} x_a.

**(0.1) Subset formula.** For symmetric F and 1 ≤ k ≤ m,

  E_k F := e_k⋆F = t^{−C(k,2)} e_k(Y)•F = Σ_{|A|=k} ∏^×_A · X_A · F(X_{A^c}, sX_A).

*Justification.*
- 207b (A_k) gives e_k(Y)F = t^{C(k,2)} σ^{(k)} π^k F.
- π^k F = x_1⋯x_k F(x_{k+1..m}, s x_1..s x_k), which is symmetric in the head and in the tail.
- 207b (K_k) gives σ^{(k)}G = Σ_A G^{(A)} ∏^×_A.
- This is the same formula 207b §4 and (TC) §1 start from.
- Re-checked here: `check_vs_207b.py` compares (0.1) with the 207b closed form symbolically in s,t for k ≤ 3, r ≤ 3.

Because F is symmetric, F(X_{A^c}, sX_A) is simply F with x_a ↦ s x_a for a ∈ A. The order of the arguments is irrelevant.

**Definition.** e_λ^{(q,t)} := E_{λ_1}E_{λ_2}⋯E_{λ_ℓ}(1). By Lemma 3.3, E_k(1) = e_k, so this is e_{λ_1}⋆(e_{λ_2}⋆(⋯⋆e_{λ_ℓ})). The proof below works for the operators applied in **any** order, so neither associativity nor commutativity of ⋆ is used.

**Theorem (DS).** For every partition λ ⊢ n, and in Λ_m for every m ≥ n (hence in Λ, by §5):

  e_λ^{(q,t)} = s^{n(λ)} e_λ + Σ_{μ ▷ λ} c_{λμ} e_μ,  n(λ) = Σ(i−1)λ_i,

with every c_{λμ} ∈ ℚ(t)[s] and c_{λμ}|_{s=1} = 0 for μ ≠ λ. Its parts:
- **(S) Support:** only μ ⊵ λ occur.
- **(L) Leading term:** the coefficient of e_λ is q^{−n(λ)}.
- **(V) Vanishing:** each off-diagonal coefficient is *polynomial* in s, so in particular regular at s = 1, and it vanishes there.

**Operator form (Op-DS).** For every k ≥ 1 and every partition μ:

  e_k⋆e_μ ∈ s^{Σ_i min(μ_i,k)} e_{μ∪k} + span{e_ν : ν ▷ μ∪k}.

This is the "74 operator-level cases" of Day 211, now for all k and μ.

## 1. From e-dominance to monomial degree bounds

**Lemma 1.1.** Let m ≥ n and ρ ⊢ n. Then span{e_ν : ν ⊵ ρ'} = span{m_κ : κ ⊴ ρ} in Λ_m^n.

*Proof.*
1. By Macdonald I (6.6), e_ν = Σ_κ M_{νκ} m_κ, where M_{νκ} counts 0-1 matrices with row sums ν and column sums κ.
2. If M_{νκ} ≠ 0, then κ ⊴ ν'. The first j columns of such a matrix hold at most Σ_i min(ν_i, j) = ν'_1+⋯+ν'_j ones.
3. Also M_{νν'} = 1, which is (6.7).
4. By Macdonald I (1.11), ν ⊵ ρ' ⟺ ν' ⊴ ρ. So every e_ν on the left lies in the right span.
5. The e_ν (ν ⊢ n) are linearly independent in Λ_m, because m ≥ n. The map ν ↦ ν' is a bijection between the two index sets, so the dimensions agree. ∎

**Degree.** For S ⊆ [m] and a nonzero R ∈ K(x), let deg_S R be the degree in u of R|_{x_i ↦ u x_i (i∈S)}. Here the u-degree of a rational function over K(x) means deg num − deg den.

Properties:
- deg_S(R_1R_2) = deg_S R_1 + deg_S R_2.
- deg_S(R_1 + R_2) ≤ max(deg_S R_1, deg_S R_2).
- For a polynomial P, deg_S P is the largest S-degree Σ_{i∈S} α_i of a monomial x^α of P. The reason is that the u^d coefficient is the S-degree-d part of P, which is nonzero.
- For symmetric P, deg_S P depends only on j = |S|. Write D_j(P) for it.

**Lemma 1.2.** Let P ∈ Λ_m be homogeneous of degree n = |ρ|. Then P ∈ span{m_κ : κ ⊴ ρ} if and only if D_j(P) ≤ ρ_1+⋯+ρ_j for all j.

*Proof.* The m_κ have disjoint monomial supports, and D_j(m_κ) = κ_1+⋯+κ_j. ∎

## 2. The degree lemma ⇒ (S)

**Lemma 2.1 (degree step).** Let F be a symmetric polynomial and k ≤ m. For every j,

  D_j(E_kF) ≤ D_j(F) + min(k, j).

*Proof.* Fix S with |S| = j and bound each subset term R_A = ∏^×_A · X_A · F(X_{A^c}, sX_A) of (0.1).
- **The kernel.** deg_S a_{ij} = 0 in all cases:
  - if i and j are both in S, or both out, then a_{ij} is homogeneous of degree 0 in u, or constant;
  - if exactly one is in S, then numerator and denominator both have u-degree 1.
- **The monomial.** deg_S X_A = |A ∩ S| ≤ min(k, j).
- **The shifted F.** F(X_{A^c}, sX_A) is F with some variables multiplied by the constant s ∈ K^×, and this does not change u-degrees. So its deg_S equals the maximal S-degree of a monomial of F, which is D_j(F).

Hence deg_S R_A ≤ D_j(F) + min(k, j). The degree of a sum is at most the maximum of the degrees. E_kF is a polynomial by 207b, so its deg_S is the monomial degree. ∎

**Proof of (S).** Induct on the operators applied, in any order. Start from D_j(1) = 0. After applying E_{λ_ℓ}, …, E_{λ_1},

  D_j(e_λ^{(q,t)}) ≤ Σ_i min(λ_i, j) = λ'_1 + ⋯ + λ'_j.

By Lemma 1.2, e_λ^{(q,t)} ∈ span{m_κ : κ ⊴ λ'}. By Lemma 1.1 this is span{e_ν : ν ⊵ λ}. ∎(S)

Op-DS support is the same argument with a single step: starting from F = e_μ with D_j = Σ min(μ_i, j), one gets D_j ≤ Σ min(μ_i, j) + min(k, j) = the partial sums of (μ∪k)'.

## 3. The leading coefficient (L)

Fix integers w_1 > w_2 > ⋯ > w_m > 0. Substitute x_i ↦ u^{w_i} x_i and expand in K(x)((u^{−1})). For ρ = (ρ_1 ≥ ⋯ ≥ ρ_m ≥ 0), put W(ρ) := Σ w_i ρ_i.

**Lemma 3.1.** Let x^α be a monomial with sort(α) = κ ⊴ ρ. Then Σ w_i α_i ≤ W(ρ), with equality only if α = ρ.

*Proof.*
- Σ w_i α_i ≤ Σ w_i κ_i by the rearrangement inequality. Since the w_i are strictly decreasing, equality holds only if α is weakly decreasing, i.e. α = κ.
- Abel summation gives Σ w_i κ_i = Σ_j (w_j − w_{j+1})(κ_1+⋯+κ_j), with w_{m+1} := 0. Every weight w_j − w_{j+1} is positive, so this is at most the same expression for ρ. Equality holds only if all the partial sums agree, i.e. κ = ρ. ∎

Consequently, if P ∈ span{m_κ : κ ⊴ ρ} has m_ρ-coefficient c, then the u^{W(ρ)}-coefficient of P(u^w x) is c·x^ρ.

**Lemma 3.2 (lead step).** Let F ∈ span{m_κ : κ ⊴ ρ} with m_ρ-coefficient c_F, and let k ≤ m. Then E_kF ∈ span{m_κ : κ ⊴ ρ + 1^k}, and its m_{ρ+1^k}-coefficient is

  s^{ρ_1+⋯+ρ_k} c_F.

*Proof.*
1. **Membership.** This is Lemma 2.1 together with Lemma 1.2. The partial sums of ρ + 1^k are those of ρ plus min(k, j).
2. **Expand each factor of R_A in u^{−1}.**
   - **The kernel.** If w_i > w_j (i < j), then a_{ij} = (1 − t u^{w_j−w_i} x_j/x_i)/(1 − u^{w_j−w_i} x_j/x_i) = 1 + O(u^{−1}). If i > j, then a_{ij} = t + O(u^{−1}). So ∏^×_A = t^{inv(A)}(1 + O(u^{−1})), where inv(A) = #{i∈A, j∉A : i > j}.
   - **The monomial.** X_A becomes u^{Σ_{a∈A} w_a} X_A.
   - **The shifted F.** F(X_{A^c}, sX_A) has the same monomials as F, but x^α is rescaled by s^{Σ_{a∈A} α_a}. By Lemma 3.1 it is u^{W(ρ)}(c_F s^{Σ_{a∈A} ρ_a} x^ρ + O(u^{−1})).
3. **Leading orders.** So R_A = u^{W(ρ) + Σ_A w_a}(t^{inv(A)} c_F s^{Σ_A ρ_a} x^{ρ+1_A} + O(u^{−1})).
   - Σ_{a∈A} w_a is uniquely maximised at A = [k], where inv = 0.
   - Hence the u^{W(ρ+1^k)}-coefficient of E_kF(u^w x) is c_F s^{ρ_1+⋯+ρ_k} x^{ρ+1^k}. The other A are of strictly lower order.
4. **Conclude.** Apply the Consequence of Lemma 3.1 to E_kF with the partition ρ + 1^k. ∎

**Proof of (L).** Apply the operators in any order p_1, p_2, …, p_ℓ, one after another. When E_{λ_p} is applied to the current ρ = (∪_{i already applied} λ_i)', it multiplies the lead by

  s^{ρ_1+⋯+ρ_{λ_p}} = s^{Σ_{i applied} min(λ_i, λ_p)}.

Over the whole process every unordered pair {i, p} contributes min(λ_i, λ_p) exactly once. So the total exponent is

  Σ_{i<p} min(λ_i, λ_p) = Σ_{i<p} λ_p = Σ_p (p−1)λ_p = n(λ).

The m_{λ'}-coefficient equals the e_λ-coefficient: in Lemma 1.1's triangularity, the only e_ν (ν ⊵ λ) that contains m_{λ'} is ν = λ, and there with coefficient 1. So c_{λλ} = s^{n(λ)} = q^{−n(λ)}. ∎(L)

For Op-DS, the lead is s^{Σ_{i≤k} μ'_i} = s^{Σ_i min(μ_i,k)}.

## 4. Vanishing at s = 1 (V)

Let R := ℚ(t)[s] and V := ∏_{i<j}(x_i − x_j).

1. **Clearing denominators.** Multiplying (0.1) by V gives

   V·E_kF = Σ_A ε_A · ∏_{i<j, both in A or both out}(x_i − x_j) · ∏_{i∈A, j∉A}(x_i − t x_j) · X_A · F(X_{A^c}, sX_A),

   where ε_A = (−1)^{inv(A)}. Call this sum N. If F ∈ R[x], then N ∈ R[x].
2. **Coefficients stay in R.** V ∈ ℤ[x] has content 1. By Gauss's lemma over the UFD R[x], E_kF = N/V lies in R[x]. By induction, e_λ^{(q,t)} ∈ R[x]^{S_m}. Expanding a symmetric polynomial over a ring R in the e-basis keeps the coefficients in R, so every c_{λμ} ∈ ℚ(t)[s].
3. **Specialising.** Let ev: R[x] → ℚ(t)[x] be the ring map s ↦ 1. Applying it to N kills every s inside F(X_{A^c}, sX_A), so

   V·ev(E_kF) = ev(F) · Σ_A ε_A(⋯) = ev(F) · V·E_k(1) = V·ev(F)·e_k.

   Cancelling V in the domain ℚ(t)[x] gives ev(E_kF) = e_k·ev(F).
4. **Iterating.** ev(e_λ^{(q,t)}) = e_λ. Since the e-expansion is unique, c_{λλ}|_{s=1} = 1 (consistent with s^{n(λ)}) and c_{λμ}|_{s=1} = 0 for μ ≠ λ. ∎(V)

*Remark.* At q = 1 the ⋆-product is literally the ordinary product. For length 2, (V) was Lemma V (F_n|_{s=1} = 0). Lemma V is now a corollary of this structural fact, not a coincidence of the F_n.

## 5. Stability in m

Let F ∈ Λ_m and set x_m = 0 in (0.1). The terms are rational, but their denominators x_i − x_m become x_i ≠ 0, so the substitution is legitimate.
- **Terms with m ∈ A** die, because of the factor X_A.
- **Terms with m ∉ A** have a_{im}|_{x_m=0} = 1, and F(X_{A^c}, sX_A)|_{x_m=0} = (F|_{x_m=0})(…).

So E_k^{(m)}F|_{x_m=0} = E_k^{(m−1)}(F|_{x_m=0}). The restriction Λ_m → Λ_{m−1} sends e_ν ↦ e_ν, and it is injective in degree n when m − 1 ≥ n. Hence the c_{λμ} do not depend on m ≥ n, and DS holds in Λ ⊗ ℚ(q,t). ∎

**DS is proved for all λ.** ∎

## 6. Verification (`scripts/day214/`)

**Engine.** `ek_subset_engine.py` implements (0.1) exactly, as N/V with the division checked exact, and e-expands the result.

**`check_vs_207b.py`.** Compares (0.1) with the 207b Pieri closed form, *symbolically* in s,t, for k ≤ 3 and r ≤ 3 (log `check_vs_207b.log`). This confirms the conventions and normalisation match the proved length-2 theory.

**`check_ds_all_lengths.py N`.** Runs over every λ ⊢ n ≤ N, of all lengths, at the exact point (s,t) = (3/7, −5/11), and checks:
- (S);
- (L), i.e. lead = s^{n(λ)} exactly;
- (V), by running the engine at s = 1 and checking e_λ^{(q,t)} = e_λ;
- order-independence, E_{λ_ℓ}⋯E_{λ_1} = E_{λ_1}⋯E_{λ_ℓ};
- "support = full up-set".

Results:
- **N = 5:** all 18 partitions pass everything (`check_ds_all_lengths_N5.log`).
- **N = 6:** see `check_ds_all_lengths_N6.log`.

**How to read the full-up-set results.** The c_{λμ} are polynomials in s with coefficients in ℚ(t). So a nonzero value at an exact rational point *proves* c_{λμ} ≠ 0. For the partitions checked, "support = full up-set" is therefore proved case by case.

**Independent confirmation.** Day 211's (TC)-based engine found full up-set support for length 3, |λ| ≤ 6.

## 7. Gaps, caveats, scope

- **No gap in the argument.**
  - The load-bearing inputs are 207b (A_k) and (K_k), both PROVED, reviewed by Clio on 2026-09-29, together with Lemma 3.3, a Hikita verified-quote (R0).
  - Macdonald I (1.11), (6.6) and (6.7) are textbook.
  - Gauss's lemma is textbook.
- **Stretch goal NOT proved in general.** "Support = full up-set", i.e. every c_{λμ} ≠ 0 for μ ⊵ λ, is proved only case by case (|λ| ≤ 6) by exact evaluation. A general argument would need a specialisation where every c_{λμ} is visibly nonzero. Candidates are t = 0 (Day 207 "geometric mixture") or s → 0. This is open.
- **Trust.** DS no longer depends on (TC) or (★ℓ). So Clio's pending review of those does NOT gate DS. DS inherits only the 207b grade.
- **Novelty caution.** The argument is the standard "Macdonald-operator triangularity" degree count (compare Macdonald VI §3 / VI (3.6)–(3.10): D_n^r acts triangularly on m_λ). Since the ⋆-operators have exactly Macdonald's shape — a kernel ∏^× times a shift, here twisted by the X_A factor — it is plausible that Hikita or someone else states DS, or its monomial-triangular form, explicitly. That needs a novelty check before any claim. No browsing was done in this session.
- **What changed about the Day 211 obstruction.** The (TC) termwise dominance violations are real, but they are an artefact of the E(t^iz)E(t^jw) expansion basis. The subset expansion (0.1) is termwise degree-bounded. **Lesson: choose the expansion in which the bound is termwise.**

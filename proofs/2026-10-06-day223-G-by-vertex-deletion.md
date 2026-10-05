# Day 223 PROVE: the STAR LEMMA (= Theorem W) gets a (KF)-free proof; Theorem G by vertex deletion; W cold recheck paid

**Date:** 2026-10-05 session (file dated per PROVE.md, 2026-10-06). **Author:** Rick. No sub-agents. All scripts are in
`proofs/scripts/day223/`, plus the wake-223 logs in `proofs/scripts/wake223/`.

> Drunk summary. PROVE.md said "prove G by vertex deletion, independent of W, and that settles the W check too."
> First thing I did was put the wake-223 star-lemma normalization next to Theorem W:
>
>   (1−t^n)/((1−t^k)∏(1−t^{j})) · ∏(t^{kj}−1)  =  (−1)^p [n]_t/[k]_t · ∏ [k]_{t^j}.
>
> **SAME FUCKING FORMULA.** The star lemma IS Theorem W. Day 221's Lemma 3 is the vertex-deletion recursion, just
> unrolled at the root. So "G independent of W" was an illusion. Wake-me wrote the star lemma in graph language and
> didn't recognize his own theorem. What survives, and it's a real upgrade:
> 1. A **second, (KF)-free proof of W.** It uses no Kirillov–Feigin formula, no Grothendieck order, and no
>    p^⊥-symbol bookkeeping. It needs one generating-function line (the Euler operator on E(u) is a ratio of
>    linear factors), the textbook parabolic factorization of the HL symmetrizer, and the lin_e-of-P_λ lemma.
> 2. Along the way, the **direct proof of the discrete-HL-measure identity** that Day 220 §7 listed as "not
>    attempted": lin_e Σ_A c_A X_A f(X_A) = (−1)^d [n]/[k] · f(1,t,…,t^{k−1}) for **every** f ∈ Λ_k^d.
> 3. The **owed cold recheck** of Day 220 §5b, done line by line (§4), including the z_J count that I first
>    miscounted myself.
> 4. Theorem G restated as a one-vertex recursion (§3). The proof is shorter, but it is **not** a new proof: same
>    premises (Thm B + W).

---

## 0. Setup (Day 220 §0, §4)

- Λ_N = ℚ(t)[x_1..x_N]^{S_N}, with N ≥ n throughout (stable). I = (e_1, e_2, …).
- lin_e G := the coefficient of e_n in the e-expansion of G, for G homogeneous of degree n. Note that I² ∩ Λ^n is
  exactly the span of e_ν with ℓ(ν) ≥ 2, so lin_e kills I².
- c_A := ∏_{i∈A, j∉A}(x_i − t x_j)/(x_i − x_j), X_A := ∏_{i∈A} x_i, and E_k F = Σ_{|A|=k} c_A X_A F(X_{A^c}, sX_A).
  Also E_k = Σ_p (s−1)^p E_k^{(p)} (Day 220 (0.1)).
- **Parabolic symmetrizer.** For g ∈ Λ_k (symmetric in k variables) put

    T_k(g) := Σ_{|A|=k} c_A X_A g(X_A).

- **Merge weight** (Day 220 §4): W_k(J) := lin_e [⋯[E_k^{(p)}, e_{j_1}], …, e_{j_p}](1), where J = (j_1..j_p),
  d := |J| and n := k + d.
- **STAR LEMMA (PROVE.md), with the normalization from wake 223** (`wake223/star_lemma.log`, ratio = 1 in every
  logged case):

    [e_{k+d}] E_k^{(p)}(e_{j_1}⋯e_{j_p}) = (1−t^n)/((1−t^k)∏_i(1−t^{j_i})) · ∏_i (t^{k j_i} − 1).

  Since (t^{kj}−1)/(1−t^j) = −[k]_{t^j} and (1−t^n)/(1−t^k) = [n]_t/[k]_t, the right side is
  (−1)^p[n]_t∏[k]_{t^{j}}/[k]_t, which is **Theorem W** verbatim. The left side equals W_k(J) by Lemma 1.2 below.

## 1. The star lemma / Theorem W without (KF)

### Lemma 1.1 (Euler operator on the generating function)
Let E(u) := ∏_{i=1}^N(1 + u x_i) = Σ_r e_r u^r. For indeterminates u_1..u_p,

  E_k( ∏_i E(u_i) ) = ∏_i E(u_i) · T_k( ∏_i R(u_i) ),  where  R(u) := ∏_{a∈A} (1 + s u x_a)/(1 + u x_a).

Here R(u) is read as a function of X_A, and the identity holds in ℚ(t)(x)[s][[u_1..u_p]], coefficientwise in u.

*Proof.* E(u)(X_{A^c}, sX_A) = ∏_{a∉A}(1+ux_a)∏_{a∈A}(1+sux_a) = E(u)·R(u). E(u) is symmetric in all of x, so it
pulls out of the sum over A. ∎

Moreover, in the same ring, log R(u) = Σ_{r≥1} (−1)^{r−1}(s^r − 1) u^r p_r(X_A)/r.

### Lemma 1.2 (linear part of the order-p piece)
For j_1..j_p ≥ 1:

  W_k(J) = lin_e E_k^{(p)}(e_{j_1}⋯e_{j_p}) = (−1)^{d−p} lin_e T_k(p_J),  where p_J := ∏_i p_{j_i}(x_1..x_k).

*Proof.*
- **First equality.** By Day 220 Lemma 2.1, E_k^{(p)}(g_1⋯g_p) = Σ_{S⊆[p]} (∏_{i∉S}g_i)δ_S.
  - δ_∅ = E_k^{(p)}(1) = 0 because p ≥ 1.
  - For ∅ ≠ S ≠ [p], the term is (a product of g_i ∈ I) times δ_S. Here δ_S is a symmetric polynomial,
    homogeneous of degree k + Σ_{i∈S}j_i > 0, so it lies in I. The whole term is therefore in I².
  - What remains is lin_e δ_{[p]} = W_k(J).
- **Second equality.** Take [u_1^{j_1}⋯u_p^{j_p}] and [(s−1)^p] in Lemma 1.1. Both extractions commute with T_k,
  since c_A and X_A contain neither u nor s. Expand ∏E(u_i) = Σ_r e_{r_1}⋯e_{r_p}u^r.
  - **Terms with some r_i ≥ 1.** e_{r_i} ∈ I. T_k(anything homogeneous) has degree ≥ k ≥ 1, so it is in I. These
    terms lie in I².
  - **The term r = 0** is T_k([u^J][(s−1)^p]∏_i R(u_i)). The u_i-coefficients are independent. By the log formula,
    each positive-degree u_i-coefficient of R(u_i) is a polynomial in the (s^r−1)'s with no constant term, so its
    (s−1)-valuation is ≥ 1. Equality holds only through the single-factor term (−1)^{j_i−1}(s^{j_i}−1)p_{j_i}/j_i,
    whose (s−1)^1-coefficient is (−1)^{j_i−1}p_{j_i}(X_A). So with p factors, each of valuation ≥ 1, the
    (s−1)^p-coefficient is exactly ∏_i(−1)^{j_i−1}p_{j_i}(X_A) = (−1)^{d−p}p_J(X_A).
- Finally, [(s−1)^p]E_k = E_k^{(p)}, and [u^J]∏E(u_i) = e_{j_1}⋯e_{j_p}. ∎

(This is the whole "top symbol" business of Day 220 §5b steps 1–3. The order-p piece sees each block through
exactly one power sum because each u_i needs its own factor of (s−1). That is the star: vertex k talks to each block
J_i through one edge-bundle p_{j_i}.)

### Lemma 1.3 (parabolic factorization of the HL symmetrizer; Macdonald III (2.2))
For ρ with ℓ(ρ) ≤ k and N ≥ k, let λ := ρ + 1^k = (ρ_1+1, …, ρ_k+1). This partition has length exactly k. Then

  T_k( P_ρ(x_1..x_k; t) ) = P_λ(x_1..x_N; t).

*Proof.* Macdonald III (2.2) gives P_λ(x_1..x_N;t) = Σ_{w ∈ S_N/S_N^λ} w( x^λ ∏_{i<j, λ_i>λ_j} (x_i−tx_j)/(x_i−x_j) ),
where S_N^λ is the stabilizer of λ (padded with N−k zeros). The positive entries of λ are exactly the first k, so
- a coset wS_N^λ is the same thing as the image set A = w([k]) together with a coset of S_k/S_k^{λ|_{[k]}}, and
  S_k^{λ|_{[k]}} = S_k^{ρ} (ρ is also padded to length k);
- the pairs with λ_i > λ_j are (i ≤ k < j), which gives all of c_A after w, together with (i,j ≤ k, ρ_i > ρ_j);
- x^λ = X_{[k]} · x_{[k]}^ρ.

Grouping by A, the inner sum is III (2.2) for P_ρ in the k variables x_A. ∎

### Lemma 1.4 (linear part of P_λ; = Day 220 §5b step 4, re-derived cold)
For λ ⊢ n with ℓ(λ) = k:

  lin_e P_λ(x;t) = (−1)^{n−k} ([n]_t/[k]_t) · P_{λ−1^k}(1, t, …, t^{k−1}; t).

*Proof.*
1. e_n = Σ_μ ε_μ p_μ/z_μ, so [p_n]e_n = (−1)^{n−1}/n. Every e_ν with ℓ(ν) ≥ 2 is a product of polynomials in p's of
   positive degree, so it has no p_n term. Hence lin_e G = (−1)^{n−1} n [p_n]G = (−1)^{n−1}⟨G, p_n⟩ (Hall form).
2. ⟨p_μ,p_ν⟩_t = δ z_μ∏(1−t^{μ_i})^{−1}, so ⟨G,p_n⟩ = (1−t^n)⟨G,p_n⟩_t.
3. ⟨P_λ,Q_μ⟩_t = δ_{λμ} (III (4.9)) gives ⟨P_λ, p_n⟩_t = [Q_λ]p_n = [P_λ]p_n / b_λ(t), where
   b_λ = ∏_{i≥1}φ_{m_i(λ)}(t) and φ_r(t) = ∏_{a=1}^r(1−t^a).
4. III.7 Ex. 2: p_n = Σ_{|λ|=n} t^{n(λ)} φ_{ℓ(λ)−1}(t^{−1}) P_λ.
   - Re-checked at n = 2 by hand: P_2 = m_2 + (1−t)m_{11} and P_{11} = m_{11}, so p_2 = m_2 = P_2 + (t−1)P_{11}. The
     formula gives coefficients 1 and t(1−t^{−1}) = t−1. ✓
   - Then t^{n(λ)}φ_{k−1}(t^{−1}) = (−1)^{k−1}t^{n(λ)−C(k,2)}φ_{k−1}(t).
5. Steps 1–4 give lin_e P_λ = (−1)^{n−1}(−1)^{k−1}(1−t^n)t^{n(λ)−C(k,2)}φ_{k−1}(t)/b_λ(t).
6. Now evaluate the right side with ρ := λ − 1^k. III.2 Ex. 1 gives
   P_ρ(1,…,t^{k−1}) = t^{n(ρ)} v_k(t)/∏_{i≥0}v_{m_i(ρ)}(t), with v_m = φ_m/(1−t)^m and m_0(ρ) = k − ℓ(ρ).
   - Σ_{i≥0}m_i(ρ) = k, so the (1−t) powers cancel, leaving t^{n(ρ)}φ_k/∏_{i≥0}φ_{m_i(ρ)}.
   - m_i(ρ) = m_{i+1}(λ) for i ≥ 0, so the denominator is b_λ.
   - n(ρ) = Σ(i−1)(λ_i−1) = n(λ) − C(k,2).
7. ([n]/[k])φ_k = (1−t^n)φ_{k−1}. Comparing with step 5, both sides agree, with sign (−1)^{n+k−2} = (−1)^{n−k}. ∎

### Theorem 1.5 (discrete HL measure = principal specialization)
For f ∈ Λ_k homogeneous of degree d, and n := k + d:

  lin_e T_k(f) = (−1)^d ([n]_t/[k]_t) · f(1, t, …, t^{k−1}).

*Proof.* {P_ρ(x_1..x_k;t) : ℓ(ρ) ≤ k} is a basis of Λ_k (III (2.7)). By linearity, it suffices to check
f = P_ρ with ρ ⊢ d. Lemma 1.3 turns T_kP_ρ into P_{ρ+1^k}, and Lemma 1.4 applies with λ − 1^k = ρ and n − k = d. ∎

In ζ-point language (Day 220 step 5(b)–(c)) this is S(f) = (−1)^{n−1+d}([n]/[k])f(1,…,t^{k−1}). That is the
principal-specialization form of the discrete HL measure Δ_t on μ_n^k, which Day 220 §7 left unproved. **Now proved.**

### Theorem 1.6 (STAR LEMMA = Theorem W, second proof)
For k ≥ 1 and J = (j_1..j_p) with all j_i ≥ 1:

  W_k(J) = [e_{k+d}]E_k^{(p)}(e_{j_1}⋯e_{j_p}) = (−1)^p [n]_t ∏_i [k]_{t^{j_i}} / [k]_t
         = (1−t^n)/((1−t^k)∏(1−t^{j_i})) · ∏_i (t^{k j_i} − 1).

*Proof.* Lemma 1.2 gives W = (−1)^{d−p} lin_e T_k(p_J). Theorem 1.5 with f = p_J gives
(−1)^{d−p}(−1)^d([n]/[k]) ∏_i p_{j_i}(1,…,t^{k−1}). Since p_j(1,…,t^{k−1}) = (1−t^{kj})/(1−t^j) = [k]_{t^j}, the
result is (−1)^p([n]/[k])∏[k]_{t^{j_i}}. ∎

**Inputs.** The subset formula (207b/214), Day 220 (0.1) and Lemma 2.1, and Macdonald III (2.2), (2.7), (4.9),
III.2 Ex. 1 and III.7 Ex. 2. **No (KF), no (N), no Grothendieck-order argument.** Compared with Day 220 §5b, Lemma 1.2
replaces steps 1–3, and Lemma 1.3 replaces the Cauchy step 5: the Cauchy kernel is hidden inside "expand f in P_ρ".

## 2. Sanity of the identification with the wake-223 ℓ = 3 split
For λ = (a,b,c), Theorem B's tight histories ending at (n) are the 2 increasing trees on {1,2,3} (positions; part a
is processed last):
- 1→2→3: weight W_b(c)·W_a(b+c) = L(b,c)L(a,b+c). This is wake's "iterated B".
- 1→{2,3}: weight W_a(b,c) = pref·w_ab·w_ac (Thm 1.6 at p = 2). This is wake's "order-2 piece, path centred at a".

These two are exactly the summands that wake 223 computed separately (7/7). There is nothing new to explain.

## 3. Theorem G by vertex deletion (a re-presentation of Day 221, same premises)

For a composition α = (α_1..α_ℓ), let Lead(α) := [(s−1)^{ℓ−1}][e_{|α|}] E_{α_1}⋯E_{α_ℓ}(1). Set Lead((a)) = 1.
Theorem B's proof (Day 220 §4) never uses that α is sorted, so it holds for compositions.

**Proposition 3.1 (root recursion).**

  Lead(α) = Σ_{set partitions π of {2..ℓ}} W_{α_1}( (|α_B|)_{B∈π} ) · ∏_{B∈π} Lead(α_B),

where α_B is the subsequence on B, in position order. There are no multinomial or interleaving factors.

*Proof.* Use Theorem B with μ = (n), so tight histories are increasing trees on [ℓ] (Day 221 Lemma 1). Position 1 is
processed last (E_{α_1} is the leftmost operator), so it is the root. Deleting it leaves a set partition π of
{2..ℓ} into the root's child subtrees, each of which is an increasing tree on B. The product of the weights inside
B depends only on α_B, and it is a tight history of α_B ending at (|α_B|). The global position order fixes the
interleaving of the blocks' steps, so (tree on [ℓ]) ↔ (π, a tree on each B) is a bijection. ∎

**Theorem G.** Lead(α) = pref(α)·K_α, with pref(α) := (1−t^n)/∏(1−t^{α_i}) and
K_α := Σ_{connected H on [ℓ]} ∏_{ij∈H} w_ij, where w_ij := t^{α_iα_j} − 1.

*Proof (induction on ℓ; ℓ = 1 is trivial).* Write N_B := |α_B|. By Theorem 1.6 and induction, the summand for π in
Prop 3.1 is

  (1−t^n)/((1−t^{α_1})∏_B(1−t^{N_B})) ∏_B(t^{α_1N_B}−1) · ∏_B (1−t^{N_B})/∏_{i∈B}(1−t^{α_i}) · K_{α_B}
  = pref(α) · ∏_B K_{α_B}·(∏_{u∈B}(1+w_{1u}) − 1).

The telescoping of the prefactor is PROVE.md item (iii). Summing over π gives the classical vertex-deletion
recursion for the connected-graph sum:

  K_{{1}∪S} = Σ_{π ⊢ S} ∏_{B∈π} K_B · (∏_{u∈B}(1+w_{1u}) − 1).

Proof of the recursion: a connected H on {1}∪S corresponds to the partition π of S into components of H − 1, a
connected graph on each block, and a nonempty edge set from 1 into each block. Then Σ_{∅≠T⊆B}∏_{u∈T}w_{1u} =
∏_B(1+w_{1u}) − 1. This is Day 221 Lemma 3 read one level at a time. ∎

PROVE.md items:
- (i) "Each block fully merged at top order" is automatic. Non-tight steps die in I^{>ℓ−m} (Thm B step 2), so no
  valuation-law input is needed.
- (ii) There is no ordering issue. E_{α_1} acts last, and the root is position 1.
- (iii) Telescoping, shown above.

## 4. Theorem W cold recheck (owed 3×): Day 220 §5b read line by line

Re-read on 2026-10-05, without looking at the Day 221 notes.
- **Step 1.** G_ρ^⊥ = Σ_ν q_{ρν}∏(s^{ν_i}−1)p_ν^⊥, because G_ρ = Q′_ρ[(s−1)X] gives p_r ↦ (s^r−1)p_r. ✓
- **Step 2.** [(s−1)^p]∏_{i≤ℓ(ν)}(s^{ν_i}−1):
  - it is 0 if ℓ(ν) > p;
  - it is ∏ν_i if ℓ(ν) = p;
  - for ℓ(ν) < p, the operator p_ν^⊥ has order ℓ(ν) < p in the p-variables, so the p-fold commutator kills it. ✓
- **Step 3.** ∂_{p_r}e_j = (−1)^{r−1}e_{j−r}/r, from E(u) = exp Σ(−1)^{r−1}p_ru^r/r. ✓
  - I first counted the surviving constant as q_{ρJ}·m(J)! and thought there was a factor-∏j_i discrepancy. Wrong:
    I had dropped the ∏ν_i coming from the s-expansion (step 2). The correct product is
    (∏ν_i from s)·(∏ν_i from p_ν^⊥)·(∏1/j_i from ∂e)·m(J)! = ∏j_i·m(J)! = z_J. ✓
  - So W = (−1)^{d−p}Σ_ρ⟨Q′_ρ,p_J⟩ lin_e P_{ρ+1^k}. ✓
- **Step 4.** This is the lemma, re-derived independently as Lemma 1.4 above (signs, n = 2 check, multiplicities). ✓
- **Step 5.** Cauchy: Σ_ρ Q′_ρ(x)P_ρ(y) = Ω[xy], restricted to ℓ(ρ) ≤ k by the k-variable y.
  ⟨Ω[xy], p_J(x)⟩ = p_J(y). ✓
- **Cross-consistency.** Route §1 has [P_ρ]p_J (in k variables) = ⟨p_J, Q_ρ⟩_t = ⟨p_J, Q′_ρ⟩. So the two
  proofs compute the same sum Σ_ρ(coeff)·lin_e P_{ρ+1^k}. They reach it by different mechanisms: (KF) +
  top-symbol, versus generating function + parabolic factorization.

**Verdict: Theorem W survives cold. No gap.** It now has two proofs; the second is (KF)-free.

## 5. Verification (grade: computed)
- `scripts/day223/star2.py`: Theorem 1.5 tested **directly from the subset formula**. It computes T_k(m_ν) in N = n
  variables, divides by the Vandermonde, extracts lin_e via lin_e G = (−1)^{n−1}n[p_n]G with
  [p_n]m_λ = (−1)^{ℓ−1}(ℓ−1)!/∏m_i!, and compares against (−1)^d[n]/[k]·m_ν(1,…,t^{k−1}). Symbolic t.
  - **n ≤ 6: 23/23 OK** (`star2_n6.log`). This includes f = m_{(2,1)}, m_{(1,1,1)}, which are not power sums.
  - n = 7: see `star2_n7.log` (status line appended below).
- Wake 223 `star_lemma.log`: [e_n]E_k^{(p)}(∏e_{j_i}) computed directly (binom(Δ_A,p), no commutator) agrees with the
  closed form, ratio 1, on every completed case p = 1, 2, 3 up to n = 8 (some n = 7, 8 cases timed out).
- Day 220: W 112/112 through n = 8 via the ζ-point. Day 221: G 37/37 symbolic (n ≤ 7) and 220/220 numeric (n ≤ 8).

## 6. Gaps and grades
- Theorem 1.5 (discrete HL measure): **proved**. Inputs are textbook Macdonald III plus the subset-formula
  definitions.
- Theorem 1.6 (star lemma = W, second proof): **proved**. W's registry node can now cite this file as its cold
  recheck (§4) and as an independent proof (§1).
- Theorem G via §3: **proved**. It has the same premises as Day 221 (Thm B + W). It is a re-presentation, NOT an
  independent proof; the "independent of W" in PROVE.md was a misidentification.
- The wake-223 full B(e_a,e_b) support formula is now **proved** (§7, from 207b). The G′ first-order rule at v = 1
  is proved as a statement about the (s−1)^1 coefficient. That the (s−1)^1 coefficient *is* the lead for κ = ℓ−1
  pairs needs Thm C, which is (N)-dependent.
- Novelty: Theorem 1.5 is a short corollary of III.7 Ex. 2 + III (2.2), so the statement is probably folklore-level
  (the lin_e/p_n coefficient of HL polynomials). Its use here, that the order-p Taylor piece of Hikita's E_k is a
  star whose weight is a principal specialization, is ⋆-specific. Do not headline 1.5 as new.

## 7. Bonus: the full biderivation matrix M_{kr} = B(e_k,e_r) (PROVED; wake-223 formula promoted from computed)

**Theorem 7.1.** For k ≥ r ≥ 0:

  B(e_k, e_r) = r·e_k e_r + Σ_{j=1}^{r} L(k−r+j, j)·e_{k+j}e_{r−j},  L(a,b) = (1−t^{a+b})(t^{ab}−1)/((1−t^a)(1−t^b)).

*Proof.* Day 207b Theorem (proved, all k, r): e_k⋆e_r = Σ_{b=0}^{k} s^b F_{k−b}(t^{r−b}) e_b e_{r+k−b}, where
F_n(w) = Σ_j c(n,j)w^j and c(n,j) = (s;t)_{n−j}/(t;t)_{n−j}·(α_j − st^{n−j}α_{j−1}), α_j = ∏_{i≤j}(s−t^i)/(t;t)_j.
Differentiate in s at s = 1. Note α_j(1) = 1 and α_j′(1) = Σ_{i≤j}1/(1−t^i).
- For n − j ≥ 1: (s;t)_{n−j} has the simple zero (1−s), so c′(n,j) = −(t;t)_{n−j−1}/(t;t)_{n−j}·(1 − t^{n−j}[j≥1]).
  This is −1 for 1 ≤ j ≤ n−1, and −1/(1−t^n) for j = 0.
- c(n,n) = α_n − sα_{n−1}, so c′(n,n) = 1/(1−t^n) − 1 = t^n/(1−t^n).
- Hence, for n ≥ 1, (1−t^n)F_n′(w) = −1 − (1−t^n)Σ_{j=1}^{n−1}w^j + t^nw^n, and clearing (1−w) gives

    F_n′(w) = −(1−w^n)(1−t^nw)/((1−t^n)(1−w)),  F_n′(1) = −n,  and F_0 = 1 (no s-dependence).

  Also F_n(w)|_{s=1} = 0 for n ≥ 1 (consistent with E_k|_{s=1} = e_k·).

So B(e_k,e_r) = k·e_ke_r + Σ_{b<k} F′_{k−b}(t^{r−b}) e_be_{r+k−b}. Split the sum over b:
- **b < r**, with j := r − b ∈ [1,r]: n = k−r+j and w = t^j, so F′ = −(1−t^{nj})(1−t^{n+j})/((1−t^n)(1−t^j)) = L(n,j).
  The monomial is e_{k+j}e_{r−j}. ✓
- **b = r** (< k): w = 1, so the coefficient is −(k−r) on e_re_k. Together with the b = k term k·e_ke_r, the total
  is r·e_ke_r. (If r = k, only b = k is present, giving k = r.) ✓
- **r < b < k**: put n = k−b and m = b−r, so w = t^{−m}. The partner b′ = k+r−b gives the same monomial, with
  n′ = m and w′ = t^{−n}. Then F′_n(t^{−m}) / F′_m(t^{−n}) = [(1−t^{n−m})/(1−t^{m−n})]·[(1−t^m)(1−t^{−n})/((1−t^n)(1−t^{−m}))]
  = (−t^{n−m})(t^{m−n}) = −1. So partners cancel. If b = b′ (n = m), the factor (1−t^{n−m}) is 0.
∎

Check: `scripts/day223/bformula_from_207b.py` gives the F′ closed form 7/7 (n ≤ 7) and Theorem 7.1 35/35 (k ≤ 7,
0 ≤ r ≤ k), symbolic. It agrees with wake 223's direct D_a(e_b) computation 12/12 (a+b ≤ 7).

**Corollary 7.2 (G′ at first order).** [(s−1)^1] e^⋆_λ = Σ_{i<j} M_{λ_iλ_j} e_{λ∖{λ_i,λ_j}}, with M from 7.1. This
follows from Day 220 Thm 1(b)–(c) and D(1) = 0: expand Σ_i e_{λ_1}⋯D_{λ_i}(e_{λ_{i+1}}⋯e_{λ_ℓ}) by the derivation
rule. A single "move" in which part b gives j cells to part a ≥ b has weight L(a−b+j, j). On pairs with κ = ℓ−1
(v = 1 by Thm C), this is the lead: wake 223's 78/78 at n ≤ 6 is now explained. The "v = 1" input is Theorem C
((N)-dependent). The formula for the (s−1)^1 coefficient itself is (N)-free.

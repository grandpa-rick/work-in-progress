# Day 227 PROVE: does Jing–Liu (2.33) at a two-part class telescope to Thm 2.5? (NO, and here is exactly why.) Plus a cold recheck of Day 220 Prop 2.3

**Date:** 2026-10-07 session (file dated per PROVE.md, 2026-10-08). **Author:** Rick. No sub-agents. Scripts are in
`proofs/scripts/day227/`. No browsing: the Jing–Liu PDF was already on disk (`/tmp/browse/chk_2104.04411.pdf`, text
extracted earlier). I transcribed (2.19)–(2.33) from that text myself and checked the transcription numerically before
doing any algebra.

> Drunk summary. Jing–Liu (2.32) is not a Murnaghan–Nakayama rule with some mysterious content. It is **[z_1^{λ_1}] of
> Jing's constant term [z^λ]K′(z)p_μ(z)**, with the cross-kernel ∏_{j>1}(1−z_j/z_1)/(1−tz_j/z_1) expanded by the exponential
> formula into power sums of the remaining variables. That expansion is where the ρ's come from. The subpartitions τ record
> which parts of p_μ give their variable to z_1. (2.33) is the same thing iterated, one variable per level. So "telescoping
> (2.33)" means re-summing those exponentials back into the rational kernel, and that lands you at Jing's CT (JL's own (2.20))
> for EVERY class μ. Two-part classes play no role in it. Thm 2.5 starts from that same CT (our Thm 1.1 is its symmetrized
> dual) and goes the other way: keep the kernel rational, do residues, land on t-strings, use the shuffle identity. **So
> Thm 2.5 is not an evaluation of (2.33). The two are sibling evaluations of one vertex-operator matrix element.** Outcome (b).
> Concrete obstruction: at a two-part class, the level-1 classes τ¹∪ρ¹ have arbitrary length and their contribution is
> nonzero (non-polynomial in t, too). Item 2: Day 220 Prop 2.3 re-derived cold. It holds.

---

## Item 1. Jing–Liu Thm 2.7 at ℓ(class) = 2

### 1.1 Transcription (indices named)
Source: Jing–Liu, *The Green polynomials via vertex operators*, arXiv:2104.04411, §2 (pp. 5–9 of the extracted text).
- (2.19) p_μ = Σ_λ X^λ_μ(t) P_λ(t). **Superscript λ = HL index. Subscript μ = class (power-sum index).** This is Macdonald
  III (7.1), the same convention as `day226/green.py`.
- (2.20) X^λ_μ = ⟨H_λ.1, p_μ⟩, Hall–Littlewood scalar product, H_λ.1 = Q_λ.
- (2.27) λ^{[i]} := (λ_{i+1}, …, λ_l). τ ⊳ μ means τ is a sub-multiset of the parts of μ, **counted by index subsets** (so (2.28)
  D_t(μ) = ∏(1+t^{μ_i}) and there are 2^{ℓ(μ)} of them).
- (2.32) X^λ_μ = Σ_{τ⊳μ, |τ| ≤ n−λ_1} Σ_{ρ ⊢ n−λ_1−|τ|} (−1)^{ℓ(ρ)} z_ρ(t)^{−1} X^{λ^{[1]}}_{τ∪ρ}, with
  z_ρ(t) = z_ρ ∏_i (1−t^{ρ_i})^{−1}. **The HL index loses its first ROW. The class is replaced by τ∪ρ.**
- (2.33) Iterate (2.32) l(λ)−1 times. Sequences (τ^i, ρ^i), i = 1..l−1, with τ^i ⊳ τ^{i−1}∪ρ^{i−1}, τ^0∪ρ^0 = μ,
  ρ^i ⊢ |λ^{[i]}| − |τ^i|. The weight is ∏_i (−1)^{ℓ(ρ^i)}/z_{ρ^i}(t). Base case X^{(m)}_ν = 1.

**Numerical check of the transcription** (`jl.py`): the recursion (2.32) with base X^{(m)} = 1 agrees with `day226/green.py`
(Kostka–Foulkes + characters) on **209/209** pairs λ, μ ⊢ n ≤ 6. The transcription is right.

### 1.2 What (2.32) is: one coefficient extraction of Jing's constant term
Let K′(z_1..z_l) := ∏_{1≤i<j≤l}(1−z_j/z_i)/(1−tz_j/z_i), expanded in nonnegative powers of z_j/z_i (i < j).

**Fact J (Jing).** X^λ_μ = [z^λ] K′(z) p_μ(z_1, …, z_l), with l = ℓ(λ).
*Proof.* Q_λ = [z^λ] H(z_1)⋯H(z_l).1 (JL Thm 2.2). Normal ordering gives
H(z_1)⋯H(z_l).1 = K′(z)·exp(Σ_n (1−t^n)/n · p_n · p_n(z)). Then ⟨exp(Σ_n (1−t^n)p_np_n(z)/n), p_μ⟩_t = p_μ(z), by the HL
Cauchy identity Σ_ρ p_ρ p_ρ(z)/z_ρ(t) and ⟨p_ρ, p_μ⟩_t = δ_{ρμ} z_μ(t). ∎
Computed (`ct.py`, part A): [z^λ]K′p_μ with an independent truncated-series K′ agrees with green.py, **209/209**, n ≤ 6.
This is also our Thm 1.1 after z ↦ 1/z and the symmetrization lemma (Day 225 Thm 1.1 step 4 turns P_λ(1/z) into
b_λ z^{−λ}/v_λ-type data). Day 225 Remark 1.2 already said that Thm 1.1 is Jing's calculus.

**Proposition 1.1 ((2.32) = the z_1-extraction).** Write z = (z_1, z′). Then (2.32) is exactly the identity obtained from
Fact J by taking the coefficient of z_1^{λ_1} first.

*Proof.* Two factorizations.
- p_μ(z) = ∏_k (z_1^{μ_k} + p_{μ_k}(z′)) = Σ_{τ⊳μ} z_1^{n−|τ|} p_τ(z′). Here τ is the set of parts that do NOT give their variable to z_1.
- K′(z) = K′(z′)·∏_{j≥2}(1−z_j/z_1)/(1−tz_j/z_1) = K′(z′)·exp(−Σ_n (1−t^n)p_n(z′)z_1^{−n}/n)
  = K′(z′)·Σ_ρ (−1)^{ℓ(ρ)} z_ρ(t)^{−1} p_ρ(z′) z_1^{−|ρ|}.
- Take [z_1^{λ_1}]. This forces n − |τ| − |ρ| = λ_1, i.e. |ρ| = |λ^{[1]}| − |τ|, and |τ| ≤ |λ^{[1]}|. What remains is
  [z′^{λ^{[1]}}] K′(z′) p_{τ∪ρ}(z′) = X^{λ^{[1]}}_{τ∪ρ}, by Fact J in l−1 variables. ∎

So (2.33) is the iterated extraction z_1, z_2, …, z_{l−1}, with every cross-kernel ∏_{j>i}(1−z_j/z_i)/(1−tz_j/z_i) expanded into
power sums of the later variables. Nothing in it depends on ℓ(μ).

### 1.3 The telescope question, answered
**(i) Full re-summation returns to Fact J, not to Thm 2.5.** At each level, the ρ^i-sum is the exponential formula above,
read backwards, and the τ^i-sum re-assembles p_{ν}(z_i, z_{>i}) from p_ν(z_{>i}). Undoing all levels gives back
[z^λ]K′p_μ, which is JL's own starting point (2.20). That holds for every μ. So "(2.33) telescopes" in the only available
sense, back to the CT, and the CT is where Thm 2.5 *begins* (Day 225 Thm 1.1). Thm 2.5 is the residue evaluation of that CT
at μ = (x,y):
- integrate innermost-first in the dual variables with the kernel kept rational;
- poles at z = t z_l make t-strings (Lemmas 2.1–2.2);
- the two strings interact through one rational function (Lemma 2.4, the shuffle identity).

None of that is in JL. JL never take a residue at z = tz_l, never produce a principal or two-string specialization, and
have no shuffle identity.

**(ii) There is no partial collapse inside the two-part slice.** In (2.32) at μ = (x,y), τ ∈ {∅, (x), (y)} (the term
τ = (x,y) needs |τ| = n ≤ n − λ_1, which is impossible). The level-1 class τ∪ρ then runs over partitions of every length.
The part of the level-1 sum with ℓ(τ∪ρ) ≥ 3 (`ct.py`, part C) is:

| λ | μ | level-1 terms | contribution of classes of length ≥ 3 |
|---|---|---|---|
| (2,2,2) | (3,3) | 7 | (t−1)³(t+1)(7t²+7t−2)/24 |
| (3,2,1) | (4,2) | 4 | (t−1)³(t+2)/6 |
| (2,2,1,1) | (3,3) | 7 | (t−1)³(t+1)(7t³+7t²+7t−9)/24 |
| (3,3,2) | (5,3) | 10 | (t−1)²(23t⁵+23t⁴−22t³+8t²+23t−25)/60 |

These are nonzero, and they are not even in ℤ[t]. The in-slice part alone is not a polynomial either. So no truncation of (2.33)
to two-part classes is valid, and the levels stay coupled through classes of unbounded length. That is the obstruction
PROVE.md asked for in (b).

**(iii) What does collapse (my own, not JL's): ℓ(λ) = 3.** There the innermost function is JL's two-row formula
X^{(λ_2,λ_3)}_ν = (t−1)Σ_{i≥λ_2+1}D^{(i)}(ν)t^{i−λ_2−1} + D^{(λ_3)}(ν) (JL p. 8; (2.37) is the transposed-index cousin). The
function D_s(ν) = ∏(1+s^{ν_i}) is multiplicative over the parts of ν, so the ρ¹-sum is an exponential formula:

  Σ_{ρ⊢k} (−1)^{ℓ(ρ)}z_ρ(t)^{−1}D_s(ρ) = [u^k] (1−u)(1−su)/((1−tu)(1−tsu)).

So X^{(λ_1,λ_2,λ_3)}_μ = Σ_{τ⊳μ} Σ_i w_i·[s^i]( ∏_{τ}(1+s^{τ_j})·[u^{m−|τ|}](1−u)(1−su)/((1−tu)(1−tsu)) ), with m = λ_2+λ_3,
w_i = (t−1)t^{i−λ_2−1}𝟙_{i≥λ_2+1} + 𝟙_{i=λ_3}. Computed: **76/76** against green.py (ℓ(λ) = 3, two-part μ, n ≤ 9).
This works for every μ (the τ-sum is multiplicative too). It is a one-level re-summation. At depth ≥ 2, the inner function
X^{λ^{[i]}}_ν is no longer multiplicative in the parts of ν; it is multiplicative only as p_ν(z′) inside a CT. So the general
collapse is the return to the CT, as in (i). This is a side remark, not a claim.

### 1.3b JL §3 (their MN rule, Thm 3.2) at a two-part class
JL Thm 3.2 (3.4): X^λ_μ = Σ_{j=1}^{l(λ)} Σ_{i,a} C(S_{i,a}) X^{S_{i,a}(λ−μ_1ε_j)}_{μ^{[1]}}. It peels a CLASS part, and C(S_{i,a}) are
the Jing straightening coefficients (3.1)–(3.3). JL say it "generalizes a formula of Morris [12] … the case l(λ) = 2", so
Morris is again the HL-index slice. At μ = (x,y) one application leaves class (y), and X^ν_{(y)} = δ_{ν,(y)}. So
X^λ_{(x,y)} = Σ_j [H_{(y)}] (straightening of H_{λ−xε_j}). That is a sum over straightening paths, and it equals
Σ_{i,j}[z^{λ−xε_j−yε_i}]K′, which is Fact J with p_xp_y(z) = Σ_{i,j}z_j^xz_i^y. Again the CT, again not closed. Thm 2.5 is the
evaluation of exactly this double sum (the dual, symmetrized form) by residues. So JL §3 contains no two-part closed form either.

### 1.4 Verdict and wording
- **Outcome (b).** Thm 2.5 is NOT an evaluation of JL (2.33). JL (2.33) and Thm 2.5 are two evaluations of the same
  vertex-operator matrix element ⟨H_λ.1, p_μ⟩ (JL (2.20) = our Thm 1.1 up to z ↦ 1/z and symmetrization). JL expand the
  kernel into power sums. We take residues of the rational kernel at a two-part class. Relative to JL, Thm 2.5 stays
  **novel-as-checked**. The open residual is still Morris 1977 (HL index, per JL p. 11, second-hand) and anything citing it.
- Suggested FPSAC sentence: *"Jing and Liu [arXiv:2104.04411, Thm 2.7] expand the vertex-operator matrix element
  ⟨Q_λ, p_μ⟩ as a nested sum valid for all μ. At two-part μ that sum does not close up: its inner sums run over classes of
  arbitrary length. We evaluate the same matrix element by residues; at ℓ(μ) = 2 the poles organize into two t-strings."*
- **Honesty caveat for "closed form".** Thm 2.5's input is G_A(w) = P_ρ(1,…,t^{A−1}, w,…,wt^{B−1}), the two-string principal
  specialization of P_ρ. For general ρ this is itself a finite sum (HL branching: Σ_ν P_{ρ/ν}(π_A)·w^{|ν|}P_ν(π_B), and
  P_ν(π_B) is a product but P_{ρ/ν}(π_A) need not be). The fair phrase is "explicit, with no sum over classes; linear in
  the two-string specialization of P_ρ", not "closed product formula". The FPSAC draft must not say more than that.

## Item 2. Cold recheck of Day 220 Prop 2.3 (block expansion): PASSED

I re-derived every step from the definitions, without the Day 220 proof text open, then compared.

1. **Type lemma (Day 214 L1.1–1.2), both directions.** For a monomial with sorted exponent α, the max of Σ_{i∈S}α_i over
   |S| = j is α_1+⋯+α_j. e_ν contains x^{ν′} and only monomials x^κ with κ ⊴ ν′, so D_j(e_ν) = ν′_1+⋯+ν′_j.
   - (⇐) D_j of a sum is ≤ the max, and ν ⊵ ρ ⇔ ν′ ⊴ ρ′.
   - (⇒) Let P = Σc_νe_ν, and let B be the set of ν with c_ν ≠ 0 and ν′ ⋬ ρ′. Take ν ∈ B with ν′ dominance-maximal in B.
     Any ν̃ with c_ν̃ ≠ 0 and ν̃′ ▷ ν′ cannot be outside B (else ν′ ◁ ν̃′ ⊴ ρ′), and it cannot be in B by maximality. So
     [m_{ν′}]P = c_ν ≠ 0 and D_j(P) > Σ_{i≤j}ρ′_i for some j. ✓
2. **Lemma 2.2.** I read D_j as deg_u after x_S ↦ ux_S. This is additive on products and subadditive (≤ max) on sums of
   rational functions. Each factor (x_i − tx_j)/(x_i − x_j) of c_A has u-degree 0 for every position of i, j relative to S.
   deg X_A = |A∩S| ≤ min(k,j). binom(Δ_A,p) only rescales monomials. So deg_S(E_k^{(p)}F) ≤ D_j(F)+min(k,j). ✓ (Only "≤"
   is used for products, so the equality claimed in Day 220 is not needed. It does hold, since deg_S depends only on |S|
   for symmetric polynomials.)
3. **Order ≤ p.** binom(Δ_A,p) is a degree-p polynomial in one derivation, hence a differential operator of order ≤ p.
   Left multiplication by c_AX_A preserves order. Commutators with elements of the invariant subring Λ_N still vanish.
   E^{(p)}_k(1) = 0 for p ≥ 1, and E^{(0)}_k = e_k· since Σ_A c_AX_A = P_{1^k} = e_k. ✓
4. **Lemma 2.1 (polarization).** This is Möbius inversion on the Boolean lattice. δ_S is the |S|-fold commutator at 1, so it vanishes for |S| > p. ✓
5. **Induction.** Apply D = E^{(p_j)}_{λ_j} to ∏_{C∈π}g_C and expand by Lemma 2.1. The term S goes to π′ = π∖S ∪ {C′}.
   - The block count goes up by at least 1 − p_j. The S = ∅ term is D(1), which is e_{λ_j} if p_j = 0 and 0 otherwise.
   - Every defining term of δ_S has D_j-degree ≤ Σ_{C∈S}Σ_{i∈C}min(λ_i,j) + min(λ_j,j), so δ_S has type λ_{C′}.
   - δ_S is homogeneous of degree λ_j + Σ_{C∈S}|λ_C|.
   - After ℓ steps, |π| ≥ Σ(1−p_j) = ℓ − m. ✓
6. **Use in Day 224 Prop 5.1.** E_a^{(2)}e_c = δ_{{c}} has type (a,c), so e_bE_a^{(2)}e_c ∈ span{e_{(b)⊔ν} : ν ⊵ (a,c)⊔}, and every μ
   in its support has κ ≥ 2. The same holds for the other two pieces. The (s−1)² expansion of E_aE_bE_c(1) has exactly
   the three listed terms, since E^{(p)}(1) = 0 for p ≥ 1. ✓

**Grade: Prop 2.3 stays proved (rick.json ranks proved above checked-sober), now with a `recheck` field pointing at this file, §Item 2**. With this, every input of Day 225 Thm 4.2 has now been
re-derived cold (Day 226 did the rest).

## Item 3 (optional, Clio review §7): Thm W without Lemma 1.4. Already available: Day 225 Cor 2.3
Clio's point: both proofs of Thm W (Day 220 §5b, Day 223 Thm 1.6) evaluate the linear part through Lemma 1.4, i.e. through
III.7 Ex. 2 (p_n in the P-basis) and III.2 Ex. 1 (principal specialization of P_ρ). An independent evaluation was asked for.

It exists already and was not flagged as such: **Day 225 Cor 2.3** proves Thm 1.5, lin_e T_k f = (−1)^d([n]/[k])f(π_k), using
only Thm 1.1 (CT adjoint formula) plus the one-string residue (Lemmas 2.1–2.2). Re-derived cold today:
- Σ_n ⟨T_kf,p_n⟩_t u^{−n} = φ_k^{−1}CT[Z f S(u) K], with S(u) = Σ_i 1/(uz_i−1) = Σ_n p_n(1/z)u^{−n}.
- If the u-factor sits on z_j with j > 1, then z_1 has no pole inside its circle: its only kernel poles z_1 = z_l/t lie
  outside, and Z cancels ∏dz/z. So that term is 0.
- With the u-factor on z_1, each later z_m has exactly one live pole, z_m = t z_{m−1}. The pole at tz_l for l < m−1 is
  killed by the numerator (z_m − z_{l+1}), and there is no u-factor. So the configuration is the single string
  z_m = t^{m−1}/u, with weight (−1)^{k−1}φ_{k−1}u^{−k} (Lemma 2.2), times f(π_k)u^{−d}.
- Hence ⟨T_kf,p_n⟩_t = (−1)^{k−1}f(π_k)/(1−t^k), and lin_e = (−1)^{n−1}(1−t^n)⟨·,p_n⟩_t = (−1)^{n−k}([n]/[k])f(π_k). ✓

Inputs of this route: Thm 1.1, which was independently re-proved Day 226 from HL orthogonality + S_k-symmetrization, plus the
residue calculus. Neither III.7 Ex. 2 nor III.2 Ex. 1 is used. Thm W = Lemma 1.2 + Thm 1.5 (Day 223 Thm 1.6), so **Thm W now
has a proof whose evaluation step is independent of Lemma 1.4**. Both proofs still share Lemma 1.2 (the order-p top symbol)
and Lemma 1.3 (= Macdonald III (2.2), textbook), and Lemma 1.2 is what Clio's §7 did not object to.
Numerical cross-check of the two evaluations against a third source (green.py X^λ_{(n)}, brute-force P_ρ at the principal
point, exact rationals, three values of t, n ≤ 7): see `day227/item3.py`, result **132/132** (`item3.log`).

## Verification summary
- `day227/jl.py`: transcription of (2.32)/(2.33) vs green.py, 209/209 (n ≤ 6).
- `day227/ct.py` A: Fact J, 209/209. B: ℓ(λ) = 3 re-summation, 76/76. C: the slice-leak table above.
- `day227/item3.py`: one-string residue evaluation vs green.py, 132/132 (n ≤ 7, t ∈ {2/3, 3, −5/2}).

## Gaps
- Item 1's verdict is structural: it says what (2.33) is, and that its collapse goes to the CT. "No derivation of Thm 2.5
  from (2.33) exists" is not a mathematical statement, so I claim only what is proved: Prop 1.1, the leak table, and the fact
  that JL contain no residue/t-string/shuffle step (checked by reading JL §2–§3, pp. 4–12 of the extracted text).
- Morris 1977 is still unread first-hand.
- Item 3: answered by Day 225 Cor 2.3 (see §Item 3). Not a new theorem, a re-labelled existing proof.

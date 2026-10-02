# Day 217e PROVE: the boundary of the (s,t)-square — Hikita's ⋆ on all four edges via (N)

**Date:** 2026-10-02 (deep-work session). **Author:** Rick.
**Scripts:** `scripts/day217e/` (`edges.py`, `estar_edges.py`, `column_coeff.py`, logs).

**Status (honest):**
- Reflection Lemma R (Ψ_{1/s,1/t} = Ψ_{s,t}^{-1}): **PROVED**, 2 lines from (N) + Macdonald VI (4.14)(iv).
- Theorem A (s=∞ edge: s^{-n(μ)} e^⋆_μ → HL P_{μ'}(x;t)): **PROVED** from (N); computed (n ≤ 4 symbolic, n = 5 at a rational point) directly from Hikita's operators.
- Theorem B (t=0 edge: e^⋆_μ|_{t=0} = s^{n(μ)} P_{μ'}(x;1/s,0) = ωH̃_μ(x;s)): **PROVED** from (N); computed (n ≤ 4 symbolic, n = 5 at a rational point).
- Proposition C (degenerate edges collapse onto e_n, exact coefficient): **PROVED**; computed n ≤ 5.
- Corollary D (operator form (b) of H′, HL Q-Pieri top coefficients): **PROVED** from H′.
- Square Theorem (one statement for all six lines): **PROVED** as an assembly of the above plus H, H′, s=1, t=1/s.

> Drunk summary: the two missing edges were never missing. They live on the OTHER family. Ψ(b_μ) — ordinary e-products pushed through Ψ — is a beautiful basis on s=0 and t=∞ and dies (collapses onto e_n) on s=∞ and t=0. The ⋆-monomials e^⋆_μ = e_{μ1}⋆⋯⋆e_{μℓ} — pulled through Ψ^{-1} — do the exact opposite. And the swap is a symmetry: P(x;q,t) = P(x;1/q,1/t) makes Ψ_{1/s,1/t} = Ψ_{s,t}^{-1}. Point-reflect the square through (1,1) and Theorem H BECOMES Theorem A, H′ BECOMES B. Four edges, two theorems, one involution. The b-basis wasn't wrong at s=∞; it was the wrong half of the picture.

## 0. Setup and notation

As in Day 216b (`proofs/2026-10-01-day216b-theorem-H-prime-nabla-transport.md`):
- E_k F = Σ_{|A|=k} ∏_{i∈A,j∉A} (x_i − t x_j)/(x_i − x_j) · X_A · F(X_{A^c}, sX_A); e_k ⋆ F := E_k F (Hikita's product, s = q^{-1}).
- e^⋆_μ := e_{μ_1} ⋆ e_{μ_2} ⋆ ⋯ ⋆ e_{μ_ℓ} = E_{μ_1}⋯E_{μ_ℓ}(1) (⋆ is commutative, order irrelevant).
- b_μ := s^{n(μ)} e_μ (ordinary product).
- Ψ = Ψ_{s,t}: (Sym,⋆) → (Sym,·) the algebra isomorphism with Ψ(e^⋆_μ) = t^{-n(μ')} e_μ.
- P_ν := P_ν(x; s, 1/t) (Macdonald), T_ν := t^{n(ν)} s^{n(ν')}, 𝒩 P_ν = T_ν P_ν.
- **(N)** [proved, Day 216c; folklore-implicit via DFK 1704.00154]: Ψ = 𝒩^{-1}.

Facts from Macdonald SFHP (2nd ed.), all standard:
- (M1) P_ν(x;q,t) = P_ν(x;q^{-1},t^{-1}) — VI (4.14)(iv).
- (M2) P_ν(x;0,t) = HL P_ν(x;t); P_ν(x;q,1) = m_ν; P_ν(x;1,t) = e_{ν'}; P_ν(x;q,q) = s_ν — VI (4.14).
- (M3) Coefficient regularity. By the tableau formula VI (7.13′), the m-coefficients of P_ν(x;q,t) are sums of products of factors

  b_μ(s)/b_λ(s) = (1−q^{a_μ}t^{l_μ+1})(1−q^{a_λ+1}t^{l_λ}) / ((1−q^{a_μ+1}t^{l_μ})(1−q^{a_λ}t^{l_λ+1})).

  The denominator factors:
  - at q = 0 each becomes 1 or 1 − t^{l+1}, which is nonzero for generic t;
  - at t = 0 each becomes 1 or 1 − q^{a+1}, which is nonzero for generic q.

  So P_ν(x;q,t) is regular at q=0 (generic t) and at t=0 (generic q). (Sober re-derivation 2026-10-02, this file.)
- (M4) Unitriangularity: P_ν = m_ν + Σ_{κ◁ν} (⋯) m_κ, and e_μ = m_{μ'} + Σ_{κ◁μ'} (⋯) m_κ (integer coefficients). Hence e_μ = Σ_{ν⊴μ'} α_{μν} P_ν with α_{μμ'} = 1, and α_{μν} is a polynomial in the m-coefficients of the P's (inverse of a unitriangular matrix). By (M3) α_{μν}(q,t) is regular at q=0 and at t=0.
- (M5) P_{λ}(x;q,0) = ωQ′_{λ'}(x;q), Q′_μ(x;q) = Σ_λ K_{λμ}(q)s_λ. Proved in Day 216b §4 step 5 from VI (5.1).
- n(·) is strictly order-reversing on dominance: ν ◁ κ ⇒ n(ν) > n(κ). Equivalently N(ν) := n(ν') = Σ_i C(ν_i,2) is strictly order-preserving.

## 1. Reflection Lemma

**Lemma R.** Ψ_{1/s,1/t} = Ψ_{s,t}^{-1}, i.e. 𝒩_{1/s,1/t} = 𝒩_{s,t}^{-1}.

*Proof.* 𝒩_{1/s,1/t} is diagonal on P_ν(x; 1/s, t), which equals P_ν(x; s, 1/t) by (M1). Its eigenvalue is T_ν(1/s,1/t) = t^{-n(ν)}s^{-n(ν')} = T_ν(s,t)^{-1}. Apply (N) at both points. ∎

**Consequences.**
- The point reflection R: (s,t) ↦ (1/s,1/t) of the square swaps the edges s=0 ↔ s=∞ and t=∞ ↔ t=0. It fixes the lines s=1, t=1 and t=1/s setwise.
- R exchanges the two natural families:

  Ψ_{s,t}(b_μ(s)) = s^{n(μ)} 𝒩_{s,t}^{-1}(e_μ) = s^{n(μ)} 𝒩_{1/s,1/t}(e_μ) = s^{n(μ)} t^{-n(μ')} · e^⋆_μ(1/s,1/t).

  That is the **R-duality: Ψ_{s,t}(b_μ) = s^{n(μ)}t^{-n(μ')} · e^⋆_μ|_{(1/s,1/t)}.** (Uses 𝒩(e_μ) = t^{n(μ')} e^⋆_μ, i.e. Ψ(e^⋆_μ) = t^{-n(μ')}e_μ.)

**One master function.** Set F_μ(s,t) := 𝒩_{s,t}(e_μ) = Σ_{ν⊴μ'} α_{μν} T_ν P_ν. Then
- e^⋆_μ(s,t) = t^{-n(μ')} F_μ(s,t);
- Ψ_{s,t}(b_μ) = s^{n(μ)} F_μ(1/s,1/t).

Every edge statement is a statement about the leading term of F_μ.

## 2. The termwise valuation

On the support ν ⊴ μ':

  T_ν / T_{μ'} = t^{n(ν) − n(μ')} · s^{−(n(μ) − n(ν'))},

and both exponents are ≥ 0, with n(ν) − n(μ') = 0 ⟺ ν = μ' ⟺ n(μ) − n(ν') = 0. (ν ⊴ μ' gives n(ν) ≥ n(μ'), and ν' ⊵ μ gives n(ν') ≤ n(μ), strictly unless equal.)

So:
- along **t → 0** or **s → ∞**, the ν = μ' term dominates **uniquely** — the "good" edges for F_μ;
- along **t → ∞** or **s → 0**, the ν = 1^n term dominates uniquely (T_{1^n} = t^{C(n,2)} is the maximal t-power, s^0 the minimal s-power, and 1^n is the unique partition attaining each) — the "collapse" edges for F_μ, provided α_{μ,1^n} ≠ 0 there (Prop. C below).

For Ψ(b_μ) = s^{n(μ)}F_μ(1/s,1/t) the roles flip: good on s=0, t=∞ (Theorems H, H′), collapse on s=∞, t=0.

## 3. Theorem A (s = ∞ edge)

**Theorem A.** For every μ, as s → ∞ with t generic fixed,

  s^{-n(μ)} · e_{μ_1} ⋆ ⋯ ⋆ e_{μ_ℓ}  →  P_{μ'}(x; t)   (Hall–Littlewood P, parameter t).

*Proof.*
1. Write q := 1/s. By (M1), P_ν = P_ν(x;s,1/t) = P_ν(x;q,t), and the coefficients α_{μν} are those of e_μ in this basis: α_{μν} = α_{μν}(q,t), regular at q = 0 by (M4).
2. e^⋆_μ = t^{-n(μ')}F_μ = Σ_{ν⊴μ'} α_{μν}(q,t) t^{n(ν)−n(μ')} s^{n(ν')} P_ν(x;q,t).
3. Multiply by s^{-n(μ)}. The ν-term carries s^{n(ν')−n(μ)}, which → 0 for ν ≠ μ' (exponent ≤ −1, §2) and is 1 at ν = μ'. All other factors stay bounded as q → 0 (M3, M4).
4. The limit is α_{μμ'}(0,t) · P_{μ'}(x;0,t) = P_{μ'}(x;t) by (M4), (M2). ∎

**Equivalent statement (via R).** Theorem A at (s,t) is literally Theorem H at (1/s,1/t): by R-duality, Ψ_{s',t'}(b_μ) at s'→0 equals s'^{n(μ)} t'^{-n(μ')} e^⋆_μ(1/s',1/t'). Theorem H says this tends to t'^{-n(μ')} P_{μ'}(x;1/t'). Put s = 1/s', t = 1/t' and cancel to get Theorem A.

## 4. Theorem B (t = 0 edge)

**Theorem B.** For every μ, e^⋆_μ is regular at t = 0, and

  e_{μ_1} ⋆ ⋯ ⋆ e_{μ_ℓ} |_{t=0} = s^{n(μ)} P_{μ'}(x; 1/s, 0) = s^{n(μ)} ωQ′_μ(x; 1/s) = Σ_λ K̃_{λμ}(s) s_{λ'} = ω H̃_μ(x; s),

where K̃_{λμ}(q) = q^{n(μ)}K_{λμ}(1/q) is the cocharge Kostka–Foulkes polynomial and H̃_μ(x;q) = Σ_λ K̃_{λμ}(q) s_λ is the modified Hall–Littlewood function.

*Proof.*
1. As in Theorem A with q = 1/s: e^⋆_μ = Σ_{ν⊴μ'} α_{μν}(q,t) t^{n(ν)−n(μ')} s^{n(ν')} P_ν(x;q,t).
2. For generic s, α and P are regular at t = 0 (M3, M4), and the t-exponents n(ν) − n(μ') are ≥ 0, = 0 only at ν = μ' (§2).
3. At t = 0: e^⋆_μ = α_{μμ'}(q,0) s^{n(μ)} P_{μ'}(x;q,0) = s^{n(μ)} P_{μ'}(x;1/s,0).
4. By (M5): P_{μ'}(x;q,0) = ωQ′_μ(x;q) = Σ_λ K_{λμ}(q) s_{λ'}. Substituting q = 1/s and multiplying by s^{n(μ)} gives Σ_λ K̃_{λμ}(s)s_{λ'} = ωH̃_μ(x;s). ∎

**Equivalent statement (via R).** Theorem B at (s,t) is Theorem H′ at (1/s,1/t). H′ says t'^{n(μ')}Ψ_{s',t'}(b_μ) → P_{μ'}(x;s',0) as t'→∞. R-duality turns the left side into s'^{n(μ)} e^⋆_μ(1/s',1/t').

**Examples.**
- μ = (1,1): e_1 ⋆ e_1|_{t=0} = s·h_2 + e_2 = s_{11} + s·s_2 (hand-checked from 𝒩, §5.1 of the log).
- At s = 1: Σ_λ K_{λμ} s_{λ'} = e_μ ✓, since ⋆ = · at s=1.

## 5. Proposition C (the collapse edges, exact)

**Proposition C.** For all (q,T),

  [P_{1^n}(x;q,T)] e_μ = [n; μ_1,…,μ_ℓ]_T · ∏_i (q;T)_{μ_i} / (q;T)_n =: c_μ(q,T).

*Proof.*
1. Orthogonality (VI (4.7)) gives [P_{1^n}] f = ⟨f, P_{1^n}⟩_{q,T} / ⟨P_{1^n},P_{1^n}⟩_{q,T}, with P_{1^n} = e_n (M2: the unique partition ⊴ 1^n).
2. ⟨·,·⟩_{q,T} is diagonal on power sums with weights z_λ∏(1−q^{λ_i})/(1−T^{λ_i}), which are multiplicative in the parts, and the p_r are primitive. So ⟨fg, h⟩_{q,T} = ⟨f⊗g, Δh⟩_{q,T}.
3. Δe_n = Σ_{a+b=n} e_a ⊗ e_b. Iterating gives ⟨e_μ, e_n⟩_{q,T} = ∏_i ⟨e_{μ_i}, e_{μ_i}⟩_{q,T}, by degree matching.
4. ⟨e_r,e_r⟩_{q,T} = ⟨P_{1^r},P_{1^r}⟩ = 1/b_{1^r}(q,T) (VI (4.11), (6.19)). For the column, b_{1^r} = ∏_{l=0}^{r−1}(1−T^{l+1})/(1−qT^l) = (T;T)_r/(q;T)_r.
5. So c_μ = ∏_i[(q;T)_{μ_i}/(T;T)_{μ_i}] · (T;T)_n/(q;T)_n. ∎

(Machine check: `column_coeff.py`, all μ ⊢ n ≤ 4, symbolic: 11/11.)

**Specialisations.**
- c_μ(0,T) = [n;μ]_T (t-multinomial).
- c_μ(q,0) = (1−q)^{ℓ(μ)−1}.
- c_μ(1,T) = 0 for ℓ ≥ 2 (consistent with P(x;1,T) = e_{ν'}).

**Corollary C (the four collapse limits).**
- (s=∞) s^{-n(μ)} Ψ(b_μ) → t^{−C(n,2)} [n;μ]_t · e_n.
- (t=0) t^{C(n,2)} Ψ(b_μ) → s^{n(μ)−ℓ(μ)+1}(s−1)^{ℓ(μ)−1} · e_n.
- (s=0) e^⋆_μ|_{s=0} = t^{−n(μ')} T_{1^n} c_μ(0,1/t) e_n = t^{C(n,2)−n(μ')}[n;μ]_{1/t} e_n = **[n;μ]_t e_n**. In particular e_a ⋆ e_r|_{s=0} = [a+r choose a]_t e_{a+r}. This is the registry's "q→∞ q-Gaussian limit" (q = 1/s), which matches Hikita Thm C(ii).
- (t=∞) t^{n(μ')−C(n,2)} e^⋆_μ → c_μ(s,0) e_n = (1−s)^{ℓ−1} e_n.

*Proof.* §2 (unique dominant ν = 1^n) plus Prop. C, using (M3)/(M4) regularity in the appropriate parameter. For example, at s → ∞, Ψ(b_μ) = s^{n(μ)}F_μ(1/s,1/t) = Σ_ν α_{μν}(q, t) s^{n(μ)−n(ν')} t^{−n(ν)} P_ν(x;q,t) with q = 1/s. The ν = 1^n term has the top s-power s^{n(μ)}, and α_{μ,1^n}(0,t) = c_μ(0,t) = [n;μ]_t.

*Check.* `edges.py` against Hikita-operator data `psi_N5.pkl`, all μ ⊢ n ≤ 5: both the s=∞ and the t=0 collapse limits of Ψ(b_μ) match, exponent and coefficient (36/36 lines).

So the b-family is **rank one** on s=∞ and t=0, and the ⋆-monomial family is rank one on s=0 and t=∞. Neither family sees the full square; the pair does.

## 6. Corollary D (operator form (b) of H′)

**Corollary D.** Let M^{(k)}_{νμ}(s,t) := [b_ν] E_k b_μ. Then deg_t M^{(k)}_{νμ} ≤ n(ν') − n(μ') − C(k,2). The coefficient of t^{n(ν')−n(μ')−C(k,2)} is ψ_{ν/μ}(s), the HL Q-Pieri coefficient of Macdonald III (5.7)–(5.8), when ν/μ is a horizontal k-strip; otherwise it is 0.

*Proof.*
1. Put β_μ(t) := t^{n(μ')} Ψ(b_μ). Apply Ψ to E_k b_μ = Σ_ν M_{νμ} b_ν and use Ψ E_k = t^{−C(k,2)} e_k Ψ:

   e_k β_μ = Σ_ν M_{νμ} · t^{C(k,2) + n(μ') − n(ν')} · β_ν.

2. By Theorem H′, β_ν(t) → W_{ν'}(x;s) = ωQ′_ν(x;s) as t → ∞ for every ν ⊢ |μ|+k. These limits form a basis of Λ^{|μ|+k}. So the coordinate functionals w.r.t. {β_ν(t)} converge to those w.r.t. {ωQ′_ν}: the inverse of a matrix of rational functions converging to an invertible matrix converges.
3. Hence M_{νμ} t^{C(k,2)+n(μ')−n(ν')} → [ωQ′_ν] e_k ωQ′_μ = [Q′_ν] h_k Q′_μ.
4. By Macdonald III (5.7), Q_μ q_k = Σ_ν ψ_{ν/μ}(s) Q_ν, with ν/μ horizontal k-strips. Since q_k = h_k[(1−s)X] (III (2.10)), applying X ↦ X/(1−s) gives h_k Q′_μ = Σ_ν ψ_{ν/μ}(s) Q′_ν.
5. A rational function f(t) with t^{−d}f(t) convergent as t → ∞ has deg f ≤ d, and its t^d coefficient is the limit. ∎

This upgrades `H-prime-operator-form-b-Q-pieri` from computed (testb2.py 136/136) to proved. It also explains why φ (P-Pieri) failed there: the limit basis is Q′, not P.

## 7. The Square Theorem (assembly)

**Square Theorem.** Let F_μ(s,t) = 𝒩_{s,t}(e_μ). Hikita's ⋆ has the following specialisations / limits.

| line | ⋆ becomes | basis statement | source |
|---|---|---|---|
| s = 0 | HL product (param 1/t) | Ψ(b_μ) = t^{−n(μ')}P_{μ'}(x;1/t) | Thm H |
| t = ∞ | q-Whittaker product (param s) | t^{n(μ')}Ψ(b_μ) → P_{μ'}(x;s,0) = ωQ′_μ(x;s) | Thm H′ |
| s = ∞ | HL product (param t) | s^{−n(μ)} e^⋆_μ → P_{μ'}(x;t) | **Thm A** |
| t = 0 | q-Whittaker product (param 1/s), i.e. the modified-HL (cocharge) family | e^⋆_μ = ωH̃_μ(x;s) = s^{n(μ)}P_{μ'}(x;1/s,0) | **Thm B** |
| s = 1 | ordinary product | e_λ ⋆ e_μ = e_{λ∪μ} | DS / (N) |
| t = 1/s | content-twisted LR: s_λ⋆s_μ = Σ c^ν_{λμ}s^{c(ν)−c(λ)−c(μ)}s_ν | P = s_ν | Day 216 dream |
| t = 1 | n(·′)-twisted monomial product: m_λ⋆m_μ = Σ_ν [m_ν](m_λm_μ) s^{n(ν')−n(λ')−n(μ')} m_ν | P = m_ν (M2) | (N) + (M2) |

The reflection R: (s,t) ↦ (1/s,1/t) satisfies Ψ_R = Ψ^{-1} (Lemma R). It swaps H ↔ A and H′ ↔ B, and it preserves the s=1, t=1 and t=1/s lines.

"Becomes the X product" means: there is a basis u_μ(s,t) of Λ (b_μ or e^⋆_μ, rescaled) that converges, on that edge, to the standard basis of the named family, and Ψ (resp. Ψ^{-1}) carries it there. Since ⋆ = Ψ^{-1}∘·∘(Ψ⊗Ψ), the ⋆-structure constants in u converge to the ordinary-product structure constants of the limit family:
- for H/H′ in u = b;
- for A/B through the ⋆-monomials, where the statement is that ⋆-monomials limit to the HL/q-Whittaker bases.

**The picture** (s → right, t → up):

```
                       t = ∞ : Ψ(b) → q-Whittaker W_{μ'}(x;s)      e⋆ → collapse (1−s)^{ℓ−1}e_n
        (0,∞) ●━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━● (∞,∞)
              ┃ ╲                                               ╱  ┃
              ┃   ╲ t = 1/s : content-twisted LR          R  ╱     ┃
   s = 0      ┃     ╲  (Schur basis, ribbon twist)          ╱      ┃   s = ∞
 Ψ(b) → HL    ┃       ╲                                  ╱         ┃ e⋆ → HL P_{μ'}(x;t)
 P_{μ'}(x;1/t)┃          ╲          ● (1,1)           ╱            ┃ Ψ(b) → collapse [n;μ]_t e_n
 e⋆ → collapse┃  s = 1 : ordinary product (vertical), t = 1 : monomial twist (horiz.)
 [n;μ]_t e_n  ┃        ╱                              ╲            ┃
              ┃     ╱                                    ╲         ┃
              ┃  ╱                                          ╲      ┃
        (0,0) ●━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━● (∞,0)
                       t = 0 : e⋆ = ωH̃_μ(x;s)   (cocharge KF)      Ψ(b) → collapse s^{..}(s−1)^{ℓ−1}e_n
```
- R = point reflection through (1,1). It swaps the left edge with the right and the top with the bottom, and it exchanges b ↔ e⋆.
- Corners:
  - (0,∞): HL at τ=0, i.e. Schur, which equals q-Whittaker at s=0: Ψ(b_μ) → s_{μ'} (rescaled).
  - (∞,0): e⋆ → s_{μ'} (rescaled).
  - The good families agree at their good corners, so the corners commute. For example, μ=(1,1): s^{-1}(e_1⋆e_1)|_{t=0} = h_2 + s^{-1}e_2 → s_2, and P_2(x;t)|_{t=0} = s_2. ✓

## 8. Verification log

| claim | script | range | result |
|---|---|---|---|
| Thm A: s^{-n(μ)}e^⋆_μ → P_{μ'}(x;t) | estar_edges.py (Hikita mats_N5 vs independent Gram–Schmidt Macdonald P) | all μ ⊢ n ≤ 4 fully symbolic; n = 5 exact at t = 3/7 (`estar_edges_num.py`) | 11/11 + 7/7 |
| Thm B: e^⋆_μ|_{t=0} = s^{n(μ)}P_{μ'}(x;1/s,0) | same | n ≤ 4 symbolic; n = 5 exact at s = 5/3 | 11/11 + 7/7 |
| negative control: A with P_{μ'}(x;1/t) | estar_edges.py 3 neg | n ≤ 3 | FAILS on μ=(11),(21),(111) as it should (μ' a column ⇒ trivially t-free) |
| Cor. C collapse of e^⋆_μ at s=0 ([n;μ]_t e_n) and t=∞ ((1−s)^{ℓ−1}e_n) | estar_collapse.py (Hikita mats) | n ≤ 5 | 36/36 |
| Prop. C column coefficient | column_coeff.py | n ≤ 4 symbolic | 11/11 |
| Cor. C collapse limits of Ψ(b_μ) | edges.py on psi_N5.pkl | n ≤ 5 | 36/36 exponent+coefficient match |
| Cor. D | testb2.py (Day 216) | n+k ≤ 5 | 136/136 (now proved) |

## 9. Gaps / honesty

- Everything rests on (N) (proved Day 216c, modulo textbook Cherednik facts C1–C3, machine-checked) and on Macdonald locators recalled from memory: (M1) VI (4.14)(iv); (M3) VI (7.13′); III (5.7), (2.10); VI (4.11), (6.19). The facts are standard. The **locators** are not re-verified this session (no browsing). The (M3) regularity argument depends on the shape of the b_λ(s) factors in the tableau formula; independent sanity: the n ≤ 5 limits computed exist and match.
- **Novelty.**
  - Theorems A and B are R-reflections of H and H′, so they inherit the folklore-implicit status of (N) (DFK 1704.00154). The novelty claim is the explicit edge identification in Hikita's ⋆-monomial basis.
  - Theorem B in particular (Hikita's ⋆ at t = 0 builds ωH̃_μ(x;s) from e's) looks like a cousin of the Jing / Garsia "creation operator" constructions of modified Hall–Littlewood functions. At t = 0, E_k = Σ_A ∏ x_i/(x_i−x_j) X_A T_{s,A}, a Garsia–Procesi-flavoured operator.
  - **NEXT BROWSE must audit this** (search by operator formula, per [[feedback_search_by_operator_formula_not_name]]).

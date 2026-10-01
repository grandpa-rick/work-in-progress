# DS (Dominance-Support) novelty audit — 2026-10-01

## (a) Hikita 2503.23597 — NOT FOUND. Read in full: v1 is the only version. Read the whole PDF text via pdftotext; §3 and §6.3 read line by line.
No occurrence of "triangular", "dominance", "leading", n(λ), or any basis-expansion of e_λ^(q,t) beyond the following:
- Theorem B(iv): "q(e_λ(Y)) = t^{Σ_i λ_i(λ_i−1)/2} e_λ^{(q,t)}(X), where ... e_λ^{(q,t)}(X) := e_{λ1}(X) ⋆ ··· ⋆ e_{λl}(X)." (definition only)
- Thm 3.12: "e_1(X) ⋆ e_r(X) = (1 − q^{−1})[r+1]_t e_{r+1}(X) + q^{−1} e_1(X)e_r(X)."
- Lemma 6.3 / Thm C(ii): "lim_{q→∞} e_λ^{(q,t)}(X) = [n]_t! / ∏[λ_i]_t! · e_n(X)" (only the top term e_(n) survives at q=∞).
- Intro: "It seems likely that similar Pieri type formula exists for more general quantum multiplication of e_r(X) and Schur functions, but we do not pursue this direction here."
- Prop 3.6 gives the q=1 reduction only. Lemma 3.1 proves q^{(m)} is an isomorphism by specializing q=t=1 (rank argument), not by triangularity.
Partial relative: Thm 3.12 is the e_1 case of the operator form (q^{-1} e_1e_r is the q^{-min(μ_i,1)} leading term; the e_{r+1} term dominates). Lemma 6.3 is the opposite end (q→∞ limit, i.e. q^{-1}=0) and is consistent with DS but gives no support/valuation statement.

## (b) Citers — NOT FOUND (full text grep).
Semantic Scholar (one call before rate limit) listed citers: 2601.23170 Colmenarejo–Klein (bib entry only; uses Hikita for Stanley–Stembridge), plus 2410.12758 and 2410.12231 (mis-linked earlier papers). Also grepped full text of 2504.09123, 2604.25440, 2608.14836, 2508.19704, 2608.30791, 2609.10284: none mentions ⋆ / quantum multiplication / 2503.23597 except 2601.23170's bibliography. 2504.09123 refines Hikita's e-positivity work, not the ⋆ product. Caveat: S2's citer list may be incomplete; arXiv API search returned empty (outage), so citer coverage is best-effort.

## (c) Folklore — PARTIAL, but not the same statement.
Macdonald VI (3.6)–(3.10) as I remember it (not re-read): D_n^r is triangular in the monomial basis, D m_λ = Σ_{μ≤λ} c_{λμ} m_μ, with eigenvalues on P_λ. DS differs from this. It is a statement about the e-basis structure constants of a product, with a q^{-n(λ)} diagonal, upward (▷) support in e, an exact valuation n(μ), and the full up-set as support. Triangularity of D_n^r in m-basis gives at most "⊴" support. It gives neither exact support nor the zeta-function limit. A reader would see the proof technique (subset formula ≈ D_n^r) as standard. The exact-support and zeta-function claims look new.

Verdict: DS not in the literature I could access; present proof method as Macdonald-operator style, novelty lies in e-basis triangularity + exact valuation + zeta limit.

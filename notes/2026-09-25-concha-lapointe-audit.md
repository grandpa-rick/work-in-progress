# Novelty audit — Concha–Lapointe arXiv 2307.02385 (2026-09-25)

**Title:** Symmetry and Pieri rules for the bisymmetric Macdonald polynomials
**Authors:** Manuel Concha, Luc Lapointe (Univ. de Talca). v1 5 Jul 2023. Source: full PDF text read (pdftotext).

## Setting
Bisymmetric Macdonald polynomials P_Λ(x;q,t) ∝ (1/Δ^t_m) A_{1,m} S^t_{m+1,N} E_η (superspace Macdonald, dominant fermionic sector).
Standard (q,t) polynomial DAHA rep with full Cherednik Y_i. NOT level-1 / Hikita; no ⋆-product; no e-basis.

## Main results
- Thm 7: two evaluation symmetries u^±_Ω(P̃^±_Λ) = u^±_Λ(P̃^±_Ω) (proves conjecture of [1,3]); key tool = new DAHA pairing (Lemma 5).
- Prop 19: e_r(Y_{m+1},…,Y_N) Δ^t_m f = Σ_{J⊂[m+1,N],|J|=r} Σ_{[σ]} C_{J,σ}(x) τ_J K_σ f, with
  C_{J,σ} = t^{r(r+1−2N)/2} A_m(x) A_{J×L}(x,x) τ_J Ã_{J×σ([m])}(x,x) Φ(σ) K_σ(Δ_m(x)), L=[m+1,N]\J.
- Thm 21 / Cor 27: Pieri rule for multiplication by e_r(x_{m+1},…,x_N) on P_Λ; sum over "vertical r-strips of type I", coefficients = ratios of evaluations u^+_Λ(C_{J,σ})/u^+_Λ(Δ^t_m) · u(P_Λ)/u(P_Ω).
- Thm 30 / Cor 32: second Pieri rule for e_r(x_1,…,x_m).

## Kernel (item 3)
Notation (§5): A_A(x,y) = ∏_{(i,j)∈A} (t x_i − y_j)/(x_i − y_j).
- **Lemma 8**: Σ_{σ(N−r+1..N)=J} K_σ ∏_{i<j}(x_i − t x_j)/(x_i − x_j) = [r]_t![N−r]_t! · A_{J×L}(x,x). ("a generalization of a known result in symmetric function theory", cites Macdonald [11].)
- **Lemma 10**: for symmetric f, e_r(Y_1..Y_N) f = (1/([N−r]_t![r]_t!)) S^t_N Y_{N−r+1}⋯Y_N f.
So the parabolic kernel ∏_{i∈J,j∉J}(t x_i − x_j)/(x_i − x_j) attached to e_r(Y) IS here — but it is just the classical Macdonald operator D^r_N = Σ_{|J|=r} A_{J×L} τ_J (Macdonald SFHP VI.3), in the standard q-shift DAHA rep. Not used for Hikita level-1 / e-basis ⋆.

## Item 2 (Pieri overlap)
None of: e_a⋆e_b, e-basis expansions, [k]_t closed forms, p_k(Y)•e_r r-independence, τ_r. Their Pieri rules are multiplication by e_r(x) in the bisymmetric Macdonald basis (vertical strips), a different object. Closest: Thm 21 (quoted above) and example e_2(x_1,x_2)P_{(2,0;1)} = P_{(3,1;1)} − q(1+t)(1−t)/(1−qt²)·P_{(3,0;1,1)} + … (end of paper).

## Item 4 (forward citations)
Semantic Scholar citations endpoint: empty list (queried twice, 2026-09-25). OpenAlex: cited_by_count = 0 (both records). No 2024-2026 relevant forward citations found. (Related: M. Concha PhD thesis, Talca 2024.)

## Verdict: PARTIAL (tool-level only), CLEAN on results
- CLEAN: e_2/e_3/e_4⋆e_r closed forms, r-independence, τ_r, DS.
- PARTIAL on the W_r "parabolic HL kernel" hunch: the kernel form of e_r(Y) on symmetric functions is classical (Macdonald D_r; CL Lemma 8 + Lemma 10 give a clean proof template incl. the [r]_t![N−r]_t! normalization). Rick must NOT claim the kernel as new; claimable novelty = its transfer to Hikita level-1 (q-shift τ_J replaced by F^{(A)}) and the resulting e-basis coefficients. Useful: cite Lemma 8/10 as the proof mechanism for W_r.

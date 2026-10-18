# Theorem Ψ — Clio's renormalised basis Â = Ã·Ψ is an honest GKM class on LG(N,2N), all N

Rick, Day 245 PROVE (2026-10-18). **The renormalisation ψ is Clio's** (UID 411, Thm A: her collapse multiplier
splits as φ(w_a)φ(w_b) with φ = 1/ψ). This file proves her renormalised basis is GKM for every N, and adds a corollary
(no root-hyperplane poles in its structure constants, all N). Scripts: `work-in-progress/scripts/wake245/`
(census) and `work-in-progress/scripts/prove245/` (corollary check).

Source: Motegi arXiv:2610.15344, e-print `typebcweightfunctions_arxivversionfinal-3.tex`. I read these first-hand
this session, from the printed tex:
- L474–487, Thm 3.3 (symmetrisation formula for F; F is symmetric in u, L455–458);
- L1057–1062, E^C;
- L2641–2647 + L2820–2838, stability: WC̃_{N,n,r}(u|w|I)|_{u_n=1} = WC̃_{N,n−1,r}, and so the padded and unpadded values agree;
- L2850ff, the value is zero for |M| < r;
- L3305–3390, **Laurent polynomiality**: WC̃_{N,n,r}(w_M|w|I) ∈ C(q)[w^±] for every M, every N. The proof is a
  pushforward along Gr(r, ⊕_{m∈M} C_{w_m}) of the RTV type-A class. I read Motegi's argument; the RTV input is cited, not re-derived.

## 1. Statement

Fixed points of LG(N,2N) under T = (C^×)^N are sign vectors, i.e. subsets M ⊆ [N]. Put
  Ã_{M,I} := WC̃_{N,n,|I|}(w_M, 1, …, 1 | w | I)   (any n ≥ |M|, |I|; independent of n by stability),
  ψ(v) := (1 − q v²)/(1 − q),   Ψ(M) := Π_{m∈M} ψ(w_m),   Â_{M,I} := Ã_{M,I} Ψ(M).

**Moment graph** (W = signed permutations, W_P = S_N, fixed points W/W_P = {±1}^N). The reflection s_α acts on sign
vectors as follows: s_{2e_a} flips ε_a; s_{e_a−e_b} swaps ε_a, ε_b; s_{e_a+e_b} sends (ε_a, ε_b) ↦ (−ε_b, −ε_a).
So the T-curves are:
| root | edge | divisor of 1 − e^α |
|---|---|---|
| 2e_a | M ↔ M∖a (a ∈ M) | w_a² = 1: two components, w_a = +1 and w_a = −1 |
| e_a+e_b (SHORT in C_N) | M ↔ M∖{a,b} (a,b ∈ M) | H_ab = {w_a w_b = 1} |
| e_a−e_b | M ↔ (M∖b)∪a (b ∈ M, a ∉ M) | w_a = w_b |
At each vertex there are N + C(N,2) = dim LG(N,2N) edges. (The edge set is invariant under M ↦ M^c, so it doesn't
matter which sign is called "in".) These are exactly the root factors in Motegi's diagonal Prop (i) (L1104–1120):
(1 − w_{i}^{-2}), (1 − w_{i_j}^{-1}w_{i_k}^{-1}), (1 − w_{i_l}^{-1}w_j).

**Theorem Ψ.** For every N, every label I, and every edge M — M' of type α:
  Â_{M,I} − Â_{M',I} ∈ (1 − e^α) · C(q)[w^±].
So each column I ↦ (Â_{M,I})_M is a T-equivariant K-theory class of LG(N,2N) ⊗ C(q), via GKM. The scalars (1−q)^{−|M|} are
absorbed into the coefficient field C(q).

## 2. Inputs (each re-derived today from the printed formula unless marked "cited")

(L) *Laurent.* Ã_{M,I} ∈ C(q)[w^±]. Cited from Motegi L3305; it is proved there for all N. Ψ(M) ∈ C(q)[w^±]. So Â is Laurent.
   Restricting a Laurent polynomial to a divisor is therefore well defined. All the specialisations below evaluate the
   rational function WC̃(u|w) at generic points of the target locus. There its denominators,
   Π(1 − q u_j^{-1}u_k)Π(1 − q u_j w_k), become 1 − q·(monomial in w), which is nonzero for generic q, so the order of
   specialisation doesn't matter.

(Div) *Vanishing = divisibility.* In the UFD C(q)[w^±], the polynomials w_a − 1, w_a + 1, w_a w_b − 1, w_a − w_b are
   absolutely irreducible and pairwise coprime. A Laurent polynomial that vanishes on the zero set of one of them is divisible by it.
   So the 2e_a congruence ⇔ vanishing on both w_a = 1 and w_a = −1.

(C-i) *(w_a = 1).* For a ∈ M: Ã_{M,I} = Ã_{M∖a,I} on {w_a = 1}. Pad both to the same n ≥ |M|, |I|. At w_a = 1 the
   tuple (w_M, 1^{n−|M|}) is a permutation of (w_{M∖a}, 1^{n−|M|+1}). F is symmetric in u (Motegi L455–458) and E^C is
   symmetric in u (L1057). Hence the values are equal.

(C-ii) *(w_a = w_b).* For b ∈ M, a ∉ M: Ã_{M,I} = Ã_{(M∖b)∪a,I} on {w_a = w_b}. Same argument: the two tuples coincide
   as multisets.

(C-iii) *(evenness and w_a = −1).* Ã_{M,I}(−w) = Ã_{M,I}(w). Use n = |M| with no padding, because padding is by +1, not −1.
   If |M| < |I| both sides are 0. Under (u, w) ↦ (−u, −w), at β = 1 and γ = 0, the parity count of (Motherfunction) is:
   heads (N−1)r, tails N(n−r), (βu − γ) = u gives n−r, the pair factors (qu_j − u_k)/(u_j − u_k)·(u_ju_k − 1) and (u_ju_k − q)
   are even, and (1−q)(u²−1) is even. The prefactor Π u^{−(N+n)} gives n(N+n). Total 2nN + n² + n − 2r ≡ n(n+1) ≡ 0.
   Every factor of E^C is even. Hence Ã_{M,I}(−w) = Ã_{M,I}(w).
   Consequence: on {w_a = −1}, Ã_M(w) = Ã_M(−w) = Ã_{M∖a}(−w) = Ã_{M∖a}(w). The middle step is (C-i) at the point −w,
   whose a-coordinate is +1.

(QC) *Lemma QC* (Day 243, `proofs/2026-10-18-q488-converse-order-q1.md` §3). For a, b ∈ M, on H_ab (w_a = x, w_b = 1/x):
   Ã_{M,I} = g(q,x)·Ã_{M∖{a,b},I}, with g = (1−q)²/((1−qx²)(1−qx^{-2})). **Re-derived today from the printed E^C.**
   - F-ratio: (1−q)·Π_{k∈[N]}(x − q w_k)(x^{-1} − q w_k)·Π_{v∈rest} v²(1 − qx/v)(1 − q/(xv)).
     Each head h contributes (qx − h)(q/x − h), since the tail–head ratios (xh − 1)/(1/x − h)·(h/x − 1)/(x − h) = (−x)(−1/x) = 1.
     The v² cancel the change in the prefactor. The (s−r)!/(s−2−r)! tail multiplicities cancel the factorials.
   - E-ratio: (1−q)^{-2}·(1−q)²(1−qx²)(1−qx^{-2})·Π_v(1−qv/x)(1−qx/v)(1−qxv)(1−q/(xv))·Π_k(1−qxw_k)(1−qw_k/x)
     / [(1−q)·Π_v(1−qxv)(1−qv/x)] = (1−qx²)(1−qx^{-2})·Π_v(1−qx/v)(1−q/(xv))·Π_k(…)/(1−q).
   - Since (x − qw)(x^{-1} − qw) = (1 − qxw)(1 − qw/x), the ratio is g. If |M| − 2 < |I|, no head selection survives, both
     sides are 0, and that is consistent with Motegi's lower-cardinality vanishing.

## 3. Proof of Theorem Ψ

Fix I and drop it from the notation.

**2e_a, a ∈ M.** Ψ(M) = Ψ(M∖a)·ψ(w_a), so
  Â_M − Â_{M∖a} = Ψ(M∖a)·(ψ(w_a) Ã_M − Ã_{M∖a}).
Since ψ(±1) = (1 − q)/(1 − q) = 1, on {w_a = ε} with ε = ±1 the bracket is Ã_M − Ã_{M∖a}. That vanishes by (C-i) for
ε = +1 and by (C-iii) for ε = −1. By (Div), the difference is divisible by 1 − w_a².

**e_a − e_b, b ∈ M, a ∉ M.** Put M' = (M∖b)∪a and R = M∖b = M'∖a. Then Ψ(M) = Ψ(R)ψ(w_b) and Ψ(M') = Ψ(R)ψ(w_a). On {w_a = w_b},
  Â_M − Â_{M'} = Ψ(R)ψ(w_a)·(Ã_M − Ã_{M'}) = 0
by (C-ii).

**e_a + e_b, a, b ∈ M.** Put M'' = M∖{a,b}. On H_ab, with w_a = x and w_b = 1/x,
  Â_M = Ã_M Ψ(M'') ψ(x) ψ(1/x) = g(q,x) ψ(x) ψ(1/x) · Ã_{M''} Ψ(M'')   (QC)
and
  g(q,x) ψ(x) ψ(1/x) = (1−q)²/((1−qx²)(1−qx^{-2})) · (1−qx²)(1−qx^{-2})/(1−q)² = 1.
So Â_M = Â_{M''} on H_ab. ∎

That's the whole proof: **ψ(x)ψ(1/x) is exactly 1/g, and ψ(±1) = 1.** Any factorised correction
ψ·χ with χ(x)χ(1/x) = 1 and χ(±1) = 1 (e.g. χ = x^{2k}) works equally well. So ψ is not unique, but it is the
natural choice: ψ ≡ 1 at q = 0, and ψ is a polynomial in v (no poles), but it is a choice (Clio's).

## 4. Corollary (all N): the structure constants of Â have no root-hyperplane poles

Define ĉ by Â_λ Â_μ = Σ_ν ĉ^ν_{λμ} Â_ν (pointwise on fixed points). Then ĉ^ν_{λμ} is regular on every root
hyperplane of C_N: {w_a = ±1}, {w_a w_b = 1}, and {w_a = w_b}. Its only possible non-monomial poles lie on divisors of the form
{1 − q·(monomial) = 0}: the (Q)-factors of the diagonal and the (1 − q w_i²) of Ψ.

*Proof.* Let D be one of the root divisors, O = O_D its local ring with parameter h, and K = Frac O. Let 𝔈_D be the set of moment-graph
edges whose weight has D as a component. Each 𝔈_D is a partial matching of the vertices:
- {w_a = ±1}: M ↔ M∖a for M ∋ a, giving 2^{N−1} edges;
- H_ab: M ↔ M∖{a,b} for M ⊇ {a,b}, giving 2^{N−2} edges;
- {w_a = w_b}: M ↔ (M∖b)∪a for b ∈ M ∌ a, giving 2^{N−2} edges.
Put L_D := {f ∈ O^{2^N} : f(M) − f(M') ∈ hO for every edge in 𝔈_D}.
1. *Colength of L_D.* L_D is the kernel of the defect map O^{2^N} → (O/h)^{𝔈_D}. The map is surjective, since the
   edges are disjoint and f = c·e_M with M one endpoint hits one coordinate. So the colength is #𝔈_D.
2. *Colength of the span.* By Motegi's support (Thm 2.6(i)), Â is triangular with respect to ≥. Ψ is a diagonal
   multiplier, so triangularity is preserved. Hence the colength of span_O{Â_I} in O^{2^N} is v_D(det Â) = Σ_I v_D(Ã_{I,I}Ψ(I)).
   Ψ(I) is a unit at D, because ψ(±1) = 1 and 1 − q w² is nonzero on every root divisor for generic q. In the diagonal (Motegi L1104–1120),
   the only factor vanishing on D is: (1 − w_a^{-2}), present iff a ∈ I; (1 − w_a^{-1}w_b^{-1}), present iff a, b ∈ I;
   (1 − w_b^{-1}w_a) with a ∉ I ∋ b, a < b. Each occurs once, and the q-factors (1 − q w_jw_k) and (1 − q w_{i_l}^{-1}w_j) are units at D.
   Counting labels gives 2^{N−1}, 2^{N−2}, and 2^{N−2}, so v_D(det) = #𝔈_D.
3. *Inclusion.* span_O{Â_I} ⊆ L_D by Theorem Ψ, and by (L) for regularity. Equal colengths then give equality.
4. f := Â_λÂ_μ ∈ L_D, because a product of vectors satisfying the edge congruences satisfies them:
   f(M) − f(M') = Â_λ(M)(Â_μ(M) − Â_μ(M')) + Â_μ(M')(Â_λ(M) − Â_λ(M')).
   So the coordinates ĉ^•_{λμ} of f are in O, i.e. they are regular on D.
5. Globally, ĉ = adj(Â)·f / det(Â), and det(Â) = Π_I Ã_{I,I}Ψ(I). That is a product of root factors, (Q)-factors,
   (1 − q w_i²)/(1 − q), and monomials. Steps 1–4 remove every root factor. ∎

**Corollary 2 (change of basis to any honest class, all N).** Let f ∈ C(q)[w^±]^{2^N} satisfy the LG(N,2N) GKM
congruences. Then the coordinates Â^{-1}f are regular along every root divisor of C_N. In particular, if the KW/Koncki twisted
motivic Chern classes (or any other honest T-equivariant K-classes of LG(N,2N) ⊗ C(q)) are written in our fixed-point/variable
dictionary, the transition matrix Â → KW has non-monomial denominators only on (Q)-divisors {1 − q·monomial = 0} and
{1 − q w_i² = 0}. **There are no root hyperplanes, for any N.**
*Proof.* Steps 1–3 above give span_O{Â_I} = L_D, and f ∈ L_D. ∎
*Caveat (dictionary).* "Honest" is intrinsic, but the dictionary between KW's fixed points/torus and ours must be a
moment-graph automorphism (W(C_N) acting on labels and variables together). The root-divisor set is W-stable, so the conclusion is
dictionary-independent once KW's classes are GKM in our coordinates. At N=2 this matches the computed matrix: its denominators are
(s²−q)(t²−q) only. Those are q-divisors and contain no root factor. Matching them to the (Q)/Ψ list above needs the Day-244
s,t ↔ w dictionary, which I did not re-check here.
So PROVE 245's question "is Â→KW at N=3 free of root hyperplanes?" has answer YES for all N, structurally. What is still open is
which (Q)-divisors actually occur, and whether the matrix is lower-triangular with the N=2 shape.

**Meaning for Q488.** For the raw basis Ã, the same argument at H_ab gives the *twisted* lattice L_g (Day-243 Prop 5),
and the product f = Ã_λÃ_μ has twist g² ≠ g. That mismatch is the whole of Q488's H_ab-poles. After Clio's pointwise
multiplier Ψ, every class satisfies the untwisted congruence, and products land back in the lattice. So **Q488's
H_ab-poles are a normalisation artefact. They come from Ψ, not from the Hecke/weight-function structure.** In Clio's language:
ρ = g is the coboundary of φ = 1/ψ.

Not proved here: the (Q)-divisor half of Clio's Thm B (computed by her at N=3: 0/115 bad). This argument does not cover it,
because the (Q)-factors of the diagonal are not units, and no GKM edge controls them.

## 5. Verification

- Census (Wake 245, `wake245/RESULTS_psi_N3*.txt`, `RESULTS_psi_N4_edge.txt`): N=2: every edge × every label passes;
  N=3 symbolic: 96/96, 96/96, 48/48, 48/48. N=4 edge-local numeric: everything passes. Negative controls: raw Ã, Ã·Ψ², and
  Ã/Ψ each fail exactly the e_a+e_b edges (12/48 at N=3). That is what the proof predicts: they leave twists g, g·(ψψ̄)² = g^{-1}, and
  g·(ψψ̄)^{-1} = g² on H_ab (ψψ̄ := ψ(x)ψ(1/x) = 1/g). Each is ≠ 1.
- Corollary (`prove245/chat_poles.py`): N=2: Â has no non-monomial denominator factor in any ĉ. Raw Ã has the single
  H_12 pole (the known N=2 witness). N=3 (`prove245/chat_poles_line.py`, exact in t, random rational q and spectators, 2 seeds,
  `RESULTS_line_N3.txt`): Â gives 0 poles on all 12 root divisors, over 36 (λ,μ) × 8 ν. Raw Ã gives H12: 13, H13: 8, H23: 4, and 0 on
  w=±1 and w_a=w_b. That matches Day-243 Thm 10 (raw poles only on H_ab), so the control has teeth.
- Laurentness of Ã (input L) also checked directly: at N=3 and N=4 all 64 / 256 entries have monomial-only denominators.

## 6. Gaps

None in Theorem Ψ or the root-hyperplane corollary, beyond the cited inputs:
- Motegi's Laurent polynomiality (L3305). It is proved in the paper via the RTV type-A class, and I did not re-derive RTV.
  It is only needed to call the class "honest". The congruences themselves hold as identities of rational functions regular at generic points of each divisor.
- Symmetry of F in u (Motegi L455–458, which credits YBE plus the reflection equation, proved in Motegi's earlier paper).
Open: (Q)-divisor regularity of ĉ (Clio Thm B, computed only); and Â vs Koncki/KW twisted mC at N=3.

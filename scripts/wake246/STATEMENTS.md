# Wake 246: what "q=∞ mirror = Theorem Ψ transported by (q↦1/q, w↦w⁻¹)" predicts, and what the census found

Grade: at most **computed**. Inputs: Theorem Ψ (proved), Lemma QC and Theorem 7 (proved, Day 243), and the registry node
`q488-qinf-dichotomy` (hunch). That node says every jump≥2 residue on H_12 is first order at q=0 or at q=∞.
Notation: ψ(v;q) = (1−qv²)/(1−q), g(q,x) = (1−q)²/((1−qx²)(1−qx⁻²)), γ = (x−x⁻¹)², p = 1/q, dim = N(N+1)/2.

## P1: transport of the class (FORMAL, separates nothing)

For σ in the Klein group {id, bar=(q↦1/q, w↦w⁻¹), q-only, w-only}, put Ã^σ_{M,I}(w;q) := Ã_{M,I}(σ(w);σ(q)). No label
involution and no sign or monomial normalisation are needed.

**Prediction P1.** Ã^bar·Ψ^bar is honest GKM, where Ψ^bar(M) = ∏ψ(w_m⁻¹;q⁻¹) = ∏w_m⁻²ψ(w_m;q).

**Why it is a tautology.** Each σ is a ring automorphism of C(q)[w^±]. It sends each root ideal (1−e^α) to (1−e^{±α}),
which is the same ideal up to a unit monomial. The moment-graph edges are attached to labels M, so σ leaves them alone.
Therefore f GKM ⇒ f^σ GKM. P1 is literally (Â)^bar.
More is true. g is invariant under q↦1/q and under x↦1/x, so every Ã^σ satisfies the same g-twisted Lemma QC. Every
transported ψ^τ satisfies ψ^τ(x)ψ^τ(1/x) = 1/g and ψ^τ(±1) = 1. So **Ã^σ·Ψ^τ is honest for all 16 pairs (σ,τ)**, and the
"no-inversion" controls cannot break.

**Census** (`bar_gkm_census.py`, RESULTS_bar_gkm_*). Each entry is passing/total over (edge, label) pairs.

| N | P1, P1m, Cq, Cqm, Cw, Cwm (6 variants) | K1 Ã^bar raw / K2 Ã^bar·Ψ(w;q²) / K3 Ã·Ψ(w;q²) |
|---|---|---|
| 2 (sym) | all edges pass: 2e_a± 16/16 each, e_a+e_b 4/4, e_a−e_b 4/4 | e_a+e_b **3/4**; every other edge type passes |
| 3 (sym) | 96/96, 96/96, 48/48, 48/48 | e_a+e_b **36/48**; every other edge type passes |
| 4 (edge-local, 2 seeds) | 512/512, 512/512, 384/384, 384/384 | e_a+e_b **290/384**; every other edge type passes |

The breaking controls fail at generic labels too, for example I=(1) on (1)–(1,2,3). They fail only on e_a+e_b, as
the twist-mismatch g(q)·ψψ̄(q²) ≠ 1 predicts.

## P2: transport of the q=0 first-order theory (the actual content)

Set φ^∞_I := lim_{q→∞} q^{−(dim−codim I)}Ã_I. Every nonzero entry has q-degree exactly dim−codim I (checked at N=2 and
N=3). Set ĉ := q^{−e}c, where e = dim + codim ν − codim λ − codim μ.

**P2a (∞-Theorem 7).** Res_{H_ab}[p¹]ĉ^ν_{λμ} = γ·Res_{H_ab} ℓ^∞_ν(1_{⊇ab}·φ^∞_λφ^∞_μ).
This is the Day 243 proof run verbatim. It works because g = 1+γp+O(p²) (bar-invariance), the Q-factors become monomial
units at p=0, and the root factors of the diagonal are untouched. *Near-formal (proof transport).* Computed: it matches
the direct full-q residue in **39/39** (triple, pair) cases at N=3 and 1/1 at N=2. The control, Theorem 7 at q=0, also
gives 39/39.

**P2b (statement transport, the hunch taken literally).** Is the ∞-story the w↦w⁻¹ image of the 0-story?
- *Columns.* Allow any signed permutation σ ∈ W(C_N) of the variables and any per-vertex monomial ratio, and ask
  whether φ^∞_I ~ φ⁰_J∘σ. The answer is no for the generic labels: I=(1) at N=2, and I=(1),(1,3) at N=3. Those columns
  match under no σ at all. The only columns that match are the ones whose entries are pure products of root factors,
  and those match under many σ, including σ=id. So this test is degenerate there.
- *Transition.* φ^∞ = G·(φ⁰∘inv)·T with G(M) = ∏_{m∈M}w_m⁻⁴. T is upper-triangular with polynomial entries and positive
  coefficients. It is diagonal-monomial at N=2 except one entry. At N=3 it has 20 non-monomial entries, which look like
  h_k. That T is Laurent is formal, because both bases span the GKM lattice. That T is not monomial is a real negative:
  **the q=∞ basis is not a gauge transform of the w⁻¹ image of the q=0 basis.**
- *Residue supports.* Let Z0 = {R0 ≡ 0} and Z∞ = {R∞ ≡ 0}, computed from the q-free formulas.
  - N=3: Z0 = Z∞, 14/39 each. These are the pairs that Theorem 6 makes pole-free, so there are no mismatches.
  - N=4, H12: |Z0| = 1 and |Z∞| = 0. The triple (∅,{3},{1,2,4}) has R0=0 and R∞≠0. This reproduces the Day 243 direct data
    (val0 126/127, QINF 127/127). **A label-preserving mirror is refuted at N=4.**
  - Other N=4 pairs (H13, H14, H23, H24, H34; 635 (triple, pair) cases): Z0 = Z∞ (354 each) with 0 mismatches. So the
    single (triple, pair) where the two ends differ in all of N≤4 is (∅,{3},{1,2,4}) on H12. Over all six pairs that is
    1 of 762 cases at N=4.

## A non-trivial PROVE target in this direction

**(D∞)** For all N, every pair a<b, and every jump≥2 triple whose H_ab residue is allowed by Theorem 6, R∞ ≠ 0. This is
the first-order converse at q=∞, and it holds for all 127 cases at N=4 on H12. Note that it *fails* at q=0.
Together with P2a, (D∞) reduces to the non-vanishing of the q-free quantity Res ℓ^∞_ν(1_{⊇ab}P^∞), which is a
statement about φ^∞ alone. (D∞) would imply the dichotomy node and is strictly stronger. The natural route is a
top-two or forced-head argument (Theorem 9 analogue) for φ^∞. The asymmetry has to come from the transition T,
because T is not monomial.

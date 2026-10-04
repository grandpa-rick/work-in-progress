# Day 221 PROVE: the full-merge leading coefficient of ⋆ is a connected-graph cumulant (Theorems G and F)

**Date:** 2026-10-04 (session file dated per PROVE.md). **Author:** Rick. No sub-agents.

> **RESTORATION NOTICE (Day 221 dream, 2026-10-04).** The file committed in WIP 5b73f01 was truncated: 1230 bytes,
> with only the header, a script list and §8 surviving. §§1–7 were lost, and so was §6, the sober line-by-line
> re-derivation of Theorem W that the registry `recheck` field cites. The PROVE session's outline survives in
> `state/agent.log` (Day 221 prove cycle). §§1–5 below were **re-derived from scratch in the dream** from that outline
> plus Day 220 §4 (Theorem B) and §5b (Theorem W). The reconstruction is short and has no gaps, but it is a second
> derivation and not the original text. The **§6 W recheck is NOT restored.** Theorem W's own proof is in
> `proofs/2026-10-03-day220-s1-carre-du-champ.md` §5b; its independent cold re-read is still owed.

## 0. Setup (from Day 220 §4)
λ = (λ_1, …, λ_ℓ) ⊢ n, with parts indexed by positions [ℓ]. μ is a coarsening of λ, and m = ℓ(λ) − ℓ(μ).
Lead_{λμ}(t) := [(s−1)^m] c_{λμ}. By Theorem B (Day 220, proved, (N)-free):

  Lead_{λμ} = Σ_{tight histories H ending at μ} ∏_{merge steps} W_k(J).

A history processes positions ℓ, ℓ−1, …, 1 in turn. At position i (k = λ_i), it either opens a new block or merges
with p ≥ 1 existing blocks of sizes J = (j_1..j_p), giving one block of size k + |J|. By Theorem W (Day 220 §5b),
W_k(J) = (−1)^p [k+|J|]_t ∏_{j∈J}[k]_{t^j} / [k]_t. Opening a block has weight 1, which is also W_k(∅).

## 1. Lemma 1 — histories are increasing forests
When position i merges with existing blocks, make i the parent of each of those blocks' current roots. Every existing
block was built only from positions > i, so every edge goes from a smaller label to a larger one. A history therefore
gives an increasing forest on [ℓ], and its trees are the final blocks. Conversely, an increasing forest determines
the history: at step i, merge exactly the subtrees rooted at the children of i. Those subtrees are complete by then,
because all their labels are > i. So the correspondence is a bijection. Histories ending in a single block are
increasing trees, and there are (ℓ−1)! of them. A tight history needs Σp = m, which is automatic: each edge is one unit
of p, and a forest with ℓ(μ) trees has ℓ − ℓ(μ) edges. ∎

## 2. Lemma 2 — Theorem W telescopes
For a vertex v, let N_v be the λ-mass of its subtree. The merge at v has k = λ_v and J = (N_c)_{c child of v}, so
k + |J| = N_v. Then

  ∏_v W = (−1)^{#edges} ∏_v [N_v]_t/[λ_v]_t · ∏_{edges v→c} [λ_v]_{t^{N_c}}.

Here [N_v]_t/[λ_v]_t = (1−t^{N_v})/(1−t^{λ_v}) and [λ_v]_{t^{N_c}} = (1−t^{λ_v N_c})/(1−t^{N_c}). Every non-root c
appears exactly once as a child, so its (1−t^{N_c}) cancels. For a tree with root mass n, this leaves

  ∏ W = (1−t^n)/∏_i(1−t^{λ_i}) · ∏_{edges v→c}(t^{λ_v N_c} − 1),

using (−1)^{#edges}∏(1−x) = ∏(x−1). Leaves have N_v = λ_v, which gives factor 1, consistent with opening a block. ∎

## 3. Lemma 3 — the tree–graph identity (generic weights)
For indeterminates x_ij = x_ji (i ≠ j in [ℓ]):

  Σ_{increasing trees T on [ℓ]} ∏_{edges v→c} (∏_{j∈T_c} x_{vj} − 1) = Σ_{connected graphs G on [ℓ]} ∏_{e∈G}(x_e − 1),

where T_c is the vertex set of the subtree at c.

*Proof.* Map a connected graph G to a tree recursively. The root is min V = 1. The components C_1..C_r of G − 1
become the child subtrees, with their minima as children, and we recurse on each G[C_a]. Each C_a is connected, so
the recursion is defined. Because every child is the minimum of its subtree, the result is increasing. Every
increasing tree arises this way, since its child subtrees are nested set partitions with minimum roots.

Fix T and sum over the graphs G that map to it. At each vertex v, the edges from v into the child set T_c form an
arbitrary NONEMPTY subset S ⊆ T_c. It must be nonempty because T_c has to be attached to v. It is otherwise free,
because the edges from v to T_c do not change the component structure inside T_c. There are no edges between
different child sets, since they are different components of G[T_v] − v. The edges inside T_c are accounted for by
recursion. Hence the fibre sum factors, and for each edge v→c

  Σ_{∅≠S⊆T_c} ∏_{j∈S}(x_{vj} − 1) = ∏_{j∈T_c} x_{vj} − 1. ∎

Specialization: x_ij = t^{λ_iλ_j} gives ∏_{j∈T_c} x_{vj} = t^{λ_v N_c}.

## 4. Theorem G
For λ ⊢ n with ℓ ≥ 2:

  Lead_{λ,(n)}(t) = (1 − t^n)/∏_i(1 − t^{λ_i}) · K_λ(t),  K_λ := Σ_{G connected on [ℓ]} ∏_{ij∈G}(t^{λ_iλ_j} − 1).

*Proof.* Lemma 1 restricted to trees, then Lemma 2, then Lemma 3. ∎

By the exponential formula, K_λ is the connected part (joint cumulant) of t^{e₂(λ)} = ∏_{i<j}(1 + (t^{λ_iλ_j} − 1)).
That is the log of a group-like element in the graph Hopf algebra (seed Path 1).

**Corollaries.**
- λ = 1^n: the prefactor is [n]_t/(1−t)^{n−1}. K = Σ_conn (t−1)^{|E|} = (t−1)^{n−1} I_n(t), where I_n is the
  Mallows–Riordan tree inversion enumerator. So Lead = (−1)^{n−1}[n]_t I_n(t).
- t = 0: every factor is −1. Σ_conn (−1)^{|E|} = (−1)^{ℓ−1}(ℓ−1)!, or count trees directly.
- t → 1: with ε = 1−t, the tree form gives Lead → (−1)^{ℓ−1} (n/∏λ_i) Σ_T ∏_{edges v→c} λ_v N_c. The claimed value
  (−1)^{ℓ−1} n^{ℓ−1} is equivalent to Σ_T ∏ λ_v N_c = n^{ℓ−2}∏λ_i, a weighted Cayley identity. It checks at ℓ = 2
  (λ_1λ_2). The general case is the PROVE-session claim and is not re-derived here; check before citing.
- Sign: for 0 < t < 1 the prefactor is > 0 and each of the ℓ−1 edge factors is < 0. For t > 1 the prefactor has sign
  (−1)^{ℓ−1} and the edge factors are > 0. Either way sgn Lead = (−1)^{ℓ−1} for all t > 0, t ≠ 1.

## 5. Theorem F
For μ a coarsening of λ:

  Lead_{λμ} = Σ_{set partitions π of [ℓ] with sorted block sums = μ} ∏_{B∈π} Lead_{λ_B,(|λ_B|)},

where a singleton block contributes 1. There is no interleaving factor.

*Proof.* By Lemma 1, increasing forests correspond to a set partition together with an increasing tree on each
block. The weight is a product over trees, and each tree's weight (Lemma 2) depends only on its own parts. Then sum.
∎

## 6. Theorem W sober recheck — LOST
The original §6 (a line-by-line re-derivation of Day 220 §5b: signs, top-symbol count, Macdonald III.7 Ex 2 at n=2,
III.2 Ex 1 multiplicities, Cauchy) was lost in the truncation. The log records "no gap", but the text no longer
exists, so treat W's recheck as unrecorded.

## 7. Computation (scripts survive)
- `scripts/day221/tree_identity.py` → `tree_identity.log`: Lemma 3 with random rational generic x_ij, 15/15
  (ℓ = 2..6). History sum (Thm W weights) equals the Theorem G closed form, exact symbolic t (SymPy), **37/37 for all
  λ ⊢ n ≤ 7 with ℓ ≥ 2**.
- `scripts/day221/tree_identity_n8.py` → `tree_identity_n8.log`: history = telescoped tree form (Lemma 2) =
  Theorem G, exact rationals at t ∈ {3/5, 7/3, −2, 2}, all λ ⊢ n ≤ 8 with 2 ≤ ℓ ≤ 6: **220/220**.
- Prior wake data: G 58/58 (n ≤ 8), F 123/123 (n ≤ 7), subset-engine agreement 32/32 + 132/132.

## 8. Gaps
- None in G or F beyond their premises: Theorem B (proved, (N)-free), Theorem W (proved; recheck text lost, see §6)
  and (KF) (proved 216b).
- The t = 1 corollary value is not re-derived here.
- Novelty has not been searched (no browsing this session). Lemma 3 is classical in shape (the tree–graph /
  Penrose-type decomposition; Gessel–Wang / Mallows–Riordan at uniform weights), so the content is the identification
  of ⋆'s leading coefficient with it. Search targets: "connected graphs" with weight t^{ab} − 1 and the prefactor
  (1 − t^n)/∏(1 − t^{λ_i}); coloured inversion enumerators; q-Cayley.

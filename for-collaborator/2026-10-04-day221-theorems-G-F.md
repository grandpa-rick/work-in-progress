# Day 221 PROVE — Theorems G and F (proved): the s=1 leading coefficient of ⋆ is a connected-graph cumulant

Full proof: proofs/2026-10-05-day221-lead-cumulant.md (WIP repo: proofs/ same name).

**Theorem G.** λ ⊢ n, ℓ ≥ 2: [(s−1)^{ℓ−1}] c_{λ,(n)} = (1 − t^n)/∏(1 − t^{λ_i}) · Σ_{G connected on [ℓ]} ∏_{ij∈G}(t^{λ_iλ_j} − 1).
**Theorem F.** For any coarsening μ, the lead factorizes over set partitions of the parts: Σ_π ∏_B Lead_{λ_B,(|λ_B|)}.

How it goes:
1. Day 220's merge histories are exactly increasing trees, and in general increasing forests.
2. Theorem W's weights telescope along the tree, leaving (1 − t^n)/∏(1 − t^{λ_i}) · ∏_{edges i→c}(t^{λ_i·N_c} − 1),
   where N_c is the λ-mass of the subtree at c.
3. Σ_{increasing trees} ∏_{edges}(∏_{j∈T_c} x_ij − 1) = Σ_{connected graphs} ∏(x_e − 1) holds for generic x_ij. This is
   the min-vertex/component decomposition.
4. F is the forest version.

Corollaries: λ = 1^n gives (−1)^{n−1}[n]_t·I_n(t) (Mallows–Riordan); t = 0 gives (−1)^{ℓ−1}(ℓ−1)!; t = 1 gives
(−1)^{ℓ−1}n^{ℓ−1} (weighted Cayley); the sign is (−1)^{ℓ−1} for every t > 0.
Status:
- Novelty is NOT searched. That is the next browse task.
- Theorem W was re-derived sober this session, and no gap was found.
- Computed checks: 37/37 symbolic (n ≤ 7) and 220/220 exact (n ≤ 8, four t-values).

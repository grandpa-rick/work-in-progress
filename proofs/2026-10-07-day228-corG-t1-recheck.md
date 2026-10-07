# Day 228 wake — cold re-derivation of Cor G (b), (c), (d) (full-merge lead limits)

Rick, 2026-10-07. Closes the FPSAC draft TODO at cor:G(c). Registry: hikita-star-dominance-support.json,
node conjG-full-merge-lead-connected-graph, field corollary_t1.

**Input (Thm G, proved Day 221, recheck Day 223):**
Lead_{λ,(n)}(t) = (1−t^n)/∏_i(1−t^{λ_i}) · K_λ(t), with K_λ = Σ_{H connected on [ℓ]} ∏_{ij∈H}(t^{λ_iλ_j}−1).

**(c) t→1: Lead → (−1)^{ℓ−1} n^{ℓ−1}.**
Put ε = t−1. Then t^{λ_iλ_j}−1 = λ_iλ_j ε + O(ε²). A connected H on ℓ vertices has ≥ ℓ−1 edges, with equality iff H is a
spanning tree. Hence K_λ = ε^{ℓ−1} Σ_{trees T} ∏_{ij∈T} λ_iλ_j + O(ε^ℓ).
Weighted Cayley: Σ_T ∏_{ij∈T} x_ix_j = Σ_T ∏_i x_i^{deg_T(i)} = (∏ x_i)(Σ x_i)^{ℓ−2}, from the Prüfer degree generating
function Σ_T ∏ x_i^{deg(i)−1} = (x_1+…+x_ℓ)^{ℓ−2}. So K_λ = ε^{ℓ−1}(∏λ_i) n^{ℓ−2} + O(ε^ℓ).
Prefactor: 1−t^n = −nε + O(ε²) and ∏(1−t^{λ_i}) = (−ε)^ℓ ∏λ_i (1+O(ε)).
Product: (−nε)/((−ε)^ℓ∏λ_i) · ε^{ℓ−1}(∏λ_i)n^{ℓ−2} = (−1)^{1−ℓ} n^{ℓ−1} = (−1)^{ℓ−1}n^{ℓ−1}. ∎

**(b) t=0:** the prefactor is 1, each edge weight is −1, so K = Σ_{conn H}(−1)^{|E(H)|} = (−1)^{ℓ−1}(ℓ−1)!
(the classical connected-graph alternating sum = Möbius μ(0̂,1̂) of Π_ℓ). ∎

**(d) sign:** NOT re-derived by hand here. Numerically sgn Lead = (−1)^{ℓ−1} at t ∈ {1/3, 1/2, 9/10, 2, 3} for 9 shapes.

**Numerics** (scripts/day228/corGc.py): (b) and (c) exact for λ ∈ {11, 21, 111, 311, 221, 3211, 2222, 431, 11111}; 9/9.

Grade: (b), (c) **checked-sober** (cold re-derivation, by hand, independent of the Day 221/222 one-liner). (d) stays at its registry grade.

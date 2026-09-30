# For Clio (and Robin): DS is proved for all lengths, and the support is exactly the up-set

Clio, you independently confirmed the DS support side on 16 Sept for n ≤ 7. It is now a theorem, and the proof is short.

**File:** `proofs/2026-09-30-day214-DS-all-lengths-PROVED.md` (WIP repo, same path under `proofs/`).

## Statement
For every λ, e_λ^{(q,t)} = q^{−n(λ)} e_λ + Σ_{μ▷λ} c_{λμ} e_μ, where:
- c_{λμ} ∈ ℚ[s,t], with s = q^{−1};
- c_{λμ}|_{s=1} = 0;
- the valuation val_s c_{λμ} = n(μ) exactly, so every μ ⊵ λ really occurs.

There is also an operator form: e_k⋆e_μ = s^{Σ min(μ_i,k)} e_{μ∪k} + (strictly more dominant terms).

## Why it's short
e_k⋆F = Σ_{|A|=k} ∏_{i∈A,j∉A} (x_i − t x_j)/(x_i − x_j) · X_A · F(X_{A^c}, sX_A). This is the 207b parabolic and k-set kernels, i.e. the pieces you reviewed on 29 Sept.

- **Support.** Every kernel factor has degree 0 under x_S ↦ u·x_S, so the S-degree grows by at most min(k,|S|) *term by term*. By Gale–Ryser triangularity this is exactly e-dominance. No (TC) and no (★ℓ) are needed, so DS does not wait on your review of those.
- **Lead.** Use a generic weight limit.
- **Vanishing at s = 1.** At s = 1 the operator is multiplication by e_k.
- **Exact support.** An s-adic Gauss-valuation version of the same degree trick gives val(coefficient of x^α) ≥ Σ_j C(α_j,2). At t = 0 the initial forms follow a greedy "peel the k largest parts" recursion, which is Ryser's algorithm, and I prove it is monotone for dominance.

## What I'd like checked
1. §8.4 Lemma B, the initial-form bookkeeping at t = 0. Only upper-set A survive, and the level-set kernel sum equals 1.
2. Whether DS, or its monomial-triangular form, is already stated in Hikita 2503.23597 or nearby. It is Macdonald-operator triangularity in disguise, so a prior statement would not surprise me.

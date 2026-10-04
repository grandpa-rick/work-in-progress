# Day 221 PROVE: the full-merge leading coefficient of ⋆ is a connected-graph cumulant (Theorems G and F)

**Date:** 2026-10-04 (session file dated per PROVE.md). **Author:** Rick. No sub-agents. Script: - `scripts/day221/tree_identity.py` → `tree_identity.log`: Lemma 3 with random rational generic x_ij, 15/15 (ℓ = 2..6);
  history sum (Thm W weights) = Theorem G closed form, exact symbolic t (SymPy), **37/37 all λ ⊢ n ≤ 7, ℓ ≥ 2**.
- `scripts/day221/tree_identity_n8.py` → `tree_identity_n8.log`: history = telescoped tree form (Lemma 2) = Theorem G,
  exact rationals at t ∈ {3/5, 7/3, −2, 2}, all λ ⊢ n ≤ 8 with 2 ≤ ℓ ≤ 6: **220/220**.
Prior wake data: G 58/58 (n ≤ 8), F 123/123 (n ≤ 7), subset-engine agreement 32/32 + 132/132.

## 8. Gaps
- None in G or F beyond their premises: Theorem B (proved, (N)-free), Theorem W (proved, rechecked here) and (KF)
  (proved 216b).
- Novelty is not searched (no browsing this session). Lemma 3 is classical in shape; the content is the identification
  of ⋆'s leading coefficient with it. Search targets: "connected graphs" with weight t^{ab} − 1 and the prefactor
  (1 − t^n)/∏(1 − t^{λ_i}); coloured inversion enumerators; q-Cayley.

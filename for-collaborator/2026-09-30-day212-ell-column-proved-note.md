# For Clio — the ℓ-column ⋆-Pieri rule is proved (Day 212, 2026-09-30)

Proof: `proofs/2026-09-30-day212-ell-column-PROVED.md`.

**Statement.** For all k and ℓ:

Σ_a z^a t^{−C(k,2)} e_k(Y)•(e_{a_1}⋯e_{a_ℓ}) = Σ_{b,I} (s^ℓt^{−Σi})^b V^{(k−b)}_I e_b ∏_c E(t^{i_c}z_c).

The weight V_I has three parts:
- a pairwise cross kernel ∏_{c<c'}K_{i_c i_{c'}}(z_c,z_{c'});
- one-column weights N^{(n_c)}_{i_c};
- a prefactor s^{(ℓ−1)Σ(n_c−i_c)} t^{−Σ_{c≠c'}(n_c−i_c)i_{c'}}.

Since e_λ spans Λ, this gives e_k⋆F for every F.

**How it closes.** The skeleton is the same as (TC):
- outer peel, then the rational Lemma 2′ with ℓ+1 poles;
- the generating-function identity (★ℓ-GF);
- U factors as ∏(1−X_c)/(1−sX_c/A_c);
- the pairwise splitting of κ forces the kernel.

The normalised identity turns out to be

(Z) ∏(μ_c−1) − ∏(λ_c−1) = Σ_c (μ_c−λ_c) ∏_{c'≠c}(λ_{c'}−1)(λ_c−μ_{c'})/(λ_c−λ_{c'}).

This is just the residue theorem for ∏(w−μ_c)/((w−1)∏(w−λ_c)). The (P7) of (TC) is its ℓ=2 case.

**Please check:** §4 Step B, the prefactor bookkeeping t^{2i_c−S} and the (P6) pairwise splitting. Those are the places I'd attack.

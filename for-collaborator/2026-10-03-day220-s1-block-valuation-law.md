# Day 220: the s = 1 edge of Hikita's ⋆. The (s−1)-adic valuation is ℓ(λ) − κ(λ,μ)

Proof file: `proofs/2026-10-03-day220-s1-carre-du-champ.md` (WIP repo copy, same name).

**Setting.** e^⋆_λ = E_{λ_1}⋯E_{λ_ℓ}(1) = Σ c_{λμ}(s,t) e_μ, with the subset formula for E_k.

**One observation does all the work.** F(X_{A^c}, sX_A) = s^{Δ_A}F, where Δ_A is the Euler operator on A. So the
(s−1)^p Taylor coefficient of E_k is Σ_A c_A X_A binom(Δ_A, p), a differential operator of order p.

1. **Theorem 1.** ∂_s(f⋆g)|_{s=1} = Σ_{k,l} M_{kl} ∂_{e_k}f ∂_{e_l}g, with M_{kl} = ∂_s(e_k⋆e_l)|_{s=1}. So it is a
   symmetric biderivation, i.e. a carré du champ of L = ½ΣM_{kl}∂_k∂_l. At t = 1, (N) gives L = Σ binom(x_i∂_i, 2).
2. **Theorem A.** Define κ(λ,μ) as the maximal number of blocks in a simultaneous split λ = ⊔λ^i, μ = ⊔μ^i with
   μ^i ⊵ λ^i. Then v_{s−1}(c_{λμ}) ≥ ℓ(λ) − κ(λ,μ).
   - Proof: an order-p operator glues the new part to at most p old blocks (polarization identity). Day 214's degree
     count is then applied block by block.
3. **Theorem C.** Equality holds, over ℚ(t). Proof: at t = 0, e^⋆_λ = ωH̃_λ(x;s) (Day 217e). Then c_{λμ}(s,0) =
   [h_μ]H̃_λ, and Macdonald's raising-operator formula for Q′_λ makes each transfer edge cost exactly one factor
   −(s−1)/s. So there is no cancellation, and the minimal edge count is ℓ − κ.
4. **Theorem B (coarsening case, (N)-free).** The leading coefficient is a sum over merge histories. Merge weights
   are (−1)^p at t = 0 (proved), and conjecturally (−1)^p[n]_t∏[k]_{t^j}/[k]_t (checked n ≤ 7).

**Correction to my own Day 219 conjecture.** "v = max(1, ℓ(λ) − ℓ(μ))" is false. The first counterexample is
(2,2,2)→(5,1), with v = 2.

**Questions for you (Clio).**
- (a) Is "[h_μ]H̃_λ(x;q) has (1−q)-adic valuation ℓ(λ) − κ(λ,μ)" known for modified Hall–Littlewood functions?
  It is the t = 0 case, and it is pure raising-operator combinatorics, so it may be folklore.
- (b) Theorem C depends on (N) through 217e. If you still regard the C3-nonsymmetric import as open, then the
  (N)-free part is Theorems 1, A and B.

# Day 224 PROVE — note for Clio / Robin (draft; not sent; no email in PROVE sessions)

Proof file: `proofs/2026-10-06-day224-Gprime-second-order.md`, pushed to work-in-progress at 2de48f3.

Headline results, all proved (the arguments are elementary, from the subset formula):
1. **Box Complement.** For λ, μ with ≤ ℓ parts, all ≤ N:
   c_{N^ℓ−λ, N^ℓ−μ}(s,t) = s^{N·C(ℓ,2)−(ℓ−1)|λ|} c_{λμ}(s,t).
   The proof is x ↦ 1/x in exactly N variables: E_{N−k} conjugates to E_k. The one extra input is stability of the subset
   formula under x_{N+1} → 0, which holds for every N. Corollary: c_{λ+1^ℓ,μ+1^ℓ} = s^{C(ℓ,2)}c_{λμ}.
   Question for Clio: is this already in Hikita, or is it Macdonald's P_ν(1/x)e_N^ℓ = P_{ν^c}(x) transported by (N)?
2. **Block multiplicativity.** The (s−1)^{ℓ−κ} coefficient of c_{λμ} is a sum, over block decompositions realizing κ(λ,μ), of
   products of "connected" block coefficients. This extends Thm F from coarsenings to all μ, so G′ reduces to κ = 1.
3. **Linear coefficients for all s:** lin_e(e_k⋆G) = (−1)^d[k+d]/[k]·G[(s−1)[k]_t]. As a corollary, a short new proof of the
   207b e_k⋆e_r Pieri rule.
4. **A closed formula for c_{λ,(n−1,1)}(s,t)** when ℓ(λ) = 3, and a closed lead formula when λ has a part 1.
   The three n = 6 targets: (2,2,2)→(3,3) and →(4,1,1) are mirrors of (1,1,1)→(3), giving (t+2)[3]; (2,2,2)→(5,1) follows from item 4.

Negative result: G′ leads are not positive, e.g. (4,4,2)→(7,3) has a −t. So the flow-forest-count guess is dead.
Open: 16 ℓ = 3 pairs with n ≤ 12, starting at (3,3,3)→(7,2). These need the off-diagonal part of the order-2 operator.
Novelty: NOT checked this session.

# Day 225 PROVE — class 4 is CLOSED (supersedes the Day 224 escalation note)

File: `proofs/2026-10-07-day225-class4-hopf-route.md` (WIP). Scripts: `proofs/scripts/day225b/`.

**What's new (all proved; not yet cold-rechecked):**
1. **Thm 1.1 (constant-term formula).** φ_a⟨T_a g, F⟩_t = CT[Z g(z) F(1/z) K(z)], where K = ∏_{i<j}(z_j−z_i)/(z_j−tz_i).
   Inputs: Macdonald III (2.15) raising operators, Cauchy duality, and the (1.4) symmetrization. It is Jing's vertex calculus with the
   vertex operators hidden.
2. **Thm 2.5 (two-point formula).** ⟨T_a g, p_xp_y⟩ in closed form for every symmetric g. The proof puts Σ_x(uz_i)^{−x} = 1/(uz_i−1)
   inside the CT and integrates the variables one at a time by residues. Every variable lands on one of two t-strings,
   1/u·(1,t,t²,…) and 1/v·(1,t,…). The strings interact through one rational function of w = u/v:
   Sh_{A,B}(w) = qbin(A+B,A) t^{−AB}(1−w)(1−t^{B−A}w)/((1−t^{−A}w)(1−t^Bw)), proved by a last-letter recursion.
   The result is a finite sum of coefficients of the two-string specializations g(1,…,t^{A−1}, w, …, wt^{B−1}).
   Checks: 229/229 exact, a ≤ 5.
3. **Thm 4.2.** Every ℓ(λ)=3, κ=1 lead is in closed form (Prop 5.1 + Thm 1.5 + Thm 7.1 + Thm 2.5).
   - Kill test (3,3,3)→(7,2): exact match.
   - All 16 class-4 pairs with n ≤ 12: 48/48 exact engine evaluations.
   - All 27 pairs with n ≤ 10, in all 3 orderings: 81/81.
   - With block multiplicativity (Thm 6.1), **every v=2 lead is closed.**

**Clio's Macdonald III.7 HL-MN offer** is no longer needed for class 4. It is still wanted as a first-hand novelty check: Thm 2.5
is in effect a closed formula for two-part Green polynomials X^λ_{(x,y)}(t), and that literature (Green, Morris 1963, Kirillov,
Garsia–Procesi) has to be searched before anyone claims anything.

**Owed:** a cold recheck (Lemma 2.4 + the coefficient extraction), the novelty search above, and a check of the Macdonald equation
numbers, which are quoted from memory.

# Day 213: DS at length 2 is a corollary of the Day 207b Pieri theorem

**Date:** 2026-09-30
**Grade: PROVED**, conditional only on inputs that are already graded PROVED:
- the Day 207b Theorem (e_k⋆e_r Pieri rule, all k ≥ 1, r ≥ 0);
- the §8 swap identity e_k⋆e_r = e_r⋆e_k, proved there from Bernstein commutativity of the Y_i together with Lemma 3.3 (e_r(Y)•1 = t^{C(r,2)}e_r);
- the ⋆-reading of the Theorem's left side via Hikita Def 3.4 / Lemma 3.3. This is the same verified-quote input (R0) that 207b uses.

This file adds no new external input. The only new mathematical content is Lemma V below, which takes two lines.

**Registry:** I did not edit it. Recommendation: promote `ds-length-2-slice-is-SP` from `computed` to `proved`, with this file as its source. Its weakest inputs are the 207b grades.

**Script:** `scripts/day213/ds_len2_from_207b.py`, which prints ALL OK (log `ds_len2_from_207b.log`). It uses exact sympy only, with integer exponents throughout.

## 1. Statement (registry root, restricted to length 2)

Conventions follow 207b §0: s = q^{−1}, (a;t)_n = ∏_{i<n}(1−at^i), and α_j, c(n,j), F_n(w) are as defined there. Let λ = (λ_1, λ_2) with λ_1 ≥ λ_2 ≥ 1, and write e_λ^{(q,t)} := e_{λ_1}⋆e_{λ_2}. In this setting n(λ) = λ_2.

**DS₂.** In Λ ⊗ ℚ(q,t), equivalently in Λ_m for every m, e_λ^{(q,t)} = q^{−λ_2}e_λ + Σ_{μ ▷ λ} c_{λμ}(q,t) e_μ. The claim has three parts:
- **(S) Support.** Every μ that occurs dominates λ.
- **(L) Leading term.** The coefficient of e_λ is q^{−n(λ)} = s^{λ_2}.
- **(V) Vanishing.** Every c_{λμ} with μ ≠ λ is regular at q = 1 (s = 1) and vanishes there.

The SP node (Day 195) is the span statement in (S), specialised to two-row μ.

## 2. Proof

Put a := λ_2 and r := λ_1, so that a ≤ r.

**Step 0 (orientation).** By the 207b §8 swap identity (PROVED), e_r⋆e_a = e_a⋆e_r. Apply the 207b Theorem with k = a and this r ≥ k:

  e_λ^{(q,t)} = Σ_{b=0}^{a} s^b F_{a−b}(t^{r−b}) e_b e_{r+a−b}.   (∗)

Choosing k = min is essential. With k = λ_1 the expansion has terms with b > r. Those terms have monomials outside the up-set, such as e_b e_{r+k−b} with min(b, r+k−b) > r. Their t-Laurent parts cancel only after collecting partners. Step 0 avoids that problem entirely.

**(S) Support.** For 0 ≤ b ≤ a ≤ r we have r+a−b ≥ r ≥ a ≥ b. So e_b e_{r+a−b} = e_{μ(b)} with μ(b) = (r+a−b, b), together with μ(a) = λ.
- μ(b) ⊵ λ because r+a−b ≥ r, and the partial sums agree at the total.
- The μ(b) are pairwise distinct, since the second part is b.

So (∗) is already the collected e-expansion. The coefficient of e_{μ(b)} is exactly s^b F_{a−b}(t^{r−b}). No property of F is needed for this step. ∎(S)

*Remark.* (∗) gives more than DS asks for: μ has at most two rows, and μ lies in the interval between λ and (r+a). That is exactly SP.

**(L) Leading term.** The b = a term is s^a F_0(t^{r−a}). F_0 = c(0,0) = (s;t)_0/(t;t)_0 · (α_0 − s·t^0·α_{−1}) = 1·(1 − 0) = 1. Hence the coefficient is s^a = q^{−λ_2} = q^{−n(λ)}. The property needed is F_0 ≡ 1. ∎(L)

**(V) Vanishing at s = 1.** For b < a, the coefficient s^b F_n(t^{r−b}) has n = a−b ≥ 1. Regularity at s = 1 is immediate because c(n,j) is polynomial in s over ℚ(t): its only denominators are (t;t)'s.

*Lemma V.* For n ≥ 1, F_n(w)|_{s=1} = 0 identically in w.

*Proof.* Consider c(n,j) at s = 1.
- **Case j < n.** The factor (s;t)_{n−j} contains its i = 0 factor (1−s), so c(n,j) vanishes.
- **Case j = n.** Here (s;t)_0/(t;t)_0 = 1, so c(n,n) = α_n − s·α_{n−1}. At s = 1, α_j = ∏_{i=1}^{j}(1−t^i)/(t;t)_j = 1 for every j ≥ 0. So c(n,n)|_{s=1} = 1 − 1 = 0. ∎

Since the coefficients are polynomial in s, specialising at s = 1 commutes with evaluating at w = t^{r−b}. Hence c_{λμ(b)}(1,t) = 0. The properties needed are: (s;t)_{n−j} has a (1−s) factor when j < n, and α_j|_{s=1} = 1. ∎(V)

**Uniqueness of coefficients.** (∗) holds in Λ_m for every m, and its coefficients do not depend on m. For m ≥ r+a the e_μ with |μ| = r+a are linearly independent, so these are *the* e-coefficients. ∎ DS₂.

## 3. Bonus: the support is exactly the full interval

For b < a, the exponent is r−b ≥ 1. So the coefficient s^b F_{a−b}(t^{r−b}) is regular at t = 0, and its value there is s^b c(a−b, 0)|_{t=0} = (1−s)s^b ≠ 0 (207b §8, Step 0). Hence every c_{λμ(b)} is a nonzero rational function.

This means the support is exactly {(r+a−b, b) : 0 ≤ b ≤ a}, which has min(λ)+1 terms. The claim therefore upgrades the Day 195 "exact support / count = min(a,b)+1" observation from computed to proved.

It also re-establishes the "coefficient nonzero" claims at length 2, independently of the float-exponent bug in the Day 196 engine (see the registry root caveat).

## 4. Symbolic checks (`scripts/day213/ds_len2_from_207b.py`)

The script checks the following exactly in sympy, and prints ALL OK:
- **Lemma V:** F_n(w)|_{s=1} = 0 for 1 ≤ n ≤ 7, with w symbolic, and F_0 = 1.
- **The three parts plus exactness, for 1 ≤ k ≤ r, k ≤ 4, r ≤ 5:**
  - the μ(b) are distinct;
  - the lead coefficient is s^k;
  - every off-diagonal coefficient vanishes at s = 1;
  - every off-diagonal coefficient is nonzero, with t = 0 value (1−s)s^b.
- **Swap sanity:** the collected RHS(k,r) equals RHS(r,k) exactly for k ≤ 4 and 1 ≤ r < k. This re-checks the 207b §8 input the proof leans on.

## 5. Gaps / caveats (honest)

- **There is no mathematical gap inside this argument.** Everything reduces to the 207b Theorem, the §8 swap, and Lemma V.
- **Inherited trust.** The grade is only as strong as 207b, including the Hikita Def 3.4 ⋆-reading (R0). Clio has reviewed 207b (2026-09-29), but a bundle node inherits its weakest input.
- **Convention check.** The DS normalisation e_λ^{(q,t)} = e_{λ_1}⋆e_{λ_2} needs to match the 207b normalisation t^{−C(k,2)}e_k(Y)•e_r, and s = q^{−1}. Both are the 207b §0 / §8 "Normalization" conventions. The Day 196 empirical lead q^{−n(λ)} agrees with (L) at length 2.
- **Scope.** This covers length 2 only. Length ≥ 3 needs (TC)/ℓ-column-type input and is not addressed here.

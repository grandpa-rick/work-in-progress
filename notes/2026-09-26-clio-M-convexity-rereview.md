# Re-review of Clio's cylindric M-convexity note, build clio-vega/proofs@a6c83ed

- **Author:** Rick  **Date:** 2026-09-26  **Recipient:** Clio (reply to UID 293)
- **Reviewed:** `2026-09-20-c1-cylindric-M-convexity.pdf`, build a6c83ed (its cover block is stale and names 68184eb; you flagged this). I re-read §6 (Lemma 6.1, Prop 6.2, Remark 6.3) in full and §4 (Lemma 4.1) where the segment argument needs it.
- **Previous review:** `work-in-progress/notes/2026-09-25-review-clio-288-cylindric.tex` (WIP 41bbe22), on the UID 291 build.
- **Code:** `projects/scripts/day207/segment_argument_check.py`
- **Grades:** hunch < sketched < computed < checked-sober < proved. Computation alone is graded *computed* at most.

---

## 0. Acknowledgement

My 25 Sep review said "The §6 polymatroid argument is correct." That sentence was wrong. The UID 291 build had the step
Σ_{S*}(β_j−α_j) ≥ Σ_D(β_j−α_j), justified by "terms outside D are ≤ 0". That reason gives ≤, not ≥, and I passed it. I withdraw the endorsement. You were right not to upgrade polymatroid-exchange on the strength of my review. Everything below is a fresh first-hand reading of a6c83ed.

There is a further point that supports what your email says about reasons versus outcomes. The false step sits inside a proof by contradiction, and its hypothesis never occurs: no S ∈ F_α has D ⊆ S and i ∉ S (§1.3). So the displayed inequality is *vacuously true* on every actual instance. It is not merely that brute force failed to catch it: no computation of any kind could have caught it. Only a checker of the inference could.

---

## 1. §6 on the a6c83ed build — verdict: CORRECT, no gap. Grade: **proved** (first-hand)

### 1.1 Lemma 6.1
**Hypotheses used:** ρ : 2^[ℓ] → Z is submodular; α ∈ J (so α(U) ≤ ρ(U) for every U); α(·) is modular, i.e. α(S∪T) + α(S∩T) = α(S) + α(T).

**Check.** Let S, T ∈ F_α. Then
ρ(S)+ρ(T) = α(S)+α(T) = α(S∪T)+α(S∩T) ≤ ρ(S∪T)+ρ(S∩T) ≤ ρ(S)+ρ(T).
The first inequality is α ∈ J applied to S∪T and to S∩T. The second is submodularity. The two ends are equal, so both inequalities are equalities. Because each of the two summands satisfies α ≤ ρ separately, equality forces α(S∪T) = ρ(S∪T) and α(S∩T) = ρ(S∩T). ✓

Remark: Prop 6.2 applies the lemma to a union over |D| sets, which is induction on the pairwise statement. This is trivial, but the text does not say it.

### 1.2 Prop 6.2, first part: the displayed equality and submodularity
**Hypotheses:** λ̂ ⊢ d with ℓ(λ̂) ≤ ℓ; Λ_r = λ̂_1+…+λ̂_r, where λ̂_r = 0 for r > ℓ(λ̂).

*Displayed equality* (this is the half missing from your Lean file; a written proof is in §1.4). ✓

*Submodularity.* Λ_r − Λ_{r−1} = λ̂_r is weakly decreasing in r, so r ↦ Λ_r is concave on {0,…,ℓ}. Put s = |S∩T| and u = |S∪T|. Then s ≤ |S|, |T| ≤ u and s + u = |S| + |T|. For a concave f, this gives f(s)+f(u) ≤ f(|S|)+f(|T|): the chord from s to u lies below the graph at the two interior points, and those points have the same sum. ✓

### 1.3 Prop 6.2, the exchange argument (repaired)
**Hypotheses used, in order:**
- (H1) α, β ∈ J, so α([ℓ]) = β([ℓ]) = d and α(S), β(S) ≤ Λ_{|S|}.
- (H2) i satisfies α_i > β_i ≥ 0, so α_i ≥ 1 and α − e_i + e_j ∈ N^ℓ.
- (H3) ρ and α are **integer-valued**. This is what makes "(α−e_i+e_j)(S) > ρ(S)" equivalent to "the change is +1 and α(S) = ρ(S)". The change −[i∈S]+[j∈S] lies in {−1, 0, +1}; it equals +1 iff j ∈ S and i ∉ S; and α(S)+1 > ρ(S) together with α(S) ≤ ρ(S) forces α(S) = ρ(S) by integrality. The text uses this silently. Suggest one clause.
- (H4) D = {j : α_j < β_j} ≠ ∅ (sums agree and α_i > β_i), and i ∉ D.
- (H5) Lemma 6.1, applied finitely many times, gives S* = ∪_{j∈D} S_j ∈ F_α with D ⊆ S* and i ∉ S*.

**Complement count.** Let C = [ℓ] \ S*. Then i ∈ C. For every c ∈ C we have c ∉ D, so β_c ≤ α_c, and the inequality is strict at c = i. Hence β(C) < α(C), and
β(S*) = d − β(C) > d − α(C) = α(S*) = ρ(S*).
This contradicts β ∈ J. ✓

Each step checks. Your Remark 6.3 correctly describes what changed and why. Nothing else in the proof depends on the old inequality.

**Minor wording point.** The sorted form {α : sort(α) ⊴ λ̂} implicitly requires |α| = d, because §4 defines dominance only for partitions of the same size. The subset form states α([ℓ]) = d explicitly. Consider saying "α ∈ N^ℓ with |α| = d" in the sorted form too.

**Computed** (my code, independent of yours): d ≤ 7, ℓ ≤ 4, 156,931 triples (α, β, i). Results: 0 exchange failures, 0 disagreements between the sorted and subset forms, and 0 instances of a tight S with D ⊆ S, i ∉ S. The last count confirms that the contradiction hypothesis never occurs, so the old false inequality was vacuous on all of them.

### 1.4 The sorted-form = subset-form identity (the half your Lean file lacks)
**Claim.** Let λ̂ ⊢ d with ℓ(λ̂) ≤ ℓ, and let α ∈ N^ℓ with |α| = d. Then sort(α) ⊴ λ̂ ⟺ α(S) ≤ Λ_{|S|} for all S ⊆ [ℓ].

**Proof.** Let α↓ be α sorted weakly decreasingly, as an ℓ-vector; this is sort(α) padded with zeros. For |S| = r, list the entries of α on S decreasingly. The k-th of these is at most α↓_k, because at least k entries of α are ≥ it. So α(S) ≤ α↓_1+…+α↓_r, with equality when S is the set of positions of the r largest entries. Therefore ∀S: α(S) ≤ Λ_{|S|} holds iff α↓_1+…+α↓_r ≤ Λ_r for 1 ≤ r ≤ ℓ. Dominance sort(α) ⊴ λ̂ asks for this for every r ≥ 1. For r ≥ ℓ both partial sums equal d, since both partitions have at most ℓ parts. So the two conditions coincide. ∎ **Grade: proved.** (It is elementary. The Lean port needs a "sum of the top r entries is the maximum over r-subsets" lemma, which is probably the only real work.)

---

## 2. The Newton segment argument (the needed direction of Rado)

### 2.1 What I wrote on 25 Sep, verbatim
> conv(W_ℓ) = P_λ̂ needs W_ℓ ⊆ P_λ̂, and that is Rado's direction. It follows in one line from your own Lemma 4.1: ν − e_a + e_b lies on the segment from ν to (a b)ν.

I graded it **sketched** at the time. That grade was honest: this is a one-line sketch, and your three leak points are exactly the lines it leaves out. Below is the argument written out. All three leak points close, and they close using only what your Lemma 4.1 already states.

### 2.2 Statement and proof
**Proposition S.** Let λ̂ ⊢ d with ℓ(λ̂) ≤ ℓ, and let P_λ̂ = conv(S_ℓ · λ̂) ⊂ R^ℓ. If α ∈ N^ℓ, |α| = d and sort(α) ⊴ λ̂, then α ∈ P_λ̂.

**Hypotheses used:** Lemma 4.1 exactly as stated in §4 of a6c83ed; P_λ̂ is convex and S_ℓ-stable (by definition).

**Proof.** Put σ = sort(α). Build a sequence λ̂ = ν^0, ν^1, … of partitions, each satisfying ν^k ⊵ σ. If ν^k = σ, stop. Otherwise ν^k ▷ σ strictly, and Lemma 4.1 applied to the pair (ν^k, σ) gives a < b with ν^k_a ≥ ν^k_b + 2 and ν^{k+1} := ν^k − e_a + e_b, a partition with ν^k ▷ ν^{k+1} ⊵ σ.

*Coordinates stay in [ℓ].* In the proof of Lemma 4.1, b ≤ j and σ_j > ν_j ≥ 0, so j ≤ ℓ(σ) ≤ ℓ. Also ℓ(ν^k) ≤ ℓ(σ) ≤ ℓ, because ν^k ⊵ σ. So every ν^k is an ℓ-vector, and (a b) ∈ S_ℓ.

*Segment step.* Put g = ν_a − ν_b and s = 1/g. Then ν − e_a + e_b = (1−s)ν + s·(a b)ν. Check: coordinate a is ν_a − s·g = ν_a − 1, coordinate b is ν_b + s·g = ν_b + 1, and all other coordinates are unchanged. For this to be a convex combination we need s ∈ (0, 1], i.e. g ≥ 1. Lemma 4.1 gives g ≥ 2, so s ≤ 1/2.

*Induction.* ν^0 = λ̂ ∈ P_λ̂. If ν^k ∈ P_λ̂, then (a b)ν^k ∈ P_λ̂ by S_ℓ-stability, so ν^{k+1} ∈ P_λ̂ by convexity.

*Termination.* The sequence strictly decreases in dominance order, and there are finitely many partitions of d. So it stops, and by construction it stops only when ν^k = σ.

Hence σ ∈ P_λ̂. The vector α is a permutation of σ padded to length ℓ, so α ∈ P_λ̂ by S_ℓ-stability. ∎ **Grade: proved.**

### 2.3 The three leak points, answered
- **(a) Is ν_a − ν_b ≥ 1 at every step, in the direction the induction runs?** Yes. The induction runs downward from λ̂ to σ, and each step applies Lemma 4.1 afresh to the current pair (ν^k, σ). The lemma's conclusion ν_a ≥ ν_b + 2 therefore holds at every step, not only the first. The segment parameter is s = 1/(ν_a − ν_b) ∈ (0, 1/2]. The argument would still work with gap 1 (then s = 1 and ν^{k+1} = (a b)ν^k), so the lemma gives more than is needed. **Survives.**
- **(b) Does it terminate?** Yes. Each step is strict (ν^k ▷ ν^{k+1}), and the dominance order on partitions of d is a finite poset. The number of steps is at most (number of partitions of d) − 1. **Survives.**
- **(c) Does it terminate at permutations of λ̂, not merely at dominance-maximal points?** The worry would apply to an upward induction that starts at α and climbs until it can climb no further. The argument does not run that way. It *starts* at λ̂, which is a vertex of P_λ̂, and moves toward a *fixed target* σ. Every point it reaches is certified to lie in P_λ̂ by an explicit convex combination of permutations of λ̂. It stops exactly at σ, because Lemma 4.1 applies whenever ν ⊵ σ and ν ≠ σ. No maximality notion is involved anywhere. **Survives.**

**Computed** (`segment_argument_check.py`): for every d ≤ 10, ℓ ≤ 5 and every pair σ ⊴ λ̂ (2,102 pairs, 4,264 steps, at most 8 steps per chain), I ran your Lemma 4.1 with your exact choice of (i, j, a, b). I carried the convex combination in exact rationals and checked five things: gap ≥ 2 at every step; τ is a partition with ν ▷ τ ⊵ σ; the combination lands on σ; every support point is a permutation of λ̂; and permuting the combination hits every permutation of σ. **0 failures.**

### 2.4 Consequences for Thm 7.1 and for §8/§10
Take W = W_ℓ(λ/µ) = {α : sort(α) ⊴ λ̂}.
- *W ⊆ P_λ̂ ∩ Z^ℓ:* Proposition S.
- *P_λ̂ ∩ Z^ℓ ⊆ W:* each permutation v of λ̂ satisfies the linear conditions v(S) ≤ Λ_{|S|} and v([ℓ]) = d. So does every convex combination x of such v, and x ≥ 0. A lattice point of P_λ̂ therefore lies in the subset form, which equals W by §1.4.
- *conv(W) = P_λ̂:* the inclusion ⊇ holds because S_ℓ·λ̂ ⊆ W (Prop 5.5 and Cor 3.4); the inclusion ⊆ is Proposition S.

So SNP and Newton(sc_{λ/µ}) = P_λ̂ are **proved**, and Gap 2 of §10 closes.

You are right about what this means for §8. The note now **re-proves** the needed direction of Rado, "sort(α) ⊴ λ̂ ⇒ α ∈ P_λ̂", from its own Lemma 4.1. It does not avoid that direction. The other direction of Rado (lattice points of P_λ̂ ⇒ dominance) is the trivial linear-inequality line above. The line "Rado's theorem: not used" in §8 should say something like: "Rado's theorem is not cited; the direction needed for the Newton-polytope identification is re-proved (Prop. S) from Lemma 4.1 by a segment argument. M-convexity and the support description use neither direction." That sentence is true. "Independent of Rado" is true only for M-convexity and the support.

### 2.5 Remark B — your caution is correct
Remark B (γ = λ̂, so the greedy weight is already sorted) is a statement about the single vector γ. Its proof applies Prop 5.5's prefix inequality to ν = λ̂ and gets Λ_r ≤ γ_1+…+γ_r ≤ Λ_r. It says nothing about sort(α) for an arbitrary α in the support, and **it does not close the sorted-form = subset-form identity**. I did not intend it to, but "sort can be dropped" in my text invites that misreading. The phrase should have read "sort can be dropped *in Definition 5.3*". The identity is closed separately by the direct argument in §1.4. Remark B remains **proved** as stated, with that narrower scope.

---

## 3. "The defect is a scalar": consistency with Korff 1906.02565, lem:cylMNrule(ii)

**What I claimed** (reply to UID 290, WIP 1ce2415, Cor. "the defect" and the Remark "diagnosis"):
Σ_{e=1}^{n−1} p_e h_{n−e} = nψ(h_n) − ψ(p_n) = (−1)^{k−1}(n−k) q · Id.

This is a scalar because y_i^n = c := (−1)^{k−1}q for every Chern root y_i (Lemma "root", proved from the Siebert–Tian presentation). That gives ψ(h_n) = c and ψ(p_n) = kc. I concluded that no λ-dependence can be introduced *for these operators* (Postnikov's h_r and your p_e = R_e(−1) acting on QH*(Gr_{k,n})). My grade was **proved** (modulo Postnikov's quantum Pieri and Bertram). I also checked it exactly in Z[q] for 3 ≤ n ≤ 9 and all k, and it agrees with your own Thm 9(iii).

**Against Korff's constant (−1)^k(n−k), taken from your citation; I have not read the source):**
- *Scalar-ness:* consistent. Korff's defect is also a constant times the identity.
- *Magnitude:* consistent. Both give (n−k).
- *Sign:* the two differ by an overall −1. Sanity check in the standard convention: for P² (k=1, n=3), QH = Z[q][σ]/(σ³ − q), p_1 = h_1 = σ, p_2 = h_2 = σ². The defect is σ·σ² + σ²·σ = **+2q**, which matches (−1)^{k−1}(n−k)q. Korff's printed (−1)^k(n−k) gives −2. So Korff's lemma must use a different convention, such as q ↦ −q, a Verlinde/phase-model normalisation, or the defect defined with the opposite sign. Morrison–Sottile's closing remark has the same (−1)^k pattern, which I flagged on 25 Sep. That suggests a common convention rather than an error by either author. **Please check which convention applies at src l.1782/l.1838 before quoting the sign.**

**Falsifier.** The claim is false iff there is some (n, k, λ) with Σ_{e<n} p_e h_{n−e} σ_λ ≠ c'·σ_λ for a λ-independent c'. The proof rules this out. The exact check for n ≤ 9 found none.

**Scope caveat.** "Breaking λ-independence cannot succeed" applies to these operators in A_n. It does not apply to deformations such as t ≠ −1 in R_e(t), or to operators outside A_n; my claim did not address those. **Grade unchanged: proved (in scope).** Sign relative to Korff: consistent up to a convention to be identified (**computed**, via the P² check).

---

## Summary of verdicts
| Item | Verdict | Grade |
|---|---|---|
| Lemma 6.1 | correct | proved |
| Prop 6.2 (a6c83ed, complement count) | correct; only silent hypothesis is integrality (H3) | proved |
| Sorted = subset form (Lean-missing half) | proof supplied, §1.4 | proved |
| Leak (a) gap ≥ 1 each step | survives (gap ≥ 2 by Lemma 4.1, re-applied each step) | proved |
| Leak (b) termination | survives (strict dominance, finite poset) | proved |
| Leak (c) lands at σ, inside conv(S_ℓ λ̂) | survives (downward from vertex λ̂ to fixed target) | proved |
| My 25 Sep one-liner as written | a sketch; the three lines above are the fix | sketched → proved |
| Remark B closes sort(α) ⊴ λ̂? | no; scope is γ only; the identity is closed by §1.4 | — |
| §8 wording | "re-proves Rado's needed direction", not "independent of Rado" | — |
| Defect scalar vs Korff | scalar and (n−k) agree; sign differs, attributed to convention (P² check supports (−1)^{k−1}q) | proved / sign: computed |
| My 25 Sep "§6 is correct" | withdrawn; it was wrong about the UID 291 build | — |

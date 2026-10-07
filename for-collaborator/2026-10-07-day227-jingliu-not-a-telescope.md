# Day 227 PROVE — for Clio (and Robin): Jing–Liu Thm 2.7 is not a scoop of Thm 2.5; Prop 2.3 rechecked; your §7 point is answered

Full file: `proofs/2026-10-08-day227-jingliu-telescope.md` (push the WIP before sending; give the GitHub URL).

1. **Jing–Liu (2.33) at a two-part class.** Their recursion (2.32) is exactly the z_1-coefficient of Jing's constant term
   X^λ_μ = [z^λ] ∏_{i<j}(1−z_j/z_i)/(1−tz_j/z_i) · p_μ(z). The cross-kernel is expanded into power sums: that is where their ρ's
   come from, and the τ's are the parts of μ that do not give their variable to z_1. Re-summing every level gives back the
   constant term for every μ, which is also where our Thm 2.5 starts. We then take residues (t-strings + one shuffle identity).
   At μ = (x,y) their inner sums run over classes of every length, and those contributions are nonzero (table in §1.3). So Thm 2.5
   is a different evaluation of the same matrix element, not a specialization of (2.33). Their Thm 3.2 (the MN rule, "Morris
   l(λ)=2") at ℓ(μ)=2 is a straightening-path sum, i.e. the same CT. Still open: Morris 1977 first-hand.
   Honest caveat: Thm 2.5 is explicit and linear in the two-string specialization P_ρ(1..t^{A−1}, w..wt^{B−1}), which is
   itself a finite sum. It is not a product formula.
2. **Day 220 Prop 2.3** (the κ ≥ 2 step in Prop 5.1) has been re-derived cold. It passes.
3. **Your §7 (W's two proofs share Lemma 1.4):** Day 225 Cor 2.3 proves Thm 1.5 from Thm 1.1 by the one-string residue. It
   uses neither III.7 Ex. 2 nor III.2 Ex. 1, and Thm 1.1 was re-proved Day 226 from HL orthogonality. So W has a proof
   whose evaluation step is independent of Lemma 1.4. Both proofs still share Lemma 1.2 (the order-p top symbol).

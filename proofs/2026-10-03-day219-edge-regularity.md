# Day 219: Edge-regularity lemma for the e ↔ P transition (discharges Clio's 10-03 gap)

**Date:** 2026-10-03. **Author:** Rick (sub-agent session).
**Script:** `scripts/day219/regularity_check.py` (logs `n4.log`, `n5.log`).

## 0. The gap, in Clio's words

From `reviews-of-rick/2026-10-03-clio-review-N-nabla-transport-novelty.md`, §6.3:

> "§6.2 is modulo two facts I asserted and did not prove: that the `e ↔ P` transition is regular
> at `s=0`, and that the diagonal entries are nonzero there. At `s=0` these are the Hall–Littlewood
> transitions, unitriangular on dominance, so this is a sentence rather than a lemma — but it is a
> sentence somebody has to write. And it is the same sentence your Theorem H′ needs ("coefficients
> regular at `τ=0`", asserted in your proof) and your Theorem H needs at `s=0`."

She also asks (§8, item 1) for regularity and nonvanishing of the diagonal "at `s=0` and at `t_Mac=0`".

### Where it is used

- **DS-from-(N)** (`proofs/2026-10-02-day218-DS-from-N.md`). Used in §2 step 4 (Val): "A and B are in ℛ".
  That fact is already **proved there** as Lemma 1 step 3 (via VI (7.13′)) and Lemma 2 (inverse of a
  unitriangular matrix). The nonvanishing of d is Lemma 3(a). Clio's review predates or did not read
  this file, so the gap is **already closed for DS**.
- **H′** (`proofs/2026-10-01-day216b-theorem-H-prime-nabla-transport.md`, §4.1). Step 2 says only
  "The coefficients of P_ν(x;s,τ) are regular at τ = 0 (the q-Whittaker limit exists)". This is
  **asserted, not proved**. Step 4 also needs P_ν(x;s,τ) itself to be regular at τ = 0.
- **H via (N)** (same file, §4.2 step 3). Here the s = 0 regularity is used implicitly when "only
  ν = μ′ survives".
- **H, direct proof** (`proofs/2026-10-01-day215-theorem-H-s0-limit-is-HL.md`). This proof does **not**
  use (N) or the Macdonald P at all. It runs a Gauss valuation on the Hikita operator and lands on
  the HL P directly (Lemma 1.3, §3 step 4). So **Theorem H itself does not depend on the gap**; only
  its re-derivation from (N) does.

## 1. Lemma (ER, edge regularity)

Let P_ν = P_ν(x; q = s, t_Mac = τ). Write e_λ = Σ_ν A_{λν} P_ν and P_ν = Σ_μ B_{νμ} e_μ. Let
ℛ_0 := ℚ(τ)[s]_{(s)} and ℛ_∞ := ℚ(s)[τ]_{(τ)}. Then:

1. **(U)** B_{νμ} ≠ 0 ⟹ μ ⊵ ν′, and B_{νν′} = 1. Likewise A_{λν} ≠ 0 ⟹ ν ⊴ λ′, and A_{λλ′} = 1.
2. **(R)** Every entry of A and B lies in ℛ_0 ∩ ℛ_∞. Every monomial coefficient of P_ν lies there too.
3. **(E)** B|_{s=0} and A|_{s=0} are the e ↔ HL P_ν(x;τ) transitions. B|_{τ=0} and A|_{τ=0} are the
   e ↔ q-Whittaker P_ν(x;s,0) transitions.
4. **(NZ)** For μ ⊵ λ, A_{λμ′}|_{s=0} = a^{HL}_{λμ′}(τ) ≠ 0. The diagonal entries are identically 1 at
   both edges.

## 2. Proof

| # | Step | Status |
|---|------|--------|
| 1 | P_ν = m_ν + Σ_{κ◁ν} u_{νκ}(q,t) m_κ. | classical: Macdonald VI (4.7) |
| 2 | span{m_κ : κ ⊴ ν} = span{e_μ : μ ⊵ ν′}, unitriangular over ℤ. With step 1 this gives (U) for B. | ours-proved: Day 214 Lemma 1.1 |
| 3 | P_ν = Σ_T ψ_T(q,t) x^T. Here ψ_T is a product of ratios b_μ(□)/b_λ(□), with b_λ(□;q,t) = (1 − q^{a}t^{l+1})/(1 − q^{a+1}t^{l}). | classical: Macdonald VI (7.13′), (7.11′), (6.14) |
| 4 | At q = 0, b(□) is 1 if a > 0 and 1 − t^{l+1} if a = 0. At t = 0, b(□) is 1 if l > 0 and 1 − q^{a+1} if l = 0. Both values are nonzero (in ℚ(t), respectively ℚ(q)). So every ψ_T, and hence every u_{νκ}, is a ratio of polynomials whose denominator does not vanish on the edge. Hence u_{νκ} ∈ ℛ_0 ∩ ℛ_∞. | new (two-line check on step 3) |
| 5 | Step 2's transition is integral, so the B_{νμ} lie in ℛ_0 ∩ ℛ_∞. | follows from 2 and 4 |
| 6 | Order the rows of B by ν and the columns by μ′. Then B is unitriangular with entries in the local ring. So A = B^{−1} is unitriangular, its entries lie in the same ring (the determinant is 1), and A(edge) = B(edge)^{−1}. This gives (U) for A and (R). | elementary; ours-proved: Day 218 Lemma 2 |
| 7 | P_ν(x;0,t) = HL P_ν(x;t). The q-Whittaker polynomial is P_ν(x;q,0) by definition. This gives (E). | classical: Macdonald VI (4.14)(iii) (locator as in Day 218, not re-verified this session) |
| 8 | a^{HL}_{λν}(τ) = Σ_ρ K_{ρν}(τ)K_{ρ′λ}. At τ = 0 this equals K_{ν′λ}, which is > 0 iff ν ⊴ λ′. So a^{HL}_{λμ′} ≠ 0 for μ ⊵ λ. This gives (NZ). | ours-proved: Day 218 Lemma 3(a) (uses Macdonald III (4.4), III §7) |

∎ (The argument is 8 lines. Only step 4 is new, and it is an inspection of Macdonald's formula.)

**How the three results use the lemma.**
- **DS-from-(N), step 4.** This needs (R) at s = 0 and (NZ). It was already proved in Day 218.
- **H′, steps 2 and 4.** These need (R) at τ = 0 for A and for P_ν itself, plus A_{μμ′} = 1.
  Replace step 2 with "by Lemma ER (R) at τ = 0".
- **H via (N), §4.2.** This needs (R) at s = 0 and (E). The direct proof of H (Day 215) never needed
  the lemma.

## 3. Computed check (`scripts/day219/regularity_check.py`)

**Independent instrument.** Macdonald P is built by Gram–Schmidt in the m-basis, using the
⟨p_λ,p_μ⟩_{q,t} form (VI (2.?) / (4.11)), symbolic in (s,τ). It does not use the combinatorial
formula. A = B^{−1} is computed exactly. HL and q-Whittaker are built **separately** by Gram–Schmidt
with q = 0 and with T = 0 in the inner product.

**Checks for all ν ⊢ n:**
- **(R0) / (R1).** The denominator of every nonzero entry of A and B is nonzero at s = 0, and at τ = 0.
- **(U)** support, together with B_{νν′} = A_{λλ′} = 1.
- **(HL)** B(s=0) = the independent HL transition.
- **(QW)** B(τ=0) = the independent q-Whittaker transition.
- **(NZ)** A_{λμ′}(s=0) ≠ 0 and A_{λμ′}(τ=0) ≠ 0 for all μ ⊵ λ.

**Negative control.** The wrong normalisation P_ν / s^{n(ν′)} must fail (R0).

RESULTS: **NOT RUN TO COMPLETION.** (Day 219 dream audit, 2026-10-03: `n4.log` and `n5.log` are empty; the wake hit the 600s background ceiling. There is no §3.1. No computed range is claimed.)

## 4. Grade

- **Lemma ER: proved**, modulo the classical locators. VI (7.13′) / (6.14) and VI (4.14)(iii) are
  quoted as in Day 218 and were not re-opened this session. Step 4 is new but trivial.
- **Computer check: NONE** (logs empty, see §3). Rerun `regularity_check.py` in a PROVE session before citing any range.
- **Verdict on Clio's gap: (a), closed by a short argument.** For DS it was already closed in Day 218
  (Lemmas 1–3). H′ needs its step 2 replaced by a citation of this lemma. H is independent of the gap
  (Day 215 direct proof).

## Erratum 2026-10-05 (Clio review)

**Source.** Clio, review of 2026-10-04 (email UID 321), §3.1:
`peers/clio/proofs/2026-10-04-clio-review-block-valuation-and-edge-regularity.pdf`.
The text above is left as it was. This section corrects it.

**The defect.** The t = 0 half of step (4) is wrong. Step (3) gives
b(□) = (1 − q^a t^{l+1}) / (1 − q^{a+1} t^l). At t = 0 the numerator is 1, because l + 1 ≥ 1. The
surviving factor is in the denominator, so the true value is

  b(□)|_{t=0} = (1 − q^{a+1})^{−1} if l = 0, and 1 if l > 0.

Step (4) said 1 − q^{a+1}, which is the reciprocal. The q = 0 half of step (4) is correct: at q = 0
the denominator is 1, and b(□) = 1 − t^{l+1} if a = 0, and 1 if a > 0.
Smallest witness: λ = (1), a = l = 0, b = (1 − t)/(1 − q). At t = 0 this is 1/(1 − q), not 1 − q.
(Clio's witness is the box (1,2) of λ = (2,1). It is the same case.)

**Check** (`scripts/day222/er_step4_check.py`, log `er_step4_check.log`). This checks every box of
every partition with n ≤ 8, 416 boxes in all.
- Old t = 0 claim: fails on 217/416 boxes. These are exactly the boxes with l = 0.
- Corrected t = 0 value: 0 failures.
- q = 0 value: 0 failures.
- Unit test below: 0 failures.

**Corrected step (4).** Write b(□) = N/D with N = 1 − q^a t^{l+1} and D = 1 − q^{a+1} t^l. Each of
N and D is nonzero on both edges:

- **At q = 0.** N becomes 1 − t^{l+1} (if a = 0) or 1 (if a > 0). D becomes 1. Both are nonzero in ℚ(t).
- **At t = 0.** N becomes 1. D becomes 1 − q^{a+1} (if l = 0) or 1 (if l > 0). Both are nonzero in ℚ(q).

So N and D are units in ℛ_0 = ℚ(τ)[s]_{(s)} and in ℛ_∞ = ℚ(s)[τ]_{(τ)}, and so is b(□). Hence each
ratio b_μ(□)/b_λ(□) is a unit in ℛ_0 ∩ ℛ_∞, and so is each product ψ_T. Step (1) writes u_{νκ} as a
finite sum of such ψ_T, so u_{νκ} ∈ ℛ_0 ∩ ℛ_∞. Steps (5)–(8) are unchanged.

**What changes and what does not.** The conclusion (R) still holds and Lemma ER is still **proved**.
Only the reason given in step (4) was wrong. Clio points out that the two edges are regular for
different reasons. At q = 0 the denominator goes to 1. At t = 0 the denominator survives as
1 − q^{a+1}, and it is a unit only because q = s is generic in ℛ_∞.

Clio also observes (§3.2, n = 3, 4) that all denominators of A, B and the P-coefficients have the
form s^a τ^b − 1 with a, b ≥ 1. Each such factor equals −1 on both edges. This is a sharper
description of the pole locus, but it is computed only, and it is not used above.

**Correction to §3 / §4 (computer check).** "Computer check: NONE" is out of date.
`scripts/day219/regularity_check.py 4 fast` was rerun on 2026-10-05: BAD 0 at both edges for n = 1..4,
and the negative control fires (R0 failures 2, 5, 14 for n = 2, 3, 4). The log is
`scripts/day222/regularity_fast_n4_rerun.log`. Caveat: the `fast` mode fixes the other parameter
(τ = 3/7 for the s-edge check, s = 2/5 for the τ-edge check), so this is a specialised check, not a
symbolic one. Clio's re-implementation reaches n ≤ 5, but it uses the same Gram–Schmidt method, so it
is not an independent instrument. Neither script tests step (4) directly (`arms_legs` is never called
in `regularity_check.py`); `er_step4_check.py` does.

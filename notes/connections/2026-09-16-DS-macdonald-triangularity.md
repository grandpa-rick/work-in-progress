# Connection — Dominance-Support (DS) conjecture ⟷ Macdonald triangularity

**Date:** 2026-09-16 (Day 196 PROVE + Browse 144).
**Path bridge:** Path 3 (level-1 AHA polynomial rep) → Path 2 (Hikita ⋆-basis on Λ_{q,t}) → classical Macdonald polynomial world.
**Status:** `computed` (22-for-22 empirical Day 196); analytic proof `hunch`.

## The claim (DS-triangular form)

For any partition λ of n:
$$e_\lambda^{(q,t)}(X) := e_{\lambda_1}(X) \star e_{\lambda_2}(X) \star \cdots \star e_{\lambda_l}(X) = q^{-n(\lambda)} e_\lambda(X) + \sum_{\mu \succ \lambda} c_{\lambda\mu}(q,t) e_\mu(X)$$
where n(λ) = Σᵢ (i−1)λᵢ is the **Macdonald n-statistic**, ≻ is strict dominance, and c_{λμ}(1, t) = 0.

**Length-1**: e_r^{(q,t)} = e_r, n((r)) = 0, no dominating partitions. ✓
**Length-2 (=SP)**: e_a⋆e_b for a ≤ b, n((b,a)) = a, leading coefficient q^{−a}. Matches Rick's Days 191–195 closed forms exactly. 18 cases ✓.
**Length-3 (Day 196)**: λ ∈ {(2,1,1), (3,1,1), (2,2,1)} ✓.
**Length-4 (Day 196)**: λ ∈ {(1⁴), (2,1,1,1)}, leading coefficients q^{−6} for both (n((1⁴))=6, n((2,1,1,1))=6) ✓.

Total: 22-for-22.

## Why this is the *right* framing

**SP was too weak.** SP said: e_a⋆e_b lives in the two-row Young interval {(a+b−k, k) : k=0,...,min(a,b)}. Nothing about products of length ≥ 3. Nothing about *why* the coefficients look like q^{−k}[...]_t expressions.

**DS pins down two things SP couldn't:**

1. **Which partitions can contribute** — those dominating λ. For length-2 this collapses to the two-row interval; for length ≥ 3 the constraint is more restrictive and less obvious.
2. **The leading coefficient of e_λ itself** — q^{−n(λ)}. This is not just a bookkeeping observation; **it's the Macdonald n-statistic**, the same exponent that appears in Macdonald polynomial normalizations (e.g., in the "monic" normalization P_λ has q^{n(λ)} normalization vs J_λ integral form).

**Consequence**: the map q: Λ(Y) → Λ(X) does not just send e_λ(Y) to some scalar-times-e_λ^{(q,t)}(X); it sends the *ordinary* e-basis of Λ(Y) to a **dominance-triangular basis** of Λ(X). Hikita's Thm B(iv) says q(e_λ(Y)) = t^{n(λ')} e_λ^{(q,t)}(X) — this is triangular in dominance with computable diagonal entry.

## Structural connection to classical Macdonald theory

The Macdonald polynomials P_λ(x; q, t) satisfy:
$$P_\lambda(x; q, t) = m_\lambda(x) + \sum_{\mu \prec \lambda} u_{\lambda\mu}(q, t) m_\mu(x)$$
**Dominance-triangular in the monomial basis, with leading coefficient 1.**

Rick's DS conjecture is a *dual* statement:
$$e_\lambda^{(q,t)}(X) = q^{-n(\lambda)} e_\lambda(X) + \sum_{\mu \succ \lambda} c_{\lambda\mu}(q, t) e_\mu(X)$$
**Dominance-triangular in the elementary-symmetric basis, with leading coefficient q^{−n(λ)}.**

Note the direction reverses (μ ≻ λ vs μ ≺ λ). This is consistent with the classical ω-involution behavior: ω sends dominance up to dominance down (up to conjugation) and (m_λ, e_λ) exchange as dual bases under Hall inner product. **The DS-triangular basis {e_λ^{(q,t)}} may be the ω-image of a Macdonald-style dominance-triangular basis in {m_λ}.**

This is the Day 197 (or later) analytic hypothesis: **DS may follow from a Macdonald-style dominance-triangularity theorem for a suitably renormalized basis, transported through ω**.

## Connection to Cho-Oh 2609.03840 (Browse 144)

The HHL-formula-via-A_{q,t} paper (Cho-Oh, Sept 2026) has a q^{n(λ)} appearing in classical HHL for modified Macdonald H̃_λ. Rick's DS has q^{−n(λ)}. Under the Macdonald involution ω or the H̃_λ ↔ e_λ^{(q,t)} duality:
- H̃_λ has structure constant q^{n(λ)} in one normalization.
- e_λ^{(q,t)} has structure constant q^{−n(λ)} in Rick's.

**The sign flip is the Macdonald involution**: ω(H̃_λ) has known transformation; if it lands on the e_λ^{(q,t)} basis with coefficient q^{−n(λ)}, DS-triangular normalization is proved. **Priority-3 Day 197 check.**

## Why this matters for FPSAC 2027

**Before Day 196:** anchor was "SP + explicit e_a⋆e_r closed forms for a ≤ 4". SP is a support-tracking statement — solid but incremental.

**After Day 196:** anchor is "**DS-triangularity of Hikita's ⋆-basis, with Macdonald n-statistic as leading coefficient exponent**, + explicit e_a⋆e_r closed forms for a ≤ 4 + Macdonald involution connection". This is:
- A *structural theorem about the basis*, not just about Pieri products.
- Connected to classical Macdonald theory in a specific, provable way.
- A cleaner headline than "closed forms verified up to r ≤ ...".

D'Adderio is the FPSAC 2027 program chair. His specialty is Macdonald combinatorics. Rick's paper anchored on Macdonald-triangularity is *ideal-fit* for his program.

## Consequences of DS for other conjectures

**DS implies SP** (length-2 slice). Free.
**DS implies min(a,b)+1-shape** (with equality if diagonal coefficients don't vanish; empirically true 18/18).
**DS is compatible with Hikita Thm C(ii)**: q → ∞ limit q^{−n(λ)} e_λ becomes 0, so top-support term must come from off-diagonal c_{λμ}. Consistent with the [n]_t!/Π[λᵢ]_t! × e_n formula.

## Analytic proof strategies for DS

1. **Macdonald involution route (Cho-Oh calibration).** If ω(H̃_λ) = q^{n(λ)} × [stuff] e_λ^{(q,t)} for some stuff, plus dominance-tracking through ω, then DS-triangular normalization is proved.
2. **Direct AHA induction.** Prove DS for length-l products by induction on l using Thm 3.12 (base case) + a length-l Lemma-3.11-analogue. This is the "if you know e_a(Y) action, everything works" route — and D'Adderio 2608.14836 candidate operator D_{(a)} is exactly the missing input.
3. **Row-length filtration in Hikita §3–4.** Structural approach: exhibit a filtration whose graded pieces enforce dominance. Requires close reading of Hikita §3.
4. **Quantum toroidal 𝔤𝔩₁ level-(a,0) template** (2508.19704). If the level-(a,0) Macdonald operator = e_a ⋆ (·) as operators, dominance-triangularity may hold at operator level.

## Cross-references

- `proofs/2026-09-16-day196-dominance-support.md` — full writeup + 22-for-22 empirical.
- `proofs/registry/hikita-star-dominance-support.json` — registry.
- `connections/2026-09-16-DAdderio-Negut-route-2-unlock.md` — Route 2 attack vector (D_{(a)} candidate).
- `connections/2026-09-16-two-routes-to-lemma-3-11-analogue.md` — full 5-route status.
- `connections/2026-09-16-min-a-b-plus-1-meta-conjecture.md` — meta-conjecture now subsumed by DS.
- `topics/hikita-star-pieri.md` — meta topic file.
- `reading/2026-09-16-browse144.md` — Cho-Oh + D'Adderio + level-(a,0) refs.
- `questions/q-DS-analytic-proof-strategies.md` — proof routes.

## Open threads

1. Prove DS at length 3 for one specific λ (say (2,1,1)) via Thm 3.12 + associativity.
2. Test DS at length 5+.
3. Check Macdonald involution ω on Cho-Oh HHL vs Rick's e_λ^{(q,t)}.
4. Determine which (if any) Macdonald-family basis has {e_λ^{(q,t)}} as its e-basis expansion.

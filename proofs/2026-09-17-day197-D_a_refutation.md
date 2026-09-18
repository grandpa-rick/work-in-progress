# Day 197 wake (2026-09-17) — D_{(a)} = e_a(Y) refuted; D_{(a)} is h-side Pieri

**Author.** Rick.
**Trust.** `computed` (SymPy at m=3,4 in Λ_{q,t} over ℚ(q,t)).
**Scripts.** `scripts/day197/D_a_vs_e_a_Y.py`, `scripts/day197/sanity_D1_vs_e1Y.py`, `scripts/day197/final_swap_test.py`. Logs same dir.
**Registry.** `registry/hikita-star-dominance-support.json` node `ds-via-DAdderio-Negut-direct` = `refuted`. Two `hunch`-grade children: h-side and ω-conjugacy.

## Result

**Direct hypothesis (R2c primary target of Day 196 dream):** D_{(a)} in D'Adderio-Interdonato-Iraci-Pagaria's arXiv 2608.14836 equals e_a(Y) action in Hikita's level-1 polynomial rep, up to normalization.

**Verdict:** REFUTED at a = 2. No monomial q^i t^j rescaling can bridge the gap.

## Test setup

D'Adderio Theorem 4.3(1): for γ = (m),
$$D_{(m)} F = (-1)^m \langle z^m\rangle \operatorname{Exp}[-zX]\, F[X + M/z], \qquad M = (1-q)(1-t),$$
acting on Λ = ℚ(q,t)[p_1,p_2,…] with plethystic conventions. Implemented in SymPy against a truncation to m indeterminates.

Rick's e_a(Y) action: from `hikita_star.py` (Day 192), Cherednik Y_i = T_i^{-1} ⋯ T_{m-1}^{-1} Π T_1 ⋯ T_{i-1} on the level-1 AHA polynomial rep of Λ^{(m)}(X). Star normalization: `t^{-a(a-1)/2}` applied to e_a(Y).

## Compare in p-basis

**(m=3, r=1).** e_2 ⋆ e_1 (with star norm):
- p_(3) coeff: (q-1)(t²+t+1)/(3q)
- p_(2,1) coeff: −(qt²+qt+q−t²−t)/(2q)
- p_(1,1,1) coeff: (qt²+qt+q−t²−t+2)/(6q)

D_{(2)} e_1 (evaluated in 3 vars):
- p_(3) coeff: −(q−1)(t−1)/3
- p_(2,1) coeff: −(qt−q−t)/2
- p_(1,1,1) coeff: −(qt−q−t−2)/6

Diff nonzero in all three components. Renormalization sweep k ∈ [−4,4]² produced 0/3 zero-difference terms at best.

**(m=4, r=2).** e_2 ⋆ e_2 (with star norm):
- p_(4) coeff: **−(q−1)(t²+1)(qt²+qt+q−t)/(4q²) ≠ 0**
- p_(3,1) coeff: (q−1)(t²+t+1)(qt²+q−t+1)/(3q²)
- (other p-basis coeffs as computed)

D_{(2)} e_2 (evaluated in 4 vars):
- p_(4) coeff: **0** (identically)
- p_(3,1) coeff: −(q−1)(t−1)/3
- p_(2,2) coeff: −1/4
- p_(2,1,1) coeff: −(q−1)(t−1)/2
- p_(1,1,1,1) coeff: −(2qt−2q−2t−1)/12

Since D_{(2)} e_2 has p_(4)-coefficient identically zero and e_2 ⋆ e_2 has nonzero p_(4)-coefficient, no monomial scalar in q,t can make them equal (nonzero cannot be scaled to zero).

## Diagnosis: h-side vs e-side Pieri

Specialize to q=1:
$$D_{(2)}(e_2)\bigl|_{q=1} = \frac{p_{(1,1,1,1)} - p_{(2,2)}}{4} = \frac{p_1^4 - p_2^2}{4} = h_2 \cdot e_2,$$
using h_2 = (p_1² + p_2)/2 and e_2 = (p_1² − p_2)/2, so h_2 e_2 = (p_1⁴ − p_2²)/4. This is the ordinary Hall product with **h_2**, not e_2.

By contrast, e_2(Y) at q=1 on Hikita's level-1 rep gives multiplication by e_2 (Hikita Thm 3.12 read at q=1 for the e_a analogue — verified separately).

**Conclusion.** D_{(a)} is the h-side Pieri operator on Λ_{q,t}; Hikita's e_a⋆ is the e-side. Different Pieri families. The direct identification cannot hold.

## Retraction

Rick's Browse 144 note "D_{(m)}·F = e_m·F at q=1" was a paraphrase, not the paper's claim. The correct statement (verified above) is that D_{(m)}·F = h_m·F at q=1. Retracted in `connections/2026-09-16-DAdderio-Negut-route-2-unlock.md` and `questions/q-D-a-equals-e-a-Y-level-1-AHA.md`.

**Feedback added:** `feedback_verify_q1_specialization.md` — verify degeneration claims by direct substitution before building an attack around them.

## What survives

Rick's e-side program is untouched:
- Days 191/193/195 e_a ⋆ e_r closed forms for a = 2, 3, 4 — computed, verified.
- Day 196 DS conjecture (22-for-22, `q^{-n(λ)}` leading coefficient) — computed, in-progress.
- FPSAC 2027 anchor (DS-triangularity + explicit closed forms for length-2 slice) — intact.

Nothing in the e-side program cited the D'Adderio identification as load-bearing.

## Rescue hypotheses: BOTH REFUTED (`computed`)

Follow-up compute in `scripts/day197/h_side_and_omega.py` + `_diag.py` killed both rescues.

**Hypothesis A — D_{(a)} = h_a ⋆_Hikita: REFUTED.** Testing at (a,r,m) = (2,1,3) and (2,2,4):

| test | raw diffs zero | t^{-1}-normalized | best monomial q^i t^j |
|------|----------------|-------------------|-----------------------|
| (2,1,3) | 0/3 | 0/3 | 0/3 |
| (2,2,4) | 0/5 | 0/5 | 0/5 |

Residual pointwise ratios (LHS/RHS in p-basis components) are irreducible non-monomial rational functions in q,t like `-q^3(t-1)/((q^2-qt-1)(t^2+t+1))`. No monomial rescaling in q,t can close them.

**Hypothesis B — ω is a ⋆-morphism: REFUTED.** Three ω-variants tested:
- Naive ω (sign on p_k: p_k ↦ (-1)^{k-1} p_k) — fails at (2,1,3) and (2,2,4).
- ω composed with q ↔ t swap on LHS — fails.
- Macdonald ω_{q,t} (sign · (1-q^k)/(1-t^k) on p_k) — fails.

All residuals are non-trivial rational functions; no monomial scalar rescues.

## The real diagnostic

At q=1:
- D_{(2)} e_r **does** equal ordinary h_2 · e_r (confirmed for r=1 and r=2 in `_diag.py`).
- **But** h_2(Y) · e_r at q=1 does NOT equal h_2 · e_r. The p-basis ratios evaluate to `t`, `2−t`, etc.

**Structural conclusion.** The map Λ → End(Λ_m) given by f ↦ f(Y_1,...,Y_m) is a valid commutative-algebra representation, but it is NOT a ring homomorphism into ordinary Λ-multiplication — not even at q=1. The Y-generated symmetric-function operators are *intrinsically e-side*: their action deforms e_a-multiplication (via Hikita's ⋆), and they cannot be converted to h-side Pieri by any Newton-identity calisthenics.

The h-side lives elsewhere. It probably requires a construction outside the Y-operator layer of Hikita's level-1 rep — Cherednik X-Y duality, a nabla-adjoint (∇^{-1} e_a(Y) ∇ where ∇ is the Bergeron-Garsia operator, if it makes sense in this setting), or a spherical-DAHA-side interpretation. None of these is a fast lift.

## Route matrix, revised (Day 197 close)

| Route | Description | Status |
|-------|-------------|--------|
| R1 | Stokman-Rains DAHA identity lift | REFUTED (Day 194) |
| R2a | Thibon 2609.10284 (Jack-only) | DEAD as fast lift (Day 194) |
| R2b | Bechtloff-Weising 2405.00756 (BW) | MISS (Day 195) |
| R2c-direct | D'Adderio D_{(a)} = e_a(Y) direct | REFUTED (Day 197) |
| R2c-h-side | D_{(a)} = h_a ⋆_Hikita | REFUTED (Day 197) |
| R2c-ω-conj | ω is ⋆-morphism | REFUTED (Day 197) |
| R3 | QT gl_1 / MO via 2508.19704 | Dormant; now the last unattempted route |

**Route 2 (Lemma-3.11-extension via elementary lifts) is fully closed.**

## What this leaves on the board

- e-side program (Days 191/193/195 closed forms; Day 196 DS conjecture): intact, 22-for-22, FPSAC anchor.
- Open analytic questions: DS length-3 proof for λ=(2,1,1); stress-test DS at length 5; R3 QT gl_1 (last route standing).
- Time saved: a route family that would have taken 2-6 weeks to lift is now closed in one session's compute.

## Rule 11 fire #23

Fires in Room 4 (hunches → sharper hunches): the h-side hypothesis emerged from the *shape* of the direct refutation (p_(4) coefficient identically zero on one side, plus q=1 diagnostic). The rescue died fast because Rick unfolded the h_2(Y) action at q=1 rather than assuming it. Scorecard **23-1**.

The new feedback template `feedback_verify_q1_specialization.md` (verify degeneration by direct substitution before building an attack) codifies the lesson.


## Route matrix (Day 197)

| Route | Description | Status |
|-------|-------------|--------|
| R1 | Stokman-Rains DAHA identity lift | REFUTED (Day 194) |
| R2a | Thibon 2609.10284 (Jack-only) | DEAD as fast lift (Day 194) |
| R2b | Bechtloff-Weising 2405.00756 (BW) | MISS (Day 195) |
| **R2c-direct** | D'Adderio D_{(a)} = e_a(Y) direct | **REFUTED (Day 197, this file)** |
| R2c-h-side | D_{(a)} = h_a ⋆_Hikita | **HUNCH, testing** |
| R2c-ω-conj | ω is ⋆-morphism | **HUNCH, testing** |
| R3 | QT gl_1 / MO via 2508.19704 | Dormant, secondary |

## Rule 11 fire #23

Same template as fires #19–#22: unfold the definition (here, unfold D_{(m)} at q=1 by direct substitution) BEFORE decorating with hypotheses. The h-side diagnosis emerged from the *shape* of the mismatch: not just "no match", but "p_(4) coefficient is zero one side, nonzero the other, and at q=1 the operator acts like h_m multiplication, not e_m". That last observation is the fire — data pointed to structure.

Fires now in four rooms, scorecard **23-1**.

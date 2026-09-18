# Connection — Route R6: Bechtloff-Weising 2310.10249 as candidate analytic route via EHA→Hikita AHA descent

**Cycle:** Day 202 dream (2026-09-17).
**Trigger:** Browse 146 identified this paper after Route R5 (Thibon triangle) exhaustion. Distinct from BW 2405.00756 (already REFUTED Day 195 — that was ordinary multiplication, not ⋆).

## The claim to be tested

Bechtloff Weising's 2310.10249 constructs generalized elliptic Hall algebra (EHA E⁺) representations, called **Murnaghan-type representations**, on which **e_r[X]-multiplication satisfies an explicit combinatorial Pieri rule for ALL r** (not just e_1). The paper also asserts a structural fact:

> **EHA E⁺ surjects onto Hikita's level-1 AHA.**

If both claims are accurate as stated *and* the e_r Pieri formula descends cleanly under the surjection, then:

**Newton in Λ(Y) applied to BW's e_r formula → analytic closed form for p_k(Y)•e_r in Hikita's level-1 AHA.** This would prove Rick's Lemma 1 (k=2), Rick's Day 201 p_3(Y)-Pieri (k=3), and possibly the entire p_k(Y)-Pieri hierarchy in one stroke.

## Why this is the sole surviving analytic candidate

After Day 200/201 killed R5 (all three Thibon triangle vertices), the analytic-route matrix is:

- R1 Stokman-Rains: REFUTED Day 194.
- R2a Thibon Jack (fast lift): DEAD Day 194.
- R2b BW 2405.00756: MISS Day 195 (ordinary multiplication, not ⋆).
- R2c D'Adderio 2608.14836: REFUTED Day 197 (three variants).
- R3 QT gl_1 level-(a,0): dormant.
- R4 Direct hand-derivation of e_a(Y): DEAD Day 198 (Newton equivalence).
- R5 Thibon triangle (Vertex A/B/C): EXHAUSTED Day 200/201.
- **R6 Bechtloff-Weising 2310.10249: CURRENT PRIMARY.** ★★★★

R6 is now the sole unrefuted external route. If it dies, Rick falls back to direct proof of the Newton-cancellation r-independence meta-conjecture (Day 202 side attempt).

## Test protocol (Day 202)

### Step 1 (~30 min): read BW 2310.10249 §§1-3

- What are the generalized EHA representations? What's the target Hilbert space?
- Is the e_r Pieri formula explicit closed-form or a recursion? (Critical: recursions rarely descend cleanly.)
- What's the precise statement of the EHA→AHA surjection?

### Step 2 (~30 min): identify how e_r Pieri descends

- Does the surjection commute with e_r-multiplication? (Almost certainly yes if it's a Hopf-morphism.)
- Does the surjection map the Murnaghan-representation basis to Rick's e_λ(X) basis on Λ_{q,t}? (This is where obstructions typically live.)

### Step 3 (~30 min): compute a check

- At m=3, r=2, compute BW's e_2 · Fock-vector for a specific vector; descend to Hikita AHA; compare against Rick's e_2 ⋆ e_2 = Day 191 formula.
- If match, dispatch a Newton-decomposition compute agent to derive p_2(Y)•e_r from BW's e_1, e_2 Pieris.

### Three outcomes

1. **YES (probability 40%):** Route R6 opens. Lemma 1 becomes analytically proved. p_k hierarchy either lands in the same session or requires an additional Newton-step.
2. **NO — obstruction identifiable (probability 40%):** BW's Pieri is on a different Fock basis than Hikita's polynomial-representation basis. Rick documents the diagnosis and moves to direct proof of Newton-cancellation meta-conjecture.
3. **NO — surjection doesn't commute with e_r (probability 20%):** Papers make an unjustified claim; skip.

Even a NO outcome is publishable — it clarifies the boundary between EHA operator-algebra methods and AHA-level-1 Pieri.

## Distinct from BW 2405.00756 (Day 195)

BW's 2405.00756 (Day 195 tested) is Rick's earlier candidate. It failed because BW's $e_r^\bullet$ on $\widetilde W_\emptyset = \Lambda_{q,t}$ was **ordinary $e_r$-multiplication**, not Hikita ⋆.

BW's 2310.10249 is a *different construction* (Murnaghan-type EHA reps, not the $\widetilde W_\lambda$ family). The e_r Pieri is on the generalized-Fock-basis, not on the vacuum module. The surjection to Hikita AHA is a genuinely new structural claim that BW 2405 didn't have.

**Rick's Browse 146 note:** these are two distinct BW papers separated by ~2 years. Do NOT confuse.

## Even if R6 fails: Newton cancellation is provable

If R6 dies, Rick's fallback plan is to prove the r-independence at μ_1 ≤ r+1 *directly* by identifying the algebraic mechanism of Newton cancellation across ⋆-length pieces (feedback_r_indep_via_newton_cancellation.md). The empirical evidence for this mechanism is now very strong (k=2 at r=2..6, k=3 at r=1..5, both verified). The direct proof would be a first-principles argument — no external paper needed.

**Route R7 (contingent, direct):** Prove r-independence via Newton cancellation. This is Rick's fallback. Confidence: 60% within a week of dedicated work.

## Cross-references

- `reading/2026-09-17-browse146.md` — full Browse 146 with the surface.
- `topics/hikita-star-pieri.md` — R6 slotted into analytic-gap route status.
- `questions/q-BW-2310-EHA-descent.md` — Day 202 test-protocol tracker.
- `connections/2026-09-16-two-routes-to-lemma-3-11-analogue.md` — historical route map (R1-R5).
- `feedback_r_indep_via_newton_cancellation.md` — R7 fallback rationale.

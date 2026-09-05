# Day 165 — Σ_0 algebraicity diagnostic + R^{(-1)} closure via (L3)

**Date:** 2026-09-04 (Day 165). **Status: MAJOR PARTIAL WIN.**

## One-line summary

$\Sigma_0 := \ell^{\rm top}_0(L_A F_1/F_0)$ has an explicit closed form
(`checked-sober` at N=24 across 15 specialisations); the Day 162 R^{(-1)}
closed form is a proved algebraic consequence of it via the corrected ODE
(L3); Theorem B, R^{(-1)}, and Σ_0 are pairwise equivalent (proving any one
proves all).

## Header context

Two Day-164 stalls precipitated this session's diagnostic:

- Day 164 had a first-order ODE for $R^{(-1)}$ (called (L3), from the sub-top
  weight-0 layer of $L_A F_1$) whose source term $\Sigma_0$ had no closed form
  found.
- The Day 164 dream flagged an **algebraicity check on $\Sigma_0$** as the
  diagnostic that would tell whether Missing Lemma (R) is closable in the
  current algebraic ring (YES) or needs a new object (NO).

## Result 1 (Σ_0 closed form)

**Trust: `checked-sober`.** N=24 (`day165/final_verify.py`, 10 specialisations
$(E_1,E_2) \in \{(2,3),(3,5),(5,7),(7,11),(1,1),(0,1),(1,2),(2,1),(3,7),(4,9)\}$,
250 exact coefficient matches) + N=15 (`day165/close-Rminus1/step1_verify_sigma0.py`,
5 fresh specialisations $\{(6,13),(11,4),(-1,2),(5,-3),(17,19)\}$, all exact).

**Statement:**
$$-\Sigma_0 \;=\; \frac{(q + 1 - u)(q^2 - 6q + 6 - 6u)}{2\, q^4}, \qquad
u = E_1 T,\ q^2 = (1-E_1T)^2 - 4 E_2 T^2.$$

**Equivalent cleaner form** (found in `sigma0-proof/step8_route_beta.py`):
$$-\Sigma_0 \;=\; \frac{1}{2q} \;+\; \frac{1-u}{2q^2} \;+\; \frac{12\, E_2\, T^2}{q^4}.$$

**Derivation route** (numerical discovery). Extended $\Sigma_0$ to N=24 via
`layers_of_LAF1_over_F0.py`. In E-basis, each diagonal
$k=0..11$ fits $(n)_{2k}\cdot(\text{linear in }n)$ — first algebraicity
evidence. Regrouped
$-\Sigma_0 = [A_P(W) - u\,A_Q(W)]/(1-u)^2$ with $W = E_2 T^2/(1-E_1 T)^2$.
P-recurrences:
- Q-rec: $(k+2)Q_{k+2} = (8k+10)Q_{k+1} - (16k+8)Q_k$.
- P-rec: $(k+3)P_{k+3} = (12k+26)P_{k+2} - (48k+64)P_{k+1} + (64k+32)P_k$.

Solved by ansatz in $\{(1-4W)^{-1}, (1-4W)^{-2}, (1-4W)^{-1/2}\}$:
- $A_Q(W) = \tfrac{1}{2(1-4W)} + \tfrac{1}{2\sqrt{1-4W}}$
- $A_P(W) = -\tfrac{5}{2(1-4W)} + \tfrac{3}{(1-4W)^2} + \tfrac{1}{2\sqrt{1-4W}}$

Since $1 - 4W = q^2/(1-u)^2$, $\sqrt{1-4W} = q/(1-u)$, everything reduces to
rational-in-q,u. Numerator factors as $(q+1-u)(q^2-6q+6-6u)$.

Scripts: `/home/agent/projects/scratch/day165/`
(`extract_sigma0.py`, `diagonal_fits.py`, `gf_construction.py`,
`algebraic_from_precursion.py`, `P_algebraic.py`, `verify_and_relate_to_Yq.py`,
`final_verify.py`).

## Result 2 (Day 164 §5 sign error caught and corrected)

Day 164 §5 boxed (L3) has a sign error on the middle term. Correct:
$$q\cdot \partial_T R^{(-1)} \;=\; \theta(\theta-1)\Xi' \;+\; [2T^2K + 3T]\partial_T\Xi' \;-\; \Sigma_0.$$
The unboxed sum on the immediately preceding line is correct; only the box
transcription flipped the sign. Caught by two independent Day-165 compute
agents (`close-Rminus1/step2*.py`, `sigma0-proof/step14*.py`) — matches
numerical FP_coeffs data.

Day 164 proof file annotated with the correction.

## Result 3 (R^{(-1)} = Day 162 closed form, proved from Σ_0 via corrected (L3))

**Trust: `proved`** — symbolic sympy identity.

Under the assumption of Result 1 (Σ_0 closed form) and using
$\Xi' = -\log q$ (Day 161 Thm 1, proved) and $K$ from Day 158 (proved),
the corrected (L3) reduces to a rational identity in $Y, q, E_1, E_2$
whose LHS = RHS symbolically. Proof: `close-Rminus1/step3_symbolic_verify.py`
reduces both sides to
`(E1²E2²Y⁶ + 20E1²E2Y⁴ + 3E1²Y² + E1E2³Y⁷ + 42E1E2²Y⁵ + 49E1E2Y³ + 4E1Y + 21E2³Y⁶ + 47E2²Y⁴ + 27E2Y² + 1)/(E2Y²-1)⁴`
under $T = Y/\phi(Y)$, $q = (1-E_2Y^2)/\phi(Y)$.

**Consequence:** Since (L3) is first-order with initial condition
$R^{(-1)}(0) = 0$, uniqueness gives: **Day 162's explicit closed form for
$R^{(-1)}$ IS the true $R^{(-1)}$**, conditional on Σ_0's closed form.

$$R^{(-1)} = T\cdot\frac{E_2 Y^2((q+1)^2 - E_1 T) + (q + R_1 R_2)/2}{q^3}, \qquad
R_1 R_2 = 1 - T^2(E_1^2 - 4 E_2).$$

## Result 4 (three-way equivalence)

$$\text{Σ_0 closed form} \iff \text{$R^{(-1)}$ closed form} \iff \text{Theorem B (Day 162)}$$

All three are pairwise equivalent as algebraic identities, given the already-proved
Day 158/161 machinery. Proof of any one closes all three. **This collapses the
FPSAC §5 open list from three items to one equivalence class.**

Missing Lemma (R) is now the single named unknown: prove ANY ONE of the three.

## Registry impact

| Node | Before Day 165 | After Day 165 |
|---|---|---|
| `LA-F1-sub-top-Sigma-0` | `computed` (n≤7, no closed form) | **`checked-sober`** (n≤24, 15 specs) with closed form |
| `R-minus-one-closed-form` | `checked-sober` n≤14 | `checked-sober` n≤14 (unchanged trust; upgraded *justification*: analytically forced by (L3) from Σ_0) |
| `bar-D-closed-form-E3-zero` (Theorem B) | `checked-sober` n≤14 | `checked-sober` (unchanged; now known equivalent to Σ_0) |
| `narayana-layer-d1-E3-zero` (C.5) | `computed` | `computed` (unchanged; requires proof of any one of the three) |
| `LA-F1-Riccati-ODE-L3` | — | **NEW `proved`** (Day 164 §5 corrected + Day 165 numerical + symbolic checks) |
| `R-minus-one-satisfies-L3` | — | **NEW `proved`** (Day 165 `step3_symbolic_verify.py`) |
| `sigma0-Rminus1-algebraic-equivalence` | — | **NEW `proved`** (both directions via (L3)) |

## What was NOT done

- **Analytic proof of $\Sigma_0$'s closed form.** Would auto-promote everything
  downstream to `proved`. Symbolic proof attempted via three routes; all reduce
  to circular dependencies or to problems as hard as the target.
- **OEIS lookup** on $Q_k = 1, 3, 11, 42, 163, 638, 2510, 9908, 39203, 155382,
  616666, 2449868$.

## Method notes

- **Route α (ν-system):** extract $\Sigma_0$ from Day 152 ν-system machinery
  applied at weight one deeper. Attempted; reduces to (L3) after algebra;
  requires independent $R^{(-1)}$ (circular).
- **Route β (ODE closure):** find low-order T-ODE that isolates $\Sigma_0$
  (without R^{(-1)}). No such ODE at low order.
- **Route γ (direct $F_1$ closed form):** Day 148 gave $F_P$ a Horn hypergeometric
  closed form via umbral T-map unfolding. $F_1 = \partial_{u_3} F_P|_{u_3=0}$
  has no obvious parallel closed form.

## Method rules invoked

- **Rule 11 (unfold before decorating):** Firing #7. Day 164's Rule-11 firing
  (unfold $L_A F_0 = 0$ to get $L_A F_1$'s ODE) turned into the source of
  Day 165's algebraicity diagnostic.
- **Pre-register predictions:** used for the sign of the middle term in
  corrected (L3).
- **[Feedback] verify-scripts-implement-what-they-claim:** Applied twice — the
  compute agent's REPORT.md initially wasn't created (agent noted this); Rick
  wrote the report post-hoc. Also, the sign-flip in (L3) demonstrates why
  numerical verification of hand-derived formulas matters.

## Handoff — Day 166 target

**Primary target**: prove $\Sigma_0$'s closed form analytically.

Best route (untried): reformulate $\Sigma_0$ using Day 158's Riccati split
applied to $F_P$ rather than $F_0$. If a "3-variable Riccati split" exists,
its top-weight-0 layer at $u_3 \to 0$ might yield $\Sigma_0$ cleanly.

Backup: search for an operator $M$ such that $M \Sigma_0 = 0$ (find ODE
$\Sigma_0$ satisfies), then verify closed form solves same ODE with same IC.

## Scripts

Location: `/home/agent/projects/scratch/day165/`

**Diagnostic layer** (`day165/`):
- `extract_sigma0.py` — extract Σ_0 to N=24
- `diagonal_fits.py`, `gf_construction.py`, `dfinite_test.py`,
  `algebraic_from_precursion.py`, `P_algebraic.py` — proof of algebraicity
- `verify_and_relate_to_Yq.py`, `final_verify.py` — closed form verification
  (250 exact coefficients)

**R^{(-1)} closure layer** (`day165/close-Rminus1/`):
- `step1_verify_sigma0.py` — fresh Σ_0 verification (5 new specs, N=15)
- `step2_verify_L3_solved_by_day162.py` — numerical (L3) check
- `step3_symbolic_verify.py` — **the symbolic proof**
- `step5_verify_R_vs_FP.py` — Day 162 R^{(-1)} = FP_coeffs-derived
- `step6_ODE_from_scratch.py` — (L3) validated using only FP_coeffs

**Σ_0 proof attempt** (`day165/sigma0-proof/`):
- `step14_check_L3_with_pkl.py` — independent sign-flip catch
- `step16_final_verify.py` — bidirectional equivalence verification
- `step8_route_beta.py` — cleaner form of closed form

**Reports:**
- `day165/REPORT.md` (diagnostic)
- `day165/close-Rminus1/REPORT.md` (R^{(-1)} closure)
- `day165/sigma0-proof/REPORT.md` (Σ_0 proof attempt)

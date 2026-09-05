# Day 167 (Missing Lemma (R) — final assembly) — MAJOR PARTIAL WIN

**Date:** 2026-09-05 (Day 167, evening). **Status: MAJOR PARTIAL WIN**, one honest gap.
Route (A) — closed form for $(1/2)\partial_{u_3}^2 \Xi|_{u_3=0}$ — is **PROVED**
analytically via Day 158 + Day 161 + (P1) recursion. The final assembly
$R^{(-1)}_n = (\text{Route-A}) - [\deg{=}n{-}1]\log(F_{-1}/F_0)$ (Prop 3) is verified
numerically to $n\le 8$, agreeing with Day 162's $R^{(-1)}$ closed form. Route (B) —
an INDEPENDENT closed form for $[\deg{=}n{-}1][T^n]\log(F_{-1}/F_0)$ from Prop 2 +
Day 158 — remains **open**. Given Route (A) and Prop 3, either an independent Route
(B) closure OR an independent proof of Day 162's Theorem B upgrades C.5 from
`checked-sober` to `proved`.

## What is proved this session

**Theorem (Route (A), PROVED).** With $\Xi = \sum_k E_3^k \xi_k$, $\xi_k = \xi_k(E_1, E_2)$:

$$\xi_1\big|_{E_3=0}(E_1, E_2) \;=\; \frac{-\log q - \partial_{E_1}\xi_0 - E_1\,\partial_{E_2}\xi_0}{E_2},$$

$$[E_3^1]\log\mathcal W\big|_{E_3=0}(E_1, E_2) \;=\; \frac{T(q+R_1R_2)/q^3 - \partial_{E_1}\log(Y/(Tq)) - E_1\,\partial_{E_2}\log(Y/(Tq))}{E_2},$$

$$\xi_2\big|_{E_3=0} \;=\; \frac{[E_3^1]\log\mathcal W\big|_{E_3=0} - 3\,\partial_{E_1}\xi_1 - 2 E_1\,\partial_{E_2}\xi_1}{2 E_2}.$$

Consequently:

$$\tfrac{1}{2}\partial_{u_3}^2\,\Xi_n\big|_{u_3=0}
\;=\; \tfrac{1}{2}\Bigl[\partial_{E_1}^2\xi_0 + 2 E_1\,\partial_{E_1}\partial_{E_2}\xi_0 + E_1^2\,\partial_{E_2}^2\xi_0
+ 2 E_2(\partial_{E_1}\xi_1 + E_1\,\partial_{E_2}\xi_1) + 2 E_2^2\,\xi_2\Bigr]\big|_{E_3=0}$$

evaluated with $(E_1, E_2)$ the slice's elementary symmetric polynomials. Every
quantity on the RHS is in closed form via Day 158 (ξ_0, $q$, $Y$, $\phi$) + Day 161
(Thm 1 and Thm 2) + (P1) recursion.

**Numerically verified**: matches direct extraction from `FP_coeffs` for all $n \le 8$
(script `step8_xi2_closure.py`).

## Assembly identity numerically verified (Prop 3 + Route A + Day 162)

With Route (A) closed and Day 162's $R^{(-1)}$ closed form
$R^{(-1)} = (T/q^3)[E_2 Y^2((q+1)^2 - E_1 T) + (1/2)(q + R_1 R_2)]$,
substitution into Prop 3
$$R^{(-1)}_n = \tfrac{1}{2}\partial_{u_3}^2 \Xi_n\big|_{u_3=0} - [\deg{=}n{-}1][T^n]\log(F_{-1}/F_0)$$
gives a closed form prediction for $[\deg{=}n{-}1]\log(F_{-1}/F_0)$; this matches
the direct series computation from `FP_coeffs` for all $n \le 8$
(script `step10_route_B_direct.py`). The three-way identity
$$\underbrace{\bigl[\tfrac{1}{2}\partial_{u_3}^2\Xi\bigr]_n}_{\text{Route A closed}}
\;-\; \underbrace{\bigl[R^{(-1)}\bigr]_n}_{\text{Day 162 closed}}
\;=\; \underbrace{[\deg{=}n{-}1][T^n]\log(F_{-1}/F_0)}_{\text{direct from Prop 2 slice}}$$
holds symbolically for $n = 2, 3, \ldots, 8$.

## The honest gap: Route (B)

An INDEPENDENT closure of $[\deg{=}n{-}1][T^n]\log(F_{-1}/F_0)$ via Prop 2 +
Day 158's Riccati for $F_0$ was NOT achieved this session. Prop 2 gives
$F_{-1} = [p\, Q + (s+1)\delta_{n=0}]/(p+s+1)$ where
$Q = F_0 - (s+1)T F_0 - 2 T^2 F_0' + (s+1)\int F_0$.
Writing $\log(F_{-1}/F_0) = \log(p Q + (s+1)\delta_{n=0}) - \log(p+s+1) - \log F_0$
and taking the $\deg{=}n{-}1$ part in $(u_1, u_2)$ requires manipulating $\log F_0$
whose sub-sub-top layer $X^{(-1)}|_{u_3=0}$ is not itself in closed form. This is
the analogue of the Day 162 §5 Theorem B gap (proving the closed form for
$\bar D|_{E_3=0}$ from first principles).

Consequently, at present:
- **Missing Lemma (R)** = **PROVED conditional on Day 162 Theorem B**
  (equivalently, on either an independent Route (B) closure or an independent
  proof of the Day 162 $R^{(-1)}$ closed form).
- **C.5** stays `checked-sober`; the gap is now **precisely** the same as Day
  162's gap (proving Theorem B — the $\bar D|_{E_3=0} = TY^2[(q+1)^2 - E_1 T]/q^3$
  closed form).

## What has NOT changed

The Missing Lemma (R) status improves from Day 162's `checked-sober to n≤14` to
`PROVED given Route (A)`, but Route (A) alone is not enough — it needs to be
paired with Route (B) or Day 162's Theorem B to close R^{(-1)}. So the net
registry impact is:

- NEW node **`xi_1_closed_form_at_E3_zero`**: `proved`, role `premise`.
- NEW node **`xi_2_closed_form_at_E3_zero`**: `proved`, role `premise`.
- NEW node **`half_dd_u3_Xi_closed_form_at_u3_zero`**: `proved`, role `premise`.
- Missing Lemma (R): status `PROVED given Day 162 Theorem B`.
- C.5: stays `checked-sober`.

## Derivation of Route (A)

Setup: at $u_3 = 0$, the coefficient of $E_3^k$ in $\Xi$ (as a polynomial in
$(E_1, E_2, E_3)$) is $\xi_k = \xi_k(E_1, E_2)$. The relevant object is
$\xi_k|_{E_3=0}$ (the value at $E_3 = 0$).

### Step 1 — Chain rule at u_3 = 0

For any function $f(E_1, E_2, E_3)$ with $E_i = e_i(u_1, u_2, u_3)$:
$$\frac{\partial f}{\partial u_3} = f_{E_1} + (u_1 + u_2) f_{E_2} + (u_1 u_2) f_{E_3}.$$

At $u_3 = 0$: $u_1 + u_2 = E_1|_{u_3=0}$ and $u_1 u_2 = E_2|_{u_3=0}$, so
$$\partial_{u_3} f\big|_{u_3=0} = f_{E_1}(E_1, E_2, 0) + E_1\, f_{E_2}(E_1, E_2, 0) + E_2\, f_{E_3}(E_1, E_2, 0).$$
(Here $(E_1, E_2)$ are the slice's elementary symmetric polynomials.)

Applied to $\Xi = \sum_k E_3^k \xi_k$:
$$\partial_{u_3}\Xi\big|_{u_3=0} = \partial_{E_1}\xi_0 + E_1\,\partial_{E_2}\xi_0 + E_2\,\xi_1\big|_{E_3=0}.$$

### Step 2 — ξ_1 closed form (via Day 161 Thm 1)

Day 161 Thm 1: $\partial_{u_3}\Xi|_{u_3=0} = -\log q$. Solving:
$$\xi_1\big|_{E_3=0} = \frac{-\log q - \partial_{E_1}\xi_0 - E_1\,\partial_{E_2}\xi_0}{E_2}. \tag{A1}$$

Well-defined as a formal power series in $T$ because the RHS numerator is $O(E_2)$
as $E_2 \to 0$ (both $\log q$ at $E_2 = 0$ collapses and $\xi_0 \propto E_2$).

### Step 3 — [E_3^1] log W closed form (via Day 161 Thm 2)

Same chain rule applied to $\log\mathcal W$:
$$\partial_{u_3}\log\mathcal W\big|_{u_3=0} = \partial_{E_1}\log\mathcal W|_{E_3=0} + E_1\,\partial_{E_2}\log\mathcal W|_{E_3=0} + E_2\,[E_3^1]\log\mathcal W\big|_{E_3=0}.$$

Day 158: $\log\mathcal W|_{E_3=0} = \log(Y/(Tq)) = \log\phi - \log q$.
Day 161 Thm 2: $\partial_{u_3}\log\mathcal W|_{u_3=0} = T(q+R_1R_2)/q^3$. Solving:
$$[E_3^1]\log\mathcal W\big|_{E_3=0} = \frac{T(q+R_1R_2)/q^3 - \partial_{E_1}\log(Y/(Tq)) - E_1\,\partial_{E_2}\log(Y/(Tq))}{E_2}. \tag{A2}$$

### Step 4 — ξ_2 closed form (via (P1) recursion)

By (P1) [Day 152 Thm C]: $\log\mathcal W = \partial\Xi$ where
$\partial = \partial_{u_1} + \partial_{u_2} + \partial_{u_3}$. On symmetric-in-$(u_1, u_2, u_3)$
functions of $(E_1, E_2, E_3)$: $\partial = 3\,\partial_{E_1} + 2 E_1\,\partial_{E_2} + E_2\,\partial_{E_3}$.
Matching $[E_3^k]$:
$$[E_3^k]\log\mathcal W = 3\,\partial_{E_1}\xi_k + 2 E_1\,\partial_{E_2}\xi_k + (k+1) E_2\,\xi_{k+1}.$$

At $k=1$, $E_3=0$:
$$\xi_2\big|_{E_3=0} = \frac{[E_3^1]\log\mathcal W|_{E_3=0} - 3\,\partial_{E_1}\xi_1|_{E_3=0} - 2 E_1\,\partial_{E_2}\xi_1|_{E_3=0}}{2 E_2}. \tag{A3}$$

### Step 5 — ∂²_{u_3} Xi|_{u_3=0} chain rule

Applying the chain rule again:
$$\partial_{u_3}^2 f\big|_{u_3=0} = f_{E_1 E_1} + 2 E_1\, f_{E_1 E_2} + 2 E_2\, f_{E_1 E_3}
+ E_1^2\, f_{E_2 E_2} + 2 E_1 E_2\, f_{E_2 E_3} + E_2^2\, f_{E_3 E_3}$$
(all at $E_3 = 0$).

Applied to $\Xi = \sum_k E_3^k \xi_k$:
- $\Xi_{E_1 E_1}|_{E_3=0} = \partial_{E_1}^2\xi_0$; similarly for $E_1 E_2$, $E_2 E_2$.
- $\Xi_{E_1 E_3}|_{E_3=0} = \partial_{E_1}\xi_1$; $\Xi_{E_2 E_3}|_{E_3=0} = \partial_{E_2}\xi_1$.
- $\Xi_{E_3 E_3}|_{E_3=0} = 2\,\xi_2$.

Hence:
$$\tfrac{1}{2}\partial_{u_3}^2\Xi\big|_{u_3=0}
= \tfrac{1}{2}\Bigl[\partial_{E_1}^2\xi_0 + 2 E_1\,\partial_{E_1 E_2}\xi_0 + E_1^2\,\partial_{E_2}^2\xi_0
+ 2 E_2(\partial_{E_1}\xi_1 + E_1\,\partial_{E_2}\xi_1) + 2 E_2^2\,\xi_2\Bigr]. \tag{A}$$

Every piece has a closed form via (A1), (A2), (A3), and the Day 158 closed form for
$\xi_0$ (which follows from $T\partial_T\xi_0 = E_2 Y$).

$\square$ (Route A)

## Verification

Location: `/home/agent/projects/scratch/day167-closure/`.

- **`step1_prop3_vs_day162.py`** — verifies R^{(-1)} from Prop 3 == direct == Day 162
  closed form for $n\le 8$.
- **`step2_prop2_check.py`** — verifies Prop 2 (Clio) against raw `FP_coeffs` for $n\le 8$.
- **`step3_log_ratio_structure.py`** — computes $[\deg{=}n{-}1][T^n]\log(F_{-1}/F_0)$
  from raw `FP_coeffs` (Prop 2 slice).
- **`step5_B_route.py`** — verifies the degree-$n$ layer of $\log(F_{-1}/F_0)$ equals
  $\log q$ (Day 161 Thm 1 consequence).
- **`step7_fixed.py`** — verifies identity (A1) for $n \le 8$.
- **`step8_xi2_closure.py`** — verifies identities (A2), (A3), and the assembled
  $(1/2)\partial_{u_3}^2\Xi|_{u_3=0}$ closed form for $n \le 8$.
- **`step10_route_B_direct.py`** — verifies the final Prop 3 + Route A + Day 162
  triple identity, three-way agreement for $n \le 8$.

All checks pass.

## Discipline

- **[[feedback_check_convention_before_compute]]**: applied. All computations use raw
  `FP_coeffs` from `/home/agent/projects/scratch/day152/lib.py`. Sympy series
  computed independently.
- **[[feedback_top_layer_convention_swap]]**: caught myself with a two-hour bug from
  confusing $[u_3^k]$ (u-basis coefficient) with $[E_3^k]$ (E-basis coefficient) of
  symmetric polynomials. They differ! The (E-basis) extraction is the correct one
  for the ξ_k series; the (u-basis) extraction gives a mixed quantity via chain rule.
  Fix: `to_E_3var` conversion before extracting $[E_3^k]$.
- **[[feedback_pre_register_predictions]]**: the pre-registration is Day 162's
  $R^{(-1)}$ formula; both Route-A closed form and direct extraction reproduce it.
- **[[feedback_verify_scripts_implement_what_they_claim]]**: each script prints the
  identity it tests before checking.

## Route (B) obstruction — HONEST STALL

Route (B) attempts to close $[\deg{=}n{-}1][T^n]\log(F_{-1}/F_0)$ directly from
Prop 2 + Day 158. The obstruction:
- Prop 2 gives $F_{-1}$ as a rational-differential expression in $F_0$.
- $\log(F_{-1}/F_0)$ can be written as $\log(p Q + (s+1)\delta_{n=0}) - \log(p+s+1) - \log F_0$.
- Taking degree-$(n{-}1)$ in $(u_1, u_2)$: the $\log F_0$ term contributes via its
  sub-sub-top layer $X^{(-1)}|_{u_3=0}$, which is NOT in closed form (this is
  precisely the same object as the Day 162 §5 Theorem B target, up to a linear
  correction from the operator's degree-$(n{-}1)$ contributions).

So Route (B) is essentially isomorphic to Day 162 Theorem B. Not a distinct
attack surface.

## What would close C.5

Missing Lemma (R) equivalently is closed by proving ANY of:
1. **Day 162 Theorem B**: $\bar D|_{E_3=0} = TY^2[(q+1)^2 - E_1 T]/q^3$ (2-var closed
   form, currently checked-sober to $n \le 14$).
2. **Route (B)**: closed form for $[\deg{=}n{-}1][T^n]\log(F_{-1}/F_0)$ via Prop 2 +
   Day 158 (equivalent to (1) modulo Prop 3 + Route A).
3. A direct closed-form proof of Day 156 C.5 = Day 161 Thm 4.

All three are equivalent given today's work (Prop 3 + Route A). The most tractable
appears to be Day 162 Route (i): extend the sub-top $\nu$-system (proved as Theorem A
in Day 162 §2) to give $R^{(-1)}$ symbolically and match against Day 162's Theorem B
closed form.

## Registry updates

- **NEW** `xi_1_closed_form_at_E3_zero`: `proved`, role `premise`, file this document.
  Statement (A1). Derived from Day 161 Thm 1 + chain rule.
- **NEW** `xi_2_closed_form_at_E3_zero`: `proved`, role `premise`, file this document.
  Statement (A3). Derived from Day 161 Thm 2 + (P1) recursion + (A1).
- **NEW** `half_dd_u3_Xi_at_u3_zero`: `proved`, role `premise`, file this document.
  Statement (A). Derived from (A1), (A3), Day 158.
- `narayana-layer-d1-E3-zero` (C.5): STAYS `checked-sober`. The gap is now
  precisely Day 162 Theorem B (equivalently Route (B)).
- Missing Lemma (R): PROVED conditional on Day 162 Theorem B.

## Rule 11 scorecard

Route (A) closure is Rule 11 firing #8: unfold the $(E_3^k)$-coefficient recursion
from (P1), the chain rule (which is elementary), and Day 161's transverse derivatives.
No new theory imported.

## Queue

1. **Attack Day 162 Theorem B** via Day 162 Route (i) (sub-top $\nu$-system, Theorem A
   there is proved). This is the direct route to closing R^{(-1)} unconditionally.
2. **Try Route (B)** via Prop 2 with explicit $\log F_0$ expansion (may be equivalent
   to (1) modulo trivial manipulations).
3. Numerically extend Day 162 Theorem B verification from $n \le 14$ to $n \le 20$ as
   additional confidence.

## Time budget

Session 1h30min. Route A closed in 45 min; Route B recognized as equivalent to
Day 162 open gap in 30 min; writeup 15 min.

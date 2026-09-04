# Day 163 — Theorem B (closed form for $\bar D|_{E_3=0}$) — attempt

**Date:** 2026-09-04 (Day 163). **Status: HONEST STALL.** Two internal sub-lemmas are proved
cleanly (the $\theta(\theta+2)Y^2$ reduction, and the chain-rule relation between
$\partial_{u_3}$ and $\partial_{E_1}+E_1\partial_{E_2}+E_2\partial_{E_3}$ on the $u_3=0$ slice).
Route (i) (sub-top $\nu$-system extraction of $R^{(-1)}$) is *localised* to a specific missing
ingredient (see §4). Route (iv) (direct $\partial/\partial E_3$ via chain rule) is derived and
shown to be **equivalent to Day 161 Theorem 3** — hence adds no new information. Theorem B stays
**`checked-sober` to $n\le 14$**, not upgraded. The exact stall point and a recommended attack
are named in §5.

## 1. What is actually proved (unconditional, new today)

* **Lemma 1 ($\theta(\theta+2)Y^2$ reduction).** With $Y = T\phi(Y)$, $\phi = 1+E_1Y+E_2Y^2$,
  $q = 1-E_1T-2TE_2Y$, and $\theta = T\,d/dT$,
  $$\theta(\theta+2)\,Y^2 \;=\; \frac{2Y^2\bigl[(q+1)^2 - E_1T\bigr]}{q^3}. \tag{L1}$$
  Consequently, **the Lagrange form $[T^n]\bar D_n|_{E_3=0} = (n+1)[Y^{n-3}]\phi^{n-1}$ is EQUIVALENT
  to Theorem B $\bar D|_{E_3=0} = TY^2[(q+1)^2-E_1T]/q^3$** as algebraic identities in
  $\mathbb Q[E_1,E_2][[T]]$. Neither direction adds unproven content beyond the other.

* **Lemma 2 (chain-rule on the $u_3=0$ slice).** For any symmetric $F \in \mathbb Q[E_1,E_2,E_3][[T]]$
  (equivalently: a symmetric polynomial in $u_1, u_2, u_3$),
  $$\partial_{u_3}F\big|_{u_3=0} \;=\; \Bigl[\partial_{E_1}F + E_1\,\partial_{E_2}F + E_2\,\partial_{E_3}F\Bigr]_{E_3=0}, \tag{L2}$$
  where the RHS is evaluated at $E_3=0$ (which equals $u_3=0$ under symmetry). *Proof.* Directly
  from $E_1 = u_1+u_2+u_3$, $E_2 = u_1u_2+u_1u_3+u_2u_3$, $E_3 = u_1u_2u_3$: at $u_3=0$,
  $\partial E_1/\partial u_3 = 1$, $\partial E_2/\partial u_3 = u_1+u_2 = E_1$,
  $\partial E_3/\partial u_3 = u_1u_2 = E_2$. Chain rule. $\square$
  Verified for $\log\mathcal W$ to $n\le 6$ (`scratch/day163/step9_verify_chain_rule.py`).

* **Lemma 3 (chain-rule closure of Theorem B).** Combining L2 for $D = X^{(0)} - (1/2)\log\mathcal W$
  and $\bar D|_{E_3=0} = \partial D/\partial E_3|_{E_3=0}$ (since $D = E_3\bar D$),
  $$E_2\,\bar D|_{E_3=0} \;=\; R^{(-1)} - \tfrac{T}{2}\cdot\frac{q+R_1R_2}{q^3} \tag{L3}$$
  where $R^{(-1)} := \partial_{u_3}X^{(0)}|_{u_3=0}$. **This is Day 161 Theorem 3 re-derived via
  chain-rule.** No new information relative to Day 161. Route (iv) provides no independent leverage.

## 2. Proof of Lemma 1 (the reduction)

Proceed algebraically. First:
- $\theta Y = TY' = Y/q$ (Day 158 (Q2)).
- $\theta Y^2 = 2Y\theta Y = 2Y^2/q$.
- $(\theta + 2)Y^2 = 2Y^2/q + 2Y^2 = 2Y^2(q+1)/q$.
- $\theta[2Y^2(q+1)/q] = 2\theta(Y^2)\cdot(q+1)/q + 2Y^2 \cdot \theta[(q+1)/q]$.

Compute $\theta[(q+1)/q]$: since $\theta 1 = 0$,
$$\theta[(q+1)/q] = \theta[1 + 1/q] = -\theta q/q^2 = -\theta(q)/q^2.$$

Hence
$$\theta[2Y^2(q+1)/q] = \frac{4Y^2(q+1)}{q^2} - \frac{2Y^2\theta(q)}{q^2}
= \frac{2Y^2}{q^2}\bigl[2(q+1) - \theta(q)\bigr].$$

That is,
$$\theta(\theta+2)Y^2 = \frac{2Y^2}{q^2}\bigl[2(q+1) - \theta(q)\bigr]. \tag{*}$$

Now use $q^2 = (1-E_1T)^2 - 4E_2T^2$ (Day 158 (Q1)):
$$2q\,\theta(q) = T\frac{d}{dT}q^2 = -2E_1T(1-E_1T) - 8E_2T^2,$$
so
$$q\,\theta(q) = -E_1T(1-E_1T) - 4E_2T^2 = -E_1T + E_1^2T^2 - 4E_2T^2.$$

Compute $q\bigl[2(q+1) - \theta(q)\bigr] = 2q + 2q^2 - q\theta(q)$:
$$= 2q + 2q^2 - \bigl[-E_1T + E_1^2T^2 - 4E_2T^2\bigr]
= 2q + 2q^2 + E_1T - E_1^2T^2 + 4E_2T^2.$$

Substitute $2q^2 = 2[(1-E_1T)^2 - 4E_2T^2]$:
$$= 2q + 2(1-E_1T)^2 - 8E_2T^2 + E_1T - E_1^2T^2 + 4E_2T^2$$
$$= 2q + 2 - 4E_1T + 2E_1^2T^2 - 8E_2T^2 + E_1T - E_1^2T^2 + 4E_2T^2$$
$$= 2q + 2 - 3E_1T + E_1^2T^2 - 4E_2T^2.$$

Compare to $(q+1)^2 - E_1T = q^2 + 2q + 1 - E_1T$. Using $q^2 = 1 - 2E_1T + E_1^2T^2 - 4E_2T^2$:
$$(q+1)^2 - E_1T = 1 - 2E_1T + E_1^2T^2 - 4E_2T^2 + 2q + 1 - E_1T
= 2q + 2 - 3E_1T + E_1^2T^2 - 4E_2T^2.$$

**Match.** Hence $q[2(q+1) - \theta(q)] = (q+1)^2 - E_1T$, and substituting into (*):
$$\theta(\theta+2)Y^2 = \frac{2Y^2 \cdot [(q+1)^2 - E_1T]/q}{q^2} = \frac{2Y^2[(q+1)^2 - E_1T]}{q^3}.
\qquad\square$$

Verified numerically to $n\le 14$ (`scratch/day163/step2_verify_theta_id.py`; both
series-by-series match and symbolic derivation confirm zero residual).

**Corollary.** The Lagrange form $[T^n]\bar D_n|_{E_3=0} = (n+1)[Y^{n-3}]\phi^{n-1}$ is equivalent to
Theorem B. *Proof.* Given the Lagrange form, extended Lagrange–Bürmann gives
$(n+1)[Y^{n-3}]\phi^{n-1} = (n+1)(n-1)/2\cdot[T^{n-1}]Y^2$. Summing $\sum_n T^n (n+1)(n-1)/2\cdot[T^{n-1}]Y^2 = (T/2)\theta(\theta+2)Y^2$
and applying Lemma 1: $\bar D|_{E_3=0} = TY^2[(q+1)^2-E_1T]/q^3$. Conversely, given
Theorem B, expand $TY^2[(q+1)^2-E_1T]/q^3$ as a series in $T$ and match to $(n+1)[Y^{n-3}]\phi^{n-1}$
via Lemma 1 reversed. $\square$

## 3. Route (i) diagnostic — sub-top $\nu$-system at $u_3 = 0$

Day 162 Theorem A (sub-top $\nu$-system) gives $\mu_i = \ell^{\rm top}_0(\lambda_i)$ as solutions of
a linear system. At $u_3 = 0$, this specialises cleanly.

**Sub-lemma (sub-top $\nu$-system at $u_3 = 0$).** At $u_3 = 0$, on the diagonal $t_i = T$:
$$\mu_2^{(0)} := \mu_2|_{u_3=0} = 0, \qquad \mu_3^{(0)} := \mu_3|_{u_3=0} = 0,$$
$$q\,\mu_1^{(0)} = T\,D_1(L_1)|_{u_3=0}$$
where $L_1|_{u_3=0} = TE_2\phi(Y) = E_2 Y$, so $D_1(L_1)$ can be computed explicitly.

*Sketch.* From Day 162 Thm A:
\begin{align*}
\mu_2 &= t_2[\nu_1\tilde d_3 + \tilde d_1\nu_3 + D_1(L_2+L_3)]\\
\mu_3 &= t_3[\nu_2\tilde d_3 + \tilde d_2\nu_3 + D_2(L_2+L_3)]
\end{align*}
At $u_3 = 0$: $\nu_3 = 0$, $L_2 = L_3 = 0$. So $\mu_2 = t_2 \nu_1 \tilde d_3$, $\mu_3 = t_3 \nu_2 \tilde d_3$
where $\tilde d_3 = \mu_2 + \mu_3$. Adding: $\tilde d_3 = t_2 \nu_1 \tilde d_3 + t_3 \nu_2 \tilde d_3$,
so $\tilde d_3(1 - t_2\nu_1 - t_3\nu_2) = 0$. On the diagonal $t_i = T$, the factor becomes
$1 - T(\nu_1 + \nu_2) = 1 - T(E_1 + 2E_2Y) = q$, a unit. Hence $\tilde d_3 = 0$, giving
$\mu_2 = \mu_3 = 0$ at $u_3 = 0$.

Then $\mu_1 = t_1[\nu_1\tilde d_2 + \tilde d_1\nu_2 + D_1(L_1+L_3)]$ with $\tilde d_2 = \mu_1 + \mu_3 = \mu_1$,
$\tilde d_1 = \mu_1 + \mu_2 = \mu_1$, $L_3 = 0$:
$$\mu_1 = t_1(\nu_1 + \nu_2)\mu_1 + t_1 D_1(L_1),$$
so on diagonal $q \mu_1 = T D_1(L_1)$. $\square$

Verified numerically (`scratch/day163/step4_solve_subtop_at_u3_zero.py`): $\mu_2|_{u_3=0} = \mu_3|_{u_3=0} = 0$
to $|t|\le 4$; $\mu_1|_{u_3=0}$ matches raw extraction.

**However, this alone does NOT give $R^{(-1)} = \partial_{u_3}X^{(0)}|_{u_3=0}$.** From Day 152 Fact I:
$$\log F_P = \varrho(S) + \varrho\log(\mathcal R/V(u)),\qquad \mathcal R = e^{-S}V(M)e^{S}.$$
Applying $\ell^{\rm top}_0$:
$$X^{(0)} = \varrho\,\ell^{\rm top}_0(S) + \varrho\,\log\ell^{\rm top}_0(\mathcal R/V(u)) = \varrho\,\ell^{\rm top}_0(S) + \varrho\,\log\bigl[\ell^{\rm top}_3(\mathcal R)/V(u)\bigr].$$

The first term is accessible from the sub-top $\nu$-system: $\theta \ell^{\rm top}_0(S) = \sum_i \mu_i$.

The second term requires computing $\ell^{\rm top}_3(\mathcal R) = \ell^{\rm top}_3(e^{-S}V(M)e^S)$, which is a
non-trivial *residual* involving operator conjugation. Numerically
(`scratch/day163/step5_check_gap.py`), the gap between $\theta R^{(-1)}$ (extracted from raw $F_P$)
and $\varrho(\sum\partial_{u_3}\mu_i)|_{u_3=0}$ is NON-ZERO:
$$\theta R^{(-1)} - \varrho\sum_i \partial_{u_3}\mu_i|_{u_3=0} = 16T^4u_1^3 + 140T^4u_1^2u_2 + \ldots$$
This is precisely the transverse derivative of the "residual" $\varrho\log[\ell^{\rm top}_3(\mathcal R)/V(u)]$.

**Route (i) stall point.** Route (i) reduces the problem to computing
$\partial_{u_3}[\varrho\log(\ell^{\rm top}_3(\mathcal R)/V(u))]|_{u_3=0}$ in closed form. The naive
approach (unfolding $\mathcal R = e^{-S}V(M)e^S$) requires operator gymnastics with $S$ (whose closed
form on the top layer is known via $\nu$'s, but the sub-layers are not clean).

## 4. Route (iv) reduces to Day 161 Theorem 3 (no new information)

**Setup.** $D = E_3 \bar D$; $\bar D$ is a polynomial in $E_1, E_2, E_3$, so
$\bar D|_{E_3=0} = \partial D/\partial E_3|_{E_3=0}$.

By Lemma 2 (chain rule at $u_3 = 0$):
$$\partial_{u_3}D|_{u_3=0} = \bigl[\partial_{E_1}D + E_1\partial_{E_2}D + E_2\partial_{E_3}D\bigr]_{E_3=0}.$$

Now $D|_{E_3=0} = 0$ (since $D \in E_3\mathbb Q[E][[T]]$), so $\partial_{E_1}D|_{E_3=0} = 0$ and
$\partial_{E_2}D|_{E_3=0} = 0$. Hence $\partial_{u_3}D|_{u_3=0} = E_2\bar D|_{E_3=0}$.

On the other hand,
$$\partial_{u_3}D|_{u_3=0} = \partial_{u_3}X^{(0)}|_{u_3=0} - \tfrac{1}{2}\partial_{u_3}\log\mathcal W|_{u_3=0}
= R^{(-1)} - \tfrac{T}{2}\cdot\frac{q+R_1R_2}{q^3}$$
by Day 161 Theorem 2 (proved). Combining:
$$E_2\bar D|_{E_3=0} = R^{(-1)} - \tfrac{T}{2}\cdot\frac{q+R_1R_2}{q^3}. \tag{L3}$$

**This is exactly Day 161 Theorem 3.** So the chain-rule derivation reproduces Day 161's identity;
it does not give an independent proof of Theorem B. **The gap remains: prove the closed form for
$R^{(-1)}$**.

## 5. Where the argument stalls, and recommended next steps

**Precise stall point.** All routes converge to the same missing lemma:

> **Missing Lemma (equivalent to Theorem B).** As formal power series in $T$ with coefficients in
> $\mathbb Q[E_1, E_2]$:
> $$R^{(-1)} \;=\; \frac{T}{q^3}\Bigl[E_2Y^2\bigl((q+1)^2 - E_1T\bigr) + \tfrac{1}{2}(q + R_1R_2)\Bigr] \tag{R}$$
> where $R^{(-1)} := \partial_{u_3} X^{(0)}|_{u_3=0}$, $X^{(0)} := \ell^{\rm top}_0(\log F_P)$.

Both Route (i) and Route (iv) reduce to this. Route (ii) (direct ODE for $F_1$ coupled to $F_0$)
is a separate approach, but numerical exploration (`scratch/day163/step7_F1_ODE_fit.py`) does not
suggest a clean ansatz $L_A F_1 = P(T)F_0 + Q(T)F_0'$ with low-degree $P, Q$.

**Recommended next steps for Rick:**

1. **Extend Route (i)** by explicitly computing $\ell^{\rm top}_3(\mathcal R)/V(u)$ from
   $\mathcal R = e^{-S}V(M)e^S$. Since $V(M) = \prod_{i<j}(M_i - M_j)$ and $\ell^{\rm top}_1(M_i) = \nu_i$
   (Day 152 §4 Step 1), $\ell^{\rm top}_3(V(M)) = V(\nu)$ is a start. Then compute conjugation
   corrections at weight 3.

2. **Attack via generating function residue.** The Lagrange form
   $(n+1)[Y^{n-3}]\phi^{n-1} = \frac{n+1}{2\pi i}\oint \phi^{n-1}\frac{dY}{Y^{n-2}}$ suggests a
   contour integral representation. Combined with the definition of $\bar D_n$ via extraction of
   $E_3$-coefficient from $D = X^{(0)} - (1/2)\log\mathcal W$, one might get a residue-calculus
   proof.

3. **Direct check via Day 152 Theorem C.** $\log\mathcal W = -\sum_i\log\rho_i$ is a *closed form*
   in $\nu_i$'s. Compute $\partial(\log\mathcal W)/\partial E_3$ at $E_3=0$ via the $\nu$-system.
   Then attempt the same for $X^{(0)}$ using an operator-theoretic version of Day 152 Fact I.

4. **Layer $d = 2$ hunt for structural analog.** If layer $d=2$ at $E_3=0$ has a similar closed
   form $\sim T^k Y^m/(1-E_2Y^2)^\ell$, the pattern would suggest a unified generating function
   for all layers, from which layer $d = 1$ (Theorem B) could be derived.

## 6. Numerical certification (this session)

All scripts in `/home/agent/projects/scratch/day163/`:

| Claim | Script | Verdict |
|---|---|---|
| Rick's $f(Y)$ form equals $TY^2[(q+1)^2-E_1T]/q^3$ | `step1_verify_fY.py` | ✓ $n\le 14$ |
| Lemma 1: $\theta(\theta+2)Y^2 = 2Y^2[(q+1)^2-E_1T]/q^3$ | `step2_verify_theta_id.py` | ✓ series ($n\le 14$) + symbolic |
| $R^{(-1)}$ closed form matches raw | `step3_subtop_at_u3_zero.py` | ✓ $n\le 8$ |
| Sub-top $\nu$-system at $u_3 = 0$: $\mu_2 = \mu_3 = 0$ | `step4_solve_subtop_at_u3_zero.py` | ✓ |
| Non-zero gap (Route (i) residual) | `step5_check_gap.py` | ✓ diagnostics |
| Lemma 2 (chain rule) for $\log\mathcal W$ | `step9_verify_chain_rule.py` | ✓ $n\le 6$ |

## 7. Registry updates

* `bar-D-closed-form-E3-zero`: **stays `checked-sober` to $n\le 14$**. No upgrade.
* `narayana-layer-d1-E3-zero` (C.5): **stays `computed`**. Gap still localised to Theorem B.
* NEW node **`theta-theta+2-Y2-reduction`** (Lemma 1): **`proved`**, role `premise`.
  Statement: $\theta(\theta+2)Y^2 = 2Y^2[(q+1)^2-E_1T]/q^3$. File: this document, §2.
* NEW node **`u3-chain-rule-lemma`** (Lemma 2): **`proved`**, role `premise`.
  Statement: On symmetric polynomials at $u_3=0$,
  $\partial_{u_3} = \partial_{E_1} + E_1\partial_{E_2} + E_2\partial_{E_3}$. File: this document, §1.
* Route (iv) status: **derived and proved equivalent to Day 161 Thm 3** — no independent progress.

## 8. Discipline scorecard

* **[[feedback_pre_register_predictions]]**: Applied. Pre-registered numerical fingerprint
  ($R^{(-1)}$ coefficients from PROVE.md) matches computed $\bar D|_{E_3=0}$ via Lemma 3 to $n\le 8$.
* **[[feedback_check_convention_before_compute]]**: Applied. All computations use
  `FP_coeffs` from `scratch/day152/lib.py` (the true library object), not a paraphrase.
* **[[feedback_operator_respects_slice]]**: Applied. Chain-rule closure at $u_3=0$ correctly
  handles the transverse $\partial_{u_3}$ direction; Lemma 2 makes explicit which $E$-derivatives
  it involves.
* **Honest stall**: This session honestly reports failure to close Theorem B, with precise
  identification of the remaining gap (Missing Lemma (R)) and three recommended attacks.

## 9. Summary for one-line handoff

Theorem B (Day 162's closed form for $\bar D|_{E_3=0}$) is not proved this session. The full
route reduces to a single unproved series identity for $R^{(-1)}$ (Missing Lemma (R) in §5).
Two new proved sub-lemmas — the $\theta(\theta+2)Y^2$ reduction and the chain-rule at $u_3=0$ —
tighten the statement but do not close the gap. Route (i) is localised to computing the
residual $\ell^{\rm top}_3(\mathcal R)/V(u)$ (an operator-theoretic conjugation term).

## 10. Queue for Day 164+

1. **Direct attack on $\ell^{\rm top}_3(\mathcal R)/V(u)$**: expand $\mathcal R = e^{-S}V(M)e^S$
   at weight 3, using $\ell^{\rm top}_1(M_i) = \nu_i$ and $\ell^{\rm top}_1(S) = \tilde\Xi$.
2. **Alternative: Route (iii) via residue formula** for $(n+1)[Y^{n-3}]\phi^{n-1}$ directly.
3. **Layer $d=2$ E-positive form search** — if a pattern emerges, it may reveal the generating
   principle for $d=1$.

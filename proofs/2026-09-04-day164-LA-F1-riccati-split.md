# Day 164 — Riccati split of $L_A F_1$; top-layer closed form; STALL on sub-top

**Date:** 2026-09-04 (Day 164). **Status: HONEST STALL WITH TWO REAL RESULTS.**

Two positive results, one clean stall:

* **Result 1 (Rule-11 firing).** From Day 158's ODE $L_A F_0 = 0$, we derive the identity
  $$L_A F_1 \;=\; F_0\cdot\Bigl[T^2 Q'' + \bigl(2 T^2 G + (E_1+3)T - 1\bigr) Q'\Bigr] \tag{L1}$$
  where $Q := F_1/F_0$, $G := F_0'/F_0$, $L_A = T^2\partial_T^2 + [(E_1+3)T-1]\partial_T + (1+E_1+E_2)$
  (Day 158 Prop A). This is a **first-order ODE for $Q$**, once $L_A F_1/F_0$ is known independently.

* **Result 2 (Top-layer closed form).** Splitting (L1) at $u$-weight 1:
  $$\ell^{\rm top}_1(L_A F_1/F_0) \;=\; q', \qquad q^2 = (1-E_1T)^2 - 4 E_2 T^2. \tag{L2}$$
  Verified analytically (proof in §3) and numerically to $n\le 7$
  (`scratch/day164/verify_top_layer.py`). This is a NEW derivation route to Day 161 Thm 1
  (via the $F_1$ Riccati), independent of Day 152's $\nu$-system.

* **STALL:** Splitting (L1) at $u$-weight 0 gives an ODE
  $$q\cdot\partial_T R^{(-1)} \;=\; \theta(\theta-1)\Xi' - (2T^2K + 3T)\frac{q'}{q} - \Sigma_0 \tag{L3}$$
  where $\Sigma_0 := \ell^{\rm top}_0(L_A F_1/F_0)$ and $\Xi' := \partial_{u_3}\Xi|_{u_3=0} = -\log q$.
  If $\Sigma_0$ had a known closed form, (L3) would determine $R^{(-1)}$ uniquely
  (via integration; the initial condition $R^{(-1)}(0) = 0$ fixes the constant).
  **We failed to find a closed form for $\Sigma_0$**; numerically it has E-coefficients
  that don't obviously admit a rational expression in $Y, q, E_1, E_2, T$ (§5).

**Registry updates:**
* `LA-F1-Riccati-split` (NEW): **`proved`** (Result 1, §2). Role: `premise`.
* `LA-F1-top-layer-q-prime` (NEW): **`proved`** analytically (Result 2, §3).
  Role: `premise`. This is an alternative derivation of Day 161 Thm 1.
* `LA-F1-sub-top-Sigma-0` (NEW): **`computed`** to $n\le 7$; no closed form found.
  Role: `attempt` — required for closing Missing Lemma (R).
* `R-minus-one-closed-form`: stays **`checked-sober` to $n\le 14$**. Not upgraded.
* `bar-D-closed-form-E3-zero` (=Theorem B): stays **`checked-sober` to $n\le 14$**.
* `narayana-layer-d1-E3-zero` (C.5): stays **`computed`**.

**Backup 2 result:** Tested layered Lagrange pattern $\bar D_{(k)} := [E_3^k] D$ for $k=2$
(§4). $\bar D_{(2)}$'s numerical data does NOT fit a Lagrange-Bürmann form: half-integer
coefficient $\tfrac{4617}{2}\cdot E_1^2$ at $n=8$ rules out $\bar D_{(2)} = c_n [Y^{n-\alpha}]\phi^{n-\beta}$
with rational $c_n$ times integer multinomial. Layered generating function story does NOT extend
cleanly to $k \ge 2$.

---

## 1. Setup

Following Day 158/161/162. Let $F_P = \Psi^+(e^{Te_2}) \in \mathbb Q[E_1,E_2,E_3][[T]]$
(Day 148 §1). Set
$$F_0(T;u_1,u_2) := F_P|_{u_3=0}, \qquad F_1(T;u_1,u_2) := \partial_{u_3} F_P|_{u_3=0},$$
$$Q := F_1/F_0, \qquad G := F_0'/F_0.$$

Day 158 Prop A: $L_A F_0 = 0$ where
$$L_A = T^2\partial_T^2 + [(E_1+3)T - 1]\partial_T + (1 + E_1 + E_2). \tag{A}$$

Day 158 (Q1)/(Q2)/H/K:
* $q^2 = (1-E_1T)^2 - 4 E_2 T^2$.
* $Y = T\phi(Y)$, $\phi = 1 + E_1 Y + E_2 Y^2$.
* $\theta Y = T Y' = Y/q$.
* $H := \ell^{\rm top}_2(G) = E_2 Y/T$ (Day 158 Lemma 4.1).
* $K := \ell^{\rm top}_1(G) = [E_2 Y(2q+1) + E_1 q]/q^2$ (Day 158 (K)).

Day 161 Thm 1 (PROVED): $\partial_{u_3}\Xi|_{u_3=0} = -\log q$. Call this $\Xi'$.

Weight-decomposition: at $[T^n]$, degree-in-$u$ is bounded by $n+1$ for $\log F_P$
(Day 149 Fact II). So $Q_n$ has degree $\le n$: top layer $\ell^{\rm top}_0(Q) = \Xi'$,
sub-top $\ell^{\rm top}_{-1}(Q) = R^{(-1)}$.

## 2. Proof of Result 1 (Identity L1)

Since $L_A F_0 = 0$, expand $L_A(F_0 Q) = L_A(F_1)$:
$$(F_0 Q)' = F_0' Q + F_0 Q', \qquad (F_0 Q)'' = F_0'' Q + 2 F_0' Q' + F_0 Q''.$$

Substitute into $L_A(F_0 Q)$:
\begin{align*}
L_A(F_0 Q) &= T^2(F_0'' Q + 2 F_0' Q' + F_0 Q'') + [(E_1+3)T - 1](F_0' Q + F_0 Q') + (1+E_1+E_2)F_0 Q\\
&= \bigl[T^2 F_0'' + ((E_1+3)T - 1)F_0' + (1+E_1+E_2)F_0\bigr]Q + T^2(2 F_0' Q' + F_0 Q'') + ((E_1+3)T - 1)F_0 Q'\\
&= \bigl[L_A F_0\bigr]Q + F_0\bigl[T^2 Q'' + (2 T^2 G + (E_1+3)T - 1) Q'\bigr]\\
&= F_0\bigl[T^2 Q'' + (2 T^2 G + (E_1+3)T - 1) Q'\bigr].
\end{align*}
This is (L1). $\square$

## 3. Proof of Result 2 (Top-layer identity, L2)

Extract $\ell^{\rm top}_1$ of both sides of $L_A F_1/F_0 = T^2 Q'' + (2 T^2 G + (E_1+3)T - 1)Q'$.

**Weight bookkeeping** (Day 152 Lemma 1.1):
* $Q$ has wt $\le 0$; top layer $\ell^{\rm top}_0(Q) = \Xi' = -\log q$.
* Hence $Q'$ has wt $\le 1$; top layer $\ell^{\rm top}_1(Q') = \partial_T\Xi' = -q'/q$.
* $G$ has wt $\le 2$; top layer $\ell^{\rm top}_2(G) = H = E_2 Y/T$.
* $T$ has wt $-1$; $T^2$ has wt $-2$; $E_1 T$ has wt $0$; $-1$ has wt $0$; $E_1$ has wt $1$.

**Term-by-term weight-1 extraction:**

* $T^2 Q''$: $Q''$ has wt $\le 2$, but $T^2$ contributes $-2$: wt $\le 0$. So $\ell^{\rm top}_1(T^2 Q'') = 0$.

* $2 T^2 G Q'$: wt $\le -2 + 2 + 1 = 1$. Top-1 part = $T^2\cdot\ell^{\rm top}_2(G)\cdot\ell^{\rm top}_1(Q') = T^2\cdot(E_2 Y/T)\cdot(-q'/q) = -T E_2 Y q'/q$.

* $((E_1+3)T - 1)Q'$: Split into $E_1 T\cdot Q'$ (wt $\le 0+1 = 1$), $3T\cdot Q'$ (wt $\le -1+1 = 0$), and $-Q'$ (wt $\le 1$):
  - $\ell^{\rm top}_1(E_1 T\cdot Q') = E_1 T\cdot\ell^{\rm top}_1(Q') = -E_1 T q'/q$.
  - $\ell^{\rm top}_1(3T\cdot Q') = 0$ (since $3T\cdot Q'$ has wt $\le 0$).
  - $\ell^{\rm top}_1(-Q') = -\ell^{\rm top}_1(Q') = q'/q$.

**Sum:**
$$\ell^{\rm top}_1(L_A F_1/F_0) = 0 + (-2T E_2 Y q'/q) + (-E_1 T q'/q + 0 + q'/q) = \frac{q'}{q}\bigl(1 - E_1 T - 2 T E_2 Y\bigr) = \frac{q'}{q}\cdot q = q'.$$

$\square$

**Consistency check with Day 161 Thm 1.** Since $\ell^{\rm top}_1(Q') = -q'/q = \partial_T(-\log q)$
and $\ell^{\rm top}_0(Q) = -\log q$, the derivation is INTERNALLY CONSISTENT with Day 161 Thm 1:
this section provides a NEW proof route via the $F_1$ Riccati (bypasses the $\nu$-system).

**Numerical certification** (`scratch/day164/verify_top_layer.py`): match to $n\le 7$
(equality of $\ell^{\rm top}_1(L_A F_1/F_0)_n$ with $q'_n$ as $u$-polys, exact `sympy` comparison).

## 4. Result 3 (Backup 2 refutation)

Computed $\bar D_{(k)} := [E_3^k] D$ where $D = X^{(0)} - (1/2)\log\mathcal W$ (Day 156). Results:
* $\bar D_{(1)} = \bar D|_{E_3=0}$: matches Day 162 Theorem B's Lagrange form (verified $n\le 9$).
* $\bar D_{(2)}$: first nontrivial at $n=6$ (value $72$), then $540 E_1$ at $n=7$, then
  $\tfrac{4617}{2} E_1^2 + 894 E_2$ at $n=8$. The presence of $\tfrac{4617}{2}$ (half-integer)
  rules out an integer-multinomial Lagrange form
  $\bar D_{(2)} = c_n [Y^{n-\alpha}]\phi^{n-\beta}$ with $c_n \in \mathbb Q$.

**Consequence:** the "layered generating function" idea (bar_D_k-lagrange-pattern in
PROVE.md's Backup 2) is REFUTED for $k=2$. This closes off a Day 164+ backup route.

Script: `scratch/day164/bar_D_k1_test.py`.

## 5. STALL: sub-top of $L_A F_1/F_0$ has no clean closed form

Splitting (L1) at $u$-weight 0 gives, similarly to §3:

Let $\Xi' = -\log q$, and $R^{(-1)} = \ell^{\rm top}_{-1}(Q)$. Then:

* $\ell^{\rm top}_0(T^2 Q'') = T^2\ell^{\rm top}_0(Q'') = T^2\partial_T^2\ell^{\rm top}_0(Q) = T^2 \Xi''$
  (using that $\ell^{\rm top}_0(Q_n)$ is the top layer; $[T^n](T^2 Q'') = n(n-1) Q_n$ has same top as $n(n-1) \Xi'_n$).
  As an operator: $T^2\partial_T^2 = \theta(\theta-1)$. So this term = $\theta(\theta-1)\Xi'$.

* $\ell^{\rm top}_0(2T^2 G Q')$: two contributions summing to wt 0:
  - $\ell^{\rm top}_2(G)\cdot\ell^{\rm top}_{-1}(T^2 Q') = H\cdot T^2 \partial_T R^{(-1)}\cdot\text{(scale)}$;
  - $\ell^{\rm top}_1(G)\cdot\ell^{\rm top}_0(T^2 Q') = K\cdot T^2 \partial_T \Xi'$.

  After careful bookkeeping: $\ell^{\rm top}_0(2T^2 G Q') = 2 T^2 H \partial_T R^{(-1)} + 2 T^2 K \partial_T \Xi'$.

* $\ell^{\rm top}_0(((E_1+3)T - 1)Q')$: 
  - $E_1 T\cdot\ell^{\rm top}_0(Q') = E_1 T\cdot\partial_T R^{(-1)}$.
  - $3T\cdot\ell^{\rm top}_1(Q') = 3T\cdot\partial_T \Xi'$.
  - $-\ell^{\rm top}_0(Q') = -\partial_T R^{(-1)}$.
  Sum: $(E_1 T - 1)\partial_T R^{(-1)} + 3T \partial_T \Xi'$.

Total:
$$\Sigma_0 = \theta(\theta-1)\Xi' + \bigl[2T^2 H + E_1 T - 1\bigr]\partial_T R^{(-1)} + \bigl[2T^2 K + 3T\bigr]\partial_T\Xi'.$$

Now, using $2T^2 H = 2 T E_2 Y$ and $q = 1 - E_1 T - 2 T E_2 Y$:
$$2 T^2 H + E_1 T - 1 = 2 T E_2 Y + E_1 T - 1 = -(1 - E_1 T - 2 T E_2 Y) = -q.$$

Substituting and rearranging for $\partial_T R^{(-1)}$:
$$\boxed{q\cdot \partial_T R^{(-1)} \;=\; \theta(\theta-1)\Xi' + \bigl[2T^2 K + 3T\bigr]\partial_T\Xi' - \Sigma_0.}\tag{L3}$$

**[SIGN CORRECTION 2026-09-04 (Day 165)]** The boxed (L3) originally had a sign error
on the middle term (was $-[2T^2K+3T]\partial_T\Xi'$; correct is $+$). The unboxed sum on
the immediately preceding line ($\Sigma_0 = \theta(\theta-1)\Xi' + [2T^2H+E_1T-1]\partial_T R^{(-1)} + [2T^2K+3T]\partial_T\Xi'$)
is correct; the box transcription flipped the sign. Verified by 3 independent Day-165
numerical checks. Caught by two Day-165 compute agents (`close-Rminus1/step2*.py`,
`sigma0-proof/step14*.py`).

**This is the target ODE.** If we knew $\Sigma_0$ in closed form, we could integrate to get $R^{(-1)}$.

### 5.1 Numerical data for $\Sigma_0$

Extracted from `scratch/day164/layers_of_LAF1_over_F0.py` (in $E$-basis):
| $n$ | $-\Sigma_0(n)$ |
|---|---|
| 0 | $1 + E_1$ |
| 1 | $1 + E_1 + 4 E_2$ |
| 2 | $3 + 4 E_1 + E_1^2 + 15 E_2 + 6 E_1 E_2$ |
| 3 | $13 + 21 E_1 + 9 E_1^2 + E_1^3 + 81 E_2 + 57 E_1 E_2 + 8 E_1^2 E_2 + 8 E_2^2$ |
| 4 | (see log; coefficients: $E_1^4$: 1, $E_1^2 E_2$: 138, $E_2^2$: 107, plus lower-degree terms) |

The $-E_1^n$ column gives $\sum -E_1^n T^n = -1/(1-E_1T)$ — clean.

The $E_1^{n-2} E_2$-coefficient is $-n(n-1)(4n+7)/2$ (fit by cubic to n=2..7 data);
generating function: $-3 E_2 T^2 (5 - E_1 T)/(1 - E_1 T)^4$.

But higher $E_2$-power coefficients do NOT admit clean rational fits. Specifically the
$E_1^{n-4} E_2^2$ coefficients at $n=4,5,6,7$ are $-107, -631, -2181, -5761$ with
third differences non-constant, so no cubic in $n$.

**No closed form for $\Sigma_0$ found this session.**

### 5.2 Why the operator-attack analog also fails here

One might hope that computing $L_A F_P$ (in 3 vars) and then dividing by $F_P$ would
give a clean 3-variable ratio $L_A F_P/F_P = -u_3 \cdot G(u_1,u_2,u_3,T)$ with $G$ nice.
Then $L_A F_1/F_0 = -G|_{u_3=0}$. Numerical exploration
(`scratch/day164/LA_FP_over_FP.py`) confirms $L_A F_P/F_P$ vanishes at $u_3=0$, but the
$u_3^1$ coefficient is a complicated non-symmetric polynomial in $u_1, u_2$ with no
clean rational-function structure.

## 6. What was NOT accomplished

* **Missing Lemma (R)** (PROVE.md §"Target identity"): still open, still `checked-sober`.
* **Theorem B** (Day 162): still `checked-sober` to $n\le 14$.
* **C.5** (Day 156): still `computed`. FPSAC §5 cannot yet be upgraded.

**Precise stall point:** Find a closed form for
$$\Sigma_0 \;=\; \ell^{\rm top}_0(L_A F_1/F_0).$$

## 7. Recommended queue for Day 165+

1. **Direct computation of $F_1$ mod $F_0$** via umbral operators. If a closed form for
   $F_1$ (analogous to Day 158's $F_0 = \sum(T^k/k!)A_k(u_1)A_k(u_2)$) can be found,
   $L_A F_1$ follows directly.

2. **Route (i) from Day 163**: computing $\ell^{\rm top}_3(\mathcal R)/V(u)$ via BCH.
   Not tried this session (deferred; requires operator gymnastics).

3. **Route (iii)** (from Day 162 §5, sketched): direct residue proof via
   $Y = T/\phi(Y)$ change of variable in the Lagrange integral. The extended L-B
   template (feedback_lagrange_burmann_extended) may apply.

4. **Sanity check on operator-annihilating candidates**: search for a 3-variable
   operator $L_A^{(3)}$ (possibly nonlinear in $u_3$) that annihilates $F_P$; if found,
   $L_A F_1 = -(\partial_{u_3} L_A^{(3)}) F_0$ would give closed form.

## 8. Discipline scorecard

* **[[feedback_check_convention_before_compute]]**: Applied. All computations use
  `FP_coeffs` from `scratch/day152/lib.py` (the true library object).
* **[[feedback_pre_register_predictions]]**: Applied to the top-layer identity.
  Predicted from numerical pattern ($-2(n+1)E_2 Y_n$ shape), then verified analytically
  in §3.
* **[[feedback_verify_scripts_implement_what_they_claim]]**: Applied. Caught a bug in
  my initial $Y$-recursion in `compute_Q_layers.py`; corrected in `verify_top_layer.py`.
* **Rule 11 firing #5**: Unfolding $L_A F_0 = 0$ (an existing theorem) into an equation
  for $F_1$'s log-derivative gives the Riccati split. This is the same style as
  Day 158's original derivation, one order up.
* **Honest stall**: The sub-top identity requires an independent closed form for
  $\Sigma_0$ which we did not find. Reported precisely, not hidden.

## 9. One-line handoff

Day 164's Riccati split gives a first-order ODE for $R^{(-1)}$ that would close
Missing Lemma (R) if the source term $\Sigma_0 = \ell^{\rm top}_0(L_A F_1/F_0)$ had a
known closed form; verified top-layer via new route ($=q'$); found no closed form for
$\Sigma_0$; C.5 stays `computed`.

## 10. Scripts

Location: `/home/agent/projects/scratch/day164/`.

| Script | Purpose | Verdict |
|---|---|---|
| `bar_D_k1_test.py` | Backup 2 test | ✓ $\bar D_{(1)}$ matches Lagrange to $n\le 9$; $\bar D_{(2)}$ has half-integer, no clean pattern |
| `compute_Q_layers.py` | Extract Q's layers (buggy $Y$; see verify_top_layer.py) | superseded |
| `compute_LA_F1.py` | Verify identity (L1) numerically | ✓ diff = 0 to $n\le 7$ |
| `layers_of_LAF1_over_F0.py` | Extract all layers of $L_A F_1/F_0$ | ✓ up to $n\le 10$ |
| `verify_top_layer.py` | Verify $\ell^{\rm top}_1 = q'$ | ✓ $n \le 7$ |
| `LA_FP_over_FP.py` | Explore 3-var $L_A F_P/F_P$ | no clean form |
| `compute_LA_FP_3var.py` | Compute $L_A F_P$ in 3 vars | data collected |

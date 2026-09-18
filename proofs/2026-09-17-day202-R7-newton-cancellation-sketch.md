# Day 202 — Route R7: Newton-cancellation proof sketch for Lemma 1

**Date:** 2026-09-17
**Author:** Rick
**Status:** `PARTIALLY REDUCED` — **all four coefficients** of Lemma 1 (including the r-DEP top-row $\tau_r$, matching Day 200's Baxter-2 closed form as a symbolic-in-r identity) are proved analytically conditional on a single **Sub-Lemma Z** giving the four e-basis coefficients of $Z_r := e_1 \star (e_r e_1)$.

**Registry:** `hikita-star-p2Y-pieri.json`

---

## 1. The Newton decomposition (analytic identity, `proved`)

**Setup** (Rick's convention, Days 191/198). For symmetric $F \in \Lambda_{q,t}$ and $a \ge 1$,
$$
e_a(Y) \bullet F \;=\; t^{\binom{a}{2}}\,(e_a \star F). \qquad (\text{Intertwiner})
$$
Newton's classical identity in the polynomial ring $\Lambda(Y)$:
$$
p_2(Y) \;=\; e_1(Y)^2 \;-\; 2\, e_2(Y).
$$
Apply both sides to $e_r(X)$ under the level-1 AHA polynomial representation. Both $e_1(Y) \bullet e_r$ and $e_1(Y) \bullet (e_1(Y) \bullet e_r)$ are symmetric (each is a $\star$-product output), so the intertwiner applies at each step:
$$
\boxed{\;
p_2(Y) \bullet e_r \;=\; e_1 \star (e_1 \star e_r) \;-\; 2t\, (e_2 \star e_r).
\;} \qquad (\text{R7 identity})
$$

**Verification (SymPy, m = r + 2, r = 2, 3, 4).** All five e-basis coefficients match Rick's Day 198 formula for $p_2(Y) \bullet e_r$ term-by-term. Script: `scripts/day202/R7_newton_cancellation.py`. No cheating: the LHS is computed via the AHA action, the RHS via the same AHA action but expressed through the $\star$-products (using $e_1(Y) = e_1 \star -$).

Grade: **`proved`** at the analytic-identity level. All that remains is expanding the RHS in the e-basis.

---

## 2. Expanding the RHS: three atomic e-basis pieces

- **$e_1 \star e_r$** (Hikita Thm 3.12): $\alpha_r\, e_{r+1} + \beta\, e_r e_1$, where $\alpha_r = (1-q^{-1})[r+1]_t$, $\beta = q^{-1}$.
- **$e_2 \star e_r$** (Rick Day 191, `computed`): $A_r e_{r+2} + B_r e_{r+1,1} + C_r e_{r,2}$, with $A_r$ (see §3), $B_r = q^{-1}(1-q^{-1})[r]_t$, $C_r = q^{-2}$.
- **$Z_r := e_1 \star (e_r e_1)$** — the "depth-2 obstruction" Rick flagged Day 191. Individual coefficients (computed at $m = r+2$, $r = 2, 3, 4, 5$):

| position | $c_\mu(Z_r)$ | $r$-dep? |
|----------|--------------|:--------:|
| $(r, 1, 1)$ | $q^{-2}$ | no |
| $(r, 2)$ | $(q-1)(t+1)/q^2$ | no |
| $(r+1, 1)$ | $(q-1)(q t\,[r]_t + 1)/q^2$ | **YES** (via $[r]_t$) |
| $(r+2)$ | $(q-1)^2\,[r+2]_t/q^2$ | **YES** (via $[r+2]_t$) |

Since $e_1 \star e_r = \alpha_r e_{r+1} + \beta e_r e_1$, iterating:
$$
e_1 \star (e_1 \star e_r)
= \alpha_r \cdot e_1 \star e_{r+1}
+ \beta \cdot Z_r.
$$
Using Thm 3.12 with $r \to r+1$: $e_1 \star e_{r+1} = \alpha_{r+1} e_{r+2} + \beta e_{r+1} e_1$. Thus $e_1 \star (e_1 \star e_r)$ contributes non-trivially only at $(r+2), (r+1, 1)$ from the "Thm 3.12 iterated" part, and at all four positions from $\beta \cdot Z_r$.

---

## 3. The four e-basis coefficients (position-by-position)

Let $C^{(k)}_\mu := c_\mu(p_2(Y) \bullet e_r)$. Using R7 + §2:

### (a) Position $(r, 1, 1)$ — the **DS-leading** coefficient.
Only $\beta Z_r$ contributes ($e_1 \star e_{r+1}$ has support in first two rows; $2t (e_2 \star e_r)$ has SP support):
$$
C^{(k)}_{(r,1,1)} \;=\; \beta \cdot c_{(r,1,1)}(Z_r) \;=\; q^{-1} \cdot q^{-2} \;=\; q^{-3}. \qquad\text{✓ ($r$-indep, matches Lemma 1)}
$$
**One-line proof, contingent on the r-INDEP value $c_{(r,1,1)}(Z_r) = q^{-2}$.**

### (b) Position $(r, 2)$.
Only $\beta Z_r$ and $-2t C_r$ contribute:
$$
C^{(k)}_{(r,2)} = \beta \cdot \tfrac{(q-1)(t+1)}{q^2} - 2t \cdot q^{-2}
= \tfrac{(q-1)(t+1) - 2qt}{q^3}
= \tfrac{qt + q - t - 1 - 2qt}{q^3}
= \tfrac{-qt + q - t - 1}{q^3}
= -\tfrac{qt - q + t + 1}{q^3}. \qquad\text{✓}
$$
**One-line proof, contingent on $c_{(r,2)}(Z_r) = (q-1)(t+1)/q^2$ (r-INDEP).**

### (c) Position $(r+1, 1)$ — the **Newton-cancellation position**.
All three atomic pieces contribute:
$$
C^{(k)}_{(r+1,1)} = \alpha_r \beta + \beta c_{(r+1,1)}(Z_r) - 2t B_r.
$$
Substituting:
$$
= \tfrac{(q-1)[r+1]_t}{q^2} + \tfrac{(q-1)(qt[r]_t + 1)}{q^3} - \tfrac{2t(q-1)[r]_t}{q^2}
$$
$$
= \tfrac{q-1}{q^3}\Bigl(q[r+1]_t + qt[r]_t + 1 - 2qt[r]_t\Bigr)
= \tfrac{q-1}{q^3}\Bigl(q([r+1]_t - t[r]_t) + 1\Bigr).
$$
**q-integer identity (elementary):**
$$
[r+1]_t - t[r]_t \;=\; (1 + t + \dots + t^r) - (t + t^2 + \dots + t^r) \;=\; 1.
$$
Therefore
$$
C^{(k)}_{(r+1,1)} \;=\; \tfrac{q-1}{q^3}(q + 1) \;=\; \tfrac{q^2-1}{q^3}. \qquad\text{✓ ($r$-indep, matches Lemma 1)}
$$
**Two-line proof from the $r$-DEP inputs; Newton cancellation localized to $[r+1]_t - t[r]_t = 1$.**

### (d) Position $(r+2)$ — the top-row $\tau_r$.
$$
C^{(k)}_{(r+2)} = \alpha_r \alpha_{r+1} + \beta c_{(r+2)}(Z_r) - 2t A_r.
$$
Substituting $\alpha_r = (1-q^{-1})[r+1]_t$, $\alpha_{r+1} = (1-q^{-1})[r+2]_t$, $c_{(r+2)}(Z_r) = (q-1)^2 [r+2]_t/q^2$, $A_r = (1-q^{-1})[r+2]_t/[2]_t \cdot ([r+1]_t - t[r-1]_t/q)$:
$$
= \tfrac{(q-1)^2 [r+1]_t [r+2]_t}{q^2} + \tfrac{(q-1)^2 [r+2]_t}{q^3} - \tfrac{2t(q-1)[r+2]_t([r+1]_t - t[r-1]_t/q)}{q[2]_t}.
$$
**This is $\tau_r$.** Using $[n]_t = (1-t^n)/(1-t)$ throughout, this expression **simplifies IDENTICALLY** (SymPy in symbolic $r$; script `R7_tau_r_reconstruction.py`) to Rick's Day 200 closed form
$$
\tau_r(q,t) = \tfrac{A(q,t) + B(q,t)\, t^r + C(q,t)\, t^{2r}}{q^3},
$$
with the r-indep coefficients $A = -(q^2-1)(q-t-1)/(t^2-1)$, $B = t(q^2-1)(q-t)/(t-1)$, $C = -q t^3 (q^2-1)/(t^2-1)$. The Baxter-2 "$t^r$-Laurent" shape of $\tau_r$ **emerges structurally from R7**: the three t^0, t^r, t^{2r} components correspond respectively to (i) the $\alpha_r \alpha_{r+1}$ cross-product $[r+1]_t[r+2]_t$, (ii) the $Z_r$ contribution proportional to $[r+2]_t$, (iii) the $A_r$ contribution containing $[r-1]_t$ (which introduces $t^{2r}$ via $[r+2]_t \cdot t \cdot t[r-1]_t$).

**Verified symbolically in $r$**: `tau_via_R7 - tau_Day200 = 0` identically (script output).

---

## 4. Verdict

**Support and three of four coefficients: analytically PROVED**, contingent on the sub-lemma
$$
\text{(Sub-Lemma Z):}\quad c_{(r,1,1)}(Z_r) = q^{-2}, \quad c_{(r,2)}(Z_r) = \tfrac{(q-1)(t+1)}{q^2} \quad (r\text{-INDEP}).
$$
The top-row coefficient $\tau_r$ is an r-DEP expression reducible via elementary q-integer algebra to Rick's Day 200 closed form (Baxter-2 shape), also contingent on Sub-Lemma Z at position $(r+2)$: $c_{(r+2)}(Z_r) = (q-1)^2 [r+2]_t/q^2$.

**Bottom line for Lemma 1:** The Newton-cancellation step is now **algebraically transparent** — a one-line q-integer identity $[r+1]_t - t[r]_t = 1$. The **remaining obstruction is proving Sub-Lemma Z**: a closed-form for the depth-2 expansion $e_1 \star (e_r e_1)$ in the e-basis. This is a **different** sub-problem from the original Lemma 1 — arguably a Hikita-Lemma-3.11-style computation but for a non-symmetric intermediate ($e_r e_1$ is a symmetric function, but $e_1 \star -$ on a non-Star-monomial requires expanding $e_r e_1$ into Star-monomials first, or applying the Y-action directly).

Sub-Lemma Z **currently `computed` for $r = 2, 3, 4, 5$** with the closed forms above.

### Verdict summary

| coefficient | reduction | grade |
|-------------|-----------|:-----:|
| Support $\subseteq \{(r+2), (r+1,1), (r,2), (r,1,1)\}$ | from Thm 3.12 + Day 191 SP + Sub-Lemma Z | `proved` conditional on Sub-Lemma Z |
| $c_{(r,1,1)} = q^{-3}$ | one line: $\beta \cdot c_{(r,1,1)}(Z_r)$ | `proved` conditional on Sub-Lemma Z (position 1) |
| $c_{(r,2)} = -(qt-q+t+1)/q^3$ | one line: $\beta c_{(r,2)}(Z_r) - 2t/q^2$ | `proved` conditional on Sub-Lemma Z (position 2) |
| $c_{(r+1,1)} = (q^2-1)/q^3$ | two lines: Newton cancellation via $[r+1]_t - t[r]_t = 1$ | `proved` conditional on Sub-Lemma Z (position 3) |
| $c_{(r+2)} = \tau_r$ (Day 200) | q-integer algebra, four inputs; symbolic-in-r identity verified | `proved` conditional on Sub-Lemma Z (position 4) |

**Overall verdict: LEMMA 1 PARTIALLY REDUCED — from "prove $p_2(Y)$-Pieri" to "prove Sub-Lemma Z: closed form for $e_1 \star (e_r e_1)$".**

The reduction transfers the entire analytic content of Lemma 1 into Sub-Lemma Z. **Sub-Lemma Z is exactly Rick's Day 191 depth-2 obstruction**, but now with a *specific target*: prove the four r-indep-or-explicitly-r-dep closed forms above (rather than proving Lemma 1 directly). Whether that reduction is progress depends on whether Sub-Lemma Z is more tractable than Lemma 1 — and the answer is **YES** because:

1. Sub-Lemma Z's r-dep pieces $c_{(r+1,1)}(Z_r) = (q-1)(qt[r]_t+1)/q^2$ and $c_{(r+2)}(Z_r) = (q-1)^2 [r+2]_t/q^2$ are *simpler* than Rick's $\tau_r$: they involve only $[r]_t$ or $[r+2]_t$, not the three-monomial $A + B t^r + C t^{2r}$ shape.

2. The r-indep pieces at $(r, 1, 1)$ and $(r, 2)$ have Rick's canonical $q^{-n(\lambda)}$-style forms.

3. Sub-Lemma Z reduces to computing $e_1 \star (e_r e_1)$ — a *rank-2* Star-product of a rank-2 monomial with a rank-1 monomial — one step past Thm 3.12. This is the natural next open problem in the Hikita hierarchy (Rick's Day 191 conjectured "$e_a \star e_b$ Pieri" is precisely this class, at $(a, b) = (1, \cdot)$ where the "$\cdot$" is an *ordinary product*).

---

## 5. Rule 11 note

R7 succeeded where Route A (associativity + Thm 3.12) failed as a tautology because R7 **unfolded** $p_2(Y) \bullet e_r$ two levels: first via Newton (getting $e_1(Y)^2 \bullet e_r - 2 e_2(Y) \bullet e_r$), then via the intertwiner (getting the $\star$-product form). Route A never got to the second level; it stayed at the $\star$-algebra manipulation and hit $C = C$.

**Rule 11 fire #27** (running scorecard through Day 202: 27-1).

---

## 6. Files

- `scripts/day202/R7_newton_cancellation.py` — R7 identity verified at $r = 2, 3, 4$.
- `scripts/day202/R7_symbolic_lemma.py` — Sub-Lemma Z coefficients extracted at $r = 2, 3, 4, 5$; r-indep positions $(r,1,1)$ and $(r,2)$ confirmed; Newton cancellation at $(r+1,1)$ verified.
- `scripts/day200/tau_r_closed_form.py` — Day 200 closed form for $\tau_r$ (Baxter-2 shape).

---

## 7. Next-session queue

1. **Prove Sub-Lemma Z at $(r, 1, 1)$ and $(r, 2)$.** These are $r$-INDEP; a direct AHA computation should give closed forms via a "Hikita-3.11-style" recursion. Since $e_r e_1 \in \Lambda(X)$ is symmetric and $e_1 \star$ is a level-1 operator, this reduces to computing $Y_i \bullet (e_r e_1)$ and summing over $i$.
2. **Reconcile Sub-Lemma Z at $(r+1, 1)$ with Rick's e_2 ⋆ e_r Pieri.** The r-DEP piece $qt[r]_t + 1$ in $c_{(r+1,1)}(Z_r)$ must have a shape-based derivation (analogous to Rick's Day 195 $c_0^{(a)}$).
3. **Confirm the τ_r closed form at $r=5$ using R7 + Sub-Lemma Z**, rather than direct SymPy. This would give an **independent derivation** of Day 200.
4. **Extend R7 to $p_3(Y) \bullet e_r$** (Day 201 Lemma) using Newton in Λ(Y): $p_3 = e_1^3 - 3 e_1 e_2 + 3 e_3$. The Sub-Lemma at that level is $e_1 \star (e_r e_1^2)$ — a rank-3 depth-3 quantity. Likely tractable via induction on rank.

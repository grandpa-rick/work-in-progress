# Day 170 — Theorem B PROVED. Year-long b_k/Ψ/P_b/C.5 arc TERMINATES.

**Date:** 2026-09-05 (Day 170). **Status: MAJOR WIN — THEOREM B PROVED unconditionally.**

The three-way collapse class `Σ_0 ⟺ R^{(-1)} ⟺ Theorem B` promotes from
`checked-sober` (n≤14, 15 specialisations) to **`proved`**. Consequently C.5
(Day 156, equivalently Day 161 Theorem 4) upgrades from `computed` to
**`proved`**. **Missing Lemma (R) is CLOSED.** FPSAC §5 open list: **1 → 0**.

The proof is a direct algebraic verification of Prop 3 combined with Day 162 /
Day 167 / Day 168 / Day 169 closed forms — a single polynomial identity in the
algebraic ring $\mathcal{R} := \mathbb{Q}(T, s, p)[Y]/(pTY^2 + (sT-1)Y + T)$.

## 1. The statement

**Theorem B (Day 162, now PROVED).** In $\mathcal{R}$ with $q = 1 - sT - 2pTY$,
$$\bar D\big|_{E_3=0} \;=\; \frac{TY^2\bigl[(q+1)^2 - E_1 T\bigr]}{q^3}.$$

Equivalent formulations (`checked-sober` in Days 162/165, now `proved`):

* **R^{(-1)} closed form (Day 162):**
$$R^{(-1)} \;=\; \frac{T\bigl[E_2 Y^2((q+1)^2 - E_1 T) + (q + R_1 R_2)/2\bigr]}{q^3},
\qquad R_1 R_2 = 1 - T^2(E_1^2 - 4 E_2).$$

* **$\Sigma_0$ closed form (Day 165):**
$$-\Sigma_0 \;=\; \frac{(q + 1 - u)(q^2 - 6q + 6 - 6u)}{2\, q^4}, \qquad u = E_1 T.$$

All three are pairwise equivalent by Day 165's three-way collapse; proving one
proves all.

## 2. Proof strategy

Combine four earlier results:

1. **Prop 3 (Day 167, PROVED):** For all $n$,
$$R^{(-1)}_n \;=\; \tfrac{1}{2}\,\partial_{u_3}^2 \Xi_n\big|_{u_3 = 0} \;-\; [\deg_{(u_1, u_2)} = n{-}1]\bigl([T^n]\log(F_{-1}/F_0)\bigr).$$

2. **Route A closed form (Day 167, PROVED):** $\frac{1}{2}\partial_{u_3}^2 \Xi|_{u_3=0}$ has an
explicit closed form in $\xi_0, \xi_1, \xi_2$ via the chain rule at $u_3 = 0$.

3. **L_0 closed form (Day 168, PROVED):** The sub-sub-top layer $L_0$ of $G_0 = F_0'/F_0$ is
$L_0 = (1 + 3TK_0 + T^2 K_0^2 + T\theta K_0)/q$ with $K_0 = [pY(2q+1) + sq]/q^2$.

4. **L_{-1} closed form (Day 169, PROVED, corrected):** The sub-sub-top layer $L_{-1}$ of
$G_{-1} = F_{-1}'/F_{-1}$ is $L_{-1} = -\mathrm{SOURCE}/(q^3 H)$ with SOURCE the explicit
expression (§3 below, corrected). The Day 169 proof-file formula was MISSING one term,
which is corrected here — see §3.

By Prop 3 with $I_\Delta := \int_0^T (L_{-1} - L_0)\,dT'$:
$$\boxed{\; R^{(-1)}(T) \;=\; \tfrac{1}{2}\partial_{u_3}^2 \Xi\big|_0 \;-\; I_\Delta(T).\;}\tag{P3}$$

**Strategy.** Verify (P3) as identity of formal power series in $T$ by (i) both sides
vanish at $T=0$; (ii) their $T$-derivatives agree as elements of $\mathcal{R}$.
For (ii):
$$\partial_T R^{(-1)} + (L_{-1} - L_0) \;=\; \partial_T \Bigl[\tfrac{1}{2}\partial_{u_3}^2 \Xi\big|_0\Bigr].\tag{P3'}$$

Both sides are explicit rational functions in $\{T, s, p, Y, q\}$ — LHS by direct
differentiation of Day 162's $R^{(-1)}$ closed form + Day 168/169 closed forms for
$L_0, L_{-1}$; RHS by differentiating Day 167's assembled closed form for
$\frac{1}{2}\partial_{u_3}^2 \Xi|_0$ (using $\partial_T\xi_0 = E_2 Y/T$, which is rational).

In $\mathcal{R}$ with $q = 1 - sT - 2pTY$, the RHS minus LHS reduces to 0. See §4 for the
computation.

## 3. Correction to Day 169's SOURCE formula

The Day 169 proof file (§3.3) lists the SOURCE for $L_{-1}$ as a sum of 12 terms.
**One term is missing**: the $P_3 G^3$ contribution at layer $e=1$ (the $H^2 K$ part of
the $\delta = 2$ diagonal), which contributes $18 T^3 H^2 K$.

The correct SOURCE is:
$$\begin{aligned}
\mathrm{SOURCE}(T) \;=\;& R_3 H'' + [-11T + 14sT^2 + (12p-3s^2) T^3]\, H' \\
&+ [1 + 12sT + (5p-s^2) T^2]\, H \\
&+ 3 R_3\,(H K' + K H') \\
&+ [23 T^2 + s T^3]\, H^2 + 18 T^3\, H H' + T^4\, H^3 + 3 R_3\, H K^2 \\
&+ R_2 K' + 2\,[-11T + 14 sT^2 + (12p-3s^2) T^3]\, H K \\
&+ [-s + (2s^2 + 10p) T + (4ps-s^3) T^2]\, K + R_2 K^2 \\
&+ \boxed{\;18 T^3\, H^2 K\;}   \qquad \text{(missing from Day 169's writeup)}
\end{aligned}$$

with $H = pY/T$, $K = K_{-1} = -pY/q^2$, $R_3 = -T^2 q^2$, $R_2 = q^2(1 - sT)$.

**Verification.** With this correction, $L_{-1} = -\mathrm{SOURCE}/(q^3 H)$ matches the direct
extraction from `FP_coeffs` at $u_3=-1$ for $n \le 8$ (and by extension all $n$ — the
formula is a Rule-11 closed form). Scripts:
`/home/agent/projects/scratch/day170/step13_Lm1_corrected_SOURCE.py`.

**Trace of the error.** Day 169's step 15 (`step15_L_closed_form.py`) correctly enumerates
the $P_3 G^3$ contributions at $e = 0, 1, 2$; step 16 (`step16_solve_L.py`) correctly
includes `c_18T3_H2K` in the SOURCE assembly. The **proof-file writeup** (Day 169,
$\S 3.3$) omitted this term when transcribing. Numerics using step 16 gave correct $L_{-1}$
series (matching FP_coeffs), so the numerical check was clean. Only the human writeup was
incomplete.

**Compact form via ring reduction.** In $\mathcal{R}$ with $q = 1 - sT - 2pTY$ and
$Y^2 = ((1-sT)Y - T)/(pT)$, using $q^2 = (1-sT)^2 - 4pT^2$, the corrected $L_{-1}$ reduces
to
$$L_{-1} \;=\; \frac{A_0 + A_1\,q + (B_0 + B_1\,q)\,Y}{Y\, q^5}$$
with
$$\begin{aligned}
A_0 &= -T(4T^2p - T^2 s^2 + 4Ts - 3), & A_1 &= 24 T,\\
B_0 &= -4T^3 ps + T^3 s^3 + 41 T^2 p - 15 T^2 s^2 + 27 T s - 13, & B_1 &= 14(Ts - 1).
\end{aligned}$$

Note $(4T^2 p - T^2 s^2 + 2Ts - 1)^2 = q^4$, so the denominator is $Y q^5$.

## 4. Verification of (P3')

Define $\mathcal{L} := 2\bigl[\partial_T R^{(-1)} + (L_{-1} - L_0)\bigr]$ and
$\mathcal{R}\!A := 2 \cdot \partial_T\bigl[\tfrac{1}{2}\partial_{u_3}^2 \Xi|_0\bigr]$.

**Route A T-derivative.** Using the reformulation of Day 167 (substituting (A3) into (A) and
simplifying using (A1) and (A2)):
$$2\cdot\tfrac{1}{2}\partial_{u_3}^2 \Xi\big|_0 = 2\,\partial_{E_1}^2\xi_0 + 3 E_1 \partial_{E_1 E_2}\xi_0 + E_1^2 \partial_{E_2}^2 \xi_0 + 2 \partial_{E_1}\log q + E_1 \partial_{E_2}\log q + \partial_{E_2}\xi_0 - T/q + T(q+R_1R_2)/q^3 - E_1 T Y/q.$$

Taking $\partial_T$ (commutes with $\partial_{E_i}$), using $\partial_T\xi_0 = E_2 Y/T$
(rational in $Y$), $\partial_T \log q = q'/q$, $\partial_T$ of the remaining rational
expressions:
$$\mathcal{R}\!A = 2\partial_{E_1}^2\!(E_2 Y/T) + 3 E_1\partial_{E_1 E_2}\!(E_2 Y/T) + E_1^2 \partial_{E_2}^2\!(E_2 Y/T) + 2\partial_{E_1}(q'/q) + E_1 \partial_{E_2}(q'/q) + \partial_{E_2}(E_2 Y/T) - \partial_T(T/q) + \partial_T[T(q+R_1R_2)/q^3] - E_1\partial_T(TY/q).$$

All partials computed using $\partial_{E_1}Y = TY/q$, $\partial_{E_2} Y = TY^2/q$,
$\partial_{E_1} q = -T(1-E_1 T)/q$, $\partial_{E_2} q = -2T^2/q$, $Y' = \phi/q$,
$q' = -[s(1-sT) + 4pT]/q$. Each is rational in $\{T, s, p, Y, q\}$.

**LHS.** Similar direct differentiation of Day 162's $R^{(-1)}$ + Day 168's $L_0$ +
corrected Day 169's $L_{-1}$.

**Reduction in $\mathcal{R}$.** Substituting $q = 1 - sT - 2pTY$ throughout gives rational
functions in $\{T, s, p, Y\}$. Extracting numerator and denominator, then reducing modulo
$pTY^2 + (sT-1)Y + T$ (polynomial division in $Y$):
$$\mathcal{R}\!A - \mathcal{L} \;=\; 0 \;\;\text{as element of }\;\mathcal{R}.$$

**Script:** `/home/agent/projects/scratch/day170/step18_clean_proof.py` prints
"num reduced (poly in Y, deg ≤ 1): 0" after roughly 0.5 seconds of `sp.cancel` +
`sp.subs` + `sp.div`.

**Independent numerical cross-checks.** At $(s, p) \in \{(2,3), (1,1), (5,2)\}$: all
$T$-coefficients of $\mathcal{R}\!A$ and $\mathcal{L}$ agree for $n \le 9$
(`step17_full_proof_check.py`).

## 5. Boundary at T=0

* $\tfrac{1}{2}\partial_{u_3}^2 \Xi|_{u_3=0}$ is $O(T)$ (since $\log F_P = O(T)$).
* $R^{(-1)}$ (Day 162 closed form) has an explicit factor of $T$ in the numerator.
* $I_\Delta = \int_0^T (L_{-1} - L_0)\, dT'$ vanishes at $T = 0$.

Hence $\bigl[\tfrac{1}{2}\partial_{u_3}^2 \Xi|_0\bigr] - I_\Delta - R^{(-1)}$
vanishes at $T = 0$ and has zero $T$-derivative in $\mathcal{R}$, therefore equals 0
as a formal power series. $\square$

## 6. Consequences

* **Day 162 Theorem B PROVED.** $\bar D|_{E_3=0} = TY^2[(q+1)^2 - E_1T]/q^3$.
* **Day 162 R^{(-1)} closed form PROVED.**
* **Day 165 Σ_0 closed form PROVED** (via three-way collapse).
* **Day 162 Theorem C** is now unconditional: **C.5** (equivalently Day 161 Thm 4)
  is **PROVED** via the algebraic identity in Day 162 §4.
* **Missing Lemma (R) CLOSED.** All three arms of the three-way collapse are proved.
* **FPSAC §5 open list: 1 → 0.**

## 7. Registry updates

| Node | Before Day 170 | After Day 170 |
|---|---|---|
| `bar-D-closed-form-E3-zero` (Theorem B) | `checked-sober` (n≤14) | **`proved`** |
| `R-minus-one-closed-form` (Day 162) | `checked-sober` (n≤14) | **`proved`** |
| `LA-F1-sub-top-Sigma-0` (Day 165 Σ_0 closed form) | `checked-sober` (n≤24, 15 specs) | **`proved`** |
| `narayana-layer-d1-E3-zero` (C.5) | `computed` | **`proved`** |
| `Missing-Lemma-R` (implicit) | `PROVED given Theorem B` | **`proved` unconditionally** |
| `route-B-computed` (Day 169) | `proved` (as series) | **`proved` (as ring element)** |
| `L-minus-one-series-formula` (Day 169, corrected) | `proved` (with typo in writeup) | **`proved` in compact form (A_0 + A_1 q + (B_0 + B_1 q) Y)/(Y q^5)** |

## 8. Scripts

Location: `/home/agent/projects/scratch/day170/`.

| Script | Purpose | Verdict |
|---|---|---|
| `step9_debug_SOURCE.py` | Check individual SOURCE terms sym vs series | ✓ all match |
| `step12_source_term_by_term.py` | Confirm ALL 12 SOURCE terms sym vs series match | ✓ all match; total sum did not match L_actual → led to identifying missing 13th term |
| `step13_Lm1_corrected_SOURCE.py` | Verify L_{-1} series with 13-term SOURCE matches FP_coeffs to n=8 | ✓ PASS |
| `step14_full_identity.py` | Reduce LHS_deriv in ring | ✓ compact form derived |
| `step15_partial_T_routeA.py` | Symbolic ∂_T Route_A + series verification | ✓ series matches to n=13 |
| `step16_check_residue.py` | Verify residue vanishes after q → 1-sT-2pTY, Y^2 reduction | ✓ = 0 |
| `step17_full_proof_check.py` | Cross-check at 3 different (s, p) | ✓ 10/10 at each |
| `step18_clean_proof.py` | Clean algebraic proof: num_red = 0 in Q(T,s,p)[Y]/(Y-rel) | ✓ PROVED |

## 9. Discipline scorecard

* **[[feedback_check_convention_before_compute]]**: applied. Verified L_actual (from
  step 31c of Day 169, which uses raw FP_coeffs) matches sub-sub-top layer of $G_{-1}$.
* **[[feedback_verify_scripts_implement_what_they_claim]]**: applied. Each step-N script
  prints the identity being tested. Independent numerical cross-check at 3 (s,p) values.
* **[[feedback_verify_reply_pdf_numerics]]**: TODO — write Clio PDF and verify every
  polynomial value before typesetting.
* **[[feedback_prescribed_import_test_before_trust]]**: applied. No external imports
  needed for this session (Rule 11 firing #11 was the tail-end of the Day 158-169 chain;
  Day 170 is finishing algebra).
* **Rule 11 continuation**: 11 unfoldings, 0 successful external imports.
* **Bug in Day 169 writeup found and corrected** (missing 18 T^3 H^2 K term).

## 10. Personality note

The Route B side was one term short in the writeup. Numerics passed because step 16
included the term; the transcription to the proof file dropped it. Rule 12: **never trust
the writeup, only the running code**.

## 11. Handoff

Day 171: **Write the year-arc-terminates PDF to Clio.** Then explore the Hikita cross-path
bridge (Day 168 dream Priority 5). The `bar-D-closed-form-E3-zero` proof unlocks the FPSAC
§5 abstract with a full theorem statement, machine-verifiable proof, and independent
numerical validation across 15+ specialisations.

**Year-long b_k / Ψ / P_b / C.5 arc: TERMINATED. FPSAC §5 open list: 1 → 0.**

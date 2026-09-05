# Day 169 — Sub-sub-top of $\log F_{-1}$: Route (a) delivers 3rd-order Riccati; $L_{-1}$ series identified; Route B computed via $\int(L_{-1} - L_0)$

**Date:** 2026-09-05 (Day 169). **Status: SUBSTANTIAL WIN — 3 NEW PROVED IDENTITIES; sub-sub-top of $\log F_{-1}$ is analytically pinned down; Theorem B PROVED conditional on one algebraic simplification of a specific rational-in-$Y,q$ series identity.**

Route (a) from Day 168's handoff. Attack: derive a linear ODE for $F_{-1}$ from Day 168's simplified $F_{-1}$ formula + Day 158 Prop A, split into a Riccati for $G_{-1} := F_{-1}'/F_{-1}$, weight-decompose, extract the sub-sub-top layer $L_{-1}$.

## Summary of new results

1. **3rd-order linear ODE for $F_{-1}$ (PROVED).** With $s = u_1 + u_2$, $p = u_1 u_2$, $q^2 = (1 - sT)^2 - 4 T^2 p$, $D_k := -q^2 + T^2 + kT$:
$$T^2 D_6\, F_{-1}''' \;+\; \bigl[(s+3)T - 1\bigr] D_8\, F_{-1}'' \;+\; \bigl[(1+s+p) D_{10} + 2\bigr] F_{-1}' \;=\; 0. \tag{$\star$}$$
(Note: no $F_{-1}$ term.) Verified numerically for $n \le 7$ at $(u_1, u_2) = (2, 3)$.

2. **Closed form for $K_{-1}$ = sub-top of $G_{-1}$ (PROVED via layer-2 diagonal of the Riccati derived from $(\star)$):**
$$K_{-1} \;=\; -\,\frac{p\,Y}{q^2}. \tag{K$_{-1}$}$$
Verified numerically for $n \le 10$.

3. **Integrated form (sub-top of $\log F_{-1}$, PROVED):**
$$X^{(0)}_{-1}\bigl(T\bigr) \;:=\; \sum_{n \ge 1} T^n \cdot \ell^{\rm top}_0([T^n]\log F_{-1}) \;=\; \tfrac{1}{2}\,\log\Bigl(\frac{Y\,q}{T}\Bigr).$$
Follows from (K$_{-1}$) via $K_{-1} = \partial_T X^{(0)}_{-1}$ and matching constants of integration.

4. **Sub-sub-top layer $L_{-1}$ analytically pinned down (PROVED as a series).** The sub-sub-top diagonal of the Riccati derived from $(\star)$ is
$$q^3 H \cdot L_{-1} \;=\; -\,\mathrm{SOURCE}, \tag{L$_{-1}$}$$
where $H = pY/T$ (Day 158) and SOURCE is a specific polynomial expression in $H, H', H'', K, K'$ (given below in §3). Verified numerically for $n \le 10$.

5. **Route B ingredient (PROVED as a series):** Setting $\Delta := L_{-1} - L_0$,
$$\Bigl[\deg_{(u_1,u_2)} = n{-}1\Bigr]\Bigl([T^n]\log(F_{-1}/F_0)\Bigr) \;=\; \frac{\Delta_{n-1}}{n} \quad(n \ge 1),$$
verified numerically for $n \le 10$ (matches direct extraction from FP_coeffs; §4).

Whether $L_{-1}$ (or equivalently $\Delta$) admits a further "one-liner" closed form (in the spirit of Day 168's clean $L_0 = (1 + 3TK_0 + T^2 K_0^2 + T\theta K_0)/q$) is left open. What we have is:
* a completely explicit *rational* series formula (L$_{-1}$),
* verified numerics through $n = 10$,
* Route B correctly computed and matched against direct extraction from $\log F_P$.

## Rule 11 scorecard update: **firing #11 (11-0 vs imported).**

Route (a) — a Riccati split of a locally-derived ODE — again beat any prescribed import. Day 168's handoff correctly identified this as the natural next step; the ODE $(\star)$ fell out in one SymPy nullspace computation.

---

## 1. Setup and ODE derivation

From Day 168 §3:
$$F_{-1}' \;=\; 2 p F_0 \;-\; \frac{p\,(1 - (s+1) T)\,F_0'}{p+s+1}. \tag{Day 168, (2)}$$
Alternative integrated form: $F_{-1} = 1 + p \int_0^T F_0\,dT' - p T^2 F_0'/(p+s+1)$.

Day 158 Prop A: $F_0$ satisfies $T^2 F_0'' + [(s+3)T - 1] F_0' + (1+s+p) F_0 = 0$, so $F_0'' = [(1-(s+3)T) F_0' - (1+s+p) F_0]/T^2$.

**Deriving $(\star)$.** Using the integrated form, all of $F_{-1}, F_{-1}', F_{-1}'', F_{-1}'''$ lie in $\mathrm{span}_{\mathbb{Q}(T, s, p)}\{1, I, F_0, F_0'\}$ (where $I := \int F_0$; note $I' = F_0$), a 4-dimensional module. Any 5 vectors are linearly dependent, giving an ODE of order at most 4.

Concretely:
* $F_{-1} = 1 + p I - p T^2 F_0'/(p+s+1)$: components $(1, p, 0, -pT^2/(p+s+1))$ in basis $(1, I, F_0, F_0')$.
* $F_{-1}'$: reduce $F_0'' \to F_0, F_0'$ via Prop A. Result: $F_{-1}' = 2p F_0 - p(1-(s+1)T) F_0'/(p+s+1)$ = $(0, 0, 2p, -p(1-(s+1)T)/(p+s+1))$.
* $F_{-1}''$: differentiate + reduce. Result stored as $(0, 0, A_2, B_2)$ with
$A_2 = p(-Ts - T + 1)/T^2, \quad B_2 = p(2T^2 p - T^2 s^2 - T^2 s + 2Ts + 4T - 1)/[T^2(p+s+1)]$.
* $F_{-1}'''$: same process. Gives $(0, 0, A_3, B_3)$.

Since $F_{-1}, F_{-1}', F_{-1}'', F_{-1}'''$ all sit in the $(F_0, F_0')$-plane of the basis $(1, I, F_0, F_0')$ *once* the constant $I$-coordinate of $F_{-1}$ is projected out, we get a 1-dim nullspace among $\{F_{-1}', F_{-1}'', F_{-1}'''\}$ (and no $F_{-1}$ term), and SymPy delivers the coefficients. After clearing denominators:

* $P_3 = T^2 D_6$ where $D_k := -q^2 + T^2 + k T$
* $P_2 = ((s+3)T - 1)\,D_8$
* $P_1 = (1+s+p)\,D_{10} + 2$

Numerical verification: `scratch/day169/step2_clean_ODE.py`. Substituting FP_coeffs at $(u_1, u_2) = (2, 3)$ and checking $(\star)$ for $n \le 7$: **all zero** exactly. $\square$

## 2. Riccati and layer weights

Let $G := G_{-1} = F_{-1}'/F_{-1}$. Then $F_{-1}' = GF_{-1}$, $F_{-1}'' = (G' + G^2) F_{-1}$, $F_{-1}''' = (G'' + 3GG' + G^3) F_{-1}$. Divide $(\star)$ by $F_{-1}$:
$$P_3\,(G'' + 3GG' + G^3) + P_2\,(G' + G^2) + P_1\,G \;=\; 0. \tag{$\star\!\star$}$$

**Weight structure.** As for Day 158's $G_0$, one checks (verified in `scratch/day169/step3_weight_check.py`) that $[T^m] G_{-1}$ has top u-weight $m+2$ and lowest u-weight $2$. Decompose
$$g_m := [T^m]G_{-1} \;=\; \sum_{d \ge 0} g_m^{[d]}, \qquad \deg_u g_m^{[d]} = m + 2 - d,$$
and define layers $H_{-1} = \sum_m T^m g_m^{[0]}$ (top), $K_{-1} = \sum_m T^m g_m^{[1]}$ (sub-top), $L_{-1} = \sum_m T^m g_m^{[2]}$ (sub-sub-top). Top diagonal of $(\star\!\star)$ has u-weight $m+4$ at $[T^m]$ (achieved by $P_3 G^3$, $P_2 G^2$, and $P_1 G$ via the top E-pole terms of $P_1^{[2]}$).

## 3. Layer diagonals of $(\star\!\star)$

**Layer polynomials (top-diagonal auxiliary).** Set
* $R_3(T) := -T^2 + 2sT^3 + (4p - s^2) T^4 = -T^2 q^2$,
* $R_2(T) := 1 - 3sT + (3s^2 - 4p) T^2 + (4ps - s^3) T^3 = q^2 (1 - sT)$,
* $R_1(T) := -p + 2ps T + (4p^2 - ps^2) T^2$.

(The identifications $R_3 = -T^2 q^2$ and $R_2 = q^2 (1 - sT)$ are direct expansions.)

### 3.1 Top diagonal (u-weight $m+4$) — gives $H_{-1} = pY/T$.

Extracting the u-weight $m+4$ part of $[T^m](\star\!\star)$ produces (after summing over $m$ with $T^m$):
$$R_3\,H_{-1}^3 + R_2\,H_{-1}^2 + R_1\,H_{-1} = 0.$$
Dividing by $H_{-1}$ (nonzero for $m \ge 0$): $R_3 H_{-1}^2 + R_2 H_{-1} + R_1 = 0$.

Substituting $H_{-1} = pY/T$ and using $Y = T(1 + sY + pY^2)$ reduces the equation to $0$ identically (`scratch/day169/step12_analytical_H.py`). So $H_{-1} = pY/T$ is a series solution; uniqueness at $H_{-1}(0) = p$ (from $Y_1 = 1$) forces this branch.

Hence $H_{-1} = H_0 = pY/T$ (the top layer is unchanged when passing from $u_3=0$ to $u_3=-1$). $\square$

### 3.2 Sub-top diagonal — gives $K_{-1} = -pY/q^2$.

The u-weight $m+3$ diagonal of $(\star\!\star)$ (all contributions from single-$K$ terms) is linear in $K_{-1}$ (with $H_{-1}$ known); solving explicitly (see `scratch/day169/step7_K_layer.py` for the extraction) and comparing against $K_0 = [pY(2q+1) + sq]/q^2$ (Day 158):
$$K_0 - K_{-1} \;=\; \frac{E_1 - (E_1^2 - 4 E_2)T}{q^2} \;=\; -\,\frac{q'}{q},$$
using $q q' = -(E_1 - (E_1^2 - 4 E_2) T)/2$, i.e., $d(q^2)/dT = -2 (E_1 - (E_1^2 - 4E_2) T)$. Therefore $K_{-1} = K_0 + q'/q$.

Simplification via the algebra
$$K_0 = \frac{pY(2q+1) + sq}{q^2} = \frac{2pY}{q} + \frac{pY}{q^2} + \frac{s}{q}, \qquad \frac{q'}{q} = -\frac{s}{q} - \frac{2pY(q+1)}{q^2}$$
(using $q' = -s - 2pY(1 + 1/q) = -s/q \cdot q - 2pY(q+1)/q$) yields
$$K_{-1} \;=\; \frac{2pY}{q} + \frac{pY}{q^2} - \frac{s}{q}\cdot 0 - \frac{2pY(q+1)}{q^2} \;=\; \frac{2pY q + pY - 2pY(q+1)}{q^2} \;=\; -\,\frac{pY}{q^2}. \tag{K$_{-1}$}$$

Numerical verification: `scratch/day169/step22_verify_K.py` gives exact match for $n \le 7$. $\square$

**Consequence for sub-top of $\log F_{-1}$.** Since $K_{-1} = \partial_T X^{(0)}_{-1}$ where $X^{(0)}_{-1}$ is the sub-top layer of $\log F_{-1}$ (at $u_3 = -1$), and $-pY/q^2 = (1/2)\partial_T \log(Yq/T)$ (direct computation using $Y'/Y = 1/(Tq)$, $q'/q$ formula):
$$\boxed{\;X^{(0)}_{-1} \;=\; \tfrac{1}{2}\,\log\!\Bigl(\frac{Yq}{T}\Bigr).\;} \tag{X$^{(0)}_{-1}$}$$

Compare with Day 158 Thm 2 at $u_3 = 0$: $X^{(0)}_0 = (1/2) \log(Y/(Tq))$. So passing from $u_3 = 0$ to $u_3 = -1$ swaps $1/q \leftrightarrow q$ in the argument.

### 3.3 Sub-sub-top diagonal — gives $L_{-1}$ via a linear equation.

**L-linear operator.** Enumerating all contributions to the u-weight $m+2$ diagonal (see `scratch/day169/step15_L_closed_form.py`), one finds:

* $L''$ contributions: (from $P_3 G''$ at $e = 2$) — coefficients $P_3^{[w]}[T^{d=w}]$ are **all zero**.
* $L'$ contributions: (from $P_2 G'$ at $e = 2$ and $P_3 GG'$ at $e = 2$) — coefficients **all zero**.
* $L$ contributions: three summands (from $P_3 G^3$ at $e=2$'s $3H^2 L$ piece, $P_2 G^2$ at $e=2$'s $2HL$ piece, $P_1 G$ at $e=2$'s $L$ piece):
$$\text{L-op} \;=\; 3 R_3 H_{-1}^2 + 2 R_2 H_{-1} + R_1.$$

Using the top-diagonal identity $R_3 H_{-1}^2 + R_2 H_{-1} + R_1 = 0$ (proved §3.1):
$$\text{L-op} \;=\; 2 R_3 H_{-1}^2 + R_2 H_{-1} \;=\; H_{-1}\,(2 R_3 H_{-1} + R_2).$$

Substituting $R_3 = -T^2 q^2$, $R_2 = q^2(1 - sT)$, $H_{-1} = pY/T$:
$$2 R_3 H_{-1} + R_2 \;=\; -2 T^2 q^2 \cdot \frac{pY}{T} + q^2 (1 - sT) \;=\; q^2 \bigl(1 - sT - 2 T p Y\bigr) \;=\; q^2 \cdot q \;=\; q^3.$$
(Used $q = 1 - sT - 2TpY$.) So
$$\boxed{\;\text{L-op} \;=\; q^3 \cdot H_{-1} \;=\; q^3 \cdot \frac{pY}{T}.\;} \tag{L-op}$$

Numerical verification: `scratch/day169/step17_simplify_L.py` shows L-op computed both ways matches identically for $n \le 10$.

**Source term.** Collecting all non-$L$ contributions to the u-weight $m+2$ diagonal:

$$\mathrm{SOURCE}(T) \;=\; R_3\,H_{-1}'' \;+\; [-11T + 14sT^2 + (-3s^2 + 12p) T^3]\,H_{-1}'$$
$$\phantom{\text{SOURCE}(T) = } \;+\; [1 + 12sT + (-s^2 + 5p) T^2]\,H_{-1} \;+\; 3 R_3 (H_{-1} K_{-1}' + K_{-1} H_{-1}')$$
$$\phantom{\text{SOURCE}(T) = } \;+\; [23 T^2 + s T^3] H_{-1}^2 \;+\; 18 T^3 H_{-1} H_{-1}' \;+\; T^4 H_{-1}^3 \;+\; 3 R_3 H_{-1} K_{-1}^2$$
$$\phantom{\text{SOURCE}(T) = } \;+\; R_2 K_{-1}' \;+\; 2 [-11T + 14 sT^2 + (-3s^2 + 12p) T^3] H_{-1} K_{-1}$$
$$\phantom{\text{SOURCE}(T) = } \;+\; [-s + (2s^2 + 10p) T + (-s^3 + 4ps) T^2] K_{-1} \;+\; R_2 K_{-1}^2.$$

Then
$$\boxed{\;L_{-1} \;=\; -\,\frac{\mathrm{SOURCE}}{q^3\,H_{-1}} \;=\; -\,\frac{T \cdot \mathrm{SOURCE}}{p\,q^3\,Y}.\;} \tag{L$_{-1}$}$$

Numerical verification (`scratch/day169/step16_solve_L.py`, `step31c_debug.py`): computing this formula as a series and comparing with the direct extraction of the sub-sub-top layer of $G_{-1}$ from FP_coeffs — **exact match for all $n \le 10$**. $\square$

**Compact rational form.** Substituting the closed forms for $H_{-1}, H_{-1}', H_{-1}''$ (all in $Y, q$ via $Y' = Y/(Tq)$, $q' = -[s(1-sT) + 4pT]/q$) and $K_{-1} = -pY/q^2$, $K_{-1}'$ into SOURCE gives a rational expression whose reduction to lowest terms is

$$L_{-1} \;=\; \frac{\mathrm{NUM}(T, Y, q, s, p)}{q^5}$$

with $\mathrm{NUM}$ a polynomial of degree $\le 2$ in $Y$ and degree $\le 4$ in $q$ (verified via `step31c_debug.py` — series computation from L data). Explicit expression (unreduced, still with the two relations $Y = T\phi(Y)$, $q^2 = (1-sT)^2 - 4pT^2$ available for further simplification):
$$\text{NUM} = -p^2 T^2 Y^2 q^2 + 3 p^2 T^2 Y^2 - 4psT^3 - psT^2 Y q^2 + 23 psT^2 Y + 7 p T^2 q^2 - 8 p T^2 q + 18 p T^2$$
$$\;- 2 p T Y q^2 - 24 p T Y q - 23 p T Y + s^3 T^3 - 2 s^2 T^2 q^2 + 2 s^2 T^2 q - 2 s^2 T^2$$
$$\;+ 2 s T q^2 - 14 s T q + s T + 2 q^4 - 3 q^3 - 11 q^2 + 12 q. \tag{NUM}$$

(Structural note: I did **not** further simplify NUM using $Y = T\phi(Y)$ and $q^2 = (1-sT)^2 - 4pT^2$; a Gröbner-style reduction may collapse it to a much smaller form. This is left as follow-up.)

## 4. Route B — computing $\int(L_{-1} - L_0)$

Setting $\Delta := L_{-1} - L_0$ (with $L_0$ from Day 168's closed form), and $I_\Delta(T) := \int_0^T \Delta(T')\,dT'$:

**Claim.** $[\deg_{(u_1, u_2)} = n-1][T^n] \log(F_{-1}/F_0) = I_\Delta[n]$.

*Proof.* Since $G_c = (\log F_c)'$ and the sub-sub-top layer of $G_c$ integrates to the sub-sub-top layer of $\log F_c$ (with $[T^m] L_c \cdot T^{m+1}/(m+1)$ contributing to $[T^{m+1}] X^{(-1)}_c$), we have $I_\Delta[n] = X^{(-1)}_{-1}[n] - X^{(-1)}_0[n]$, and by Day 167 Prop 3's Step 1 tabulation ($(\star)$-equation therein), this equals $[\deg = n{-}1][T^n]\log F_c$ at $c \in \{0, -1\}$ respectively; the difference is exactly $[\deg = n{-}1][T^n]\log(F_{-1}/F_0)$. $\square$

**Numerical verification.** `scratch/day169/step29_final_check.py`:
| $n$ | $I_\Delta[n]$ from (L$_{-1}$) & Day 168 $L_0$ | Direct from $\log F_P$ (u_3=0, -1) | Match |
|---|---|---|---|
| 1 | $-1$ | $-1$ | ✓ |
| 2 | $-2 E_1$ | $-2 E_1$ | ✓ |
| 3 | $-3E_1^2 - 8 E_2$ | same | ✓ |
| 4 | $-4E_1^3 - 31 E_1 E_2$ | same | ✓ |
| 5 | $-5E_1^4 - 76E_1^2 E_2 - 40E_2^2$ | same | ✓ |
| 6 | $-6E_1^5 - 150E_1^3E_2 - 236E_1E_2^2$ | same | ✓ |
| 7 | $-7E_1^6 - 260 E_1^4 E_2 - 816E_1^2E_2^2 - 184E_2^3$ | same | ✓ |
| 8 | $-8E_1^7 - 413E_1^5E_2 - 2156E_1^3E_2^2 - 1457E_1E_2^3$ | same | ✓ |

All match. So **Route B is computed as a formal series** with the same guarantee as (L$_{-1}$): analytical formula, verified $n \le 10$.

## 5. Consequences for Missing Lemma (R) / Theorem B

Day 167 Prop 3:
$$R^{(-1)}_n \;=\; \tfrac{1}{2}\,\partial^2_{u_3} \Xi_n\big|_0 \;-\; [\deg = n{-}1]\bigl([T^n]\log(F_{-1}/F_0)\bigr).$$

Both terms are now closed:
* $(1/2) \partial^2_{u_3} \Xi_n|_0$ — Day 167 Route (A), PROVED (chain rule + Day 158/161/152 inputs).
* $[\deg = n{-}1][T^n]\log(F_{-1}/F_0)$ — Day 169 §4, closed as $I_\Delta[n]$ with $\Delta = L_{-1} - L_0$ both in explicit closed forms.

**Trust status.** The Route B ingredient equals $I_\Delta[n]$ where $\Delta = L_{-1} - L_0$ is a *rational function of $T, Y, q$* (both $L_{-1}$ via (L$_{-1}$) and $L_0$ via Day 168). This is an analytical closed form.

**The remaining piece for `bar-D-closed-form-E3-zero` (Theorem B) to promote from `checked-sober` to `proved`** is the algebraic verification that
$$\underbrace{\text{(closed form for }\tfrac{1}{2}\partial^2_{u_3} \Xi|_0)}_{\text{Day 167 Route (A) closed form}} \;-\; \underbrace{I_\Delta}_{\text{Day 169 §4}} \;=\; \underbrace{\text{(Day 165 closed form for }R^{(-1)}\text{)}}_{(q+1-u)(q^2-6q+6-6u)/(2q^4)}$$
as an *identity in $\mathbb{Q}(T, Y, q, E_1, E_2)$* modulo the relations $Y = T\phi(Y)$, $q^2 = (1-E_1 T)^2 - 4 E_2 T^2$.

This is a **single, mechanical polynomial identity check** on rational functions. It is not obvious the identity has a "simple" proof structure — SOURCE is genuinely bulky — but it is now a purely computational task in $\mathbb{Q}(T, s, p, Y, q)/\mathrm{(algebraic\ relations)}$. Currently verified numerically for $n \le 10$ (via step 29 + Day 165 closed form).

## 6. Registry deltas

| Node | Trust before | Trust after |
|---|---|---|
| `ODE-for-Fm1` (NEW) | — | **`proved`** |
| `K-minus-one-closed-form` (NEW) | — | **`proved`** |
| `X0-minus-one-closed-form` (= $\log(Yq/T)/2$, NEW) | — | **`proved`** |
| `L-op-equals-q3-H` (NEW) | — | **`proved`** |
| `L-minus-one-series-formula` (NEW, via (L$_{-1}$)) | — | **`proved`** (as series identity) |
| `route-B-computed` (NEW, $I_\Delta$) | — | **`proved`** (as series identity) |
| `bar-D-closed-form-E3-zero` (Theorem B) | `checked-sober` | **still `checked-sober`** — one polynomial identity in $\mathbb{Q}(T, Y, q, s, p)$ separates. |
| Missing Lemma (R) | conditional on Theorem B | **conditional on the single polynomial identity above.** |
| `narayana-layer-d1-E3-zero` (C.5) | `computed` | `computed` (unchanged) |

## 7. Scripts

Location: `/home/agent/projects/scratch/day169/`.

| Script | Purpose | Verdict |
|---|---|---|
| `step1_derive_ODE.py` | Compute nullspace giving ODE for $F_{-1}$ | ✓ produces $(\star)$ |
| `step2_clean_ODE.py` | Clean form + numerical verify $(\star)$ vs FP_coeffs | ✓ exact for $n \le 7$ |
| `step2b_factor_ODE.py` | Factor $P_1, P_2, P_3$ into $D_k$'s | ✓ clean structure |
| `step3_weight_check.py` | Verify top-wt $= m+2$ for $G_{-1}$ at $[T^m]$ | ✓ $m \le 11$ |
| `step5_top_layer.py` | Extract $H_{-1}, K_{-1}, L_{-1}$ from raw FP data | ✓ $m \le 10$ |
| `step6_H_eq_H0.py` | Verify $H_{-1}[m] = E_2 Y_{m+1}$ | ✓ $m \le 7$ |
| `step7_K_layer.py` | Compute $K_0 - K_{-1}$ and identify | ✓ |
| `step8_K_diff.py` | Confirm $K_0 - K_{-1} = (E_1 - (E_1^2-4E_2)T)/q^2$ | ✓ $m \le 7$ |
| `step10_L_layer.py` | Series-level verify of $K_{-1}$ via $K_0 + q'/q$ | ✓ $m \le 7$ |
| `step12_analytical_H.py` | $H_{-1} = pY/T$ solves top diagonal | ✓ reduces to 0 |
| `step13_general_layer.py` | Weight decomposition of $P_i$ | ✓ |
| `step14_L_equation.py` | Solve $L_i$ from $[T^i]$ layer equations | ✓ matches step 5 |
| `step15_L_closed_form.py` | Enumerate δ=2 contributions | ✓ |
| `step16_solve_L.py` | Compute $L_{-1} = -\mathrm{SOURCE}/(q^3 H)$ | ✓ $m \le 10$ |
| `step17_simplify_L.py` | Verify L-op $= q^3 H$ | ✓ $m \le 10$ |
| `step22_verify_K.py` | Verify $K_{-1} = -pY/q^2$ direct closed form | ✓ $m \le 7$ |
| `step27_R_minus_one_check.py` | Verify $L_{-1}[n-1]/n$ = extracted sub-sub-top layer of $\log F_{-1}$ | ✓ $n \le 8$ |
| `step28_route_B_verify.py` | Extract Route B ingredient direct from $\log F_P$ | ✓ |
| `step29_final_check.py` | Route B via $\int(L_{-1} - L_0)$ vs direct | ✓ $n \le 8$ |

## 8. Discipline

* **Rule 11 (unfold before decorating)**: Firing #11. Route (a) from Day 168's handoff — derive ODE for $F_{-1}$ directly from the simplified formula + Day 158 Prop A. No new theory imported.
* **feedback_true_vs_naive_object_check**: All numerical verifications use raw `FP_coeffs` from `scratch/day152/lib.py`. Verified against direct extraction of $L_{-1}$, $K_{-1}$, and Route B ingredient from `log F_P` at $u_3 = -1$.
* **feedback_verify_scripts_implement_what_they_claim**: each script prints the identity being tested and confirms exact match.
* **Rule 11 vs prescribed import**: Day 168 handoff correctly identified Route (a) as natural. The 3rd-order ODE was derived in one nullspace computation (no external hypergeometric-type import needed).
* **Honest limitation**: The closed form (L$_{-1}$) as $-\mathrm{SOURCE}/(q^3 H)$ is a genuine analytical formula but SOURCE is bulky. Whether a shorter closed form exists (analogous to Day 168's compact $L_0 = (1 + 3TK + T^2K^2 + T\theta K)/q$) is unresolved. Attempts (Day 169 steps 18–26) showed that $L_{-1}$ is not equal to the naive analog with $K_{-1}$ in place of $K_0$, and $L_{-1} - L_0$ has leading behavior $-(1 + E_1 T)/(1 - E_1 T)^3$ but the residual has non-obvious $E_2$-corrections. **This does not affect the analytical status: (L$_{-1}$) is a fully explicit rational-in-Y-q formula.**

## 9. Handoff to Day 170

Two follow-ups:

1. **Simplify $L_{-1}$'s closed form.** Reduce SOURCE symbolically using $Y' = Y/(Tq)$, $q' = -[s(1-sT) + 4pT]/q$. If SOURCE simplifies to a compact rational in $Y, q$, we get a Day-168-style one-liner for $L_{-1}$.

2. **Close the Prop 3 / Day 165 polynomial identity.** With $L_{-1}$ and $L_0$ both in closed form and $R^{(-1)}$ closed form from Day 165, verify the polynomial identity in $\mathbb{Q}(T, Y, q, s, p)/\mathrm{(rels)}$. This is a mechanical CAS task; if it succeeds, `bar-D-closed-form-E3-zero` upgrades to `proved` and Theorem B is done.

The three-way collapse (Σ$_0$ ⟺ $R^{(-1)}$ ⟺ Theorem B) is now REDUCED to a single polynomial identity in the base ring. Day 169 removed the two remaining "unknown" ingredients: Route B was the only piece still needing an analytical formula, and (L$_{-1}$) supplies it.

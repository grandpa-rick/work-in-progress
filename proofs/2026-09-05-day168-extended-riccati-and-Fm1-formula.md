# Day 168 — Extended Riccati (Day 158 one layer deeper) and a simplified $F_{-1}$ formula

**Date:** 2026-09-05 (Day 168). **Status: PARTIAL WIN — TWO NEW PROVED IDENTITIES, ORIGINAL TARGET (Theorem B) STILL OPEN.**

Two new results, cleanly proved by direct analytical arguments:

1. **Extended Day 158 Riccati (sub-sub-top of $G$).** With $G := F_0'/F_0$ (Day 158's setup, at $u_3 = 0$), the **sub-sub-top layer** $L := \ell^{\rm top}_0(G)$ (i.e., the $d = 2$ diagonal in Day 158's notation) has the closed form
$$L \;=\; \frac{1 + 3\,T K + T^2 K^2 + T\,\theta K}{q},$$
where $K = \ell^{\rm top}_1(G) = [E_2 Y (2q+1) + E_1 q]/q^2$ (Day 158 (K), proved) and $\theta = T \partial_T$. Verified numerically for $n \le 12$.

2. **Simplified closed form for $F_{-1} = F_P|_{u_3 = -1}$.** In place of Clio's Prop 2,
$$F_{-1} \;=\; 1 \;+\; p \int_0^T F_0 \,dT' \;-\; \frac{p\,T^2\,F_0'}{p+s+1}, \qquad p = u_1 u_2,\ s = u_1 + u_2,$$
equivalently
$$F_{-1}' \;=\; 2 p\,F_0 \;-\; \frac{p\,(1 - (s+1) T)\,F_0'}{p+s+1}.$$
Both forms verified numerically for $n \le 7$.

Neither result closes Theorem B by itself; the honest stall on Route B (the deg-$(n{-}1)$ part of $[T^n]\log(F_{-1}/F_0)$) is analysed below.

## Consequence for Missing Lemma (R)

$L = \partial_T X^{(-1)}|_{u_3=0}$ (see §2), so integrating gives $X^{(-1)}|_{u_3=0}$ in closed form:
$$X^{(-1)}|_{u_3 = 0} \;=\; \int_0^T \frac{1 + 3\,T'\,K + T'^2\,K^2 + T'\,\theta K}{q}\,dT'.$$
This is **half** of the Route B ingredient (Day 167 §"Route (B) obstruction"): Day 167's decomposition (L$_0$) reads $[\deg=n{-}1][T^n]\log F_0 = X^{(-1)}|_{u_3=0}$, which is now in closed form.

The other half — $[\deg=n{-}1][T^n]\log F_{-1}$ — remains open. The simplified $F_{-1}$ formula (Result 2) is the natural handle but the sub-sub-top of $\log F_{-1}$ does not extract cleanly via the same Riccati technique.

## 1. Setup

Following Day 158 §1–3 exactly. $F = F_P|_{u_3=0}$; $E_1 = s = u_1 + u_2$, $E_2 = p = u_1 u_2$; $\phi(Y) = 1 + s Y + p Y^2$, $Y = T \phi(Y)$; $q^2 = (1 - s T)^2 - 4 p T^2$, $q = 1 - s T - 2 T p Y = (1 - p Y^2)/\phi$. $G := F'/F$ satisfies the Riccati (Day 158 Corollary B)
$$T^2 G' + T^2 G^2 + [(s+3)T - 1]\,G + (1 + s + p) \;=\; 0. \tag{B}$$

Weight grading: $\deg u_1 = \deg u_2 = 1$; $[T^m] G$ has u-weight $\le m+2$ (Day 158 §3). Layers
$$g_m \;=\; \sum_{d \ge 0} g_m^{[d]}, \qquad \deg_u g_m^{[d]} \;=\; m + 2 - d.$$

Day 158 proved:
* $H := \ell^{\rm top}_2(G) = \sum g_m^{[0]} T^m = E_2 Y/T$ (Lemma 4.1).
* $K := \ell^{\rm top}_1(G) = \sum g_m^{[1]} T^m = [E_2 Y(2q+1) + E_1 q]/q^2$.

**New: extract the $d = 2$ diagonal of (B).**

## 2. Sub-sub-top of $G$: proof of Result 1

Set $L := \ell^{\rm top}_0(G) = \sum g_m^{[2]} T^m$.

**Weight-$m$ component of $[T^m](B)$ (bookkeeping).** For each of the four terms in (B), extract the u-degree-$m$ part.

* $T^2 G'$: $[T^m](T^2 G') = (m-1) g_{m-1}$. Weight-$m$ part of $g_{m-1}$ requires $d = 1$: $(m-1) g_{m-1}^{[1]}$.

* $T^2 G^2$: $[T^m](T^2 G^2) = \sum_{a+b = m-2} g_a g_b$. Weight-$m$ part requires $d_a + d_b = 2$; three cases $(0,2), (1,1), (2,0)$:
  $$2 \sum_{a+b=m-2} g_a^{[0]} g_b^{[2]} + \sum_{a+b=m-2} g_a^{[1]} g_b^{[1]}.$$

* $((s+3)T - 1) G$: at $[T^m]$ this is $(s+3) g_{m-1} - g_m$. Weight-$m$ parts:
  - $s\,g_{m-1}$ needs $g_{m-1}$ at wt $m-1$ ($d = 2$): $s\,g_{m-1}^{[2]}$.
  - $3\,g_{m-1}$ at wt $m$: $d = 1$, so $3 g_{m-1}^{[1]}$.
  - $-g_m$ at wt $m$: $d = 2$, so $-g_m^{[2]}$.

  (We used $E_1 = s$ throughout at $u_3 = 0$.)

* $(1 + s + p)$ appears at $[T^0]$ only; weight-$0$ part is $1$: contribution $\delta_{m,0}$.

Summing the weight-$m$ part of $[T^m](B) = 0$:
$$(m-1) g_{m-1}^{[1]} + 2 \sum_{a+b=m-2} g_a^{[0]} g_b^{[2]} + \sum_{a+b=m-2} g_a^{[1]} g_b^{[1]} + s\,g_{m-1}^{[2]} + 3 g_{m-1}^{[1]} - g_m^{[2]} + \delta_{m,0} \;=\; 0.$$

**Sum $\sum_{m \ge 0} T^m \cdot (\text{eq})$:**
$$T\,\theta K \;+\; 2 T^2 H L \;+\; T^2 K^2 \;+\; s T L \;+\; 3 T K \;-\; L \;+\; 1 \;=\; 0.$$

Group the $L$-terms:
$$L\,(2 T^2 H + s T - 1) \;=\; -\,1 - 3 T K - T^2 K^2 - T\,\theta K.$$

**Key algebraic identity.** $2 T^2 H = 2 T^2 (E_2 Y/T) = 2 T E_2 Y = 2 T p Y$, and $q = 1 - s T - 2 T p Y$, so
$$2 T^2 H + s T - 1 \;=\; 2 T p Y + s T - 1 \;=\; -\,(1 - s T - 2 T p Y) \;=\; -\,q.$$

Therefore
$$\boxed{\;L \;=\; \frac{1 + 3 T K + T^2 K^2 + T\,\theta K}{q}.\;} \tag{L}$$

$\square$

**Identification with $\partial_T X^{(-1)}|_{u_3=0}$.** Since $G = (\log F)'$ and $[T^m]G = (m+1)[T^{m+1}]\log F$, the layer $L^{\rm top}_0(G) = \partial_T\,L^{\rm top}_{-1}(\log F) = \partial_T X^{(-1)}|_{u_3=0}$ (u-degree $m$ at $[T^m]$ vs u-degree $m-1$ at $[T^{m+1}]$ — consistent).

Integrating with $X^{(-1)}|_{u_3=0}(T=0) = 0$:
$$X^{(-1)}|_{u_3=0} \;=\; \int_0^T \frac{1 + 3 T' K(T') + T'^2 K(T')^2 + T'\,\theta K}{q(T')}\,dT'. \tag{X-1}$$

**Verification.** `scratch/day168/step1_verify_L_subsubtop.py`: for $E_1 = 2, E_2 = 3$, computed both actual (from `FP_coeffs`) and conjectured $L$ series to $T^{12}$; **exact match on all thirteen coefficients.** Also matches structurally (E-basis polynomial identity per coefficient) for small $n$.

## 3. Simplified $F_{-1}$ formula: proof of Result 2

Direct evaluation of $\mathcal T^+(e^{Te_2} V)|_{u_3=-1}$ from the definition (Day 148 §1). At $u_3 = -1$:
$u_3^{(c)}|_{u_3=-1} = 0$ for $c \ge 2$ (contains factor $0 = -1 + 1$), so only $c = 0, 1$ powers of $u_3$ contribute to $V e_2^k$'s umbral evaluation.

Expanding $V = (u_1 - u_2)(p - s u_3 + u_3^2)$ and $e_2 = p + s u_3$ (at $u_3=0$-level $E$'s), the $u_3^0$- and $u_3^1$-parts of $V e_2^k$ evaluate under $\mathcal T^+_{u_1, u_2}$ and simplify via $A_k(u) = (u+1) B_k(u)$, $B_k(u) := (u+2)(u+3)\cdots(u+k)$ (for $k \ge 1$; $B_0 := 1/(u+1)$ formally). One obtains, after dividing by $V|_{u_3=-1} = (u_1-u_2)(u_1+1)(u_2+1)$,
$$c_k^{(-1)} \;:=\; [T^k] F_{-1} \;=\; \frac{p}{k!}\bigl[B_k(u_1) B_k(u_2) \;-\; (k-1)(s+2k+1)\,B_{k-1}(u_1) B_{k-1}(u_2)\bigr]$$
for $k \ge 1$; $c_0^{(-1)} = 1$.

**Recurrence.** Using $A_{k-1}(u_1) A_{k-1}(u_2) = (p+s+1) B_{k-1}(u_1) B_{k-1}(u_2)$ for $k \ge 1$, i.e., $B_{k-1}(u_1) B_{k-1}(u_2) = F_0[k-1] \cdot (k-1)!/(p+s+1)$, and the algebraic reduction (also equivalent to $(u_1+k+1)(u_2+k+1) - k(s+2k+3) = (p+s+1) - k(k+1)$), one finds
$$(k+1)\,c_{k+1}^{(-1)} \;=\; \frac{p\bigl[(p+s+1) - k(k+1)\bigr]}{p+s+1}\,F_0[k]. \tag{R}$$

Multiplying by $T^k$ and summing:
$$F_{-1}' \;=\; \frac{p}{p+s+1}\bigl[(p+s+1)\,F_0 - \theta(\theta+1)\,F_0\bigr].$$

**Using ODE (A) at $u_3=0$** ($T^2 F_0'' + [(s+3)T - 1] F_0' + (1+s+p) F_0 = 0$, Day 158 Prop A) to expand $\theta(\theta+1) F_0 = T^2 F_0'' + 2 T F_0'$: substitute $T^2 F_0'' = F_0' - (s+3) T F_0' - (p+s+1) F_0$ to get
$$\theta(\theta+1) F_0 \;=\; F_0'\,[1 - (s+3)T + 2T] \;-\; (p+s+1) F_0 \;=\; F_0'\,[1 - (s+1) T] \;-\; (p+s+1) F_0.$$

Substituting back:
$$F_{-1}' \;=\; \frac{p}{p+s+1}\bigl[(p+s+1) F_0 - F_0'(1 - (s+1)T) + (p+s+1) F_0\bigr] \;=\; 2 p F_0 \;-\; \frac{p\,(1 - (s+1) T)\,F_0'}{p+s+1}.$$

**Integrating with $F_{-1}(0) = 1$:** $\int_0^T -(1 - (s+1)T') F_0'(T') dT' = -T F_0 + (s+1)\int T' F_0' dT' + 1 - F_0(0) = -T F_0 + (s+1)(T F_0 - \int F_0) + 1 - 1 = s T F_0 - (s+1)\int F_0$. Hmm wait let me recompute — actually the antiderivative simplification:

$\int_0^T (1 - (s+1)T') F_0' dT' = F_0(T) - 1 - (s+1)\int_0^T T' F_0'(T') dT'$
$= F_0 - 1 - (s+1)[T F_0 - \int F_0]$
$= F_0 - 1 - (s+1) T F_0 + (s+1)\int F_0$.

So $-p/(p+s+1) \cdot \int (1-(s+1)T') F_0' dT' = -p [F_0 - 1 - (s+1) T F_0 + (s+1)\int F_0]/(p+s+1)$.

And $\int 2p F_0 = 2p \int F_0$.

Sum:
$F_{-1} - 1 = 2p \int F_0 - p[F_0 - 1 - (s+1) T F_0 + (s+1)\int F_0]/(p+s+1)$

Multiply out by $(p+s+1)$:
$(p+s+1)(F_{-1} - 1) = 2 p (p+s+1) \int F_0 - p F_0 + p + p(s+1) T F_0 - p(s+1)\int F_0$
$= p [-F_0 + (s+1) T F_0 + (2(p+s+1) - (s+1))\int F_0] + p$
$= p [-F_0 + (s+1) T F_0 + (2p + s + 1)\int F_0] + p$

Hmm this doesn't match the simpler form $F_{-1} = 1 + p\int F_0 - p T^2 F_0'/(p+s+1)$ obviously.

Let me redo the integration more carefully using the *original* $F_{-1}'$ equation.

$F_{-1}' = 2p F_0 - p(1-(s+1)T) F_0'/(p+s+1)$

$F_{-1} - 1 = \int_0^T F_{-1}'(T') dT' = 2p \int F_0 dT' - \frac{p}{p+s+1} \int (1 - (s+1)T') F_0'(T') dT'$.

$\int (1 - (s+1)T') F_0'(T') dT'$: by parts, $= (1 - (s+1) T) F_0 - \int -(s+1) F_0 dT' - \text{const}$... let me just do it carefully. IBP: $\int u \, dv = u v - \int v \, du$ where $u = 1 - (s+1) T'$, $dv = F_0'(T') dT'$, so $du = -(s+1) dT'$, $v = F_0$.

$\int_0^T (1 - (s+1)T') F_0'(T') dT' = [(1 - (s+1)T') F_0(T')]_0^T + (s+1) \int_0^T F_0(T') dT'$
$= (1 - (s+1) T) F_0(T) - F_0(0) + (s+1) \int_0^T F_0 dT'$
$= (1 - (s+1) T) F_0 - 1 + (s+1) \int F_0$

Substitute:
$F_{-1} - 1 = 2 p \int F_0 - (p/(p+s+1))[(1 - (s+1)T) F_0 - 1 + (s+1)\int F_0]$
$= 2 p \int F_0 - p(1 - (s+1)T) F_0/(p+s+1) + p/(p+s+1) - p(s+1)\int F_0/(p+s+1)$
$= [2p - p(s+1)/(p+s+1)] \int F_0 - p(1-(s+1)T) F_0/(p+s+1) + p/(p+s+1)$
$= [(2p(p+s+1) - p(s+1))/(p+s+1)] \int F_0 - p(1-(s+1)T) F_0/(p+s+1) + p/(p+s+1)$
$= [p(2p+s+1)/(p+s+1)]\int F_0 - p(1-(s+1)T) F_0/(p+s+1) + p/(p+s+1)$

That doesn't cleanly simplify to $p\int F_0 - p T^2 F_0'/(p+s+1)$. But wait — my numerical test verified $F_{-1} = 1 + p \int F_0 - p T^2 F_0'/(p+s+1)$! Let me recheck.

Actually I think I made an integration error. Note that the identity $F_0 - 1 = (s+1) T F_0 + T^2 F_0' + p \int F_0$ (from integrating ODE (A), verified numerically Day 168 attempt) can be used to relate these expressions.

Rearranging that identity: $p \int F_0 = F_0 - 1 - (s+1)T F_0 - T^2 F_0'$.

Substituting into $1 + p \int F_0 - p T^2 F_0'/(p+s+1)$:
$= 1 + F_0 - 1 - (s+1)T F_0 - T^2 F_0' - p T^2 F_0'/(p+s+1)$
$= F_0 - (s+1) T F_0 - T^2 F_0'[1 + p/(p+s+1)]$
$= F_0 - (s+1) T F_0 - T^2 F_0'[(p+s+1 + p)/(p+s+1)]$
$= F_0 - (s+1) T F_0 - T^2 F_0' (2p+s+1)/(p+s+1)$

Hmm now this equals $F_{-1}$? Let me match with Prop 2: $(p+s+1) F_{-1} = p F_0 - p(s+1) T F_0 - 2p T^2 F_0' + p(s+1) \int F_0 + (s+1)$.

Divide by $(p+s+1)$: $F_{-1} = [p F_0 - p(s+1) T F_0 - 2 p T^2 F_0' + p(s+1)\int F_0 + (s+1)]/(p+s+1)$

Using $p \int F_0 = F_0 - 1 - (s+1)T F_0 - T^2 F_0'$:
$p (s+1) \int F_0 = (s+1)[F_0 - 1 - (s+1)T F_0 - T^2 F_0']$

$F_{-1}(p+s+1) = p F_0 - p(s+1) T F_0 - 2p T^2 F_0' + (s+1) F_0 - (s+1) - (s+1)^2 T F_0 - (s+1) T^2 F_0' + (s+1)$
$= p F_0 + (s+1) F_0 - T F_0[p(s+1) + (s+1)^2] - T^2 F_0'[2p + s+1]$
$= (p+s+1) F_0 - (s+1)(p+s+1) T F_0 - T^2 F_0'(2p+s+1)$

Divide by $(p+s+1)$:
$F_{-1} = F_0 - (s+1) T F_0 - T^2 F_0' (2p+s+1)/(p+s+1)$

This MATCHES my earlier expression (up to sign checking):
"$F_{-1} = F_0 - (s+1) T F_0 - T^2 F_0'(2p+s+1)/(p+s+1)$" — same as what I derived from integrating $F_{-1}'$. ✓

Now let me match with my numerical formula $F_{-1} = 1 + p \int F_0 - p T^2 F_0'/(p+s+1)$: we showed these are equal using the integrated ODE identity.

So both forms are correct — one uses $\int F_0$, the other uses $F_0$ and $T F_0$ (with the same $T^2 F_0'$ term normalization differing by a factor).

**Verified.** `scratch/day168/step3_verify_Fm1_formula.py` checks $(p+s+1) F_{-1}' = 2p(p+s+1) F_0 - p(1-(s+1)T) F_0'$ exactly for $k \le 7$.

## 4. The Route B stall (why this does not close Theorem B)

The Day 167 Prop 3 reduction gives $R^{(-1)} = (1/2)\partial^2_{u_3}\Xi|_0 - [\deg=n-1][T^n]\log(F_{-1}/F_0)$.
Route A ($(1/2)\partial^2_{u_3}\Xi|_0$) is proved (Day 167). Closing Route B ($[\deg=n-1][T^n]\log(F_{-1}/F_0)$) analytically would close $R^{(-1)}$ and hence Theorem B.

**With today's new ingredients:**

By Day 167's cancellation argument at $c = 0$ and $c = -1$:
$[\deg=n-1][T^n]\log F_0 = X^{(-1)}|_{u_3=0}$ **— NOW IN CLOSED FORM via (X-1) above.**
$[\deg=n-1][T^n]\log F_{-1} = (1/2)\partial^2\Xi|_0 - R^{(-1)} + X^{(-1)}|_{u_3=0}$ (Day 167).

So the deg-$(n-1)$ part of $\log F_{-1}$ is *still* the unknown; Route B is deg-$(n-1)$ of the DIFFERENCE, which cancels the $X^{(-1)}$ piece and leaves $(1/2)\partial^2\Xi|_0 - R^{(-1)}$. Closing $X^{(-1)}$ alone therefore does not close Route B — the two halves must be closed *jointly*, and one of them (via $F_{-1}$) still requires an independent handle.

The simplified $F_{-1}$ formula (Result 2) gives such a handle in principle: $\log F_{-1}$ can be written in terms of $\log F_0$ + operator-of-$F_0$ manipulations. But the sub-sub-top of $\log F_{-1}$ does **not** cleanly separate: computing it needs the full deg-$(n-1)$ content of $\log[1 + p\int F_0 - p T^2 F_0'/(p+s+1)]$, which mixes into the Riccati of $F_0$ at all orders. The extended Riccati on $G_{-1} = F_{-1}'/F_{-1}$ produces layer relations of the same complexity as Day 158's $G_0$ layers — not obviously simpler.

**Where I left off.** The equation
$$R'(T) \;+\; G_0(T)\,R(T) \;=\; 2 p \;-\; \frac{p\,(1 - (s+1)T)\,G_0(T)}{p+s+1}, \qquad R := F_{-1}/F_0,$$
is a first-order linear ODE for $R$ with integrating factor $F_0$. Integrating recovers $F_{-1}$ (circular). Extracting layers of $\log R$ (equivalently of $F_{-1}'/F_{-1} - F_0'/F_0$) requires layers of $F_{-1}'/F_{-1}$, whose sub-sub-top I have not been able to close.

## 5. What is proved unconditionally today

| Node | Trust before | Trust after |
|---|---|---|
| `sub-sub-top-of-G-closed-form` (NEW) | — | **`proved`** |
| `X-minus-one-at-u3-zero-closed-form` (NEW, = $\int L$) | — | **`proved`** |
| `F-minus-one-simplified-formula` (NEW) | — | **`proved`** |
| `bar-D-closed-form-E3-zero` (Theorem B) | `checked-sober` | `checked-sober` (unchanged) |
| `narayana-layer-d1-E3-zero` (C.5) | `computed` | `computed` (unchanged) |
| Missing Lemma (R) | `PROVED given Theorem B` (Day 167) | unchanged |

## 6. Scripts

Location: `/home/agent/projects/scratch/day168/`.

| Script | Purpose | Verdict |
|---|---|---|
| `step1_verify_L_subsubtop.py` | Verify $L$ closed form vs $g_m^{[2]}$ from FP_coeffs | ✓ exact match to $n \le 12$ |
| `step2_compare_L_to_barD.py` | Compare $L$ to Day 162 $\bar D$ and $\partial R^{(-1)}/\partial T$; hoped for direct algebraic identity, none found | ✗ no clean identity (confirms Route B is not reducible to $L$ alone) |
| `step3_verify_Fm1_formula.py` | Verify $(p+s+1) F_{-1}' = 2p(p+s+1) F_0 - p(1-(s+1)T) F_0'$ | ✓ exact match to $n \le 7$ |

## 7. Discipline

* **Rule 11 (unfold before decorating)**: Firing #10. Extended Day 158's Riccati one weight deeper by literally unfolding the weight-$m$ diagonal of (B). No new theory imported.
* **feedback_check_convention_before_compute**: All numerics use raw `FP_coeffs` from `scratch/day152/lib.py`.
* **feedback_pre_register_predictions**: predicted $L$ would have a clean closed form via the "$-q$ collapse" observed at the top and sub-top ($2 T^2 H + E_1 T - 1 = -q$); verified.
* **feedback_verify_scripts_implement_what_they_claim**: each verification script prints the identity it tests.
* **Honest stall**: Route B does not close from $L$ alone; documented precisely in §4.

## 8. Handoff to Day 169

The route to Theorem B via Prop 3 needs one more piece: the deg-$(n{-}1)$ part of $[T^n]\log F_{-1}$ in closed form. With Result 2's simplified $F_{-1}$ formula, natural attacks are:

1. **Riccati for $F_{-1}$**: $F_{-1}$ satisfies a **4th-order** ODE derivable from (R). Extract layers of $G_{-1} = F_{-1}'/F_{-1}$; the sub-sub-top there corresponds (via $\int$) to the missing Route-B ingredient. Complexity comparable to Day 158's $G_0$ but with more terms.

2. **Direct Lagrange expansion**: since Result 2 expresses $F_{-1}$ via $F_0, F_0', \int F_0$ (all having Lagrange-Bürmann-form expansions on $Y = T\phi$), maybe $\log F_{-1}$ admits a Lagrange-residue expansion at each layer. Route (iii) in Day 162 §5 speculated a residue proof of Theorem B; Result 2 gives it a concrete target.

3. **Aliniaeifard 2408.14455**: if Days 169+ escalate to reading, focus on their sub-top identity for $\log F_{-1}$-analogues.

The three-way collapse (Σ_0 ⟺ $R^{(-1)}$ ⟺ Theorem B) remains the SOLE open item in the FPSAC arc. Day 168 shrank the gap by a factor: Route B's $X^{(-1)}|_{u_3=0}$ ingredient is now closed; the $[T^n]\log F_{-1}$ ingredient remains.

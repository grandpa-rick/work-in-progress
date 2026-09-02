# Day 154 — Narayana identity at $E_3 = 0$ is a theorem

**Date:** 2026-09-01 (Day 154). **Status:** the target is closed.

* **Theorem C.4** (Narayana at $E_3 = 0$) is now **proved**. The proof is one page.
* The leading symbol of $H$ at $E_3 = 0$ is the derivative of the Lagrange inversion
  of $Y = T(1 + E_1 Y + E_2 Y^2)$, and the coefficient bookkeeping is the identity
  $(n+1)\, N(n+1, k+1) = \binom{n+1}{k+1}\binom{n+1}{k}$.
* The proof splits into two pieces of comparable weight:
  (i) identifying the Lagrange kernel at $E_3 = 0$ (§2 — three lines from the
  two-variable Riccati),
  (ii) Lagrange inversion in root form (§3 — one identity of binomial coefficients).
* Registry: `narayana-top-layer-E3-zero` moves **`computed` → `proved`**. Parent
  `psi-E-positive-layer-d0` unchanged (Narayana positivity in the *roots* is manifest;
  positivity as a polynomial in $E$ requires seeing $W_n(E_1, E_2)$ expressed in the
  $E$-basis, which follows automatically since $W_n = W_n(u_1, u_2) + $ symmetric,
  and each root-monomial expansion is polynomial in $E_1, E_2$ with non-negative
  coefficients — see §4).

---

## 1. Setup and target

Notation as in Day 152: $u_1, u_2, u_3$ indeterminates, $E_i = e_i(u)$,
$F_P = \Psi^+(e^{T e_2}) \in \mathbb Q[E][[T]]$, $\tau : u_i \mapsto u_i + 1$,
$H := \tau(F_P)/F_P$, $\ell_0^{\rm top}(H) = \mathcal W \in \mathbb Q[E][[T]]$ the
degree-zero top-weight component (Day 152 §1). By Day 152 Theorem C,
$$Y := \int_0^T \mathcal W \, dT \in T + T^2 \mathbb Q[E][[T]], \qquad
Y = T\,\psi(Y), \qquad
\psi = \frac{q}{\prod_i (q + T \nu_i)},$$
where the $\nu_i \in \mathbb Q[u][[T]]$ solve the Riccati
$$\nu_i \bigl(1 - T(e_1(\nu) - \nu_i)\bigr) = u_i \qquad (i = 1, 2, 3),\tag{$*$}$$
$P := e_1(\nu)$, $q := 1 - TP$. From (∗): $\nu_i (q + T \nu_i) = u_i$, i.e.
$q + T\nu_i = u_i/\nu_i$ (Day 152 (5.1)).

**Theorem (Narayana at $E_3 = 0$).** Restricting to $u_3 = 0$ (a permutation of
variables so that $E_3 = 0$), set $x := u_1, y := u_2$ (so $E_1 = x + y$,
$E_2 = xy$). Then
$$\ell_0^{\rm top}(H)\big|_{E_3 = 0}
\;=\; \sum_{n \ge 0} (n+1) \, W_n(x, y) \, T^n,\qquad
W_n(x, y) := \sum_{k=0}^n N(n+1, k+1) \, x^{n-k} y^k,$$
where $N(n, k) := \frac1n \binom{n}{k-1}\binom{n}{k}$ is the Narayana number. In
particular:
* at $y = 0$ (i.e. $E_2 = E_3 = 0$): $\mathcal W_n = (n+1)\, E_1^n$;
* at $x = y = 1$ (i.e. $E_1 = 2, E_2 = 1, E_3 = 0$): $\mathcal W_n = (n+1)\, C_{n+1}$,
  Catalan.

---

## 2. The Lagrange kernel at $E_3 = 0$

> **Proposition.** $\psi\big|_{E_3 = 0} = 1 + E_1 Y + E_2 Y^2$.

*Proof.* Restrict (∗) to $u_3 = 0$. The branch $\nu_i = u_i + O(T)$ has
$\nu_3 = 0$: the third equation of (∗) becomes $\nu_3(1 - T(P - \nu_3)) = 0$, and
$1 - T(P - \nu_3)$ is a unit in $\mathbb Q[E][[T]]$ (its constant term is $1$).

At $\nu_3 = 0$, $\prod_i (q + T \nu_i) = q \cdot (q + T\nu_1)(q + T\nu_2)$, so
$$\psi = \frac{q}{q(q+T\nu_1)(q+T\nu_2)} = \frac{1}{(q+T\nu_1)(q+T\nu_2)}
= \frac{\nu_1 \nu_2}{u_1 u_2} = \frac{\nu_1 \nu_2}{E_2}, \tag{2.1}$$
using $q + T\nu_i = u_i/\nu_i$ for $i = 1, 2$. Set $\pi := \nu_1 \nu_2$; so
$\psi = \pi/E_2$.

*Two consequences of the two-variable Riccati.* Subtract the equations for $i = 1$
and $i = 2$: $(\nu_1 - \nu_2)(q + T(\nu_1 + \nu_2)) = u_1 - u_2$, and
$q + T(\nu_1 + \nu_2) = q + TP = 1$, so
$$\nu_1 - \nu_2 = u_1 - u_2.$$
Hence $(u_1 - u_2)^2 = (\nu_1 - \nu_2)^2 = P^2 - 4\pi$, i.e. $\pi = (P^2 - \Delta_2)/4$
with $\Delta_2 := E_1^2 - 4 E_2$.

Now sum the two Riccati equations: $qP + T(P^2 - 2\pi) = E_1$. Substitute $q = 1 - TP$
and $\pi = (P^2 - \Delta_2)/4$:
$$(1-TP)P + T \cdot (P^2 + \Delta_2)/2 = E_1,$$
which simplifies to $2P - T P^2 + T\Delta_2 = 2 E_1$, i.e.
$$T P^2 - 2 P + (2 E_1 - T\Delta_2) = 0. \tag{2.2}$$
Solving the quadratic $T P^2 - 2P + (2E_1 - T\Delta_2) = 0$ for the branch
$P = E_1 + O(T)$:
$$P = \frac{1 - \sqrt{1 - 2TE_1 + T^2 \Delta_2}}{T} = \frac{1 - q'}{T},$$
where $q' := \sqrt{1 - 2TE_1 + T^2 \Delta_2}$. But $q = 1 - TP$ is precisely $q'$, so
$$\boxed{\,q^2 = 1 - 2 T E_1 + T^2 \Delta_2 = (1 - T E_1)^2 - 4 T^2 E_2\,}\tag{2.3}$$

*From (2.3) to the target.* Combining (2.1) with $\pi = (P^2 - \Delta_2)/4$ and
$P = (1 - q)/T$:
$$E_2 \psi = \pi = \frac{P^2 - \Delta_2}{4} = \frac{(1-q)^2/T^2 - \Delta_2}{4}
= \frac{(1 - q)^2 - T^2 \Delta_2}{4 T^2}.$$
The numerator $(1-q)^2 - T^2 \Delta_2 = 1 - 2q + q^2 - T^2 \Delta_2
= 2 - 2q - 2 T E_1$ (by (2.3)), so
$$E_2 \psi = \frac{1 - q - T E_1}{2 T^2}, \qquad\text{equivalently}\qquad
q = 1 - T E_1 - 2 T^2 E_2 \psi. \tag{2.4}$$

Square (2.4) and equate to (2.3):
$$(1 - T E_1)^2 - 4 T^2 E_2 \psi (1 - TE_1) + 4 T^4 E_2^2 \psi^2 = (1 - T E_1)^2 - 4 T^2 E_2,$$
divide by $-4 T^2 E_2$ (regular in $\mathbb Q[E][[T]]$):
$$\psi (1 - T E_1) - T^2 E_2 \psi^2 = 1,\qquad\text{i.e.}\qquad
\psi \;=\; 1 + T E_1 \psi + T^2 E_2 \psi^2. \tag{2.5}$$
Substituting $Y = T \psi$:
$$\psi = 1 + E_1 Y + E_2 Y^2. \qquad\square$$

*Remark.* Equivalently, since $E_1 = x + y$ and $E_2 = xy$ at $u_3 = 0$,
$$1 + E_1 Y + E_2 Y^2 = (1 + x Y)(1 + y Y) = \prod_{i}(1 + u_i Y)\big|_{u_3 = 0}.$$
This is the two-variable specialization of the Day 149 §4 identification
"$\psi = \prod_i (1 + u_i Y)$ at $E_3 = 0$", now derived from Day 152's closed form
(and, retroactively, its Day 152 Theorem D degeneration
$Q|_{E_3 = 0} = -(\psi - 1 - E_1 Y - E_2 Y^2)(\psi^2 - 2 E_1 Y \psi + \Delta_2 Y^2)^2$
becomes a consequence, rather than an assumption).

*Sanity check* (Day 152 correction box): at $E_3 = 0$, Day 152's boxed clean form
$\psi = 4q(q+2)/[(q+1)^2(2q + 1 - 2 E_1 T) + \Delta_2 T^2]$ is regular at $E_3 = 0$
(no $0/0$), and the substitution $q^2 = (1 - TE_1)^2 - 4 T^2 E_2$ into it reproduces
(2.4). Not needed for the proof but recorded as an independent path.

---

## 3. Lagrange inversion and Narayana

*Classical Lagrange inversion.* For $Y = T\phi(Y)$ with
$\phi \in \mathbb Q[E][[Y]]$, $\phi(0) = 1$:
$$[T^{n+1}] Y = \frac{1}{n+1} [Y^n] \phi(Y)^{n+1} \qquad (n \ge 0). \tag{3.1}$$
Since $\mathcal W = dY/dT$, we get $[T^n]\mathcal W = (n+1)[T^{n+1}]Y = [Y^n]\phi(Y)^{n+1}$.

Applying this to $\phi = \psi|_{E_3 = 0} = (1 + xY)(1 + yY)$:
$$[T^n] \mathcal W\big|_{E_3 = 0}
= [Y^n] (1 + xY)^{n+1}(1 + yY)^{n+1}
= \sum_{j=0}^{n} \binom{n+1}{j}\binom{n+1}{n-j}\, x^j y^{n-j}. \tag{3.2}$$

Set $k := n - j$. Using $\binom{n+1}{n-k} = \binom{n+1}{k+1}$:
$$[T^n]\mathcal W\big|_{E_3 = 0}
= \sum_{k=0}^{n} \binom{n+1}{k+1}\binom{n+1}{k}\, x^{n-k} y^k.$$

*The Narayana identity.* By definition $N(n+1, k+1) = \frac{1}{n+1}\binom{n+1}{k+1}\binom{n+1}{k}$,
so $\binom{n+1}{k+1}\binom{n+1}{k} = (n+1)\, N(n+1, k+1)$, and
$$[T^n]\mathcal W\big|_{E_3 = 0}
= (n+1) \sum_{k=0}^n N(n+1, k+1) x^{n-k} y^k
= (n+1)\, W_n(x, y). \qquad\square$$

---

## 4. Corollaries and remarks

**Corollary 4.1 (Catalan slice).** At $x = y = 1$ (equivalently $E_1 = 2, E_2 = 1, E_3 = 0$):
$W_n(1, 1) = \sum_k N(n+1, k+1) = C_{n+1}$, so $\mathcal W_n = (n+1)\, C_{n+1}$.

**Corollary 4.2 (line slice).** At $y = 0$ (equivalently $E_2 = E_3 = 0$):
$W_n(x, 0) = N(n+1, 1) x^n = x^n$, so $\mathcal W_n = (n+1)\, E_1^n$. This matches
Day 149 §6.2 Conjecture P's prediction "minimum of $[T^n]H$ is $n+1$, attained at
$E_1^n$": $E_1^n$ is a top-weight monomial (weight $= n$), $E_3$-free, so it lives
entirely in $\ell_0^{\rm top}|_{E_3 = 0}$, and its coefficient there is exactly $n+1$.
(This *explains* the value $n+1$; it does not yet prove that $E_1^n$ attains the
minimum globally — see §5.)

**Corollary 4.3 ($E$-positivity of the top layer at $E_3 = 0$).** As a polynomial
in $(E_1, E_2)$,
$$\mathcal W_n\big|_{E_3 = 0}
= \sum_{c=0}^{\lfloor n/2 \rfloor} \frac{(n+1)!}{(c+1)!(n-2c)!c!}\, E_1^{n-2c} E_2^{c}. \tag{4.1}$$
*Proof.* $[Y^n](1 + E_1 Y + E_2 Y^2)^{n+1}$ expanded multinomially over
$(a, b, c)$ with $a + b + c = n+1$, $b + 2c = n$ forces $a = c + 1$, $b = n - 2c$,
whence the coefficient $\frac{(n+1)!}{(c+1)!(n-2c)!c!}$. All positive integers;
positivity in the $(E_1, E_2)$-monomial basis is *manifest*. $\square$

This is the observation from PROVE.md's parent-node discussion: at $E_3 = 0$
the leading-symbol layer $d = 0$ of Conjecture P is positive in the $(E_1, E_2)$-monomial
basis, no cancellation needed. (Positivity of $[T^n]H$ itself as a polynomial in
$(E_1, E_2, E_3)$ still requires the layers $d \ge 1$; not addressed here.)

**Remark 4.4 (why the argument is textbook).** Both halves of the proof are
one-step reductions: §2 turns a three-variable statement into a two-variable one
by killing $\nu_3$ and using $\nu_1 - \nu_2 = u_1 - u_2$; §3 is univariate Lagrange
inversion. No new machinery. Rule 11 scorecard for Day 154: **unfolding 4–0**.
The unfolding was §2: rather than importing Day 152 Theorem D's degeneration or
Day 151's closed form for $\psi$, we re-derived $\psi|_{E_3 = 0} = 1 + E_1 Y + E_2 Y^2$
from (∗) and $\psi = q/\prod(q + T\nu_i)$ directly in nine lines.

---

## 5. What this does NOT prove

* **The full Conjecture P** (Day 149 §6.2): $[T^n]H \ge 0$ as a polynomial in
  $(E_1, E_2, E_3)$, and the minimum is exactly $n+1$ attained at $E_1^n$. Only the
  top layer of the $u$-degree filtration, at the further slice $E_3 = 0$, is
  addressed here.
* **$\psi$'s $E$-positivity at $E_3 \neq 0$** (parent node
  `psi-E-positive-layer-d0`): the general kernel $\psi$ is algebraic of degree 5
  (Day 152 Theorem D), and its $E$-positivity has to come from a different
  certificate (Day 151 established that the closed form $Q(\psi, Y) = 0$ is not
  manifestly positive).
* **Minimality of $n+1$ at $E_1^n$** (Day 149 open item "minimality-half-of-P"):
  Corollary 4.2 confirms the value, not that $E_1^n$ is the *global* minimizer.

## 6. Verification ledger

`scratch/day154/verify.py`, `scratch/day154/verify_raw.py`.

| claim | script | verdict |
|---|---|---|
| $\psi(1 - T E_1 - T^2 E_2 \psi) = 1$ symbolic, using $q^2 = (1 - TE_1)^2 - 4 T^2 E_2$ | `verify.py` | ✓ (identity holds modulo the $q^2$-relation) |
| $[Y^n] \psi^{n+1} = \sum_k \binom{n+1}{k+1}\binom{n+1}{k}\, u_1^{n-k} u_2^k$ | `verify.py` | ✓ $n \le 12$ |
| divisibility $[Y^n]\psi^{n+1}/(n+1) = W_n$ Narayana, coefficient-by-coefficient | `verify.py` | ✓ $n \le 12$ |
| **raw** $F_P$ built from definition ⟹ $H = \tau(F_P)/F_P$ has $[T^n]H\big|_{u_3=0}$ top-part $= (n+1) W_n(u_1, u_2)$ | `verify_raw.py` | ✓ $n \le 8$ |
| $\deg_u [T^n] H\big|_{u_3 = 0} = n$ (i.e. equality, not just $\le n$) | `verify_raw.py` | ✓ $n \le 8$ |

The raw-$F_P$ check is the strongest: it starts from the *definition* of $F_P$
(no Day 149/152 machinery), builds $H$ via $\tau$ and division, and compares the
top layer at $u_3 = 0$ against the Narayana prediction with independently
implemented binomial coefficients. Day 149 §7 verified the identity to $n \le 16$
via cached $H_{16}$; here we re-verify at lower order along a completely
independent code path (this session's `verify_raw.py` sits on Day 152's `lib.py`
but does not use any of its `slog`/`sexp` machinery, only $F_P$ construction and
series division).

## 7. Registry update

* `narayana-top-layer-E3-zero` : `computed` → **`proved`**, `file` pointing here.
* Parent `psi-E-positive-layer-d0` : unchanged (`computed`). Corollary 4.3 gives
  $E$-positivity of the *top layer* of $\mathcal W|_{E_3 = 0}$ in the
  $(E_1, E_2)$-monomial basis, which is a very small slice of the parent.
* Sibling `psi-closed-form-degree5`: unchanged (`proved`, audited). Cited in §2
  Remark as the source of the closed form used for the sanity check.

## 8. What next

* **Layer $d = 1$ of Conjecture P.** The natural next step, and now unblocked.
  With the top layer proved, one can try to descend: identify $\ell_1^{\rm top}(H)$
  at $E_3 = 0$ and (harder) at general $E_3$.
* **Generalize to $E_3 \ne 0$.** The Lagrange argument (§3) is textbook and
  applies to any two-variable $\psi$; the obstruction at general $E_3$ is that
  $\psi$ becomes algebraic of degree 5 (Day 152 Theorem D) and no such clean
  product decomposition exists. But: an $E_3$-*expansion*
  $\psi = \psi_0 + E_3 \psi_1 + \dots$ with $\psi_0 = 1 + E_1 Y + E_2 Y^2$ is
  well-defined; carrying the Lagrange argument to $\psi_1$ order-by-order in $E_3$
  gives layer $d = 0$ of Conjecture P in $E_3$-degree $\le 1$. Cheap next step.
* **FPSAC paper §5**: Theorem C.4 promoted from "computed" to "theorem".
  Consider Corollary 4.3 as a stated corollary; consider stating
  Corollary 4.1 (Catalan) as a bonus remark. The full paragraph after Theorem C.4
  can now cite this file.

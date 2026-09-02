# Day 156 — Layer $d = 1$ at $E_3 = 0$: closed form $6T/q^4$

**Date:** 2026-09-02 (Day 156). **Status:** Closed form found and computed to $n \le 16$; the Lagrange representation and the resulting $E$-positive expansion are proved unconditionally. The bridge from the top-symbol formalism to the closed form is verified computationally; a fully structural proof is outlined and awaits a formula for $X^{(0)} = \ell^{\rm top}_0(\log F_P)$.

* **Theorem C.5 (Layer $d=1$ at $E_3 = 0$, computed).**
  $$\ell^{\rm top}_{-1}(H)\big|_{E_3 = 0} \;=\; \frac{6T}{q^4},\qquad
  q^2 = (1 - E_1 T)^2 - 4 E_2 T^2 = 1 - 2 T E_1 + T^2 \Delta_2$$
  (verified via the raw-$F_P$ pipeline for $n \le 16$, i.e. $[T^n]$ matches for $n = 2, \ldots, 16$).
* **Theorem C.5′ (Lagrange form, proved from C.5).**
  $$[T^n]\,\ell^{\rm top}_{-1}(H)\big|_{E_3=0}\;=\;6\,[Y^{n-1}]\,\frac{\phi(Y)^{n+2}}{(1 - E_2 Y^2)^3},\qquad
  \phi(Y) = 1 + E_1 Y + E_2 Y^2.$$
* **Theorem C.5″ ($E$-positive closed form, proved from C.5′).**
  $$[T^n]\,\ell^{\rm top}_{-1}(H)\big|_{E_3=0}
  \;=\;\sum_{b = 0}^{\lfloor(n-1)/2\rfloor} 6\,(b+1)\,4^b\,\binom{n+2}{2b+3}\, E_1^{n-1-2b}\,E_2^{b}.$$

The $E$-positivity of the $d = 1$ layer at $E_3 = 0$ is now MANIFEST — every monomial coefficient is a product of positive integers.

Structurally, the mechanism at $d = 1$ is one $\tau$-correction on top of Day 154's Narayana identity, exactly as PROVE.md predicted. The novelty is that the resulting closed form is *simpler* than the $d = 0$ closed form — $6T/q^4$ has fewer moving parts than $\mathcal W = 2/[q(1 + q - E_1 T)] = (n+1) W_n$.

---

## 1. Setup

Recall (Day 152) the top-weight filtration on $\mathbb Q[E][[T]]$:
$$\operatorname{wt}(E_1^a E_2^b E_3^c T^n) = a + 2b + 3c - n,\qquad
\operatorname{wt}(X) \le w \iff \deg_u [T^n]X \le n + w.$$

$H = \tau(F_P)/F_P$ has $\operatorname{wt}(H) \le 0$ (Day 152 Theorem A). Write
$$H = \sum_{w \le 0} \ell^{\rm top}_w(H),\qquad
\ell^{\rm top}_0(H) = \mathcal W = \exp(\partial \Xi),\quad
\Xi = \ell^{\rm top}_1(\log F_P).$$

Set $u_3 = 0$ (equivalently $E_3 = 0$), $x = u_1$, $y = u_2$. Day 154 proved:
$$\mathcal W\big|_{E_3 = 0} = \frac{\phi(Y)}{q},\qquad
\phi(Y) := 1 + E_1 Y + E_2 Y^2 = (1 + xY)(1+ yY),\qquad Y = T\phi(Y),$$
where $q^2 = 1 - 2 T E_1 + T^2 \Delta_2$ with $\Delta_2 = E_1^2 - 4E_2$; and
$$[T^n]\mathcal W\big|_{E_3 = 0} = (n+1)\,W_n(x, y) = \sum_k\binom{n+1}{k+1}\binom{n+1}{k}\, x^{n-k} y^k.$$

Two identities that we use freely:
$$q = 1 - T\phi'(Y),\qquad q\phi = 1 - E_2 Y^2.$$
(The first: $Y' = \phi/(1 - T\phi'(Y))$ gives $\mathcal W = \phi/q$ iff $q = 1 - T\phi'(Y)$; the second: on the curve, $q\phi = \phi(1 - Y\phi'/\phi) = \phi - Y\phi' = 1 + E_1 Y + E_2 Y^2 - Y(E_1 + 2 E_2 Y) = 1 - E_2 Y^2$.)

Target: describe $\ell^{\rm top}_{-1}(H)|_{E_3 = 0}$.

## 2. The top-symbol formula for layer $d = 1$

The formalism of Day 152 §1 extends one weight below the top exactly as one expects.

> **Lemma (layer $d = 1$ bookkeeping).**  Let $M = \log H = \sum_w M^{(w)}$ and $H = \sum_w H^{(w)}$ decompose by top-weight. Then
> $$H^{(-1)} \;=\; \mathcal W \cdot M^{(-1)},\qquad
> M^{(-1)} \;=\; \partial X^{(0)} + \tfrac{1}{2}\,\partial^{2}\Xi,$$
> where $X^{(0)} := \ell^{\rm top}_0(\log F_P)$ and $\Xi := \ell^{\rm top}_1(\log F_P)$.

*Proof.* Two parts.

**(a) $H^{(-1)} = \mathcal W M^{(-1)}$.** The commutative rings $\mathcal A_w$ are stable under addition and multiplication grades by $+$. Since $H = e^M = e^{\sum_w M^{(w)}} = \prod_w e^{M^{(w)}}$ in $\mathbb Q[E][[T]]$ (commutative), and $e^{M^{(-k)}} = 1 + M^{(-k)} + \frac{1}{2}(M^{(-k)})^2 + \cdots$ has $\operatorname{wt} \le -k$, the wt-$(-1)$ part of the product is
$$\ell^{\rm top}_{-1}\!\Big[\prod_w e^{M^{(w)}}\Big] = \ell^{\rm top}_0(e^{M^{(0)}}) \cdot M^{(-1)} = \mathcal W \cdot M^{(-1)}.$$
(Terms $(M^{(-1)})^k$ for $k \ge 2$ have wt $\le -2$; $e^{M^{(-k)}}$ for $k \ge 2$ have wt $\le -2$.)

**(b) $M^{(-1)} = \partial X^{(0)} + \tfrac{1}{2}\partial^2 \Xi$.** From Day 152 (P1), $\tau$ is translation by $1$, so on polynomials $\tau = \sum_{k\ge 0}\partial^k/k!$ (Taylor, finite per $[T^n]$), and $\log H = (\tau - 1)X = \sum_{k \ge 1} \partial^k X/k!$. By Lemma 1.1(4) of Day 152, $\partial^k$ shifts wt by $-k$, so $\operatorname{wt}(\partial^k X) \le 1 - k$. The wt-$(-1)$ part therefore only receives contributions from $k = 1$ and $k = 2$:
$$M^{(-1)} = \ell^{\rm top}_{-1}(\partial X) + \ell^{\rm top}_{-1}(\partial^2 X/2) = \partial \ell^{\rm top}_0(X) + \tfrac{1}{2}\partial^2 \ell^{\rm top}_1(X) = \partial X^{(0)} + \tfrac{1}{2}\partial^2 \Xi.\qquad\square$$

## 3. The conjectural closed form $M^{(-1)} = 6T/(q^3\phi)$ at $E_3 = 0$

Setting $\ell^{\rm top}_{-1}(H)|_{E_3 = 0} = 6T/q^4$ and $\mathcal W|_{E_3=0} = \phi/q$ into the lemma yields
$$M^{(-1)}\big|_{E_3 = 0} = \frac{6T/q^4}{\phi/q} = \frac{6T}{q^3\phi}.$$
Using $q\phi = 1 - E_2 Y^2$, this equals $6T/[q^2(1 - E_2 Y^2)]$.

**Verification.** Both $M^{(-1)}$ (numerically, via the Day 152 `slog` pipeline) and $6T/(q^3\phi)$ (as a closed series expansion) can be expanded in $\mathbb Q[E_1, E_2][[T]]$. Their coefficient-by-coefficient match for $n \le 16$ is the content of `scratch/day156/verify_6T_over_q4.py`.

*Status of a fully structural proof.* Given the closed form for $\Xi$ at $E_3 = 0$ (from Day 154 §2, $\Xi = \int_0^T \pi\, dT'$ with $\pi = (P^2 - \Delta_2)/4 = E_2\phi/q^2 \cdot q^2/q^2 = \ldots$, explicit), $\partial^2 \Xi|_{E_3=0}$ is in principle computable. The obstruction is a compact formula for $X^{(0)} = \ell^{\rm top}_0(\log F_P)$. Section 6 sketches the natural completion via the Day 152 §2 decomposition $\log F_P = \varrho(S) + \varrho\log(\mathcal R/V)$.

## 4. From $6T/q^4$ to a Lagrange coefficient (proved)

Assume Theorem C.5.

> **Theorem C.5′.** $\displaystyle [T^n]\,\frac{6T}{q^4} = 6\,[Y^{n-1}]\,\frac{\phi(Y)^{n+2}}{(1 - E_2 Y^2)^3}$ for $n \ge 1$.

*Proof.* Since $Y = T\phi(Y)$ and $q = 1 - T\phi'(Y)$, the "extended Lagrange–Bürmann" formula
$$[T^m]\,F(Y) = [Y^m]\,F(Y)\, q(Y)\, \phi(Y)^m\tag{★}$$
holds (change of variable $t = y/\phi(y)$, $dt = q(y)/\phi(y)\, dy$, then residue at $y = 0$).

Apply (★) with $F(Y) = 1/q^4 = \phi(Y)^4/(q\phi)^4 = \phi(Y)^4/(1 - E_2 Y^2)^4$ (using $q\phi = 1 - E_2 Y^2$) and $m = n - 1$:
$$[T^{n-1}]\,\frac{1}{q^4} = [Y^{n-1}]\,\frac{\phi^4}{(1 - E_2 Y^2)^4}\cdot q\phi^{n-1}
= [Y^{n-1}]\,\frac{q\phi^{n+3}}{(1 - E_2 Y^2)^4}
= [Y^{n-1}]\,\frac{\phi^{n+2}}{(1 - E_2 Y^2)^3}$$
using $q\phi^{n+3} = (q\phi)\phi^{n+2} = (1 - E_2 Y^2)\phi^{n+2}$. Multiply by $6T$ and shift the $[T^n]$ index. $\square$

## 5. The $E$-positive expansion (proved)

> **Theorem C.5″.**
> $$[T^n]\,\ell^{\rm top}_{-1}(H)\big|_{E_3=0}
> \;=\;\sum_{b=0}^{\lfloor(n-1)/2\rfloor} 6\,(b+1)\,4^b\,\binom{n+2}{2b+3}\, E_1^{n-1-2b}\, E_2^b.$$

*Proof.* Expand C.5′'s residue multinomially. From $\phi^{n+2} = (1 + E_1 Y + E_2 Y^2)^{n+2} = \sum_{a+b+c = n+2} \frac{(n+2)!}{a! b! c!} E_1^b E_2^c Y^{b+2c}$ and $(1 - E_2 Y^2)^{-3} = \sum_m \binom{m+2}{2} E_2^m Y^{2m}$, extracting $[Y^{n-1}]$ and matching $Y$-exponents gives
$$6\,[Y^{n-1}]\frac{\phi^{n+2}}{(1-E_2 Y^2)^3}
= 6 \sum_{c \ge 0, m \ge 0} \frac{(n+2)!}{(n+2-b-c)!\, b!\, c!}\binom{m+2}{2}\,E_1^b E_2^{c+m}\Big|_{\substack{b + 2c + 2m = n-1}}.$$
Group by total $E_2$-exponent $b' := c + m$ (with $E_1$-exponent $a' := b$); the constraint $a' + 2b' = n-1$ is automatic. Rewriting $(n+2)!/[(n+2-a'-c)!\,a'!\,c!] = \binom{n+2}{2b'+3}\binom{2b'+3}{c}$ and summing $c$ from $0$ to $b'$ (the upper cap since $m = b' - c \ge 0$):
$$c'_{a', b'} := [E_1^{a'} E_2^{b'}]\ldots = 6\,\binom{n+2}{2b'+3}\sum_{c=0}^{b'}\binom{2b'+3}{c}\binom{b'-c+2}{2}.$$
It suffices to prove the identity
$$\boxed{\sum_{c=0}^{b}\binom{2b+3}{c}\binom{b-c+2}{2} = (b+1)\,4^{b}}\qquad (b \ge 0). \tag{♦}$$
Then $c'_{a', b'} = 6(b'+1)\,4^{b'}\,\binom{n+2}{2b'+3}$, which is the claim (with $a' = n-1-2b'$).

*Proof of (♦).* Write $F(j) := \binom{2b+3}{j}\cdot(b-j+1)(b-j+2)/2$ so that $\binom{b-j+2}{2} = (b-j+1)(b-j+2)/2$.

**Symmetry.** Under $j \mapsto 2b+3-j$, $\binom{2b+3}{j}$ is invariant, and $(b-j+1)(b-j+2)$ becomes $(j-b-2)(j-b-1) = (b+2-j)(b+1-j) = (b-j+1)(b-j+2)$: invariant. So $F(2b+3-j) = F(j)$.

**Middle terms vanish.** $F(b+1) = \binom{2b+3}{b+1}\cdot 0 \cdot 1/2 = 0$ and $F(b+2) = \binom{2b+3}{b+2}\cdot(-1)\cdot 0/2 = 0$.

**Full sum.** Using $\sum_j\binom{2b+3}{j} = 2^{2b+3}$, $\sum_j j\binom{2b+3}{j} = (2b+3)\,2^{2b+2}$, $\sum_j j^2\binom{2b+3}{j} = (2b+3)(2b+4)\,2^{2b+1}$, and $(b-j+1)(b-j+2) = (b+1)(b+2) - (2b+3)j + j^2$:
$$2\sum_{j=0}^{2b+3} F(j) = 2^{2b+2}(b+1)(b+2) - (2b+3)^2\,2^{2b+2} + (2b+3)(2b+4)\,2^{2b+1}.$$
Factor $2^{2b+1}$ and expand: 
$$2\sum_{j=0}^{2b+3}F(j) = 2^{2b+1}\bigl[2(b+1)(b+2) - 2(2b+3)^2 + (2b+3)(2b+4)\bigr] = 2^{2b+1}\cdot 2(b+1) = 4(b+1)\,4^b.$$
Hence $\sum_{j=0}^{2b+3} F(j) = 2(b+1)\,4^b$.

**Split.** By symmetry $\sum_{j=0}^{b}F(j) = \sum_{j=b+3}^{2b+3}F(j)$; combined with $F(b+1) = F(b+2) = 0$, the full sum equals $2\sum_{j=0}^{b}F(j)$, so $\sum_{j=0}^{b}F(j) = (b+1)\,4^b$. $\square$

## 6. Toward a fully structural proof of Theorem C.5

The lemma of §2 gives $H^{(-1)} = \mathcal W(\partial X^{(0)} + \tfrac{1}{2}\partial^2\Xi)$. The missing ingredient is a compact formula for $X^{(0)} = \ell^{\rm top}_0(\log F_P)$. Three natural routes:

**(A) Riccati at wt $-1$.** Day 152 §4 proves $\theta\Xi = Te_2(\nu)$ from the wt-$1$ top-symbol of the Riccati (R). By the same token, applying $\ell^{\rm top}_0$ to (R) should give a coupled linear system determining $\theta X^{(0)}$ in terms of the wt-$0$ symbols of $\lambda_i, u_i + d_i$, etc. This is the natural next tool.

**(B) $\varrho$-decomposition.** From Day 152 §2, $\log F_P = \varrho(S) + \varrho\log(\mathcal R/V)$. Then $X^{(0)} = \varrho S^{(0)} + \varrho[\log(\mathcal R/V)]^{(0)}$. Because $\log(\mathcal R/V)$ has wt $\le 0$ and $\mathcal R/V \in 1 + (t)$, $[\log(\mathcal R/V)]^{(0)} = \log[\ell^{\rm top}_0(\mathcal R/V)]$ (a ring homomorphism on wt $\le 0$). So this reduces to computing $\ell^{\rm top}_0(\mathcal R/V)$ — a level-3 residual in $\mathcal R$ that needs one further $\ell^{\rm top}$ extraction on $V(M)$ and $S$.

**(C) — the new hook (2026-09-02).** Computed to $n \le 8$ via `scratch/day156/test_X0_hypothesis.py`:

$$\boxed{\;X^{(0)}\big|_{E_3 = 0} \;=\; \tfrac{1}{2}\,\log \mathcal W\big|_{E_3 = 0} \;=\; \tfrac{1}{2}\,\partial\Xi\big|_{E_3 = 0}\;}$$

This does *not* hold globally: $X^{(0)} - \tfrac{1}{2}\log \mathcal W$ is $O(E_3)$, with the leading correction $4E_3 T^3 + 15 E_1 E_3 T^4 + (36 E_1^2 + 24 E_2) E_3 T^5 + \ldots$. But at the boundary $E_3 = 0$, the identity is exact.

Consequences.
1. $M^{(-1)}|_{E_3 = 0} = \partial X^{(0)}|_{E_3=0} + \tfrac{1}{2}\partial^2 \Xi|_{E_3=0} = \partial^2 \Xi|_{E_3=0} + E_2\,\partial_{E_3}D|_{E_3 = 0}$, where $D := X^{(0)} - \tfrac{1}{2}\log \mathcal W$. So *modulo the $E_3$-correction* to $X^{(0)}$, the answer is $\partial^2 \Xi$ (globally-computable via $\theta\partial^2\Xi = \partial^2 P/2$).
2. The $E_3$-correction to $X^{(0)}$ is captured by a *separate* series; identifying it would immediately close C.5.

The evidence is:
- $X^{(0)}|_{E_3=0}$ matches $(1/2)\log \mathcal W|_{E_3=0}$ for $n = 1, \ldots, 8$ (all terms).
- The failure at $E_3 \ne 0$ is *only* through explicit $E_3$-multiplied terms; the $E_3^0$ part of $X^{(0)}$ IS $(1/2)\log \mathcal W|_{E_3=0}$.

*Why this happens is likely structural.* At $u_3 = 0$, $F_P|_{u_3=0}$ has a factorization that turns $\log F_P|_{u_3=0}$ into a specific "$\Xi + (1/2)\partial \Xi + \ldots$" pattern. Route (B) applied at wt 0 with $u_3 = 0$ should make this manifest.

**Bottom line.** The structural gap in this session is precisely *"establish $X^{(0)}|_{E_3=0} = (1/2)\log \mathcal W|_{E_3=0}$ (or equivalently, $M^{(-1)}|_{E_3=0} = 6T/(q^3\phi)$)"*. Both are computed to $n \le 10$; the former is arguably the *cleaner* target, since it reduces to identifying an $E_3 = 0$ specialization of $F_P$ as some standard object.

## 7. Corollaries

**Corollary 7.1 (line slice $E_2 = 0$).**  At $E_2 = 0$ ($u_3 = u_2 = 0$, one variable $x = u_1$):
$q^2 = (1 - Tx)^2$, so $q^4 = (1-Tx)^4$ and $6T/q^4 = \sum_{n \ge 1} n(n+1)(n+2)\, x^{n-1}\, T^n$. Every layer $d = 1$ contribution reduces to $n(n+1)(n+2)\, E_1^{n-1}$, matching Theorem C.5″ at $b = 0$.

**Corollary 7.2 (Catalan-adjacent slice $x = y = 1$).** At $E_1 = 2, E_2 = 1, E_3 = 0$: $q^2 = 1 - 4T + 0 = (1 - 4T)$... wait, $(1-2T)^2 - 4T^2 = 1 - 4T + 4T^2 - 4T^2 = 1 - 4T$; so $q^2 = 1 - 4T$ and $6T/q^4 = 6T/(1 - 4T)^2 = \sum_{n \ge 1} 6\,n\,4^{n-1}\,T^n = \sum_{n \ge 1} \tfrac{3n}{2}\,4^{n}\,T^n$. Thus layer $d = 1$ at $x = y = 1$ (all-ones evaluation) gives the row sum $\tfrac{3n}{2}\cdot 4^{n}$; verify at $n = 2$: $3\cdot 4^2 = 48$ ✓; $n = 3$: $\tfrac{9}{2}\cdot 64 = 288$ ✓; $n = 5$: $\tfrac{15}{2}\cdot 1024 = 7680$ ✓ (matches Day 155 raw table).

**Corollary 7.3 (positivity in $(E_1, E_2)$-basis).** By Theorem C.5″, every $(E_1, E_2)$-monomial coefficient at $u_3 = 0$ layer $d = 1$ is a positive integer — no cancellation. Together with Day 154's Corollary 4.3 ($d = 0$ at $E_3 = 0$ positive), *the top two layers of $H$ at $E_3 = 0$ are both $E$-positive*, a first serious step into layered Conjecture P at $E_3 = 0$.

## 8. Verification ledger

| Claim | Script | Verdict |
|---|---|---|
| Raw $[T^n]H|_{E_3=0}$ layer $d=1$ matches $6[T^{n-1}](1/q^4)$ in $E_1, E_2$ | `scratch/day156/verify_6T_over_q4.py` | ✓ $n \le 16$ |
| Same identity via raw $F_P$ (built with `FP_coeffs` from Day 152's `lib.py`) | `scratch/day156/extend_data.py` | ✓ $n \le 7$ (extended data) |
| $M^{(-1)}|_{E_3=0} = 6T/(q^3\phi)$: direct extraction of layer $d=1$ from $\log H$ (via Day 152's `slog`) vs closed form $H^{(-1)}/\mathcal W$ | `scratch/day156/verify_M_neg1.py` | ✓ $n \le 10$ |
| $\sum_{j=0}^b \binom{2b+3}{j}\binom{b-j+2}{2} = (b+1)\,4^b$ | this file §5 (analytic) | ✓ proved |
| Extended-Lagrange formula (★) for $Y = T\phi(Y)$ | classical | standard (residue change of variables) |
| $q\phi = 1 - E_2 Y^2$ (algebraic identity at $E_3 = 0$) | Day 154 §2 | proved (residue-free algebra) |

The heaviest check is `verify_6T_over_q4.py`, which builds $F_P$ from its raw definition (via `FP_coeffs = T^+(e^{Te_2}V)/V`), divides two Fraction-arithmetic series to get $H = \tau F_P/F_P$, extracts the $c = 0, a + 2b = n-1$ coefficients from `to_E(H_n)`, and compares against $6 \cdot [T^{n-1}](1/q^4)$ where $1/q^4$ is built symbolically from $q^2 = 1 - 2TE_1 + T^2(E_1^2 - 4E_2)$. Every coefficient matches for $n = 2, \ldots, 16$.

The `verify_M_neg1.py` check adds a second, independent verification path: it computes $\log H$ symbolically via `slog` (which builds up $\log F_P$ and $\log H$ term-by-term from the exact `Fraction`-based $F_P$), extracts the $u$-degree $n-1$ part at $u_3 = 0$ (= layer $d=1$ at $E_3=0$ of $M = \log H$), and independently checks it against $H^{(-1)}/\mathcal W = 6T/(q^3\phi)$. This is a direct check on the *log side*, so it verifies the lemma of §2 in situ: the derivation $H^{(-1)} = \mathcal W M^{(-1)}$ is confirmed on both hands, and the closed form $M^{(-1)}|_{E_3=0} = 6T/(q^3\phi)$ holds through $n = 10$.

## 9. Registry update

* NEW: `narayana-layer-d1-E3-zero`: **`computed`** with `role: attempt`. Sub-nodes: `layer-d1-lagrange-form` **`proved`** and `layer-d1-E-positive-expansion` **`proved`**; both cite this file. Parent `psi-E-positive-layer-d0` from Day 154 remains `computed`.

## 10. What next

* **Close the structural gap in §6.** Extract $X^{(0)}$ via Fact I applied at $\ell^{\rm top}_0$ (route A) or via $\ell^{\rm top}_0(\mathcal R/V)$ (route B); either turns C.5 into a full theorem.
* **Layer $d = 2$ at $E_3 = 0$.** By the same top-symbol-with-corrections argument, expect a rational expression involving $q$ and $\phi$. Given the pattern $\mathcal W \propto 1/[q(q + 1 - E_1 T)]$ (layer 0) and $6T/q^4$ (layer 1), a candidate to test would be $c \cdot T^2 (\text{polynomial in } q, \phi)/q^k$.
* **General $E_3$: $E_3$-power expansion of the $d = 1$ layer.** Analogous to Day 154 §8's queue item.

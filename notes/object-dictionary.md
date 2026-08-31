# Object Dictionary — every symbol, its ring, its normalisation

> Written Day 150 dream cycle 2, because Rule 6 v2 has fired ~10 times and **every single firing
> was "two names, one object" or "one name, two objects."** This is §1 of the FPSAC paper, so it
> is not overhead. Keep it current; when you introduce a symbol, add a row the same day.

## 0. The headline: there are only TWO normalisation knobs

> **Day 151 correction.** This section said THREE knobs, $2^3=8$ frames. Wrong: knobs 2 and 3
> are not independent — flipping the letter *forces* the divisor. Two knobs, $2^2=4$ frames.
> Found by the Day-151 frame check, `scratch/2026-08-31-day151-smu-frame-check.md` §4.

Every frame confusion in this project (Day 141 "STUCK on $\Psi$ vs $P$ frame confusion",
Day 150's Day-131/Day-149 double-derivation, the $M_\mu$/$s^*_\mu$/$\mathfrak s_\mu$ triangle)
is a disagreement about one of exactly two binary choices:

| # | Knob | Option A | Option B | Related by |
|---|------|----------|----------|-----------|
| 1 | Factorial direction | **falling** $(x)_k$, map $\mathcal T$ | **rising** $x^{(k)}$, map $\mathcal T^+$ | $\varphi: u_i\mapsto-u_i$ |
| 2 | **Which letter** | **shifted** $x=u+\rho$, $x_i=u_i+(n{-}i)$, $\rho=(2,1,0)$ | **plain** $u$, $\rho$ carried in $\lambda=\mu+\rho$ | translation by $\rho$ |

$2^2=4$ frames. The project has walked between them for forty days without a map. **When two
formulas look "the same up to normalisation," name which knob differs before writing anything.**

**The divisor is DETERMINED, not chosen.** The old "knob 3" (ordinary $V(u)$ vs shifted
Vandermonde) is not a free choice: because falling/rising factorials are monic,
$$\det\big[(x_i{+}n{-}i)_{n-j}\big] \;=\; V(x+\rho)\qquad\text{identically (verified, }n=3).$$
So the shifted Vandermonde in the plain letter **is** the ordinary Vandermonde in the shifted
letter. Once you pick your letter, the divisor is the ordinary Vandermonde in *that* letter, and
you are done. Writing (2A,3A) or (2B,3B) is a typo, not a frame.

**Meta, worth recording.** This correction was found *by the dictionary*: the Day-151 check only
happened because the table forced a knob setting to be written down, and the written setting was
then testable — and false. Rule 13 caught an error in the very table built to enforce Rule 13,
in under 24 hours. The table pays for itself. Keep it current.

## 1. The maps

| Symbol | Definition | Domain → codomain | Note |
|---|---|---|---|
| $\mathcal T$ | $u^\alpha\mapsto\prod_i(u_i)_{\alpha_i}$ (falling) | $\mathbb Q[u]\to\mathbb Q[u]$ | knob 1A |
| $\mathcal T^+$ | $u^\alpha\mapsto\prod_iu_i^{(\alpha_i)}$ (rising) | " | knob 1B |
| $\varphi$ | $u_i\mapsto-u_i$ | " | conjugates $\mathcal T\leftrightarrow\mathcal T^+$ |
| $\Psi$ | $f\mapsto\mathcal T(fV)/V$ | $\Lambda_3\to\Lambda_3$ | the project's central map (Day 125 operator form) |
| $\Psi^+$ | $f\mapsto\mathcal T^+(fV)/V$ | " | $=\varphi\Psi\varphi$ |
| $\phi$ *(Day 123)* | $s_\mu\mapsto s^*_\mu$, then specialise | $\Lambda_3\to\mathbb Q[s,t]$ | **= $\Psi^+$ composed with a specialisation.** See §5. **Day 151 FLAG:** with $s^*_\mu$ in its corrected frame, $\Psi(s_\mu)=s^*_\mu$ (23/23) while $\Psi^+(s_\mu)=\mathfrak s_\mu$ — so this row reads as $\Psi$, not $\Psi^+$. Unresolved: go re-read Day 123 and settle it. Do not guess |
| $\tau$ | $u_i\mapsto u_i+1$; $E_1\mapsto E_1{+}3$, $E_2\mapsto E_2{+}2E_1{+}3$, $E_3\mapsto E_3{+}E_2{+}E_1{+}1$ | $\mathbb Q[E][[T]]\to$ itself, $\tau T=T$ | $\tau=e^{\partial}$, $\partial=\sum_i\partial_{u_i}$ |
| $\mathcal B$ | $\Psi^+\circ(e_2\cdot)\circ(\Psi^+)^{-1}$ | $\Lambda_3\to\Lambda_3$ | $=$ the $e_2$-Pieri operator; $P_b=\mathcal B^b(1)$ |
| $\mathcal B_k$ | $V^{-1}e_k(\hat u)V$, $\hat u_i=u_iS_i$ | " | $e_k$-Pieri; **$\mathcal B_3=E_3\tau$ exactly** |
| $\ell_d$ | $\operatorname{ord}$-$d$ layer, $\operatorname{ord}(E_1^iE_2^jE_3^kT^n)=n-3k$ | | ring hom on $\{\operatorname{ord}\ge0\}$ |
| $\ell^{\mathrm{top}}_w$ | $\operatorname{wt}$-$w$ layer, $\operatorname{wt}=(i{+}2j{+}3k)-n$ | | $u$-degree minus $T$-degree |

> **Theorem A (Day 149 §5).** $\Psi^+(s_\mu)=\mathfrak s_\mu$. *In words: $\Psi$ is the map that
> replaces every Schur function by its factorial Schur counterpart.* Everything since Day 110 is a
> statement about this map.

## 2. The three "shifted Schur" objects — pick one and stop

| Symbol | Formula | Knobs | Where used |
|---|---|---|---|
| $\mathfrak s_\mu$ | $\det[u_i^{(\lambda_j)}]/V(u)$, $\lambda=\mu+\rho$, $\rho=(2,1,0)$ | 1B, 2B (plain) | Day 149. **Use this one.** $=$ Macdonald $s_\mu(u\vert a)$ at $a_l=1-l$ |
| $s^*_\mu$ | $\det\big[[u_i]_{\mu_j+n-j}\big]/V(u)$, $[y]_k=y(y{-}1)\cdots(y{-}k{+}1)$ | 1A, 2B (plain) | Days 108–131 |
| $M_\mu$ | $\det[(x_i)_{k_j}]$, $x=(a{+}2,b{+}1,c)$, $k=\mu+\rho$ | 1A, 2A (shifted) — *undivided* | Day 108, `hk-three-var-fit.py` |

> **Day 151: the $s^*_\mu$ row was in the wrong frame.** It said "Okounkov–Olshanski,
> $\det[(x_i{+}n{-}i)_{\mu_j+n-j}]/(\text{shifted Vdm})$, knobs 1A, 2A, 3B". That object is the
> **translate**, not the object: read in the same letter it equals $s^*_\mu(u+\rho)$, confirmed
> 22/23 (all $\mu$ but $\emptyset$, where the two coincide). The thing Rick actually used on
> Days 108/117/118/123/124/127/129, and the thing $\Psi$ produces, is the **plain-letter**
> falling bialternant above, knobs (1A, 2B). Verbatim source: `proofs/2026-08-23-day129-shifted-schur-l4.md`,
> Main Theorem. Evidence: `scratch/2026-08-31-day151-smu-frame-check.md`.

**$\mathfrak s_\mu$ and $s^*_\mu$ now differ by knob 1 ALONE.** Both are in the plain letter, so
the only disagreement left is falling vs rising, and the conjugating map is named:
$$\boxed{\;\mathfrak s_\mu(u)\;=\;(-1)^{|\mu|}\,s^*_\mu(-u)\;=\;(-1)^{|\mu|}\varphi(s^*_\mu),\qquad \varphi:u_i\mapsto-u_i.\;}$$
Verified 23/23. **The sign is $(-1)^{|\mu|}$, not $(-1)^{|\lambda|}$** — the $V$ in the denominator
eats $(-1)^{|\rho|}=(-1)^3$. This is Day 149's $\Psi=\varphi\Psi^+\varphi$ written out on basis
elements, with its sign.

**Confirmed dictionary rows (Day 151, all 23/23 symbolic over $\mathbb Q$, $n=3$, $|\mu|\le6$):**

| Claim | Status |
|---|---|
| $s^*_\mu=s_\mu(u\vert a)$, the classical **factorial Schur** at shift $a_l=l-1$, $l\ge1$ (so $a=(0,1,2,\dots)$) | EQUAL 23/23 |
| $s_\mu(u\vert a)$ at the *other* shift $a_l=1-l$ $=\mathfrak s_\mu$ | $(-1)^{|\mu|}\varphi$, 23/23 — same knob-1 flip |
| $M_\mu=V(x)\cdot s^*_\mu(x)$ at $x=u+\rho$ — i.e. $M_\mu$ is the **undivided numerator of $s^*_\mu$ in the shifted letter** | EQUAL 23/23 |
| $\Psi(s_\mu)=\mathcal T(s_\mu V)/V=s^*_\mu$ — re-derived from the *map* definition, independently | EQUAL 23/23 |

Convention for the factorial Schur, **state it every time**: Macdonald, *Symmetric Functions and
Hall Polynomials*, 2nd ed., I.3 Ex. 20, with
$$(u\mid a)^k \;=\; (u-a_1)(u-a_2)\cdots(u-a_k),$$
**1-indexed**, and $s_\mu(u\mid a)=\det\big[(u_i\mid a)^{\mu_j+n-j}\big]/V(u)$. At $a_l=l-1$ this
is $(u)(u{-}1)\cdots(u{-}k{+}1)=[u]_k$, which is why $s^*_\mu$ comes out on the nose.

All three objects are the same function in different coordinates. $\mathfrak s_\mu=s_\mu+(\text{lower degree})$;
all are a $\mathbb Z$-basis of $\Lambda_3=\mathbb Z[E_1,E_2,E_3]$.

## 3. The series

| Symbol | Definition | Ring | Anchor value |
|---|---|---|---|
| $E_i$ | $e_i(u_1,u_2,u_3)$ | $\mathbb Z[u]^{S_3}$ | — |
| $V$ | $\prod_{i<j}(u_i-u_j)$ | | $u$-homog. degree 3 |
| $P_b$ | $\varphi(\Psi(e_2^b))=\Psi^+(e_2^b)$ | $\mathbb Z[E]$ | $P_3=\mathfrak s_{222}+2\mathfrak s_{321}+\mathfrak s_{330}$ |
| $E_j$ *(Day 123)* | $\sum_\mu K_{\mu',(2^j)}s^*_\mu$ | $\mathbb Q[e_1,e_2,e_3]$ | **$=P_j$ in frame (1A, 2B).** See §5. *Day 151: label was (1A,2A,3B); the $s^*_\mu$ frame was recorded wrong — see §2* |
| $F_P$ | $\sum_b P_bT^b/b!=\mathcal T^+(e^{Te_2}V)/V$ | $\mathbb Q[E][[T]]$ | $1+O(T)$; $\deg_u[T^n]F_P=2n$ (sharp) |
| $H$ | $\tau(F_P)/F_P$ | " | $\deg_u[T^n]H\le n$ (Thm 2, sharp) |
| $\mathcal H$ | $\ell_0(H)=\dfrac{F(F-1)}{\vartheta(2F-3)}=\dfrac{2G-1}{(3G-1)^2(4G-1)}$ | $\mathbb Z[[\vartheta]]$ | integral, **unconditional** (Day 149 Cor 3) |
| $\Xi$ | $\ell^{\mathrm{top}}_1(\log F_P)$ | | $\theta\Xi=\tfrac{1-q}{2T}-\tfrac{E_1}2$ |
| $\mathcal W$ | $\ell^{\mathrm{top}}_0(H)=\exp(\partial\Xi)$ | | leading symbol of $H$; $=\sum\mathcal W_nT^n$ |
| $\mathcal R$ | $e^{-S}V(M)e^S$, $S=\log\mathcal M$ | | prefactor, $\operatorname{wt}\le3$ |
| $\psi$ | Lagrange kernel: $Y=T\psi(Y)$, $\sum_n(n{+}1)\mathcal W_nT^n=dY/dT$ | $\mathbb Z[E][[Y]]$? | $1{+}E_1Y{+}E_2Y^2{+}2E_3Y^3{+}E_1E_3Y^4{+}2E_2E_3Y^5{+}(E_1E_2E_3{+}5E_3^2)Y^6$ |
| $q$ | $1-Te_1(\nu)$; branch $q=1+O(T)$ of $\sum_i\sqrt{q^2{+}4Tu_i}=q+2$ | | **the master curve** (Day 149 Thm 4) |

## 4. The $b_k$ arc (Arc A) — all of it is $q$ specialised

Specialisation $\sigma$: $E_1=E_2=0$, $T=1$, $u_i=\omega^{i-1}z$, so $E_3=z^3=:\vartheta$.

| Symbol | Definition | Relation |
|---|---|---|
| $\vartheta$ | $=E_3$ under $\sigma$ | the Arc-A variable |
| $F$ | $\sum_{k\ge1}b_k\vartheta^k$ | $\boxed{q=1-2F}$ — **$F$ is the master curve under $\sigma$** |
| $M$ | $1-2F$ | $=q\vert_\sigma$. The free-cumulant object; $\kappa_n(M)\in6\mathbb Z$ |
| $G$ | $F/3$ | $G=\vartheta\,\Phi_{\mathrm A}(G)$ |
| $\Phi_{\mathrm A}$ | $\dfrac{(1-2G)^2}{(1-3G)^3(1-4G)}$ | Arc-A **Lagrange kernel**. $=1+9G+\cdots$ |
| $b_k$ | $3,27,417,7851,164124,3661389,85384566,2056373739$ | $b_k/3=\frac1k[G^{k-1}]\Phi_{\mathrm A}^k\in\mathbb Z$ |
| $a_k$ | $(1-2F)^2=1+4A$, $A=\sum a_k\vartheta^k$ | $a_k=-b_k+\sum b_ib_j$; $a_7=-48005802$, $a_8=-1147833720$ |

**Naming fix, adopt immediately.** Day 148 called the Arc-A kernel $\phi$; that collides with
Day 123's map $\phi$ *and* with $\varphi:u_i\mapsto-u_i$ — **three distinct objects spelled phi.**
Renamed $\Phi_{\mathrm A}$ here. Reserve $\psi$ for the Arc-B kernel, $\Psi$ for the map.

## 5. Resolved: Day 131's $E_j$ IS Day 149's $P_b$

See [[q-day131-vs-day149-same-statement]]. Verdict **SAME STATEMENT**, verified by hand on the
coefficients. Day 149 Corollary B is not new; Theorems A, C, E are.

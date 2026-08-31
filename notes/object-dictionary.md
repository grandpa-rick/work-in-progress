# Object Dictionary — every symbol, its ring, its normalisation

> Written Day 150 dream cycle 2, because Rule 6 v2 has fired ~10 times and **every single firing
> was "two names, one object" or "one name, two objects."** This is §1 of the FPSAC paper, so it
> is not overhead. Keep it current; when you introduce a symbol, add a row the same day.

## 0. The headline: there are only THREE normalisation knobs

Every frame confusion in this project (Day 141 "STUCK on $\Psi$ vs $P$ frame confusion",
Day 150's Day-131/Day-149 double-derivation, the $M_\mu$/$s^*_\mu$/$\mathfrak s_\mu$ triangle)
is a disagreement about one of exactly three binary choices:

| # | Knob | Option A | Option B | Related by |
|---|------|----------|----------|-----------|
| 1 | Factorial direction | **falling** $(x)_k$, map $\mathcal T$ | **rising** $x^{(k)}$, map $\mathcal T^+$ | $\varphi: u_i\mapsto-u_i$ |
| 2 | Variable shift | **shifted** $x_i=u_i+(n{-}i)$, i.e. $\rho=(2,1,0)$ absorbed | **plain** $u_i$, $\rho$ carried in $\lambda=\mu+\rho$ | translation by $\rho$ |
| 3 | Determinant divisor | ordinary Vandermonde $V(u)$ | shifted Vandermonde $\prod_{i<j}(x_i{-}x_j{-}(j{-}i))$ | $V(u+\rho)$ vs $V(u)$ |

$2^3=8$ frames. The project has walked between them for forty days without a map. **When two
formulas look "the same up to normalisation," name which knob differs before writing anything.**

## 1. The maps

| Symbol | Definition | Domain → codomain | Note |
|---|---|---|---|
| $\mathcal T$ | $u^\alpha\mapsto\prod_i(u_i)_{\alpha_i}$ (falling) | $\mathbb Q[u]\to\mathbb Q[u]$ | knob 1A |
| $\mathcal T^+$ | $u^\alpha\mapsto\prod_iu_i^{(\alpha_i)}$ (rising) | " | knob 1B |
| $\varphi$ | $u_i\mapsto-u_i$ | " | conjugates $\mathcal T\leftrightarrow\mathcal T^+$ |
| $\Psi$ | $f\mapsto\mathcal T(fV)/V$ | $\Lambda_3\to\Lambda_3$ | the project's central map (Day 125 operator form) |
| $\Psi^+$ | $f\mapsto\mathcal T^+(fV)/V$ | " | $=\varphi\Psi\varphi$ |
| $\phi$ *(Day 123)* | $s_\mu\mapsto s^*_\mu$, then specialise | $\Lambda_3\to\mathbb Q[s,t]$ | **= $\Psi^+$ composed with a specialisation.** See §5 |
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
| $\mathfrak s_\mu$ | $\det[u_i^{(\lambda_j)}]/V(u)$, $\lambda=\mu+\rho$, $\rho=(2,1,0)$ | 1B, 2B, 3A | Day 149. **Use this one.** $=$ Macdonald $s_\mu(u\vert a)$ at $a_l=1-l$ |
| $s^*_\mu$ | Okounkov–Olshanski, $\det[(x_i{+}n{-}i)_{\mu_j+n-j}]/(\text{shifted Vdm})$ | 1A, 2A, 3B | Days 108–131 |
| $M_\mu$ | $\det[(x_i)_{k_j}]$, $x=(a{+}2,b{+}1,c)$, $k=\mu+\rho$ | 1A, 2A, 3 — *undivided* | Day 108, `hk-three-var-fit.py` |

All three are the same function in different coordinates. $\mathfrak s_\mu=s_\mu+(\text{lower degree})$;
all are a $\mathbb Z$-basis of $\Lambda_3=\mathbb Z[E_1,E_2,E_3]$.

## 3. The series

| Symbol | Definition | Ring | Anchor value |
|---|---|---|---|
| $E_i$ | $e_i(u_1,u_2,u_3)$ | $\mathbb Z[u]^{S_3}$ | — |
| $V$ | $\prod_{i<j}(u_i-u_j)$ | | $u$-homog. degree 3 |
| $P_b$ | $\varphi(\Psi(e_2^b))=\Psi^+(e_2^b)$ | $\mathbb Z[E]$ | $P_3=\mathfrak s_{222}+2\mathfrak s_{321}+\mathfrak s_{330}$ |
| $E_j$ *(Day 123)* | $\sum_\mu K_{\mu',(2^j)}s^*_\mu$ | $\mathbb Q[e_1,e_2,e_3]$ | **$=P_j$ in frame (1A,2A,3B).** See §5 |
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

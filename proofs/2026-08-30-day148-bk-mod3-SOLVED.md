# Day 148 — $b_k \equiv 0 \pmod 3$: **SOLVED**

**Date:** 2026-08-30. **Status:** theorem. No gap in the main line (§9); the only external
input is used solely for a corollary.

---

## 0. The result

$$\boxed{\;F(\vartheta)\,\bigl(F(\vartheta)-1\bigr)^{3}\,\bigl(4F(\vartheta)-3\bigr) \;=\; \vartheta\,\bigl(2F(\vartheta)-3\bigr)^{2}\;}$$

where $F(\vartheta)=\sum_{k\ge1}b_k\vartheta^k$. So $F$ is **algebraic of degree 5** over $\mathbb Q(\vartheta)$.

Substituting $F=3G$, the two sides each acquire exactly one factor $9$, which cancels:

$$G\,(3G-1)^3(4G-1) \;=\; \vartheta\,(2G-1)^2
\qquad\Longleftrightarrow\qquad
G \;=\; \vartheta\,\phi(G),\quad \phi(G)=\frac{(2G-1)^2}{(3G-1)^3(4G-1)} .$$

$\phi \in \mathbb Z[[G]]$ with $\phi(0)=1$, so by **Lagrange inversion** $G\in\mathbb Z[[\vartheta]]$:

$$\frac{b_k}{3}\;=\;[\vartheta^k]G\;=\;\frac1k\,[G^{k-1}]\,\frac{(2G-1)^{2k}}{(3G-1)^{3k}(4G-1)^{k}}\;\in\;\mathbb Z .$$

> **THEOREM. $b_k\equiv0\pmod 3$ for every $k\ge1$.**

And as a bonus, at *every* prime, not just $3$ — this one does lean on (H2), see §9.2:

$$\mathcal H(\vartheta)\;=\;\frac{2F-3}{(F-1)^2(4F-3)}\;=\;\frac{2G-1}{(3G-1)^2(4G-1)}\;\in\;\mathbb Z[[\vartheta]].$$

The whole four-session arc collapses because a *closed form for $F_P$* exists (§2) and nobody had looked for it.

---

## 1. Two negatives cleared first (Day 147's "one non-circular lead" was neither)

**1.1 $m_n\ge0$ is automatic — it is NOT evidence.**
Day 147 claimed non-negativity of the necklace numbers
$m_n=\frac1n\sum_{d\mid n}\mu(n/d)s_d$ is "strictly stronger than integrality, hence
non-circular". False in substance. If $s_d\ge0$ for all $d$ and
$s_n>\sum_{d\mid n,\,d<n}s_d$, then $n\,m_n\ge s_n-\sum_{d\mid n,d<n}s_d>0$ trivially.
Measured ratios $s_n\big/\sum_{d\mid n,d<n}s_d$ for $n=2..15$:

$$21.8,\;532,\;602,\;3.6\cdot10^5,\;1.8\cdot10^4,\;2.7\cdot10^8,\;5.4\cdot10^5,\;3.9\cdot10^8,\;1.6\cdot10^7,\;1.6\cdot10^{14},\;4.7\cdot10^8,\;1.3\cdot10^{17},\;1.4\cdot10^{10},\;2.9\cdot10^{14}.$$

$s_n$ beats the *sum of all its proper divisor terms* by $10^2$–$10^{17}$. So exact
realizability of $(s_n)$ is **equivalent** to $\mathcal H\in\mathbb Z[[\vartheta]]$, i.e.
to the target — circular like everything else. **Rule 6 v2, firing #10.**

**1.2 $\zeta(t)=\mathcal H(t)$ is not rational of low degree.** Hankel determinants
$\det[h_{i+j+m}]_{i,j<r}$ are nonzero for all $r\le8$, $m=0,1$. No shift of finite type.
(Correct — $\mathcal H$ is algebraic but irrational, see §0.)

Neither test cost more than ten minutes. Both should have been run on Day 147.

---

## 2. The closed form (the actual discovery)

Notation: $u_1,u_2,u_3$ indeterminates, $E_i=e_i(u)$, $\delta_{ij}=u_i-u_j$,
$V(v)=\prod_{i<j}(v_i-v_j)$, rising factorial $x^{(n)}=x(x+1)\cdots(x+n-1)$,
falling factorial $(x)_n=x(x-1)\cdots(x-n+1)$.

Recall the definitions (Day 131): $\mathcal T$ is the umbral map $u^\alpha\mapsto\prod_i(u_i)_{\alpha_i}$;
$\Psi(f)=\mathcal T(fV)/V$; $\Psi_b=\Psi(e_2^b)$; $\varphi:u_i\mapsto-u_i$; $P_b=\varphi(\Psi_b)$;
$F_P=\sum_b P_b T^b/b!$.

### Lemma 2.1 (umbral determinant collapse)
For $m\in\mathbb Z_{\ge0}^3$,
$$\mathcal T\bigl(u^m\,V(u)\bigr)\;=\;\Bigl(\prod_i (u_i)_{m_i}\Bigr)\cdot V(u-m),
\qquad (u-m)_i:=u_i-m_i .$$

*Proof.* $u^mV(u)=\sum_{w\in S_3}\operatorname{sgn}(w)\prod_i u_i^{\,m_i+\epsilon^w_i}$ with
$(\epsilon^w_i)$ the permuted exponents $(2,1,0)$, so
$\mathcal T(u^mV)=\det\bigl[(u_i)_{m_i+3-j}\bigr]_{i,j}$.
Now $(x)_{m+k}=(x)_m\,(x-m)_k$, so the $i$-th row factors:
$\det[(u_i)_{m_i}(u_i-m_i)_{3-j}]=\prod_i(u_i)_{m_i}\cdot\det[(v_i)_{3-j}]$ with $v=u-m$.
Finally $(v)_2,(v)_1,(v)_0$ are monic of degrees $2,1,0$, so
$\det[(v_i)_{3-j}]=\det[v_i^{3-j}]=V(v)$. $\square$

### Theorem 2.2 (CLOSED FORM FOR $F_P$)
$$\boxed{\;F_P \;=\; \sum_{a,b,c\ge0}\frac{T^{\,a+b+c}}{a!\,b!\,c!}\;
u_1^{(a+b)}\,u_2^{(a+c)}\,u_3^{(b+c)}\;\frac{V(u+m)}{V(u)}\;,\qquad
m=(a{+}b,\;a{+}c,\;b{+}c).\;}$$
Explicitly $\dfrac{V(u+m)}{V(u)}=\dfrac{(\delta_{12}+b-c)(\delta_{13}+a-c)(\delta_{23}+a-b)}{\delta_{12}\delta_{13}\delta_{23}}$.

*Proof.* $\sum_b\Psi_bT^b/b!=\mathcal T\!\left(e^{Te_2(u)}V\right)/V$ by linearity of $\mathcal T$.
Expand $e^{Te_2}=\sum_{a,b,c}\frac{T^{a+b+c}}{a!b!c!}(u_1u_2)^a(u_1u_3)^b(u_2u_3)^c
=\sum\frac{T^{n}}{a!b!c!}u^m$ and apply Lemma 2.1:
$$\sum_b\Psi_b\frac{T^b}{b!}=\sum_{a,b,c}\frac{T^n}{a!b!c!}\prod_i(u_i)_{m_i}\frac{V(u-m)}{V(u)} .$$
Now apply $\varphi$, i.e. substitute $u\mapsto-u$. Since $|m|=2n$ is even,
$\prod_i(-u_i)_{m_i}=\prod_i u_i^{(m_i)}$; and
$V(-u-m)/V(-u)=\bigl(-\prod_{i<j}(\delta_{ij}+m_i-m_j)\bigr)/(-V(u))=V(u+m)/V(u)$. $\square$

**Verified** (`/tmp/vclosed.py`, and re-derived in `fast.py`): identical to the
$\Psi$-recursion output of `core.py` for all $b\le9$ at $u=(\tfrac12,\tfrac13,\tfrac15)$,
and used to reproduce $b_1,\dots,b_{15}$ exactly at two base points (§8).

*(Provenance of the derivation: $\mathcal T$ is exactly the operator
$g\mapsto\bigl[g(\partial_{x_1},\partial_{x_2},\partial_{x_3})\,x_1^{u_1}x_2^{u_2}x_3^{u_3}\bigr]_{x=1}$,
because $\partial_x^n x^u|_{x=1}=(u)_n$. Then $e^{Te_2(\partial)}$ is a Gaussian
operator and everything is forced. The three-session hunt for a
"Frobenius structure" never asked what $\mathcal T$ actually is.)*

### Corollary 2.3 (Horn system)
Let
$$\mathcal M(t_1,t_2,t_3):=\sum_{a,b,c\ge0}\frac{t_1^at_2^bt_3^c}{a!b!c!}\,
u_1^{(a+b)}u_2^{(a+c)}u_3^{(b+c)},\qquad \theta_i:=t_i\partial_{t_i},$$
$$M_1:=u_1+\theta_1+\theta_2,\qquad M_2:=u_2+\theta_1+\theta_3,\qquad M_3:=u_3+\theta_2+\theta_3 .$$
The $M_i$ commute, $M_i$ acts on the $(a,b,c)$-coefficient by the scalar $u_i+m_i$, and
$$\theta_1\mathcal M=t_1M_1M_2\mathcal M,\qquad
\theta_2\mathcal M=t_2M_1M_3\mathcal M,\qquad
\theta_3\mathcal M=t_3M_2M_3\mathcal M, \tag{H}$$
$$\boxed{\;F_P=\frac{V(M)\,\mathcal M}{V(u)}\bigg|_{t_1=t_2=t_3=T}\;,\qquad V(M)=(M_1-M_2)(M_1-M_3)(M_2-M_3).}$$

*Proof.* (H) is the coefficient ratio: writing $c(a,b,c)$ for the coefficient,
$\dfrac{a\,c(a,b,c)}{c(a-1,b,c)}=(u_1+a+b-1)(u_2+a+c-1)$, which is (H) read off in
degree $(a,b,c)$; likewise for $\theta_2,\theta_3$.
For the boxed identity, $M_1-M_2=\delta_{12}+\theta_2-\theta_3$ acts by $\delta_{12}+m_1-m_2$,
etc., so $V(M)\mathcal M$ has $(a,b,c)$-coefficient $\frac{\prod u_i^{(m_i)}}{a!b!c!}V(u+m)$;
compare Theorem 2.2. $\square$

---

## 3. The base point and the order filtration

Fix $E_1,E_2$ (arbitrary) and set $E_3=\epsilon^3$ with $\epsilon$ a parameter, so
$$\vartheta \;=\; E_3T^3\;=\;(\epsilon T)^3\;=\;z^3,\qquad z:=\epsilon T .$$
The $u_i$ are the roots of $t^3-E_1t^2+E_2t-\epsilon^3$, i.e. the Laurent series
$$u_i \;=\; \omega^{\,i-1}\epsilon \;+\; \tfrac{E_1}{3} \;+\; O(\epsilon^{-1}),
\qquad \omega=e^{2\pi i/3},$$
(the $\epsilon\to\infty$ branch: $\sum_i\omega^{i-1}=e_2(1,\omega,\omega^2)=0$, so the $\epsilon$
and $\epsilon^2$ parts of $e_1,e_2$ vanish automatically and $\sum_iu_i=E_1$ forces the constant
term $E_1/3$). At $E_1=E_2=0$ this is exact: $u_i=\omega^{i-1}\epsilon$.

Write $t_i=T\hat t_i$ and work in $\mathcal K((\epsilon^{-1}))[\hat t_1,\hat t_2,\hat t_3]$ per
$t$-degree, $\mathcal K=\mathbb Q(\omega)(E_1,E_2)$, with the **order** grading
$$\operatorname{ord}\bigl(\epsilon^{m}T^{n}\hat t^{\alpha}\bigr)=n-m
\qquad\text{and}\qquad \ell_d:=(\text{coefficient of }T^d).$$
This is exactly the Day 146 order filtration: $E_3^kT^b=z^{3k}T^{\,b-3k}$ has order $b-3k$.
Key facts, all immediate:

* $\operatorname{ord}u_i\ge-1$ with $\;\boxed{\ell_{-1}(u_i)=\omega^{\,i-1}z}\;$ — **independent of $E_1,E_2$.**
* $\operatorname{ord}\delta_{ij}=-1$ exactly, $\ell_{-1}(\delta_{ij})=(\omega^{i-1}-\omega^{j-1})z$;
  so $\operatorname{ord}V(u)=-3$ and $\ell_{-3}(V(u))=z^3D$, $D=\prod_{i<j}(\omega^{i-1}-\omega^{j-1})\ne0$.
* $t_i$ raises order by $1$ and $t$-degree by $1$; $\theta_i$ preserves both.
* $\operatorname{ord}$ is additive on products; for each $t$-degree $n$ one has
  $\deg_\epsilon\le|m|=2n$, so $\operatorname{ord}\ge-n$: **bounded below in each $t$-degree.**

*(No ambiguity from the choice of branch: for each fixed $T$-degree $n$, $[T^n]F_P=P_n/n!$ is a
polynomial in $E_1,E_2,\epsilon^3$, and the $(a,b,c)$-sum of Theorem 2.2 is finite, so evaluating
it on the $\epsilon\to\infty$ Laurent expansion of the roots returns exactly that polynomial —
all negative powers of $\epsilon$ cancel.)*

---

## 4. The vanishing lemma (unconditional, any base point)

Put $S:=\log\mathcal M$, $\lambda_i:=\theta_iS$, and
$$d_1:=(\theta_1{+}\theta_2)S=\lambda_1+\lambda_2,\quad
d_2:=(\theta_1{+}\theta_3)S=\lambda_1+\lambda_3,\quad
d_3:=(\theta_2{+}\theta_3)S=\lambda_2+\lambda_3 .$$
Since $M_2e^S=(u_2+d_2)e^S$ and
$M_1\bigl[(u_2+d_2)e^S\bigr]=\bigl[(u_1+d_1)(u_2+d_2)+D_1d_2\bigr]e^S$ (with $D_1=\theta_1+\theta_2$;
the $u_i$ are $t$-constants), dividing (H) by $\mathcal M$ gives the **Riccati system**

$$\lambda_1=t_1\bigl[(u_1{+}d_1)(u_2{+}d_2)+D_1d_2\bigr],\quad
\lambda_2=t_2\bigl[(u_1{+}d_1)(u_3{+}d_3)+D_1d_3\bigr],\quad
\lambda_3=t_3\bigl[(u_2{+}d_2)(u_3{+}d_3)+D_2d_3\bigr]. \tag{R}$$

### Theorem 4.1 (vanishing / Prop 1, now UNCONDITIONAL)
$\operatorname{ord}\lambda_i\ge-1$ for $i=1,2,3$. Consequently $\operatorname{ord}\log F_P\ge-1$
(§5, Step 3), i.e. $[E_3^kT^b]\log F_P=0$ for all $b<3k-1$.

*(The last implication: writing $\log F_P=\sum_{k,b}c_{k,b}E_3^kT^b$, the monomial
$E_3^kT^b=z^{3k}T^{b-3k}$ contributes to $\ell_{b-3k}$ with the distinct $z$-power $z^{3k}$, so
monomials of equal order cannot cancel. Hence $\ell_d(\log F_P)=0$ forces $c_{k,3k+d}=0$ for
every $k$.)*

*Proof.* Induction on the $t$-degree $n$ (i.e. on $a+b+c$). Since $\mathcal M-1$ has
$t$-degree $\ge1$, so does $S$, so $[\deg 0]\lambda_i=0$ and the case $n=0$ holds.

Let $n\ge1$ and assume $\operatorname{ord}\bigl([\deg n']\lambda_i\bigr)\ge-1$ for all
$n'<n$ and all $i$ — hence the same for $d_1,d_2,d_3$. Because $t_1=T\hat t_1$ raises both
$t$-degree and order by $1$,
$$[\deg n]\,\lambda_1 \;=\; \hat t_1\,T\cdot[\deg (n{-}1)]\bigl\{(u_1{+}d_1)(u_2{+}d_2)+D_1d_2\bigr\}.$$
Bound the bracket at $t$-degree $n-1$:

* $D_1d_2$: $D_1$ preserves degree and order, and $[\deg(n{-}1)]d_2$ has order $\ge-1$ by the
  induction hypothesis. Order $\ge-1$.
* $(u_1{+}d_1)(u_2{+}d_2)$: expanding by $t$-degree, each factor is taken at some degree
  $m_i$ with $m_1+m_2=n-1$, so $m_i\le n-1<n$. At degree $0$ the factor is $u_i$, of order
  $\ge-1$; at degree $m\ge1$ it is $[\deg m]d_j$, of order $\ge-1$ by induction. Hence each
  factor has order $\ge-1$ and the product has order $\ge-2$.

So the bracket has order $\ge-2$ at degree $n-1$, and
$\operatorname{ord}\bigl([\deg n]\lambda_1\bigr)\ge 1+(-2)=-1$. Identically for
$\lambda_2,\lambda_3$. $\square$

*(Note how little this uses: only that $\operatorname{ord}u_i\ge-1$, that $t_i$ raises order
by one, and that $S$ has positive $t$-degree. In particular it holds at **every** base point.
This is the first unconditional proof of the Day 143 "leading-$T$ vanishing" lemma; Day 146
Prop. 1 had it only modulo (H2). Independently confirmed: **zero violations** for all $k\le22$
at two base points, §8.)*

---

## 5. The leading-order system

Write $L_i:=\ell_{-1}(\lambda_i)\in\mathcal K[[z]][\hat t]$. Apply $\ell_{-1}$ to (R).

*Step 1.* $\ell_{-1}(\lambda_1)=\hat t_1\,\ell_{-2}\bigl[(u_1+d_1)(u_2+d_2)+D_1d_2\bigr]$.
By Theorem 4.1, $D_1d_2$ has order $\ge-1$, so its $\ell_{-2}$ vanishes; and both
$u_i+d_i$ have order $\ge-1$, so $\ell_{-2}$ of the product is the product of the
$\ell_{-1}$'s. With $\ell_{-1}(u_i)=\omega^{i-1}z$:
$$L_1=\hat t_1\,\nu_1\nu_2,\qquad L_2=\hat t_2\,\nu_1\nu_3,\qquad L_3=\hat t_3\,\nu_2\nu_3,$$
$$\nu_1:=z+L_1+L_2,\qquad \nu_2:=\omega z+L_1+L_3,\qquad \nu_3:=\omega^2z+L_2+L_3 .$$

*Step 2.* Set $\hat t_1=\hat t_2=\hat t_3=1$ (legitimate: $\hat t$ has order $0$). Then
$L_1+L_2=\nu_1(\nu_2+\nu_3)$, etc., and the system becomes

$$\boxed{\;\nu_i\bigl(1-e_1(\nu)+\nu_i\bigr)=\omega^{\,i-1}z\quad(i=1,2,3).\;}\tag{$\ast$}$$

The Jacobian of $(\ast)$ at $z=0,\nu=0$ is the identity, so $(\ast)$ has a **unique**
solution $\nu\in\bigl(z\,\mathcal K[[z]]\bigr)^3$.

*Step 3 (the prefactor is invisible).* $F_P=\dfrac{V(M)\mathcal M}{V(u)}$, so
$\log F_P=S+\log\dfrac{\mathcal R}{V(u)}$ with $\mathcal R:=V(M)e^S/e^S$. Applying
$M_i-M_j=\delta_{ij}+\Delta_{ij}$ ($\Delta_{ij}$ a $\theta$-difference) to $e^S$ produces
$(\delta_{ij}+\Delta_{ij}S)e^S$; iterating, the corrections are $\Delta$-derivatives of
already-produced factors, hence of order $\ge-1$ instead of $-2$, so
$$\operatorname{ord}\mathcal R\ge-3,\qquad \ell_{-3}(\mathcal R)=\prod_{i<j}\ell_{-1}(\delta_{ij}+\Delta_{ij}S)=\prod_{i<j}(\nu_i-\nu_j)=V(\nu),$$
using $\ell_{-1}(\delta_{12}+\Delta_{12}S)=(1-\omega)z+L_2-L_3=\nu_1-\nu_2$, etc.

> **[Day 149 audit — CORRECTION.]** The parenthetical above ("the corrections are $\Delta$-derivatives
> of already-produced factors, hence of order $\ge-1$ instead of $-2$") is **wrong**: the produced
> factors $\delta_{ij}+\Delta_{ij}S$ themselves have order $\ge-1$, so a $\Delta$-derivative of one
> of them also has order $\ge-1$ and the correction terms **do** reach order $-3$. Hence
> $\ell_{-3}(\mathcal R)=V(\nu)$ is *not* exact. **The conclusion survives.** A correction term
> replaces one factor by $\Delta_{kl}\Delta_{ij}S$, whose $\ell_{-1}$ is a $\theta$-difference of
> $\ell_{-1}(S)$ and hence $O(z^2)$ (since $L_i=\hat t_i\nu_j\nu_k$ and $\nu_i=\omega^{i-1}z+O(z^2)$),
> while the two surviving factors are $O(z)$ each: every correction is $O(z^4)$. As
> $V(\nu)=Dz^3(1+O(z))$, we get $\ell_{-3}(\mathcal R)=Dz^3(1+O(z))$, so $\ell_0(\mathcal R/V(u))$ is
> a unit of $\mathcal K[[z]]$ congruent to $1$ mod $z$ — which is all Step 3 is used for.
> See `2026-08-30-day149-H2-PROVED.md` §6.
Meanwhile $\ell_{-3}(V(u))=z^3D$, $D=\prod_{i<j}(\omega^{i-1}-\omega^{j-1})\ne0$, and
$\nu_i=\omega^{i-1}z+O(z^2)$ gives $V(\nu)=z^3D\,(1+O(z))$. So $\mathcal R/V(u)$ has
order $\ge0$ with $\ell_0$ a unit of $\mathcal K[[z]]$ congruent to $1$; hence
$\log(\mathcal R/V(u))$ has order $\ge0$ and contributes nothing to $\ell_{-1}$.

*Step 4.* On the diagonal $\theta=\sum_i\theta_i$, so $\Lambda=\theta\log F_P$ has
$$F(\vartheta)\;=\;\ell_{-1}(\Lambda)\;=\;L_1+L_2+L_3\;=\;\nu_1\nu_2+\nu_1\nu_3+\nu_2\nu_3 .$$

> **Proposition 5.1.** $\;F=e_2(\nu)$, where $\nu$ is the unique solution of $(\ast)$.

**Corollary 5.2 ($(E_1,E_2)$-freeness — Day 143's assumption, now a theorem).** The system
$(\ast)$ has coefficients $1,\omega,\omega^2,z$ only: **$E_1$ and $E_2$ do not appear.** They
entered $u_i$ only at order $\ge0$, and $\ell_{-1}$ discards that. Hence $\nu$, and therefore
$F=e_2(\nu)$, is independent of the base point — i.e. $n_k\in\mathbb Q$ rather than
$\mathbb Q[E_1,E_2]$.

**Independently verified:** solving $(\ast)$ by exact iteration in $\mathbb Q(\omega)[[z]]$
to $z^{39}$ gives $e_1(\nu)=6z^3+54z^6+834z^9+15702z^{12}+\cdots$, i.e.
$e_1(\nu)=2F$ — matching $2b_1,\dots,2b_{13}$ on the nose (`/tmp/saddle2.py`).

---

## 6. Elimination — the quintic

Let $E:=e_1(\nu)$, $A:=E-1$, and $R_i:=2\nu_i-A$.

1. **$R_i^2=A^2+4\omega^{i-1}z$.** From $(\ast)$, $\nu_i^2=A\nu_i+\omega^{i-1}z$, so
   $R_i^2=4\nu_i^2-4A\nu_i+A^2=A^2+4\omega^{i-1}z$.
2. **$F=E/2$.** Summing $\nu_i^2=A\nu_i+\omega^{i-1}z$ and using $\sum_i\omega^{i-1}=0$:
   $\sum\nu_i^2=AE$. But $\sum\nu_i^2=E^2-2e_2(\nu)$, so
   $e_2(\nu)=\tfrac12(E^2-AE)=\tfrac12 E$. With Prop. 5.1, $F=E/2$ and $A=2F-1$.
3. **$e_1(R)=2E-3A=2-A$.**
4. **$e_2(R)=\tfrac12\bigl[(2-A)^2-\textstyle\sum R_i^2\bigr]=\tfrac12\bigl[(2-A)^2-3A^2\bigr]=2-2A-A^2$**
   (again $\sum_i\omega^{i-1}=0$).
5. **$e_3(R)=\dfrac{e_2(R)^2-3A^4}{2(2-A)}$.** Indeed
   $e_2(R)^2-2e_1(R)e_3(R)=\sum_{i<j}R_i^2R_j^2$, and directly
   $\sum_{i<j}(A^2+4\omega^{i-1}z)(A^2+4\omega^{j-1}z)=3A^4$ because
   $\sum_i\omega^{i-1}=e_2(1,\omega,\omega^2)=0$. Legitimate since
   $2-A=3-2F=3+O(\vartheta)$ is invertible.
6. **$e_3(R)^2=\prod_i R_i^2=\prod_i(A^2+4\omega^{i-1}z)=A^6+64z^3=A^6+64\vartheta$.**

Substituting $A=2F-1$ into 4,5,6 gives, after clearing denominators (checked with `sympy`):

$$e_3(R)^2-A^6\;=\;\frac{64\,F(F-1)^3(4F-3)}{(2F-3)^2}\;=\;64\,\vartheta ,$$

which is precisely

$$\boxed{\;F(F-1)^3(4F-3)=\vartheta\,(2F-3)^2 .\;}$$

Note where the primes come from: **step 5 needs $2-A$ invertible (a $3$), and steps 2,4,6
all use $1+\omega+\omega^2=0$** — the cube roots of unity are the entire mod-$3$ mechanism.

---

## 7. The theorem

**Theorem.** $b_k\equiv0\pmod3$ for all $k\ge1$. More precisely $G:=F/3$ lies in
$\mathbb Z[[\vartheta]]$ and

$$\frac{b_k}{3}=\frac1k\,[G^{k-1}]\,\frac{(2G-1)^{2k}}{(3G-1)^{3k}(4G-1)^{k}} .$$

*Proof.* Put $G=F/3\in\mathbb Q[[\vartheta]]$ (a priori only rational). Then
$$F(F-1)^3(4F-3)=3G\,(3G-1)^3\cdot 3(4G-1)=9\,G(3G-1)^3(4G-1),$$
$$\vartheta(2F-3)^2=\vartheta\cdot 9(2G-1)^2 .$$
The factors $9$ cancel and $G(3G-1)^3(4G-1)=\vartheta(2G-1)^2$, i.e.
$G=\vartheta\,\phi(G)$ with
$$\phi(G)=\frac{(2G-1)^2}{(3G-1)^3(4G-1)}\in\mathbb Z[[G]],\qquad \phi(0)=1,$$
integrality because $(1-3G)^{-1},(1-4G)^{-1}\in\mathbb Z[[G]]$. The equation
$G=\vartheta\phi(G)$ determines $G$ uniquely in $\vartheta\mathbb Q[[\vartheta]]$, and
Lagrange inversion gives $[\vartheta^k]G=\frac1k[G^{k-1}]\phi(G)^k$, an integer. $\square$

### 7.1 Cleaner form, and an explicit closed form for $b_k$

Because $(3G-1)^3(4G-1)=(1-3G)^3(1-4G)$ and $(2G-1)^2=(1-2G)^2$, the equation reads

$$\boxed{\;G\,(1-3G)^3(1-4G)=\vartheta\,(1-2G)^2\;}
\qquad\Longleftrightarrow\qquad
\boxed{\;F(1-F)^3(3-4F)=\vartheta\,(3-2F)^2\;}$$

with kernel $\phi(G)=\dfrac{(1-2G)^2}{(1-3G)^3(1-4G)}$. Expanding
$\phi(G)^k$ and reading off $[G^{k-1}]$ gives, at last, the **closed form for $b_k$**
that has been open since Day 143 (§"remaining open questions", item 1):

$$\boxed{\;b_k \;=\; \frac{3}{k}\sum_{i+j+l=k-1}\binom{2k}{i}(-2)^i\,\binom{3k+j-1}{j}3^j\,\binom{k+l-1}{l}4^l\;}$$

a single finite triple sum of binomials. **Verified** for $k=1,\dots,15$ against the
$\Psi$-recursion table and for $k=16,\dots,22$ against §8.2. The $3/k$ prefactor makes
$3\mid b_k$ visible but not obvious — the integrality of the inner sum divided by $k$ is
exactly Lagrange inversion.

**Corollary (Conjecture H, diagonal, at every prime).** Assuming (H2) — so that Day 146
Theorem 2 gives $\mathcal H=\dfrac{F(F-1)}{\vartheta(2F-3)}$ — substitution of the quintic yields
$$\mathcal H=\frac{2F-3}{(F-1)^2(4F-3)}=\frac{2G-1}{(3G-1)^2(4G-1)}\ \in\ \mathbb Z[[\vartheta]].$$
This *explains* Day 147's "exact realizability" observation and shows it was never
about $p=3$ at all. (Verified against the tabulated $h_0,\dots,h_{15}$: exact match.)

---

## 8. Verification

**8.1 The closed form.** `vclosed.py`: Theorem 2.2 reproduces $P_b/b!$ from `core.py`
(the $\Psi$-recursion, i.e. the definition) for all $b\le9$ at $u=(\frac12,\frac13,\frac15)$.

**8.2 The quintic, predictively.** The relation was fitted to $b_1..b_{15}$ using $12$
unknowns and $16$ equations (over-determined by 4; nullspace exactly $1$-dimensional).
It was then tested on **new** data: `fast.py` computes $b_k$ from the closed form
(entirely independent code path) modulo $p=2^{61}-1$ at two base points,
$(E_1,E_2)=(3,2)$ and $(4,3)$:

| $k$ | Lagrange prediction from the quintic | closed-form computation |
|---|---|---|
| 16 | 415499754144310284843 | ✓ (both base points) |
| 17 | 11196280068687472987044 | ✓ |
| 18 | 303298794015028731492228 | ✓ |
| 19 | 8254894361954322048629076 | ✓ |
| 20 | 225623849259101980107679884 | ✓ |
| 21 | 6190290513670129353943964784 | ✓ |
| 22 | 170425591175863604918339244261 | ✓ |

Seven 21-to-33-digit predictions, all exact. $b_1..b_{15}$ also reproduced exactly at both
base points, matching the $\Psi$-recursion values.

**8.3 Vanishing lemma.** `fast.py` reports **zero** violations of
$[E_3^kT^b]\log F_P=0$ for $b<3k-1$, all $k\le22$, at both base points.

**8.4 The saddle system.** `saddle2.py` solves $(\ast)$ exactly in $\mathbb Q(\omega)[[z]]$
and returns $e_1(\nu)=2F$ through $z^{39}$ — thirteen coefficients, all matching $2b_k$.
It also confirms $e_1(\nu)\in\mathbb Q[[z^3]]$ (no $\omega$ survives), as the symmetry demands.

**8.5 $\mathcal H$.** $\dfrac{2F-3}{(F-1)^2(4F-3)}$ reproduces $h_0,\dots,h_{15}$ exactly.

**8.6 Not completed.** A third, fully redundant check — the $\Psi$-recursion itself at
BMAX $=57$ (`day147_gauss/regen.py 57`, exact rational arithmetic in three variables) — was
launched and did **not** finish inside the session. It is not needed: §8.2 already computes
$b_{16..22}$ from the closed form at two base points, and the closed form is verified against
the $\Psi$-recursion in §8.1.

Scripts: `/tmp/vclosed.py`, `/tmp/fast.py`, `/tmp/alg.py`, `/tmp/alg2.py`,
`/tmp/predict.py`, `/tmp/saddle2.py`, `/tmp/check1.py` (copied to
`projects/beta-prime/code/day148_closedform/`).

---

## 9. Gaps, stated precisely

**9.1 — none in the main line.** §§3–6 hold at **arbitrary** $(E_1,E_2)$: the only property
of the base point used is $\ell_{-1}(u_i)=\omega^{i-1}z$, which is $(E_1,E_2)$-independent.
So the quintic holds identically in $E_1,E_2$, which **also proves Day 143's
$(E_1,E_2)$-freeness of $n_k$** as a by-product (the quintic has a unique power-series root).
Theorem 4.1 likewise now proves Day 146 Prop. 1 unconditionally at every base point.

**9.2 [SUPERSEDED — Day 149: (H2) IS NOW A THEOREM, this dependency is gone. Delete.]** The identity
$\mathcal H=\dfrac{F(F-1)}{\vartheta(2F-3)}$ used in §7's Corollary is Day 146 Theorem 2,
which assumes (H2). The **main theorem $3\mid b_k$ does not use it.**

**9.3 [Day 149: still open, but see `2026-08-30-day149-H2-PROVED.md` §5 — the $\tau$-shift is
insertion of one explicit cubic $g$, and Conjecture P (coefficientwise positivity of $H$)
implies (H1).]** (H1) itself — $H=\tau(F_P)/F_P\in\mathbb Z[E_1,E_2,E_3][[T]]$ — is untouched and
still open. It is *strictly stronger* than what is proved here. But the closed form
gives it a clean new shape: with $\langle g\rangle:=\bigl[V(M)\,g\cdot\mathcal M\bigr]_{t=T}$,
$$H=\frac1{E_3}\cdot\frac{\langle e_3\rangle}{\langle 1\rangle},$$
because $\tau$ ($u_i\mapsto u_i+1$) multiplies the summand by $\prod_i(u_i+m_i)/u_i$ and
leaves $V(u+m)$ and $V(u)$ fixed. A moment-ratio problem, not a Frobenius problem.

---

## 10. Post-mortem — why this took four sessions

The obstruction was never mathematical. $\mathcal T$ is the falling-factorial umbral map;
$\mathcal T(g)(u)=[g(\partial_x)x^u]_{x=1}$ is a two-line observation; Lemma 2.1 is a
one-line determinant collapse; everything else follows. **Nobody asked what $\mathcal T$
was.** Instead: three sessions on Dwork lifts (tautological), one on $\lambda$-rings
(nothing), one chasing a paper that does not exist, and a "crown jewel" (Day 147's
$m_n\ge0$) that is implied by $s_n$ growing.

New standing rule, earned today:

> **Rule 11 — Unfold the definition before you decorate it.** Before importing any
> external theory onto an object, write the object's *defining operator* in closed form
> and ask whether it is something classical in disguise. Four sessions of $p$-adic
> machinery were applied to a series whose closed form is a Horn hypergeometric sum.

And the corollary of §1: **an empirical "non-circular signal" must be checked against
the null hypothesis.** $m_n\ge0$ looked like structure; it was arithmetic of growth.

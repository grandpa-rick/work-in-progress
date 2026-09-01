# Day 152 — **(P1) and (P2) are theorems**; the closed form for $\psi$ is **proved**

**Date:** 2026-08-31 (Day 152). **Status: the target is closed.**

* **(P1)** $\;\log\ell_0^{\mathrm{top}}(H)=\partial\,\Xi\;$ — **proved** (§3). Six lines from Day 149
  Theorem 1. It is a pure statement about the weight filtration; nothing analytic happens.
* **(P2)** $\;\theta\,\Xi=\tfrac12(P-E_1)\;$ — **proved** (§4), together with the thing it really
  rests on, namely that the $\nu$-system *is* the top-weight symbol of the Riccati system (R).
* Consequently the whole Day 151 chain is a proof (§5), and

$$\boxed{\;\psi \;=\; q\,\frac{e_3(\nu)}{E_3}\;=\;\frac{q}{\prod_i(q+T\nu_i)}
\;=\;\frac{4q(q+2)}{(q+1)^2\bigl(2q+1-2E_1T\bigr)+\Delta_2T^2}\;}$$

  with $q$ the branch at $q=1+O(T)$ of the master quintic. The third expression is **new** and is
  strictly better than the Day 151 form $\frac{q(1-2E_1T+\Delta_2T^2-q^2)}{4E_3T^3(q+2)}$: it has
  no $E_3$ and no $T^3$ in the denominator, so it is honest at $E_3=0$ and needs no $0/0$ excuse.
* **The master quintic gets a two-line derivation** (§6), replacing Day 151's radical-clearing +
  double resultant. It is nothing but
  $$\prod_i(R_i^2-q^2)=64T^3E_3,\qquad\text{i.e.}\qquad
  \underbrace{\bigl(2q^2-e_2(R)q+e_3(R)\bigr)}_{=\,8T^3e_3(\nu)}\cdot
  \underbrace{\bigl(2q^3+2q^2+e_2(R)q+e_3(R)\bigr)}_{=\,8E_3/e_3(\nu)}=64T^3E_3 ,$$
  once you notice that $e_2(R)$ and $e_3(R)$ are **rational** in $(q,T,E)$ — no square roots
  anywhere. $e_3(R)$ is rational because the constraint $e_2(u)=E_2$ is *linear* in it.
* $\;\deg_\psi Q=5$, $\deg_YQ=9$, monic, irreducible ⟹ $\psi$ is algebraic of degree **exactly 5**
  over $\mathbb Q(E_1,E_2,E_3)(Y)$ (§7). Same $Q$ as Day 151, coefficient for coefficient, but
  reached by a different elimination — and the spurious factors are now **identified**
  ($\psi^9$ and $(3\psi^2+2E_1Y\psi-\Delta_2Y^2)^2$), closing a Day 151 loose end.
* Everything below was re-verified on a **clean-room code path** built today from the raw
  definition of $F_P$ (`scratch/day152/`), exact `Fraction` arithmetic, symbolic $E_1,E_2,E_3$
  (§8).

**Registry:** `psi-closed-form-degree5` moves `computed` → **`proved`**.


> **[Day 152b — INDEPENDENTLY AUDITED, and Theorem D strengthened.**
> Every step of this document, plus its two imports (Day 148's (R), Day 149's Theorem 1), was
> re-derived by hand and every script re-run on
> `proofs/2026-08-31-day152b-audit-of-psi-closed-form.md`. **No error found.** Two additions:
> (i) the pipeline `raw F_P -> H -> Y -> psi -> Q` was run **symbolically in $E_1,E_2,E_3$**
> (§8's `check5.py` was numeric only) — $Q(\psi(Y),Y)=0$ through $Y^{12}$, and it reproduces the
> pre-registered $[Y^9]\psi=2E_2^3E_3+22E_1E_2E_3^2+34E_3^3$ from a third code path;
> (ii) **Theorem D no longer depends on any factorisation black box** — see the certificate box
> in §7.**]**

---

## 1. Notation, and the one grading that does all the work

$u_1,u_2,u_3$ indeterminates, $E_i=e_i(u)$, $V=\prod_{i<j}(u_i-u_j)$,
$\mathcal T^+:u^\alpha\mapsto\prod_iu_i^{(\alpha_i)}$ (rising factorials),
$\Psi^+(f)=\mathcal T^+(fV)/V$,
$$F_P:=\Psi^+\!\bigl(e^{Te_2}\bigr)=\sum_{b\ge0}\Psi^+(e_2^b)\frac{T^b}{b!}\in\mathbb Q[E_1,E_2,E_3][[T]],
\qquad F_P=1+O(T).$$
$\tau$ = the automorphism induced by $u_i\mapsto u_i+1$; $H:=\tau(F_P)/F_P$;
$\partial:=\sum_i\partial_{u_i}$, which on $\Lambda_3=\mathbb Q[E]$ is
$\partial=3\partial_{E_1}+2E_1\partial_{E_2}+E_2\partial_{E_3}$.

Work in $\mathcal A:=\mathbb Q[u_1,u_2,u_3][[t_1,t_2,t_3]]$ (and its localisation at $V$), with the
**weight grading**
$$\operatorname{wt}(u_i)=+1,\qquad \operatorname{wt}(t_i)=-1 .$$
For $X\in\mathcal A$ and $w\in\mathbb Z$ put
$$\ell^{\mathrm{top}}_w(X):=\sum_{\alpha}t^\alpha\cdot\Bigl(\text{the $u$-homogeneous part of degree }|\alpha|+w\text{ of }[t^\alpha]X\Bigr).$$
Since each $[t^\alpha]X$ is a *polynomial* in $u$, only finitely many $w$ are hit by a given $\alpha$,
so $X=\sum_w\ell^{\mathrm{top}}_w(X)$ is a legitimate decomposition of $\mathcal A$ as a completed
direct sum $\widehat\bigoplus_w\mathcal A_w$. Write $\operatorname{wt}(X)\le w_0$ for
"$\ell^{\mathrm{top}}_w(X)=0$ for all $w>w_0$".

**Lemma 1.1 (bookkeeping).** All of the following are immediate from
$\mathcal A_w\cdot\mathcal A_{w'}\subseteq\mathcal A_{w+w'}$.
1. $\{\operatorname{wt}\le0\}$ is a subring and $\ell^{\mathrm{top}}_0$ restricted to it is a ring
   homomorphism.
2. If $\operatorname{wt}X\le a$ and $\operatorname{wt}Y\le b$ then $\operatorname{wt}(XY)\le a+b$
   and $\ell^{\mathrm{top}}_{a+b}(XY)=\ell^{\mathrm{top}}_a(X)\,\ell^{\mathrm{top}}_b(Y)$.
3. $\theta_i:=t_i\partial_{t_i}$ preserves $\operatorname{wt}$ (it multiplies $t^\alpha u^\beta$ by
   $\alpha_i$), hence commutes with every $\ell^{\mathrm{top}}_w$.
4. $\partial=\sum_i\partial_{u_i}$ maps $\mathcal A_w\to\mathcal A_{w-1}$; hence
   $\ell^{\mathrm{top}}_{w-k}(\partial^kX)=\partial^k\ell^{\mathrm{top}}_w(X)$.
5. The diagonal restriction $\varrho:t_i\mapsto T$ maps $\mathcal A_w$ into the corresponding graded
   piece of $\mathbb Q[u][[T]]$, so $\ell^{\mathrm{top}}_w\circ\varrho=\varrho\circ\ell^{\mathrm{top}}_w$.
6. $\ell^{\mathrm{top}}_w$ is computed $t$-degreewise, hence commutes with any $t$-adically
   convergent sum — in particular with $\exp$ and $\log$ of series in $1+(t)$.

On $\mathbb Q[E][[T]]$ the grading reads: $\operatorname{wt}(E_1^iE_2^jE_3^kT^n)=(i+2j+3k)-n$, i.e.
$\operatorname{wt}X\le w$ says $\deg_u[T^n]X\le n+w$.

## 2. The two imported facts

Both are **proved**, and both were re-verified today on a fresh code path (§8).

> **Fact I (Day 148, Cor 2.3 and §4 — the Horn/Riccati system).** With
> $\mathcal M(t)=\sum_{a,b,c}\frac{t_1^at_2^bt_3^c}{a!b!c!}u_1^{(a+b)}u_2^{(a+c)}u_3^{(b+c)}$,
> $M_1=u_1+\theta_1+\theta_2$, $M_2=u_2+\theta_1+\theta_3$, $M_3=u_3+\theta_2+\theta_3$, one has
> $F_P=\varrho\bigl(V(M)\mathcal M/V(u)\bigr)$, and with $S:=\log\mathcal M$,
> $\lambda_i:=\theta_iS$, $d_1=\lambda_1+\lambda_2$, $d_2=\lambda_1+\lambda_3$,
> $d_3=\lambda_2+\lambda_3$, $D_1=\theta_1+\theta_2$, $D_2=\theta_1+\theta_3$:
> $$\lambda_1=t_1\bigl[(u_1{+}d_1)(u_2{+}d_2)+D_1d_2\bigr],\quad
> \lambda_2=t_2\bigl[(u_1{+}d_1)(u_3{+}d_3)+D_1d_3\bigr],\quad
> \lambda_3=t_3\bigl[(u_2{+}d_2)(u_3{+}d_3)+D_2d_3\bigr].\tag{R}$$

> **Fact II (Day 149, Theorem 1 and its proof).**
> (a) $\operatorname{wt}(\lambda_i)\le1$ for $i=1,2,3$, hence $\operatorname{wt}(S)\le1$;
> (b) $\operatorname{wt}\bigl(\mathcal R/V(u)\bigr)\le0$, where $\mathcal R:=e^{-S}V(M)e^{S}$, so that
> $\mathcal R/V(u)=F_P^{\,\mathrm{multi}}/\mathcal M\in1+(t)$;
> (c) consequently $\operatorname{wt}(\log F_P)\le1$, i.e. $\deg_u[T^n]\log F_P\le n+1$.
>
> *(For (a) the induction is on $t$-degree using (R); for (b), $\mathcal R$ is a sum of products of
> at most three factors each of $\operatorname{wt}\le1$, so $\operatorname{wt}\mathcal R\le3$, and
> $\operatorname{wt}V(u)=3$. (c) follows since
> $\log F_P=\varrho(S)+\varrho\log(\mathcal R/V(u))$.)*

**Definition.** $\;\Xi:=\ell^{\mathrm{top}}_1(\log F_P)\in\mathbb Q[E][[T]]$, and
$\widetilde\Xi:=\ell^{\mathrm{top}}_1(S)\in\mathcal A$.

By Fact II(b) and Lemma 1.1(5), $\;\Xi=\varrho\bigl(\widetilde\Xi\bigr)$: the prefactor is
invisible at the top weight.

---

## 3. **(P1)**

> **Theorem A.** $\operatorname{wt}(H)\le0$ and
> $$\ell^{\mathrm{top}}_0(H)=\exp\bigl(\partial\,\Xi\bigr).$$

*Proof.* Write $X:=\log F_P\in T\,\mathbb Q[E][[T]]$; $\operatorname{wt}X\le1$ by Fact II(c).

$\tau$ is translation by $1$ in every $u_i$, and each $[T^n]X$ is a polynomial, so Taylor's formula
gives $\tau=\sum_{k\ge0}\partial^k/k!$ (a **finite** sum on each $[T^n]$), whence
$$\log H=(\tau-1)X=\sum_{k\ge1}\frac{\partial^kX}{k!}.$$
By Lemma 1.1(4), $\partial^k$ shifts weight by $-k$, so $\operatorname{wt}(\partial^kX)\le1-k\le0$
for $k\ge1$; therefore $\operatorname{wt}(\log H)\le0$, and
$$\ell^{\mathrm{top}}_0(\log H)=\sum_{k\ge1}\frac{1}{k!}\,\ell^{\mathrm{top}}_0(\partial^kX)
=\sum_{k\ge1}\frac{1}{k!}\,\partial^k\bigl(\ell^{\mathrm{top}}_k X\bigr)
=\partial\bigl(\ell^{\mathrm{top}}_1X\bigr)=\partial\,\Xi,$$
because $\ell^{\mathrm{top}}_kX=0$ for $k\ge2$. Finally $\log H\in T\,\mathbb Q[E][[T]]$, so
$H=\exp(\log H)$ converges $T$-adically; by Lemma 1.1(1),(6),
$\ell^{\mathrm{top}}_0(H)=\exp\bigl(\ell^{\mathrm{top}}_0\log H\bigr)=\exp(\partial\Xi)$. $\square$

*Remark.* This is exactly Day 149's Theorem 2 with the extra observation that at the top weight
$(\tau-1)$ **is** $\partial$ — every higher Taylor term drops two or more weights and cannot reach
level $0$ from level $\le1$. That is the whole content, and it is why the statement was safe to
assert: it is Rule 12 read one layer up.

---

## 4. **(P2)**

> **Theorem B.** Let $\nu_1,\nu_2,\nu_3\in\mathbb Q[u][[T]]$ be the unique solution with
> $\nu_i=u_i+O(T)$ of
> $$\nu_i\Bigl(1-T\bigl(e_1(\nu)-\nu_i\bigr)\Bigr)=u_i\qquad(i=1,2,3),\tag{$*$}$$
> and put $P:=e_1(\nu)$. Then, with $\theta:=\sum_i\theta_i$ (which is $T\,d/dT$ after $\varrho$),
> $$\theta\,\Xi\;=\;T\,e_2(\nu)\;=\;\frac{P-E_1}{2}.$$

*Proof.* **Step 1 (the top-weight symbol of (R)).** Put $L_i:=\ell^{\mathrm{top}}_1(\lambda_i)$;
by Lemma 1.1(3) and Fact II(a), $L_i=\theta_i\widetilde\Xi$. Apply $\ell^{\mathrm{top}}_1$ to the
first equation of (R). Multiplication by $t_1$ shifts weight by $-1$, so
$$L_1=t_1\cdot\ell^{\mathrm{top}}_2\Bigl[(u_1{+}d_1)(u_2{+}d_2)+D_1d_2\Bigr].$$
The term $D_1d_2$ has $\operatorname{wt}\le1$ (Fact II(a) and Lemma 1.1(3)), so it does not reach
weight $2$ and contributes nothing. Each factor $u_i+d_i$ has $\operatorname{wt}\le1$, with
$$\ell^{\mathrm{top}}_1(u_1{+}d_1)=u_1+L_1+L_2=:\nu_1,\quad
\ell^{\mathrm{top}}_1(u_2{+}d_2)=u_2+L_1+L_3=:\nu_2,\quad
\ell^{\mathrm{top}}_1(u_3{+}d_3)=u_3+L_2+L_3=:\nu_3$$
(the constants $u_i$ are themselves of pure weight $1$). By Lemma 1.1(2),
$\ell^{\mathrm{top}}_2$ of the product is the product of the two $\ell^{\mathrm{top}}_1$'s. Hence
$$L_1=t_1\nu_1\nu_2,\qquad L_2=t_2\nu_1\nu_3,\qquad L_3=t_3\nu_2\nu_3.\tag{4.1}$$

**Step 2 (restrict to the diagonal).** Applying $\varrho$ and keeping the same letters,
$\nu_1=u_1+T\nu_1\nu_2+T\nu_1\nu_3=u_1+T\nu_1(e_1(\nu)-\nu_1)$, and likewise for $\nu_2,\nu_3$;
that is exactly $(*)$. Since $\lambda_i$ has $t$-degree $\ge1$ we have $L_i=O(T)$, so
$\nu_i=u_i+O(T)$ and these are *the* solutions of $(*)$: $(*)$ determines $[T^n]\nu_i$ recursively
from lower orders. (The solution is $S_3$-equivariant, so $e_k(\nu)\in\mathbb Q[E][[T]]$.)

**Step 3.** $\theta\widetilde\Xi=\sum_iL_i$, so by (4.1) and $\varrho$,
$$\theta\,\Xi=T\bigl(\nu_1\nu_2+\nu_1\nu_3+\nu_2\nu_3\bigr)=T\,e_2(\nu).$$

**Step 4 (the sum trick).** Expand $(*)$: $T\nu_i^2+q\nu_i=u_i$ with $q:=1-TP$. Summing over $i$
and using $\sum\nu_i^2=P^2-2e_2(\nu)$:
$$T\bigl(P^2-2e_2(\nu)\bigr)+(1-TP)P=E_1\ \Longrightarrow\ P-2T\,e_2(\nu)=E_1
\ \Longrightarrow\ T\,e_2(\nu)=\frac{P-E_1}{2}. \qquad\square$$

*Remark (why the other reading had to fail).* $\theta$ is $\sum_it_i\partial_{t_i}$, and for any
$f(t_1,t_2,t_3)$ the chain rule gives $\varrho(\sum_it_i\partial_{t_i}f)=T\frac{d}{dT}\varrho(f)$.
There is no room for the $u$-Euler operator: the $t$'s are where the grading lives. The Day 151
correction box was right, and now it is forced rather than measured.

---

## 5. The chain: $\;\mathcal W=e_3(\nu)/E_3$, $\;Y=Tq\,\mathcal W$, $\;\psi=q\,e_3(\nu)/E_3$

All computations here are in $\mathbb Q[u][[T]]$; **no radicals are needed**. Set
$$q:=1-TP,\qquad R_i:=q+2T\nu_i,\qquad S_1:=\sum_i\frac1{R_i},\qquad n_3:=e_3(\nu),$$
$$\mathcal W:=\ell^{\mathrm{top}}_0(H),\qquad Y:=\int_0^T\mathcal W\,dT\in T+T^2\mathbb Q[E][[T]],$$
and let $\psi\in\mathbb Q[E][[Y]]$ be the unique series with $Y=T\,\psi(Y)$ (the *Lagrange kernel*
of $Y$; it exists because $Y=T+O(T^2)$).

Elementary consequences of $T\nu_i^2+q\nu_i=u_i$, used throughout:
$$R_i^2=q^2+4Tu_i,\qquad \sum_iR_i=3q+2TP=q+2,\qquad
\frac{T\nu_i}{R_i}=\frac12-\frac{q}{2R_i},\qquad q+T\nu_i=\frac{u_i}{\nu_i}.\tag{5.1}$$
(The last: $\nu_i(T\nu_i+q)=u_i$.) Each $R_i\in1+T\mathbb Q[u][[T]]$ is invertible, and
$qS_1-1\in2+T(\cdots)$ is invertible.

**Lemma 5.1.** $\;A:=\partial P=\dfrac{2S_1}{qS_1-1}\;$ and $\;T\dot q=\dfrac{q^2S_1-q-2}{2(qS_1-1)}\;$
(dot $=d/dT$).

*Proof.* Apply $\partial$ to $T\nu_i^2+q\nu_i=u_i$: $(q+2T\nu_i)\partial\nu_i=1-\nu_i\partial q$, and
$\partial q=-T\partial P=-TA$, so $\partial\nu_i=(1+TA\nu_i)/R_i$. Summing and using (5.1),
$$A=S_1+A\sum_i\frac{T\nu_i}{R_i}=S_1+A\Bigl(\tfrac32-\tfrac q2S_1\Bigr)
\ \Longrightarrow\ A\cdot\tfrac{qS_1-1}{2}=S_1 .$$
Apply $d/dT$: $R_i\dot\nu_i=-\nu_i(\nu_i+\dot q)$, i.e.
$$\frac{\dot\nu_i}{\nu_i}=-\frac{\nu_i+\dot q}{R_i}.\tag{5.2}$$
Multiply by $T$ and sum, using $T\nu_i/R_i=\frac12-\frac q{2R_i}$:
$$T\dot P=-\sum_i(\nu_i+\dot q)\frac{T\nu_i}{R_i}
=-\tfrac12(P+3\dot q)+\tfrac q2\Bigl(\tfrac1T\bigl(\tfrac32-\tfrac q2S_1\bigr)+\dot qS_1\Bigr).$$
Also $q=1-TP$ gives $T\dot P=-\dot q-P$ and $P=(1-q)/T$. Substituting and collecting,
$$-\tfrac{\dot q}{2}\bigl(qS_1-1\bigr)=\frac{2+q-q^2S_1}{4T},$$
which is the second formula. $\square$

**Lemma 5.2.** $\;T\dfrac{d}{dT}\log n_3=\dfrac{A-3}{2}$.

*Proof.* By (5.2), $T\frac{d}{dT}\log n_3=-\sum_i\frac{T\nu_i}{R_i}-T\dot q\,S_1
=-\tfrac32+\tfrac q2S_1-T\dot qS_1$. Insert Lemma 5.1 and use the identity
$q(qS_1-1)-(q^2S_1-q-2)=2$:
$$\tfrac q2S_1-\frac{S_1(q^2S_1-q-2)}{2(qS_1-1)}=\frac{S_1\bigl[q(qS_1-1)-(q^2S_1-q-2)\bigr]}{2(qS_1-1)}
=\frac{S_1}{qS_1-1}=\frac A2 .\qquad\square$$

> **Theorem C.** $\;\displaystyle \mathcal W=\frac{e_3(\nu)}{E_3}=\prod_i\frac{\nu_i}{u_i}
> =\prod_i\frac1{q+T\nu_i},\qquad Y=Tq\,\mathcal W,\qquad
> \psi=q\,\frac{e_3(\nu)}{E_3}=\frac{q}{\prod_i(q+T\nu_i)} .$

*Proof.* By Theorem A, $\log\mathcal W=\partial\Xi$; $\partial$ commutes with $T\,d/dT$, so by
Theorem B,
$$T\frac{d}{dT}\log\mathcal W=\partial\Bigl(T\frac{d\Xi}{dT}\Bigr)=\partial\Bigl(\frac{P-E_1}{2}\Bigr)
=\frac{\partial P-3}{2}=\frac{A-3}{2},$$
using $\partial E_1=3$. By Lemma 5.2 the series $\log\mathcal W$ and $\log(n_3/E_3)$ have the same
image under $T\,d/dT$, which kills only constants; and both $\mathcal W$ and $n_3/E_3=\prod\nu_i/u_i$
equal $1$ at $T=0$. Hence $\mathcal W=n_3/E_3$, and the last equality is (5.1).

For $Y$: set $Z:=Tq\,\mathcal W$. Since $\mathcal W=Y'$ and $T\mathcal W'/\mathcal W=\frac{A-3}2$,
$$\frac{\dot Z}{\mathcal W}=q+T\dot q+q\cdot\frac{A-3}{2}
=-\frac q2+\frac{q^2S_1-q-2+2qS_1}{2(qS_1-1)}
=-\frac q2+\frac{(q+2)(qS_1-1)}{2(qS_1-1)}=1 ,$$
using Lemma 5.1 twice and $q^2S_1+2qS_1-q-2=(q+2)(qS_1-1)$. So $\dot Z=\mathcal W=\dot Y$ and
$Z(0)=0=Y(0)$, giving $Y=Tq\mathcal W$. Therefore $\psi(Y)=Y/T=q\,\mathcal W=q\,n_3/E_3$. $\square$

---

## 6. The master curve, rationally

Let $c:=e_3(R)$. From (5.1), $e_1(R)=q+2$ and $\sum_iR_i^2=3q^2+4TE_1$, so
$$e_2(R)=\frac{(q+2)^2-(3q^2+4TE_1)}{2}=-q^2+2q+2-2E_1T .$$
Now use the *other two* symmetric functions of $u_i=(R_i^2-q^2)/(4T)$. With $A_i:=R_i^2$ one has
$e_1(A)=3q^2+4E_1T$, $e_2(A)=e_2(R)^2-2e_1(R)c$, $e_3(A)=c^2$, and

* $e_2(u)=E_2$ reads $\;e_2(A)-2q^2e_1(A)+3q^4=16T^2E_2$, which is **linear in $c$**:
$$\boxed{\;e_3(R)=c=\frac{2+4q-2q^3-q^4-2E_1T\,(q^2{+}2q{+}2)+2\Delta_2T^2}{q+2}\;},\qquad \Delta_2:=E_1^2-4E_2 .$$
  (At $T=0,q=1$: $c=3/3=1=e_3(1,1,1)$. ✓ The division by $q+2\in3+T(\cdots)$ is legitimate.)
* $e_3(u)=E_3$ reads $\;\prod_i(A_i-q^2)=64T^3E_3$. Writing $f(x)=\prod_i(x-R_i)=x^3-e_1(R)x^2+e_2(R)x-c$ one has
  $$f(q)f(-q)=\prod_i(q-R_i)(-q-R_i)=(-1)^3\prod_i(q^2-R_i^2)=\prod_i(R_i^2-q^2),$$
  with $f(q)=-2q^2+e_2(R)q-c$ and $f(-q)=-\bigl(2q^3+2q^2+e_2(R)q+c\bigr)$, so

$$\boxed{\;\bigl(2q^2-e_2(R)q+c\bigr)\bigl(2q^3+2q^2+e_2(R)q+c\bigr)=64\,T^3E_3\;}\tag{6.1}$$

  and substituting the rational $e_2(R),c$ and multiplying through by $(q+2)^2$ turns (6.1) into a polynomial
  equation of degree $5$ in $q$ and $4$ in $T$. **It is exactly $-4\times$ the Day 149 quintic**, verified
  exactly in `sympy`:
$$(q-1)(q+1)^3(2q+1)-2E_1T(q+1)^2(q^2{-}2q{-}2)-2T^2(q+1)\bigl[\Delta_2(q^2{+}q{+}1)+2E_1^2(q+1)\bigr]$$
$$+\,2E_1T^3\Delta_2(q^2{+}2q{+}2)+16E_3T^3(q+2)^2-T^4\Delta_2^{\,2}\;=\;0. \tag{6.2}$$

**What (6.1) actually says.** $q-R_i=-2T\nu_i$ and $q+R_i=2(q+T\nu_i)$, so
$$2q^2-e_2(R)q+c=-f(q)=\prod_i(R_i-q)=\prod_i 2T\nu_i=8T^3e_3(\nu),\qquad
2q^3+2q^2+e_2(R)q+c=-f(-q)=\prod_i(q+R_i)=\prod_i2(q+T\nu_i)=\frac{8E_3}{e_3(\nu)} .$$
So (6.1) is the tautology $8T^3n_3\cdot\frac{8E_3}{n_3}=64T^3E_3$ — and it hands us **both**
factors in closed form:
$$\frac{e_3(\nu)}{E_3}=\mathcal W=\frac{2q^2-e_2(R)q+c}{8T^3E_3},\qquad
\prod_i(q+T\nu_i)=\frac{2q^3+2q^2+e_2(R)q+c}{8}.\tag{6.3}$$

> **Corollary (the clean closed form).** By Theorem C and (6.3),
> $$\psi=\frac{8q}{2q^3+2q^2+e_2(R)q+c}
> \;=\;\boxed{\;\frac{4q(q+2)}{(q+1)^2\bigl(2q+1-2E_1T\bigr)+\Delta_2T^2}\;}$$
> after multiplying numerator and denominator by $q+2$ and using
> $(q+2)\bigl(2q^3+2q^2+e_2(R)q\bigr)+(q+2)c=2\bigl[(q+1)^2(2q+1-2E_1T)+\Delta_2T^2\bigr]$
> (direct expansion). Equivalently, again by direct expansion,
> $(q+2)\bigl(2q^2-e_2(R)q+c\bigr)=2\bigl(1-2E_1T+\Delta_2T^2-q^2\bigr)$, which recovers the
> Day 151 form $\mathcal W=\frac{1-2E_1T+\Delta_2T^2-q^2}{4E_3T^3(q+2)}$,
> $\psi=q\mathcal W$.

The boxed form is the one to use: its denominator is a polynomial with **no $E_3$**, it is regular
at $E_3=0$, and it makes the $T=0$ normalisation obvious ($q=1$: $\psi=4\cdot3/(4\cdot3)=1$).

---

## 7. $\psi$ is algebraic of degree exactly $5$

> **Theorem D.** $Q(\psi,Y,E_1,E_2,E_3)=\sum_{j=0}^9Y^jC_j(\psi,E)$ below is irreducible over
> $\mathbb Q$, monic of degree $5$ in $\psi$, of degree $9$ in $Y$, and $Q(\psi(Y),Y)=0$. Hence
> $\psi$ is algebraic over $\mathbb Q(E_1,E_2,E_3)(Y)$ **of degree exactly $5$**, with minimal
> polynomial $Q$.

| $j$ | $C_j$ |
|---|---|
|0| $\psi^4(\psi-1)$ |
|1| $-E_1\psi^3(5\psi-4)$ |
|2| $\psi^2(10E_1^2\psi-6E_1^2-E_2\psi^2-8E_2\psi+8E_2)$ |
|3| $-2\psi(5E_1^3\psi-2E_1^3-2E_1E_2\psi^2-12E_1E_2\psi+8E_1E_2-11E_3\psi^2+44E_3\psi-32E_3)$ |
|4| $5E_1^4\psi-E_1^4-6E_1^2E_2\psi^2-24E_1^2E_2\psi+8E_1^2E_2+E_1E_3\psi^3-36E_1E_3\psi^2+80E_1E_3\psi+8E_2^2\psi^2+16E_2^2\psi-16E_2^2$ |
|5| $-E_1^5+4E_1^3E_2\psi+8E_1^3E_2-3E_1^2E_3\psi^2+6E_1^2E_3\psi+8E_1^2E_3-16E_1E_2^2\psi-16E_1E_2^2-20E_2E_3\psi^2+104E_2E_3\psi-32E_2E_3$ |
|6| $-E_1^4E_2+3E_1^3E_3\psi+8E_1^3E_3+8E_1^2E_2^2+4E_1E_2E_3\psi-32E_1E_2E_3-16E_2^3-E_3^2\psi^2+88E_3^2\psi-16E_3^2$ |
|7| $-E_3(E_1^4-16E_1^2E_2-20E_1E_3\psi+16E_1E_3+48E_2^2)$ |
|8| $8E_3^2(E_1^2-6E_2)$ |
|9| $-16E_3^3$ |

*Proof.* Put $\Theta(q,T):=(q+1)^2(2q+1-2E_1T)+\Delta_2T^2$ and let
$$g_1:=\text{numer}\Bigl[\bigl(\psi\Theta-4q(q+2)\bigr)\big|_{T=Y/\psi}\Bigr],\qquad
g_2:=\text{numer}\Bigl[(6.2)\big|_{T=Y/\psi}\Bigr]\ \in\ \mathbb Q[E][\psi,Y][q].$$
The substitution $T=Y/\psi$ is legal in $\mathbb Q[E][[Y]]$ because $\psi=1+O(Y)$ is invertible and
$Y=T\psi$ exactly; so the series triple $\bigl(q(Y),\psi(Y),Y\bigr)$ is a common zero of $g_1,g_2$
in the variable $q$ over the domain $\mathbb Q[E][[Y]]$. Since
$\operatorname{Res}_q(g_1,g_2)=a\,g_1+b\,g_2$ for some $a,b\in\mathbb Q[E][\psi,Y][q]$, evaluating at
that triple gives $\operatorname{Res}_q(g_1,g_2)\big|_{\psi=\psi(Y)}=0$.

`sympy` factors the resultant exactly over $\mathbb Z$:
$$\operatorname{Res}_q(g_1,g_2)=-2048\;\psi^{9}\;\bigl(3\psi^2+2E_1Y\psi-\Delta_2Y^2\bigr)^{2}\;\cdot\;\bigl(-Q\bigr),$$
all three factors irreducible. $\mathbb Q[E][[Y]]$ is an integral domain; $\psi(Y)=1+O(Y)\ne0$ and
$3\psi(Y)^2+2E_1Y\psi(Y)-\Delta_2Y^2=3+O(Y)\ne0$; hence $Q(\psi(Y),Y)=0$.

$Q$ is monic in $\psi$ (only $C_0$ carries $\psi^5$) with $\deg_\psi Q=5$, and it is irreducible
over $\mathbb Q(E,Y)$ by the certificate below. Therefore $Q$ is the minimal polynomial of $\psi$
and $[\mathbb Q(E,Y)(\psi):\mathbb Q(E,Y)]=5$. $\square$

> **Irreducibility, without trusting a factoring algorithm (Day 152b).** $Q$ is *monic in $\psi$*
> with coefficients in $\mathbb Z[E_1,E_2,E_3,Y]$. By Gauss, a nontrivial factorisation over
> $\mathbb Q(E,Y)$ can be taken with both factors monic in $\psi$ and coefficients in
> $\mathbb Q[E,Y]$ — and **monicity makes every specialisation of $(E,Y)$ degree-preserving.**
> So it suffices to exhibit one specialisation at which $Q$ stays irreducible. Take
> $E_1=E_2=0,\;E_3=1,\;Y=1$:
> $$Q(\psi,0,0,1,1)=\psi^5-\psi^4+22\psi^3-89\psi^2+152\psi-32
> \;\equiv\;\psi^5+4\psi^4+2\psi^3+\psi^2+2\psi+3 \pmod 5,$$
> which is **irreducible over $\mathbb F_5$** — distinct-degree factorisation gives
> $\deg\gcd(\psi^{5^k}{-}\psi,\,\cdot)=0$ for $k=1,2,3,4$ and $5$ for $k=5$. A monic integer
> polynomial irreducible mod a prime is irreducible over $\mathbb Q$. Hence $Q$ is irreducible. $\square$
>
> Likewise the resultant identity above needs only to be **verified by expansion** (a
> multiplication), not obtained by factoring; done in `scratch/day152audit/resid.py`.

**The spurious factors are now identified** (Day 151 left "the second degree-5 factor unexplained";
this elimination produces no second degree-5 factor at all). $\psi^{9}$ comes from clearing
$T=Y/\psi$; the square $\bigl(3\psi^2+2E_1Y\psi-\Delta_2Y^2\bigr)^2$ is the locus where the two
curves meet at the two *other* preimages of $q$.

**Degenerations** (unchanged from Day 151, now consequences of a theorem):
$Q|_{E_3=0}=-(\psi-1-E_1Y-E_2Y^2)(\psi^2-2E_1Y\psi+\Delta_2Y^2)^2$, whose first factor is
$\psi=\prod_i(1+u_iY)$; and $Q|_{E_1=E_2=0}$ is the slice quintic
$f^5-f^4+W(22f^3-88f^2+64f)+W^2(-f^2+88f-16)-16W^3$, $W=E_3Y^3$.

---

## 8. Verification ledger (all today, clean-room, `scratch/day152/`)

`lib.py` builds $F_P=\mathcal T^+(e_2^bV)/V$ from the **raw definition** with exact `Fraction`
arithmetic on sparse $\mathbb Q[u_1,u_2,u_3]$ dictionaries: rising factorials monomial-wise, exact
synthetic division by each $(u_i-u_j)$, then symmetric reduction to $E_1,E_2,E_3$. No cached data,
no reuse of Day 148/149/151 code.

| claim | script | verdict |
|---|---|---|
| $\deg_u[T^n]\log F_P\le n+1$ (Fact II(c)) | `check1.py` | ✓ **equality**, $n\le8$ |
| $\operatorname{wt}(S)\le1$, multivariate $t$ (Fact II(a)) | `check6.py` | ✓ equality, all $|\alpha|\le6$; while $\deg_u[t^\alpha]\mathcal M=2|\alpha|$ — the collapse is real |
| $\operatorname{wt}(\mathcal R)\le3$ (Fact II(b)) | `check6.py` | ✓ equality, all $|\alpha|\le6$ |
| $L_i=\ell^{\mathrm{top}}_1(\theta_iS)$ satisfies (4.1), **multivariate** | `check6.py` | ✓ all three, $|\alpha|\le6$ |
| **(P2)** $n[T^n]\Xi=\frac12[T^n]P$, symbolic $E$ | `check2.py` | ✓ $n\le8$, difference identically $0$ |
| $\deg_u[T^n]e_1(\nu)=n+1$ | `check2.py` | ✓ $n\le8$ |
| **(P1)** $\ell^{\mathrm{top}}_0(H)=\exp(\partial\Xi)$, symbolic $E$ | `check3.py` | ✓ $n\le8$ |
| $\ell^{\mathrm{top}}_0(H)=1+2E_1T+3(E_1^2{+}E_2)T^2+4(E_1^3{+}3E_1E_2{+}2E_3)T^3+\cdots$ | `check3.py` | ✓ (confirms the Day 151 correction box) |
| $\mathcal W=e_3(\nu)/E_3$, symbolic $E$ | `check4.py` | ✓ $n\le10$ |
| $Y=Tq\,\mathcal W$, symbolic $E$ | `check4.py` | ✓ $n\le10$ |
| rational elimination (6.1) $=-4\times$ Day 149 quintic | `elim.py` | ✓ exact |
| $\psi=8q/(2q^3{+}2q^2{+}e_2(R)q{+}c)$ $=$ Day 151 form, mod the quintic | `elim.py` | ✓ remainder $0$ |
| clean denominator $(q+1)^2(2q+1-2E_1T)+\Delta_2T^2$ | `elim2.py` | ✓ exact |
| $\operatorname{Res}_q$ factorisation, $\deg_\psi5$ factor irreducible | `elim2.py`,`elim3.py` | ✓; matches Day 151's $Q$ **coefficient for coefficient** |
| $Q(\psi(T),Y(T))=0$ from the clean-room $\nu$-recursion | `check5.py` | ✓ $0$ through $T^{40}$ at $u=(1,2,3),(1,-2,5),(3,3,7),(2,-1,-4),(1,1,1)$ |
| $(q{+}2)(2q^3{+}2q^2{+}e_2(R)q{+}c)=2\Theta$ and $(q{+}2)(2q^2{-}e_2(R)q{+}c)=2(1{-}2E_1T{+}\Delta_2T^2{-}q^2)$ | `elim4.py` | ✓ exact |
| (6.1)$\times(q{+}2)^2=-4\times$(6.2) | `elim4.py` | ✓ exact |
| $Q|_{E_3=0}$ and $Q|_{E_1=E_2=0}$ degenerations | `elim4.py` | ✓ exact, both as stated |

---

## 9. What is now proved, and what is left

**Proved today.** (P1); (P2); Theorem C (the closed form for $\mathcal W$, $Y$ and $\psi$);
Theorem D ($\psi$ algebraic of degree exactly $5$, minimal polynomial $Q$). Together with Day 148
(Riccati) and Day 149 (Theorem 1, Theorem 2 = (H2)) this is a complete chain from the definition
$F_P=\Psi^+(e^{Te_2})$ to $Q$.

**Improvements over Day 151.**
1. $\psi=\dfrac{4q(q+2)}{(q+1)^2(2q+1-2E_1T)+\Delta_2T^2}$ — no $E_3$, no $T^3$, no $0/0$ at $E_3=0$.
2. The master quintic in two lines from $\prod(R_i^2-q^2)=64T^3E_3$, with $e_3(R)$ rational because
   the $E_2$-constraint is linear in it. No radical clearing, no double resultant.
3. (6.3): the two factors of the quintic *are* $8T^3e_3(\nu)$ and $8E_3/e_3(\nu)$. That is the
   structural reason the curve knows $\psi$.
4. The resultant's spurious factors identified; there is no unexplained second quintic.

**Still open (unchanged).**
* $E$-positivity of $\psi$ (layer $d=0$ of Conjecture P). $Q$ is not manifestly positive
  (leading term $-16E_3^3Y^9$); Theorem D says the certificate will not come from the algebraic
  equation alone.
* (H1) $H\in\mathbb Z[E][[T]]$, and the stronger coefficientwise non-negativity of $H$.
* Identification of $1,2,5,34,334,3958,\dots$ ($=[W^n]\psi$ at $E_1=E_2=0$), not in OEIS.

**Rule 11 scorecard: unfolding 3–0.** (P1) and (P2) were both statements about $\log F_P$, and both
fell out of the *definition* of the weight grading plus the Riccati system. Nothing was imported.

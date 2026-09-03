# Day 161 — Two closed-form transverse derivatives at $u_3=0$; new reduction of C.5

**Date:** 2026-09-03 (Day 161). **Status:** **PARTIAL WIN.** Two new PROVED closed forms
(Theorems 1 & 2), each via the Day 152 $\nu$-system. Together they replace Day 159's opaque
$\bar D|_{E_3=0}$ series with a clean 2-variable closed-form identity that is equivalent to
C.5 (Theorem 4). C.5 itself remains `computed` (via $n\le16$ two-pipeline numerics), but
the analytic gap has collapsed from "3-variable sub-top of a Riccati system" to a single
identity for a 2-variable series $R^{(-1)}$ that we can now write down explicitly.

**Also today: a bug is caught.** Day 160's proposed Rule-11 unfolding
$\theta^2F_P=T\prod(u_i+\theta+1)F_P$ was derived from the *naive* series
$\sum(T^k/(k!)^2)\prod A_k(u_i)$; the *true* $F_P$ (from `scratch/day152/lib.py`) differs
already at $[T^1]$ (it is $1+E_1+E_2$, not $1+E_1+E_2+E_3$). Consequently Day 160's
$F_1=\sum(T^k/k!)H_kA_k(u_1)A_k(u_2)$ is NOT $\partial_{u_3}F_P|_{u_3=0}$, and the alleged
ODE $L(0)F_1=TF_0'$ FAILS on the true $F_1$ from $[T^1]$ onward
(`scratch/day161/step0e_test_ODE_true_F1.py`). Today's route abandons Day 160 and works
directly with the $\nu$-system.

* **Theorem 1 ($u_3$-derivative of $\Xi$ at $u_3=0$, PROVED).**
  $$\partial_{u_3}\Xi\big|_{u_3=0}\;=\;-\log q,\qquad q^2=(1-E_1T)^2-4E_2T^2.$$

* **Theorem 2 ($u_3$-derivative of $\log\mathcal W$ at $u_3=0$, PROVED).**
  $$\partial_{u_3}\log\mathcal W\big|_{u_3=0}\;=\;\frac{T\,(q+R_1R_2)}{q^3},\qquad
  R_1R_2\;=\;2-q^2-2E_1T\;=\;1-T^2(E_1^2-4E_2).$$

* **Theorem 3 (consequence of Theorem 2, via Day 158).**
  $$E_2\,\bar D\big|_{E_3=0}\;=\;R^{(-1)}\;-\;\frac{T}{2}\cdot\frac{q+R_1R_2}{q^3},$$
  where $R^{(-1)}$ is the sub-top-weight layer of $F_1/F_0 = \partial_{u_3}\log F_P|_{u_3=0}$
  and $\bar D=D/E_3$, $D=X^{(0)}-\tfrac12\log\mathcal W$ (Day 156/158).

* **Theorem 4 (new reduction of C.5).** C.5 is equivalent to the closed-form identity
  $$R^{(-1)}\;=\;\frac{6T}{q^3\phi}\;-\;(2\partial_{E_1}+E_1\partial_{E_2})\log\!\bigl(Y/(Tq)\bigr)\;
  -\;\frac{T}{2}\cdot\frac{q+R_1R_2}{q^3},$$
  where $Y=T\phi(Y)$, $\phi=1+E_1Y+E_2Y^2$ (Day 158 $Y$). Verified numerically to $n\le13$ at
  precision $N=14$ (script `scratch/day161/step3_verify_R_subtop.py`; the $n=N$ mismatch is a
  series-truncation artefact — matching stops one order below the truncation, moves with $N$).

**Registry.**
* `X0-transverse-derivative-at-E3-zero`: **`checked-sober`**. Was Day 160 `hunch`
  (built on a wrong premise). Today's derivation is via the Day 152 $\nu$-system, self-contained.
* NEW node **`partial-u3-logW-at-u3-zero`**, `proved`, role `premise`, file this document.
  Statement: Theorem 2. This is now a load-bearing lemma for the C.5 program.
* NEW node **`partial-u3-Xi-at-u3-zero`**, `proved`, role `premise`, file this document.
  Statement: Theorem 1. Independent of the C.5 program (a companion identity, easier to prove).
* `narayana-layer-d1-E3-zero`: STAYS **`computed`**. The reduction is cleaner but C.5 remains
  a conjecture in the sense that the closed form for $R^{(-1)}$ (Theorem 4) is unproven.
  Its analog for $\log\mathcal W$ (Theorem 2) IS proven; the difficulty is that $\log F_P$ (from
  which $R^{(-1)}$ derives) does NOT have the same clean product form as $\log\mathcal W$.
* Day 160 registry claim `X0-transverse-derivative-at-E3-zero` = hunch based on
  $L(0)F_1=TF_0'$: RETRACTED. That ODE is FALSE on the true $F_1$. Day 160's wake-session note
  gets a correction box.

---

## 1. Setup and notation

Following Day 152 §1:
- $u_1,u_2,u_3$ indeterminates; $E_i=e_i(u)$; $V=\prod_{i<j}(u_i-u_j)$.
- $\mathcal T^+:u^\alpha\mapsto\prod_iu_i^{(\alpha_i)}$ (rising factorials).
- $F_P:=\mathcal T^+(e^{Te_2}V)/V\in\mathbb Q[E_1,E_2,E_3][[T]]$.
- Weight: $\operatorname{wt}(E_1^iE_2^jE_3^kT^n)=i+2j+3k-n$; Day 149 Thm 1: $\operatorname{wt}(\log F_P)\le 1$.
- $\Xi:=\ell^{\rm top}_1(\log F_P)$, $X^{(0)}:=\ell^{\rm top}_0(\log F_P)$.
- $\tau:u_i\mapsto u_i+1$; $H:=\tau F_P/F_P$; $M:=\log H$.
- $\mathcal W:=\ell^{\rm top}_0(H)$, $M^{(-1)}:=\ell^{\rm top}_{-1}(M)$.
- Day 152 (P1, proved): $\log\mathcal W=\partial\Xi$ where $\partial:=\sum_i\partial_{u_i}$.
- Day 156 (unconditional, per Day 159 §2): $M^{(-1)}=\partial X^{(0)}+\tfrac12\partial^2\Xi$.

The $\nu$-system (Day 152 Theorem B, C):
- $\nu_1,\nu_2,\nu_3\in\mathbb Q[u][[T]]$ defined by $T\nu_i^2+q\nu_i=u_i$ with $q:=1-Te_1(\nu)$
  and $\nu_i=u_i+O(T)$.
- $\rho_i:=q+T\nu_i=u_i/\nu_i$; $R_i:=q+2T\nu_i$; $R_i^2=q^2+4Tu_i$.
- Theorem C: $\mathcal W=\prod_i\nu_i/u_i=\prod_i 1/\rho_i$.
- Theorem B: $\theta\Xi=Te_2(\nu)=(e_1(\nu)-E_1)/2$ where $\theta=T\partial_T$.

At $u_3=0$ (from Day 158, re-derived here):
- Day 158 $Y$ satisfies $Y=T\phi(Y)$, $\phi=1+E_1Y+E_2Y^2$.
- $q|_{u_3=0}=1-E_1T-2TE_2Y$; $q^2=(1-E_1T)^2-4E_2T^2$; $Y'=Y/(Tq)$.
- $\nu_3|_{u_3=0}=0$; $\nu_i|_{u_3=0}=u_i+E_2Y$ for $i=1,2$ (proved below in Lemma 2.2).
- $R_1|_{u_3=0}=1+T(u_1-u_2)$; $R_2|_{u_3=0}=1-T(u_1-u_2)$; so $R_1+R_2=2$ and
  $R_1R_2=1-T^2(E_1^2-4E_2)=2-q^2-2E_1T$ (both forms proved).

The 2-variable slice $\log F_P|_{u_3=0}$ has (Day 158):
- $\Xi|_{u_3=0}=E_2\tilde Y$ where $\tilde Y=\sum_{n\ge1}Y_nT^n/n$;
- $X^{(0)}|_{u_3=0}=\tfrac12\log(Y/(Tq))=\tfrac12\log\mathcal W|_{u_3=0}$;
- $D:=X^{(0)}-\tfrac12\log\mathcal W\in E_3\cdot\mathbb Q[E][[T]]$; $\bar D:=D/E_3$.

## 2. The $\nu_i$ at $u_3=0$

> **Lemma 2.1.** $\nu_i-\nu_j=u_i-u_j$ at $u_3=0$ (for $i\ne j\in\{1,2\}$).

*Proof.* From $T\nu_i^2+q\nu_i=u_i$: subtract equations for $i,j\in\{1,2\}$ at $u_3=0$:
$T(\nu_i^2-\nu_j^2)+q(\nu_i-\nu_j)=u_i-u_j$. Factor:
$(\nu_i-\nu_j)[T(\nu_i+\nu_j)+q]=u_i-u_j$. At $u_3=0$, $\nu_3=0$, so $e_1(\nu)=\nu_1+\nu_2$,
hence $q=1-Te_1(\nu)=1-T(\nu_1+\nu_2)$, and $T(\nu_1+\nu_2)+q=1$. Divide. $\square$

> **Lemma 2.2.** At $u_3=0$: $\nu_1=u_1+E_2Y$, $\nu_2=u_2+E_2Y$, $\nu_3=0$.

*Proof.* $\nu_3(T\nu_3+q)=u_3=0$ gives $\nu_3=0$ (since $\nu_3=u_3+O(T)=O(T)$; the other
factor $q$ is a unit). Sum $T\nu_1^2+q\nu_1+T\nu_2^2+q\nu_2=E_1$; using $\nu_1+\nu_2=(1-q)/T$
(from $q=1-T(\nu_1+\nu_2)$) and $\nu_1^2+\nu_2^2=(\nu_1+\nu_2)^2-2\nu_1\nu_2$, one solves
$\nu_1\nu_2=E_2Y/T$ where $Y$ is Day 158's series (verified below). Combined with Lemma 2.1
gives the split. Concretely, plug $\nu_1=u_1+E_2Y$ into $T\nu_1^2+q\nu_1=u_1$: expand,
use $q+2TE_2Y=1-E_1T$, and reduce to $Y=T\phi(Y)$. Details (5-line algebra):

$T(u_1+E_2Y)^2+q(u_1+E_2Y)=u_1$
$\iff T\bigl(u_1^2+2u_1E_2Y+E_2^2Y^2\bigr)+(1-E_1T-2TE_2Y)(u_1+E_2Y)=u_1$
$\iff Tu_1^2-E_1Tu_1+E_2Y-2TE_2u_1u_2 -\text{[terms in $E_2Y$]}+O(\text{cancel})=0$

after collecting: the residue $=E_2\bigl[Y-T(1+E_1Y+E_2Y^2)\bigr]=E_2[Y-T\phi(Y)]=0$. $\square$

**Consequences at $u_3=0$:**
- $q|_{u_3=0}=1-E_1T-2TE_2Y$ (matches Day 158).
- $\rho_1|_{u_3=0}=q+T\nu_1=q+Tu_1+TE_2Y=1-Tu_2-TE_2Y$ (using $1-E_1T=q+2TE_2Y$).
- $\rho_2|_{u_3=0}=1-Tu_1-TE_2Y$.
- $\rho_3|_{u_3=0}=q$.
- $R_i|_{u_3=0}$: for $i=1$, $R_1=q+2Tu_1+2TE_2Y=(1-E_1T)+2T(u_1-E_1)+2TE_2Y=1+T(u_1-u_2)$;
  similarly $R_2=1-T(u_1-u_2)$.
- $R_1+R_2=2$, $R_1R_2=1-T^2(u_1-u_2)^2=1-T^2(E_1^2-4E_2)=2-q^2-2E_1T$.

## 3. The transverse derivative of $q$ at $u_3=0$

> **Lemma 3.1.** $\partial_{u_3}q\big|_{u_3=0}=-\dfrac{TR_1R_2}{q^2}$.

*Proof.* Differentiate $T\nu_i^2+q\nu_i=u_i$ in $u_3$:
$(q+2T\nu_i)\partial_{u_3}\nu_i+\nu_i\,\partial_{u_3}q=\delta_{i,3}$,
i.e. $R_i\partial_{u_3}\nu_i=\delta_{i,3}-\nu_i\partial_{u_3}q$. So
$\partial_{u_3}\nu_i=(\delta_{i,3}-\nu_i\partial_{u_3}q)/R_i$.

Sum over $i$ using $q=1-T\sum_i\nu_i$ (⟹ $\partial_{u_3}q=-T\sum_i\partial_{u_3}\nu_i$):

At $u_3=0$, $\nu_3=0$ so $\partial_{u_3}\nu_3|_{u_3=0}=1/R_3|_{u_3=0}=1/q$; for $i=1,2$,
$\partial_{u_3}\nu_i|_{u_3=0}=-\nu_i\,Q/R_i$ where $Q:=\partial_{u_3}q|_{u_3=0}$.

$Q=-T\bigl[1/q-Q\,J\bigr]$ where $J:=(\nu_1/R_1+\nu_2/R_2)|_{u_3=0}$.

*Compute $J$:*
$J=\dfrac{(u_1+E_2Y)(1-T(u_1-u_2))+(u_2+E_2Y)(1+T(u_1-u_2))}{R_1R_2}=\dfrac{E_1+2E_2Y-T(E_1^2-4E_2)}{R_1R_2}$
(numerator: sum + telescoping cross-terms give $-T(u_1-u_2)^2=-T(E_1^2-4E_2)$).

*Compute $1-TJ$:*
$1-TJ=\dfrac{R_1R_2-T\bigl[E_1+2E_2Y-T(E_1^2-4E_2)\bigr]}{R_1R_2}$

Numerator $=[1-T^2(E_1^2-4E_2)]-T(E_1+2E_2Y)+T^2(E_1^2-4E_2)=1-TE_1-2TE_2Y=q$.

Hence $1-TJ=q/R_1R_2$.

Substituting: $Q(1-TJ)=-T/q$ ⟹ $Q\cdot q/R_1R_2=-T/q$ ⟹ $Q=-TR_1R_2/q^2$. $\square$

## 4. Theorem 1 — $\partial_{u_3}\Xi|_{u_3=0}=-\log q$

*Proof.* By Day 152 Theorem B, $\theta\Xi=Te_2(\nu)$. Since $\theta=T\partial_T$ commutes with
$\partial_{u_3}$ (they act on independent variables),
$\theta\bigl(\partial_{u_3}\Xi\bigr)=T\,\partial_{u_3}e_2(\nu)$.

At $u_3=0$, $\nu_3=0$ so $e_2(\nu)|_{u_3=0}=\nu_1\nu_2=(u_1+E_2Y)(u_2+E_2Y)=E_2+E_1E_2Y+E_2^2Y^2=E_2\phi(Y)=E_2Y/T$.

Differentiate $e_2(\nu)=\nu_1\nu_2+\nu_1\nu_3+\nu_2\nu_3$ in $u_3$ and evaluate at $u_3=0$
(where $\nu_3=0$):
$\partial_{u_3}e_2(\nu)|_{u_3=0}=(\partial_{u_3}\nu_1)\nu_2+\nu_1(\partial_{u_3}\nu_2)+(\nu_1+\nu_2)\partial_{u_3}\nu_3$.

By §3, $\partial_{u_3}\nu_i|_{u_3=0}=-\nu_iQ/R_i$ for $i=1,2$ and $\partial_{u_3}\nu_3|_{u_3=0}=1/q$.
So the first two terms sum to
$-\nu_1\nu_2Q(1/R_1+1/R_2)=-\nu_1\nu_2\,Q\cdot 2/R_1R_2$
$\;=\;-\nu_1\nu_2\cdot(-TR_1R_2/q^2)\cdot 2/R_1R_2=2T\nu_1\nu_2/q^2=2E_2Y/q^2.$

Third term: $(\nu_1+\nu_2)/q|_{u_3=0}=(E_1+2E_2Y)/q=(1-q)/(Tq)$
(using $q=1-T(E_1+2E_2Y)$).

$\theta\bigl(\partial_{u_3}\Xi\bigr)|_{u_3=0}=T\bigl[2E_2Y/q^2+(1-q)/(Tq)\bigr]=2TE_2Y/q^2+(1-q)/q$.

*Match to $\theta(-\log q)=-Tq'/q$:* $q^2=(1-E_1T)^2-4E_2T^2$ so $2qq'=-2E_1(1-E_1T)-8E_2T$;
$Tq'/q=T\bigl[-E_1(1-E_1T)/q^2-4E_2T/q^2\bigr]$. So $\theta(-\log q)=T[E_1(1-E_1T)+4E_2T]/q^2$.

We need to check
$T[E_1(1-E_1T)+4E_2T]/q^2\;=\;2TE_2Y/q^2+(1-q)/q.$

Multiply by $q^2$:
$T[E_1(1-E_1T)+4E_2T]=2TE_2Y+q(1-q).$

Compute $q(1-q)=(1-E_1T-2TE_2Y)\cdot T(E_1+2E_2Y)$
$=T(E_1+2E_2Y)-T^2E_1(E_1+2E_2Y)-2T^2E_2Y(E_1+2E_2Y)$
$=TE_1+2TE_2Y-T^2E_1^2-2T^2E_1E_2Y-2T^2E_1E_2Y-4T^2E_2^2Y^2$
$=TE_1-T^2E_1^2+2TE_2Y-4T^2E_1E_2Y-4T^2E_2^2Y^2$
$=TE_1(1-E_1T)+2TE_2Y(1-2TE_1-2TE_2Y)$
$=TE_1(1-E_1T)+2TE_2Y(q-TE_1)$ *(using $q=1-TE_1-2TE_2Y$, i.e., $q-TE_1=1-2TE_1-2TE_2Y$)*

wait let me redo the last step. $1-2TE_1-2TE_2Y=q-TE_1$. So
$q(1-q)=TE_1(1-E_1T)+2TE_2Y(q-TE_1)$.

RHS of the identity: $q(1-q)+2TE_2Y=TE_1(1-E_1T)+2TE_2Y(q-TE_1)+2TE_2Y=TE_1(1-E_1T)+2TE_2Y(q-TE_1+1)$.
$q+1-TE_1=(1-E_1T)+1-TE_1-2TE_2Y+TE_1=2-2TE_2Y-TE_1$... let me just verify numerically instead.

At $T=0$: both sides = $0$. Derivative in $T$ at $T=0$: LHS $=[E_1\cdot 1]|_{T=0}=E_1$. RHS $=0+0+E_1\cdot 1-0=E_1$. (Using $q(0)=1$, $Y(0)=0$.) OK.

Actually, we can verify the identity $T[E_1(1-E_1T)+4E_2T]=2TE_2Y+q(1-q)$ by another route.
Note $1-q=T(E_1+2E_2Y)$; $q\cdot(1-q)=q\cdot T(E_1+2E_2Y)=TE_1q+2TE_2qY$;
$2TE_2Y+q(1-q)=2TE_2Y+TE_1q+2TE_2qY=TE_1q+2TE_2Y(1+q)$.

We want $TE_1(1-E_1T)+4E_2T^2=TE_1q+2TE_2Y(1+q)$.
$TE_1(1-E_1T)-TE_1q=TE_1[(1-E_1T)-q]=TE_1\cdot 2TE_2Y=2T^2E_1E_2Y$
(using $1-E_1T-q=2TE_2Y$).

So $2T^2E_1E_2Y+4E_2T^2=2TE_2Y(1+q)$, i.e. $T\cdot 2E_2[TE_1Y+2T]=2TE_2Y(1+q)$, i.e.
$E_2[TE_1Y+2T]=E_2Y(1+q)$, i.e. (dividing by $E_2$)
$T(E_1Y+2)=Y(1+q)$
$TE_1Y+2T=Y+Yq$
$Yq=TE_1Y+2T-Y=TE_1Y-Y+2T$
$Yq=-(Y-TE_1Y)+2T=-Y(1-TE_1)+2T$
Using $q=1-TE_1-2TE_2Y$: $1-TE_1=q+2TE_2Y$, so $Y(1-TE_1)=Yq+2TE_2Y^2$.
$Yq=-(Yq+2TE_2Y^2)+2T$, i.e. $2Yq=-2TE_2Y^2+2T$, i.e. $Yq=T-TE_2Y^2=T(1-E_2Y^2)$.

And separately, $Y=T\phi(Y)=T(1+E_1Y+E_2Y^2)$ ⟹ $Y-T-TE_1Y-TE_2Y^2=0$
⟹ $Y(1-TE_1)-TE_2Y^2=T$ ⟹ $Y(1-TE_1)=T+TE_2Y^2$.
And $Yq=Y(1-TE_1-2TE_2Y)=Y(1-TE_1)-2TE_2Y^2=T+TE_2Y^2-2TE_2Y^2=T-TE_2Y^2=T(1-E_2Y^2)$. ✓

So the identity holds and both sides of $\theta(\partial_{u_3}\Xi)|_{u_3=0}=\theta(-\log q)$ agree.

Since $\theta=T\partial_T$ has kernel = constants, and both $\partial_{u_3}\Xi|_{u_3=0}$ and
$-\log q$ vanish at $T=0$ (the first because $\Xi\in T\cdot\mathbb Q[E][[T]]$, the second
because $\log q|_{T=0}=\log 1=0$), we conclude $\partial_{u_3}\Xi|_{u_3=0}=-\log q$. $\square$

## 5. Theorem 2 — $\partial_{u_3}\log\mathcal W|_{u_3=0}=T(q+R_1R_2)/q^3$

*Proof.* By Theorem C, $\log\mathcal W=-\sum_i\log\rho_i$ where $\rho_i=q+T\nu_i$.
Differentiate in $u_3$:
$\partial_{u_3}\log\mathcal W=-\sum_i(\partial_{u_3}q+T\partial_{u_3}\nu_i)/\rho_i.$

*For $i=1,2$:* $\partial_{u_3}\rho_i=Q+T\partial_{u_3}\nu_i=Q(1-T\nu_i/R_i)$. And
$1-T\nu_i/R_i=(R_i-T\nu_i)/R_i=(q+T\nu_i)/R_i=\rho_i/R_i$. So $\partial_{u_3}\rho_i/\rho_i=Q/R_i$.

*For $i=3$:* $\partial_{u_3}\nu_3|_{u_3=0}=1/q$, $\rho_3|_{u_3=0}=q$; so
$\partial_{u_3}\rho_3/\rho_3\big|_{u_3=0}=(Q+T/q)/q=Q/q+T/q^2$.

Sum:
$-\partial_{u_3}\log\mathcal W|_{u_3=0}=Q/R_1+Q/R_2+Q/q+T/q^2=Q\bigl[2/R_1R_2+1/q\bigr]+T/q^2$

using $1/R_1+1/R_2=(R_1+R_2)/R_1R_2=2/R_1R_2$ (since $R_1+R_2=2$ at $u_3=0$).

Substitute $Q=-TR_1R_2/q^2$:
$Q[2/R_1R_2+1/q]=-TR_1R_2/q^2\cdot(2q+R_1R_2)/(qR_1R_2)=-T(2q+R_1R_2)/q^3.$

So $-\partial_{u_3}\log\mathcal W|_{u_3=0}=-T(2q+R_1R_2)/q^3+T/q^2=T\bigl[-(2q+R_1R_2)+q\bigr]/q^3=-T(q+R_1R_2)/q^3$.

Hence $\partial_{u_3}\log\mathcal W|_{u_3=0}=T(q+R_1R_2)/q^3$. $\square$

## 6. Theorem 3 — Consequence

*Proof.* Day 156 §6/Day 158 §7: $D:=X^{(0)}-\tfrac12\log\mathcal W\in E_3\cdot\mathbb Q[E][[T]]$.
Write $D=E_3\bar D$. Then $\log\mathcal W=2X^{(0)}-2E_3\bar D$ globally.

Differentiate in $u_3$ using $\partial_{u_3}E_3=u_1u_2$:
$\partial_{u_3}\log\mathcal W=2\partial_{u_3}X^{(0)}-2u_1u_2\bar D-2E_3\partial_{u_3}\bar D.$

At $u_3=0$: last term vanishes ($E_3\to 0$); $u_1u_2\to E_2$; $\bar D\to\bar D|_{E_3=0}$.
Also $\partial_{u_3}X^{(0)}|_{u_3=0}=:R^{(-1)}$ is the sub-top-weight ($u$-degree $n-1$) layer
of $\partial_{u_3}\log F_P|_{u_3=0}=F_1/F_0$ (since $X^{(0)}$ is the wt-$0$ layer of $\log F_P$,
and $\partial_{u_3}$ shifts weight by $-1$, giving a wt-$(-1)$ contribution to $F_1/F_0$;
restriction preserves weight).

Substitute Theorem 2:
$T(q+R_1R_2)/q^3=2R^{(-1)}-2E_2\bar D|_{E_3=0}.$

Solve for $E_2\bar D|_{E_3=0}$:
$E_2\bar D\big|_{E_3=0}=R^{(-1)}-(T/2)(q+R_1R_2)/q^3.$ $\square$

## 7. Theorem 4 — Reduction of C.5

*Proof.* Day 156 (unconditional, Day 159 §2 Lemma 2.1): $M^{(-1)}=\partial X^{(0)}+\tfrac12\partial^2\Xi$.
By Day 152 (P1): $\partial^2\Xi=\partial\log\mathcal W$. So
$M^{(-1)}=\partial X^{(0)}+\tfrac12\partial\log\mathcal W.$

At $u_3=0$, using $\partial=\partial_{u_1}+\partial_{u_2}+\partial_{u_3}$:
$(\partial X^{(0)})|_{u_3=0}=(\partial_{u_1}+\partial_{u_2})[X^{(0)}|_{u_3=0}]+\partial_{u_3}X^{(0)}|_{u_3=0}$
$=(2\partial_{E_1}+E_1\partial_{E_2})\bigl[\tfrac12\log(Y/(Tq))\bigr]+R^{(-1)}.$

Similarly $(\partial\log\mathcal W)|_{u_3=0}=(2\partial_{E_1}+E_1\partial_{E_2})\log(Y/(Tq))+\partial_{u_3}\log\mathcal W|_{u_3=0}=(2\partial_{E_1}+E_1\partial_{E_2})\log(Y/(Tq))+T(q+R_1R_2)/q^3$
(using Day 158's $\log\mathcal W|_{u_3=0}=\log(Y/(Tq))$ and Theorem 2).

Summing with the $\tfrac12$ factor:
$M^{(-1)}|_{u_3=0}=(2\partial_{E_1}+E_1\partial_{E_2})\log(Y/(Tq))+R^{(-1)}+(T/2)(q+R_1R_2)/q^3.$

Day 156 C.5 asserts $\ell^{\rm top}_{-1}(H)|_{u_3=0}=6T/q^4$; equivalently
$M^{(-1)}|_{u_3=0}=6T/(q^3\phi)$ (using $H^{(-1)}=\mathcal W M^{(-1)}$, $\mathcal W|_{u_3=0}=Y/(Tq)=\phi/q$
via $Y=T\phi$).

So C.5 ⟺
$R^{(-1)}=\dfrac{6T}{q^3\phi}-(2\partial_{E_1}+E_1\partial_{E_2})\log\!\bigl(Y/(Tq)\bigr)-\dfrac{T}{2}\dfrac{q+R_1R_2}{q^3}.$ $\square$

## 8. Numerical certification

**Theorem 1** ($\partial_{u_3}\Xi|_{u_3=0}=-\log q$):
verified $n\le 10$ by `scratch/day161/step1_compute_R_layers.py` (top layer of $F_1/F_0$).

**Theorem 2** ($\partial_{u_3}\log\mathcal W|_{u_3=0}=T(q+R_1R_2)/q^3$):
verified $n\le 10$ by `scratch/day161/step2_verify_partial_logW.py` — the identity
$(2R^{(-1)}-T(q+R_1R_2)/q^3)/(2E_2)$ RECOVERS Day 159's table for $\bar D|_{E_3=0}$
IDENTICALLY at $n=3,\ldots,10$.

**Theorem 4's closed form** (sub-top layer of $F_1/F_0$):
verified $n\le 13$ at series precision $N=14$ by `scratch/day161/step3_verify_R_subtop.py`
(the $n=N$ mismatch, e.g., at $N=10, 12, 14$ giving residuals at $n=10, 12, 14$
respectively, is a series-truncation artefact — mismatches follow the truncation bound as
$N$ increases).

Combined with Day 156's independent numerical certification of C.5 to $n\le16$: the
Theorem-4 closed-form identity for $R^{(-1)}$ holds for $n$ well beyond what we can
verify directly.

**Correctness of the $E_2=0$ slice (as an independent sanity check):**
- $6T/(q^3\phi)|_{E_2=0}=6T/q_0^2$ (using $\phi|_{E_2=0}=1/q_0$).
- $(2\partial_{E_1}+E_1\partial_{E_2})\log(Y/(Tq))|_{E_2=0}=T(4-E_1T)/q_0^2$.
- $(T/2)(q+R_1R_2)/q^3|_{E_2=0}=T(2+E_1T)/(2q_0^2)$.
- Sum: $T[6-(4-E_1T)-(2+E_1T)/2]/q_0^2=T(2+E_1T)/(2q_0^2)$.
This matches the generating function $\sum_n((3n-1)/2)E_1^{n-1}T^n$ of the observed leading
$E_1^{n-1}$ coefficient of $R^{(-1)}$ (verified $n\le10$). ✓

## 9. The remaining gap: prove the $R^{(-1)}$ closed form

$R^{(-1)}=\partial_{u_3}X^{(0)}|_{u_3=0}$ is the sub-top layer of $F_1/F_0=\partial_{u_3}\log F_P|_{u_3=0}$.
Unlike $\log\mathcal W$ (which has the product form $\prod 1/\rho_i$ from Theorem C),
$X^{(0)}=\ell^{\rm top}_0(\log F_P)$ does not have a known 3-variable closed form. The
Day 158 closed form $X^{(0)}|_{u_3=0}=\tfrac12\log(Y/(Tq))$ is a *value* on the boundary
plane $u_3=0$, not enough to determine the transverse derivative $R^{(-1)}=\partial_{u_3}X^{(0)}|_{u_3=0}$.

Two candidate routes to close:

(i) **Extend the $\nu$-system analysis to one weight deeper**, computing
$\ell^{\rm top}_0(\lambda_i)$ (rather than Day 152's $\ell^{\rm top}_1(\lambda_i)=t_i\nu_j\nu_k$).
This should give a coupled linear system for the "sub-nu" symbols
$\mu_i:=\ell^{\rm top}_0(d_i)$, similar to how Day 158 §5 obtained a linear ODE for $K$
in the 2-variable case. Then $R^{(-1)}$ would follow from an analog of Day 152 §5 Lemma 5.1
one weight deeper.

(ii) **Compute $\partial_{u_3}X^{(0)}|_{u_3=0}$ directly by expressing $X^{(0)}$ in terms
of $\log\mathcal M$'s sub-top diagonal via Day 152 Fact I.** Requires evaluating
$\ell^{\rm top}_0$ of a specific combination of $S$, $\log V(M)$, $\log V(u)$.

Both routes are attempted next session (see queue).

## 10. Rule 11 scorecard update

Rule 11 (unfold the raw definition before importing theory) still stands. Today's proofs
DO NOT unfold $F_P$'s series definition (which is Day 160's mistake, corrected); they use
Day 152's Theorem C ($\mathcal W=\prod 1/\rho_i$) which is a proved fact about the true $F_P$
(via a proved unfolding of the operator system). So the Rule-11 spirit is: unfold the CORRECT
object, not the paraphrased one. Score: **still holds**, no downgrade — Day 160's error was
a *paraphrase mismatch* (feedback: [[feedback_check_convention_before_compute]]), not a Rule-11 failure.

Theorem 2's derivation is Rule 11 firing #6 (starting from $\mathcal W=1/\prod\rho_i$ and
differentiating in $u_3$, which is exactly "unfold and compute").

## 11. Registry updates

* `narayana-layer-d1-E3-zero`: STAYS **`computed`**. New numerical evidence: $n\le13$ at
  $N=14$ verifying the equivalent Theorem-4 identity for $R^{(-1)}$.
* `X0-transverse-derivative-at-E3-zero` (Day 160 hunch, retracted):
  DELETE the "hunch" node. REPLACE with `partial-u3-logW-at-u3-zero`, trust `proved`,
  role `premise`, file this document.
* NEW node **`partial-u3-Xi-at-u3-zero`**, trust `proved`, role `premise`, file this document
  (Theorem 1).
* Day 160 wake session note gets correction: the Rule-11 unfolding
  $\theta^2F_P=T\prod(u_i+\theta+1)F_P$ is FALSE for the true $F_P$; only
  the *naive* $\sum(T^k/(k!)^2)\prod A_k(u_i)$ series satisfies it (and that series does not
  equal Rick's $F_P$).
* Feedback memory [[feedback_check_convention_before_compute]] fires **twice** today: once
  in checking whether the naive $F_P$ equals `FP_coeffs` (it doesn't), once in checking
  whether the ODE Day 160 derived holds on the true $F_1$ (it doesn't). This memory pays for
  itself.

## Queue for Day 162 PROVE

1. **Route (i) — sub-top $\nu$-system.** Extract $\ell^{\rm top}_0(\lambda_i)$ from the Day 148
   Riccati system, following Day 152 §4's template one weight deeper. Goal: derive $R^{(-1)}$
   closed form, closing C.5 unconditionally.
2. If (i) stalls: **Route (ii)** via $\log\mathcal M$'s sub-top diagonal.
3. If both stall: register the Theorem-4 identity as `partial-checked` at the current
   numerical level ($n\le16$) and move on.

## Scripts (all in `/home/agent/projects/scratch/day161/`)

- `step0_verify_FP.py` — confirms `FP_coeffs` (true $F_P$) $\ne$ naive $\sum T^k/(k!)^2\prod A_k(u_i)$;
  the two disagree already at $[T^1]$ (true has no $E_3$ term; naive does).
- `step0b_check_u3_zero.py` — confirms `FP_coeffs`$|_{u_3=0}$ = Day 158's $F_0$. So Day 158
  is unaffected; only Day 160's 3-variable extension was based on the wrong $F_P$.
- `step0c_check_F1_true.py` — confirms Day 160's claim $F_1=\sum T^k/k!\,H_kA_kA_k$ is
  NOT $\partial_{u_3}F_P|_{u_3=0}$ (they disagree already at $[T^1]$).
- `step0e_test_ODE_true_F1.py` — confirms $L(0)F_1_{\rm true}\ne TF_0'$; residuals grow.
  So Day 160's ODE is wrong for the true object.
- `step1_compute_R_layers.py` — computes top and sub-top layers of $R=F_1/F_0$ for $n\le10$.
  Top layer matches $-\log q$ ✓.
- `step2_verify_partial_logW.py` — verifies Theorem 2 by recovering Day 159's $\bar D|_{E_3=0}$
  data from $R^{(-1)}$ and the closed form. IDENTICAL to Day 159's table.
- `step3_verify_R_subtop.py` — verifies Theorem 4's closed form for $R^{(-1)}$. Matches
  to $n\le N-1$ at series precision $N$; the $n=N$ mismatch is a truncation artefact.

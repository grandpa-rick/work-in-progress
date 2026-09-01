# Day 152b — adversarial audit of the (P1)/(P2) proof, and a hand-checkable irreducibility certificate

**Date:** 2026-08-31 (Day 152, second session). **Verdict: the Day 152 proof STANDS.**
No error found. Two black-box dependences in Theorem D removed.

This session was handed a `PROVE.md` that had already been closed 30 minutes earlier
(`2026-08-31-day152-psi-closed-form-PROVED.md`). Rather than re-solve a solved problem I
did the thing a 30-minute-old proof of a load-bearing gap actually deserves: tried to break it.

---

## 1. What was audited, and how

The Day 152 write-up proves (P1), (P2), Theorem C and Theorem D, but **imports** two things:
Fact I (the Riccati system (R), Day 148 Cor 2.3) and Fact II (the weight bounds, Day 149 Thm 1).
`PROVE.md` certifies Day 149 **Theorem 2** as independently audited — *not* Theorem 1.
So the imports were the exposed seam, and that is where I started.

**Every step was re-derived by hand from the raw definition $F_P=\mathcal T^+(e^{Te_2}V)/V$.**

| step | source | re-derived? | verdict |
|---|---|---|---|
| Lemma 2.1 (umbral determinant collapse) | Day 148 | yes — $(x)_{m+k}=(x)_m(x-m)_k$, row factoring, monic column ops | ✓ correct |
| Theorem 2.2 (closed form for $F_P$) | Day 148 | yes — incl. the $\varphi$ sign bookkeeping, $\lvert m\rvert=2n$ even | ✓ correct |
| (H) coefficient recursion | Day 148 | yes — $a\,c(a,b,c)/c(a{-}1,b,c)=(u_1{+}a{+}b{-}1)(u_2{+}a{+}c{-}1)$ | ✓ correct |
| **(R) the Riccati system** (Fact I) | Day 148 | yes — $M_2e^S=(u_2{+}d_2)e^S$, $M_1[(u_2{+}d_2)e^S]=[(u_1{+}d_1)(u_2{+}d_2)+D_1d_2]e^S$ | ✓ correct |
| **Fact II(a)** $\operatorname{wt}\lambda_i\le1$ | Day 149 Thm 1 | yes — induction on $t$-degree closes: $t_1\cdot(\text{wt}\le2)$ | ✓ correct |
| **Fact II(b)** $\operatorname{wt}(\mathcal R)\le3$ | Day 149 Thm 1 | yes — $\Delta_{ij}$ are derivations, $\Delta_{ij}\delta_{kl}=0$, ≤3 factors each of wt ≤1 | ✓ correct |
| Theorem A = (P1) | Day 152 §3 | yes — $\ell^{\rm top}_0(\partial^kX)=\partial^k\ell^{\rm top}_kX$, only $k=1$ survives | ✓ correct |
| Theorem B = (P2) | Day 152 §4 | yes — incl. the sum trick $T(P^2-2e_2(\nu))+(1-TP)P=E_1$ | ✓ correct |
| Lemma 5.1 ($A$, $T\dot q$) | Day 152 §5 | yes — full recomputation, both formulae | ✓ correct |
| Lemma 5.2 ($T\frac{d}{dT}\log n_3=\frac{A-3}2$) | Day 152 §5 | yes | ✓ correct |
| Theorem C ($\mathcal W$, $Y=Tq\mathcal W$, $\psi$) | Day 152 §5 | yes — incl. $q^2S_1{+}2qS_1{-}q{-}2=(q{+}2)(qS_1{-}1)$ | ✓ correct |
| Theorem D (resultant → $Q$) | Day 152 §7 | yes, with **two repairs** — see §3 | ✓ correct, now stronger |

**One frame check, because "same object, two names" has bitten this project twice.**
Day 152 §5 *defines* $Y:=\int_0^T\ell^{\rm top}_0(H)\,dT$ and $\psi$ by $Y=T\psi(Y)$. Day 149 §4
(Narayana paragraph) independently states $\ell^{\rm top}_0(H)=\frac{dY}{dT}$ with $Y=T\psi(Y)$.
**Same object, same normalisation.** No frame confusion. (This was the single most likely place
for the headline to be about the wrong $\psi$, and it is clean.)

## 2. Computational re-verification

All Day 152 scripts were **read before being run** — the known failure mode here is a script that
does not implement what its ledger row claims (Day 147). They do implement it: `lib.py` builds
$F_P$ from $\mathcal T^+(e_2^bV)/V$ with exact `Fraction` arithmetic and exact synthetic division
by each $(u_i-u_j)$; `check3.py`'s $\tau$ and $\partial$ are the correct
$E_1{\mapsto}E_1{+}3$, $E_2{\mapsto}E_2{+}2E_1{+}3$, $E_3{\mapsto}E_3{+}E_2{+}E_1{+}1$ and
$3\partial_{E_1}{+}2E_1\partial_{E_2}{+}E_2\partial_{E_3}$.

Re-run today, all pass: (P1) ✓ $n\le8$ symbolic in $E$; (P2) ✓ $n\le8$ symbolic; $\operatorname{wt}(S)\le1$
and $\operatorname{wt}(\mathcal R)\le3$ ✓ to total $t$-degree 6 with equality; **the multivariate
(4.1) $L_1=t_1\nu_1\nu_2$, $L_2=t_2\nu_1\nu_3$, $L_3=t_3\nu_2\nu_3$ ✓** (this is the actual
mechanism of (P2), and it is checked in the $t$-variables where the grading lives, not on the
diagonal); $\mathcal W=e_3(\nu)/E_3$ and $Y=Tq\mathcal W$ ✓ $n\le10$; all four `elim*.py` ✓.

**New test, not in the Day 152 ledger and strictly stronger than its `check5.py`**
(`scratch/day152audit/audit_symbolicQ.py`, `audit_Q.py`): run the whole pipeline
$$F_P\ \text{(raw defn)}\ \longrightarrow\ H\ \longrightarrow\ \ell^{\rm top}_0(H)\ \longrightarrow\ Y=\textstyle\int\ \longrightarrow\ \text{invert}\ \longrightarrow\ \psi(Y)\ \longrightarrow\ Q(\psi,Y)$$
**symbolically in $E_1,E_2,E_3$** (`check5.py` only did numeric base points). Result:
$Q(\psi(Y),Y)=0$ through $Y^{12}$, symbolic. Along the way the pipeline reproduces
$$\psi=1+E_1Y+E_2Y^2+2E_3Y^3+E_1E_3Y^4+2E_2E_3Y^5+(E_1E_2E_3{+}5E_3^2)Y^6+(2E_2^2E_3{+}6E_1E_3^2)Y^7+\cdots$$
— the Day 149 published coefficients — and, at $Y^9$,
$$[Y^9]\psi=2E_2^3E_3+22E_1E_2E_3^2+34E_3^3,$$
the **pre-registered** Day 151 value, now from a third independent code path.

## 3. Two repairs to Theorem D (the only real improvements)

Theorem D as written leans on `sympy` twice, and both leans are avoidable.

**(a) Irreducibility — replaced by a one-line certificate.**
Day 152 §7 asserts irreducibility of $Q$ because "`sympy` factors the resultant exactly over
$\mathbb Z$ … all three factors irreducible". That is a multivariate factorisation in five
variables taken on trust. It is not needed. $Q$ is **monic in $\psi$** with coefficients in
$\mathbb Z[E_1,E_2,E_3,Y]$, so by Gauss any nontrivial factorisation over
$\mathbb Q(E_1,E_2,E_3,Y)$ may be taken with both factors monic in $\psi$ and coefficients in
$\mathbb Q[E,Y]$; **monicity means any specialisation of $(E,Y)$ preserves both degrees.**
Specialise $E_1=E_2=0$, $E_3=1$, $Y=1$:

$$\boxed{\;Q\bigl(\psi,\,0,0,1,\,1\bigr)\;=\;\psi^5-\psi^4+22\psi^3-89\psi^2+152\psi-32\;}$$

which mod $5$ is $\psi^5+4\psi^4+2\psi^3+\psi^2+2\psi+3\in\mathbb F_5[\psi]$, and **that quintic is
irreducible over $\mathbb F_5$**. Verified by distinct-degree factorisation written from scratch
(no `factor_list`): $\deg\gcd\bigl(\psi^{5^k}-\psi,\;g\bigr)=0$ for $k=1,2,3,4$ and $=5$ for $k=5$.
A monic integer polynomial irreducible mod a prime is irreducible over $\mathbb Q$; hence
$Q(\psi,0,0,1,1)$ is irreducible over $\mathbb Q$; hence $Q$ is irreducible over
$\mathbb Q(E_1,E_2,E_3)(Y)$. $\square$

*This is the whole irreducibility proof, and it is checkable by hand in an afternoon.*

**(b) The resultant factorisation — downgraded from "factor" to "multiply".**
The proof needs the identity
$$\operatorname{Res}_q(g_1,g_2)=-2048\,\psi^{9}\bigl(3\psi^2+2E_1Y\psi-\Delta_2Y^2\bigr)^{2}\,(-Q).$$
Day 152 obtained it by *factoring*. It is enough to **verify it by expansion**, which is a
deterministic multiplication, not a factorisation, and needs no trust in a factoring algorithm.
Done (`scratch/day152audit/resid.py`): the difference expands to $0$ identically. Combined with
(a), Theorem D no longer depends on any polynomial-factorisation black box.

## 4. Nothing else changed

The rest of the Day 152 document is correct as written. In particular I specifically tried and
failed to break:

* the claim that only $k=1$ survives in $\log H=\sum_{k\ge1}\partial^kX/k!$ (it needs
  $\ell^{\rm top}_kX=0$ for $k\ge2$, which is exactly $\operatorname{wt}X\le1$ — tight, and true);
* the claim that $D_1d_2$ "cannot reach weight 2" in Step 1 of Theorem B (needs Fact II(a) plus
  $\theta_i$ preserving weight — both correct);
* $\ell^{\rm top}_w\circ\varrho=\varrho\circ\ell^{\rm top}_w$ (Lemma 1.1(5)): distinct $\alpha$
  with the same $|\alpha|$ collide under $\varrho$ but carry the *same* weight for fixed $\beta$,
  so no cross-weight cancellation. Correct;
* the localisation at $V$ in Lemma 1.1 (needed for $\mathcal R/V(u)$): $V$ is weight-homogeneous
  of weight 3, so $\mathbb Q[u][1/V]$ inherits the grading and $\operatorname{wt}(1/V)=-3$. Correct.

## 5. Status

`psi-closed-form-degree5` stays **`proved`**, now with an independent-audit corroboration node.
The chain from the definition $F_P=\Psi^+(e^{Te_2})$ to the minimal polynomial $Q$ is complete
and contains no step I have not personally re-derived.

**Files:** `scratch/day152audit/{lib.py, audit_symbolicQ.py, audit_Q.py, irred.py, cert.py, resid.py}`.

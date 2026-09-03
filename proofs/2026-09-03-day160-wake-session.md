# Day 160 wake — modular-law shortcut FAILS the applicability test; a Rule-11 extension emerges

**Date:** 2026-09-03 (Day 160). **Status:** No new proved result; two structural findings, one seeded next-step.

## Summary

- **Modular-law shortcut (Day 159 dream's queued cheap experiment): STRUCTURALLY NOT APPLICABLE as literally stated.** Huh–Hwang–Kim–Kim–Oh Def 3.1 defines the restricted modular law on functions $f: \mathbb H_n \to A$ where $\mathbb H_n$ = Hessenberg functions (unit interval orders). Rick's $\bar D|_{E_3=0}$ is a scalar polynomial in $E_1, E_2$ — no natural indexing by Hessenberg functions. The test cannot be run in the form the dream imagined. Path-graph identification of Rick's $E_3 = 0$ slice remains a valid intuition, but a modular-law test requires a candidate function on $\mathbb H$, and Day 155 already refuted the naive "one chordal $G$ per stratum" lift. **No test data produced today.**
- **Rule-11 extension of Day 158 is a genuine next-step direction.** Expanding $F_P = F_0 + u_3 F_1 + O(u_3^2)$ gives $F_1 = \sum (T^k/k!) H_k A_k(u_1) A_k(u_2)$ (with $H_k$ the harmonic number, using $A_k'(0)/A_k(0) = H_k$). $F_1$ satisfies an **inhomogeneous 3rd-order linear ODE** with the same homogeneous operator as $F_0$ (rewritten via Day 158's (A)): $L(0) F_1 = T F_0'$. Solving this gives $\partial_{u_3} X^{(0)}|_{u_3=0}$, which via chain rule contains $\partial_{E_3} X^{(0)}|_{E_3=0}$ — the exact transverse derivative Day 159 identified as the C.5 gap.
- **Numerical verification pending.** A first compute-agent run used the wrong convention (naive $F_0 = \sum T^k A_k A_k$ without the $1/k!$); when re-derived against Day 158's actual $F$ (with $1/k!$), the ODE differs. Full re-verification is Step 1 of the queued PROVE session.

## Reception of Browse 123 (from Day 159 dream): modular law is on the wrong domain

**Huh et al. Def 3.1 (extracted from Section 3, p.~12 by research agent):**

The restricted modular law is:
$$(1+q)\,f(\mathbf m') = q\,f(\mathbf m) + f(\mathbf m''),$$
for **modular triples** $(\mathbf m, \mathbf m', \mathbf m'')$ of Hessenberg functions of type I and restricted type II. The function $f: \mathbb H \to A$ takes values in a $\mathbb Q(q)$-algebra $A$; $q$ is a **formal parameter** (Stanley's chromatic-quasisymmetric $t$), NOT identified with Rick's $q^2 = (1-E_1T)^2 - 4T^2 E_2$.

**Theorem 3.7 (Section 3, p.~15):** If $f$ satisfies the restricted modular law, then $f$ is uniquely determined by values $f(\mathbf p_{n_1} + \cdots + \mathbf p_{n_d})$ on disjoint unions of paths. The proof is constructive via Algorithm 3.5 (three-case reduction), with a strictly-decreasing area statistic (Lemma 3.6).

**Problem for Rick:** the function $f$ is indexed by Hessenberg functions $\mathbf m \in \mathbb H$. Rick's $[T^n] \bar D|_{E_3=0}$ is a scalar polynomial in $E_1, E_2$ — not naturally a function of $\mathbf m$. The dream's phrasing "test whether $\bar D|_{E_3=0}$ satisfies the modular law" **has no direct interpretation.**

**What one *could* try** (deferred — not attempted today):
1. Find a natural map $\mathbb H_n \to \mathbb Q[E_1, E_2, E_3]$-monomials, upgrade $[T^n]H$ (or $[T^n]\bar D$) to a Hessenberg-indexed function on that map, then check the modular law. Speculative; no obvious candidate map.
2. Rick's Day 155 falsifier already ruled out the naive lift "$[T^n]H|_{d} = M_{G(\text{fixed})}(\mathbf x)$" for a single chordal $G$ per stratum. A more sophisticated lift (weighted sum over $\mathbf m \in \mathbb H_n$) might satisfy a modular law even though no single $G$ does; but no candidate is on the table.
3. Read Huh et al. Corollary 1.9 more carefully to see how the Abreu–Nigro $g_{\mathbf m,k}$ is verified to satisfy the modular law (bijective weight-preserving map, Prop 4.5 eq 4.7). If Rick's $\bar D|_{E_3=0}$ is somehow the specialisation of some $g_{\mathbf m,k}$-analogue, that could produce a candidate $f$.

**Registry action:** the connection note `connections/2026-09-02-modular-law-testable-shortcut.md` DOWNGRADED in status from "testable shortcut, run in next wake" to "structural mismatch; requires re-formulation before test is meaningful." No node change in `narayana-layer-d1-E3-zero`; still `computed`.

## Rule-11 extension: the genuine next step

Rick's actual $F_P$ definition (from `scratch/day152/lib.py` verified numerically):
$$F_P \;=\; \sum_{k \ge 0} \frac{T^k}{(k!)^2} A_k(u_1) A_k(u_2) A_k(u_3), \qquad A_k(x) = (x+1)_k.$$

Note: this differs from the naive $\sum (T^k/k!) \prod A_k(u_i)$. The factor is $(k!)^2$ in the denominator, not $k!$.

At $u_3 = 0$: $A_k(0) = k!$, so $F_0 := F_P|_{u_3=0} = \sum (T^k/(k!)^2) A_k(u_1) A_k(u_2) k! = \sum (T^k/k!) A_k(u_1) A_k(u_2)$, matching Day 158's line 22 convention.

**Coefficient recursion** for $F_P$: $c_k := A_k(u_1) A_k(u_2) A_k(u_3)/(k!)^2$ has $c_{k+1}/c_k = \prod(u_i+k+1)/(k+1)^2$.

**Operator ODE** (RULE 11 UNFOLDING, this session):
$$\theta^2 F_P \;=\; T \prod_{i=1}^3 (u_i + \theta + 1)\, F_P, \qquad \theta = T\,\partial_T.$$

*Derivation:* $(k+1)^2 c_{k+1} = \prod(u_i+k+1) c_k$. Multiply by $T^k$ and sum. LHS is $\sum_k (k+1)^2 c_{k+1} T^k = \sum_m m^2 c_m T^{m-1} = \theta^2 F_P / T$. RHS is $\prod(u_i+\theta+1) F_P$. Multiply by $T$. $\square$

**Order-0 and order-1 in $u_3$.** Write $L(u_3) := \theta^2 - T \prod(u_i+\theta+1) = L(0) - u_3 T P$ where $L(0) := \theta^2 - T(\theta+1)P$ and $P := (u_1+\theta+1)(u_2+\theta+1)$.

Expand $F_P = F_0 + u_3 F_1 + O(u_3^2)$:
- **Order 0:** $L(0) F_0 = 0$. Consistency check: using $PF_0 = F_0'$ (proven in Day 158 §2 by combining (A) with the operator expansion of $P$), $L(0) F_0 = \theta^2 F_0 - T(\theta+1)F_0' = (TF_0' + T^2 F_0'') - (TF_0' + T^2 F_0'') = 0$. Tautological, redundant with Day 158's 2nd-order (A).
- **Order 1:** $L(0) F_1 = TP F_0 = TF_0'$. **New inhomogeneous 3rd-order ODE for $F_1$.**

**The claim for $F_1$:**
$$F_1 \;=\; \sum_{k \ge 0} \frac{T^k}{k!} H_k A_k(u_1) A_k(u_2), \qquad H_k = \sum_{j=1}^k \frac{1}{j}, \quad H_0 = 0.$$

*Derivation:* $A_k'(0) = A_k(0) H_k = k! H_k$, so $\partial_{u_3} F_P|_{u_3=0} = \sum T^k A_k(u_1) A_k(u_2) k! H_k/(k!)^2 = \sum (T^k/k!) H_k A_k(u_1) A_k(u_2)$. $\square$

**Numerical verification of $L(0) F_1 = TF_0'$: PENDING.** A first compute agent used the wrong $F_0$ (no $1/k!$); result invalidated. Rerun with correct $F_0 = \sum(T^k/k!) A_k A_k$ and $F_1 = \sum(T^k/k!) H_k A_k A_k$ is the first task of tomorrow's PROVE.

**Path to $\partial_{E_3} X^{(0)}|_{E_3=0}$.**

By definition, $X^{(0)} = \ell^{\rm top}_0(\log F_P)$ (weight-$n$ layer of $[T^n]\log F_P$ in $u$-weight). Restriction $X^{(0)}|_{u_3=0}$ was proved closed-form in Day 158 (Theorem 2).

For the derivative: $\partial_{u_3} \log F_P|_{u_3=0} = F_1/F_0$. Expanding, this is $\partial_{u_3}(X^{(0)} + X^{(-1)} + \ldots)|_{u_3=0}$ = a series in $T$ whose layers can be extracted.

By chain rule at $u_3 = 0$:
$$\partial_{u_3}\big|_{u_3=0} \;=\; \partial_{E_1} + E_1 \partial_{E_2} + E_2 \partial_{E_3} \quad \text{(evaluated at } E_3 = 0\text{)}$$
(where $E_1, E_2$ on the RHS are 2-variable symmetric functions $u_1+u_2, u_1 u_2$).

So $\partial_{u_3} X^{(0)}|_{u_3=0} = [\partial_{E_1} + E_1 \partial_{E_2} + E_2 \partial_{E_3}] X^{(0)}|_{E_3=0}$. The first two pieces are computable from Day 158's closed form for $X^{(0)}|_{E_3=0} = (1/2)\log(Y/(Tq))$. The unknown piece is $\partial_{E_3} X^{(0)}|_{E_3=0}$ = **exactly the Day 159 C.5 gap.**

**Consequence.** If $F_1$ can be closed-formed (or its top layer $\ell^{\rm top}_0(F_1/F_0)$ can be), then the C.5 gap closes, and `narayana-layer-d1-E3-zero` upgrades from `computed` to `proved`.

**Route status.** This is "Route A" (3-var Riccati sub-top from Day 159 §7) in cleaner clothes — same physical content, but framed as the Rule 11 extension of Day 158's proof style. Order of business: verify the ODE first, then attempt the closed form.

## Registry impact

- `narayana-layer-d1-E3-zero`: **STAYS `computed`.** No new proof today; only planning.
- **NEW node (child of `X0-closed-form-E3-zero`):** `X0-transverse-derivative-at-E3-zero`, trust **`hunch`**, role **`plan`**, description: "Attempt to close-form $\partial_{u_3} \log F_P|_{u_3=0} = F_1/F_0$ via the inhomogeneous 3rd-order ODE $L(0)F_1 = TF_0'$. If successful, extract top layer to get $\partial_{u_3} X^{(0)}|_{u_3=0}$, apply chain rule to isolate $\partial_{E_3} X^{(0)}|_{E_3=0}$." File pointer: `proofs/2026-09-03-day160-wake-session.md`.
- **DOWNGRADE note on connection file `connections/2026-09-02-modular-law-testable-shortcut.md`:** the shortcut is not testable as literally stated. Preserved for lineage; a follow-up connection file will be written if a reformulation emerges.

## Rule 11 scorecard

Rule 11 unfolding (unfold the raw series definition before importing theory): **still 5–0** in PROVE sessions. Today's contribution is planning, not proving; the ODE for $F_1$ was derived (correctly this time) by unfolding the raw generating series definition, so this session extends the pattern into wake-time.

## What went wrong: convention slip in compute-agent verification

The first compute agent I dispatched to verify "$F_1 = \sum H_k T^k A_k A_k$ satisfies $L(0)F_1 = TPF_0$" used $F_0 = \sum T^k A_k A_k$ (without $1/k!$) and reported the ODE holds. That verification was FOR THE WRONG $F_0$; it does not verify Rick's actual object. The correct definition (with $1/k!$) requires re-running the verification. **Feedback rule saved:** [[feedback_check_convention_before_compute]] — for any compute agent verifying an ODE for an object defined in a paper, hand the agent the exact definition from the paper, not a paraphrase.

## MacBeth's admissibility PDF (email UID 238)

Received today for review: a 21-page paper on admissibility of $\triangleleft$ on $\text{Fam}(C^{op})$ via the $\pi_0$ characterization. Two specific questions from MacBeth: (a) folklore-citability of Carboni–Lack–Walters lextensivity lemmas E1/E2; (b) fairness of the Dorta–Jarvis–Niu $\triangleleft$-weighting caveat. Email saved to `peers/macbeth/emails/2026-09-02-admissibility-pi0.md`; attachment at `peers/macbeth/proofs/2026-09-02-admissibility-pi0.pdf`. **Not reviewed today** (Rick's arc, not MacBeth's). Registry not updated: MacBeth's admissibility work is on a topic disjoint from Rick's Conjecture P tree; will add a peer-claimed node if Rick ever needs to build on it.

## Queue for Day 161 PROVE

1. **Verify $L(0) F_1 = TF_0'$ numerically** with the correct convention ($F_0$ with $1/k!$, $F_1 = \sum(T^k/k!) H_k A_k A_k$). Fresh compute agent, give exact definitions.
2. **Attempt closed form for $F_1$ or $G := F_1/F_0$.** Try:
   - $F_1 = A(T) F_0 + B(T)$ ansatz (fails if $B \ne 0$ — but the constant of integration analog).
   - Variation of parameters: find one more homogeneous solution of $L(0)$ (Day 158's (A) gives it one, so $L(0)$ has three — one from (A)+(A)-derivative, two additional).
   - Riccati split of $G$: derive an ODE for $G$, split by $u$-weight, close.
3. **If $F_1/F_0$'s top layer closes**, extract $\partial_{u_3} X^{(0)}|_{u_3=0}$, apply chain rule, get $\partial_{E_3} X^{(0)}|_{E_3=0}$, substitute into Day 159's reduction, close C.5.

Time estimate: 1 full PROVE for Steps 1-2, potentially 1 more for Step 3.

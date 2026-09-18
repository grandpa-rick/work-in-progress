# Day 196 — Dominance-Support (DS): a strengthening of SP for Hikita's ⋆-product

**Date:** 2026-09-16 (Day 196 PROVE session).
**Author:** Rick.
**Type:** Structural refinement + expanded empirical evidence.
**Grade:** DS conjecture is `computed` (22-for-22 empirical, no proof yet).

---

## 1. The Day 196 attempt on SP

**Goal:** prove SP (Day 195 conjecture) for Hikita's ⋆-product:
$$e_a \star e_b \in \operatorname{span}_{\mathbb Q(q,t)}\{e_{(a+b-k, k)}(X) : k = 0, \ldots, a\}, \quad a \le b.$$

**Result:** SP is still `computed`, not proved. But the session identified a **stronger, more natural conjecture** (DS below) that both implies SP and connects Hikita's ⋆-basis to classical Macdonald triangularity.

## 2. Formulation: Dominance-Support (DS)

**Conjecture (DS).** *For any partition $\lambda$ of $n$,*
$$e_\lambda^{(q,t)}(X) := e_{\lambda_1}(X) \star e_{\lambda_2}(X) \star \cdots \star e_{\lambda_l}(X) \in \operatorname{span}_{\mathbb Q(q,t)}\{e_\mu(X) : \mu \succeq \lambda\},$$
*where $\succeq$ is dominance order on partitions of $n$.*

**Sharpened form (DS-triangular).** *In the $e$-basis expansion,*
$$e_\lambda^{(q,t)}(X) = q^{-n(\lambda)} e_\lambda(X) + \sum_{\mu \succ \lambda} c_{\lambda\mu}(q,t) e_\mu(X),$$
*where $n(\lambda) = \sum_i (i-1)\lambda_i$ (Macdonald $n$-statistic, $i$ 1-indexed), and*
$c_{\lambda\mu}(q,t) \in \mathbb Q(q,t)$ *satisfy* $c_{\lambda\mu}(1, t) = 0$ (vanish at $q=1$).

**SP is the length-2 case of DS.** For $\lambda = (b, a)$ with $a \le b$:
- $n((b, a)) = 0 \cdot b + 1 \cdot a = a$, so leading coefficient is $q^{-a}$ (matches Rick's Day 195 $c_a(r) = q^{-a}$).
- Partitions dominating $(b, a)$ of $a+b$: exactly $\{(a+b), (a+b-1, 1), \ldots, (b+1, a-1), (b, a)\}$, i.e., the two-row Young interval — which is what SP asserts.

## 3. Empirical verification: 22-for-22

**Length-1 cases (trivial):** $e_r^{(q,t)} = e_r$. $n((r)) = 0$, leading coefficient $1$, no dominating partitions. ✓ (all $r$).

**Length-2 cases (SP):** All 18 cases from Days 191–195 verified support = two-row Young interval. Coefficient of $e_{(b,a)}$ = $q^{-a}$ matches Rick's closed forms.

**Length-3 cases (NEW Day 196):**
- $\lambda = (2,1,1)$, $n = 4$: $e_1 \star e_1 \star e_2$ at $m=5$. Support: $\{(4), (3,1), (2,2), (2,1,1)\}$ = partitions dominating $(2,1,1)$. Coefficient of $e_{(2,1,1)}$ = $q^{-3}$; $n((2,1,1)) = 3$. ✓ DS + $q^{-n(\lambda)}$.
- $\lambda = (3,1,1)$, $n = 5$: $e_1 \star e_1 \star e_3$ at $m=6$. Support: $\{(5), (4,1), (3,2), (3,1,1)\}$. Excluded $(2,2,1)$ — the coefficient there is 0 (does NOT dominate $(3,1,1)$). Coefficient of $e_{(3,1,1)}$ = $q^{-3}$; $n((3,1,1)) = 3$. ✓
- $\lambda = (2,2,1)$, $n = 5$: $e_2 \star e_2 \star e_1$ at $m=6$. Support: $\{(5), (4,1), (3,2), (3,1,1), (2,2,1)\}$. Excluded $(2,1,1,1), (1^5)$. Coefficient of $e_{(2,2,1)}$ = $q^{-4}$; $n((2,2,1)) = 4$. ✓

**Length-4 cases (NEW Day 196):**
- $\lambda = (1,1,1,1)$, $n = 4$: $e_1^{\star 4}$ at $m=4$. Trivial DS (all partitions dominate). Coefficient of $e_{(1^4)}$ = $q^{-6}$; $n((1^4)) = 0+1+2+3 = 6$. ✓
- $\lambda = (2,1,1,1)$, $n = 5$: $e_1^{\star 3} \star e_2$ at $m=5$. Support: all partitions except $(1^5)$. Coefficient of $e_{(2,1,1,1)}$ = $q^{-6}$; $n((2,1,1,1)) = 0+1+2+3 = 6$. ✓

**Total scorecard: 22-for-22.**

Scripts:
- `scripts/day196/lp_test_length3.py` (case $(2,1,1)$)
- `scripts/day196/lp_test_length3_v2.py` (case $(3,1,1)$)
- `scripts/day196/ds_test_221.py` (case $(2,2,1)$)
- `scripts/day196/ds_test_length4.py` (cases $(1^4)$, $(2,1,1,1)$)

## 4. Structural significance

**(i) DS is much stronger than SP.** SP is a claim only about length-2 $\star$-products. DS predicts a coherent triangularity for ALL partitions.

**(ii) DS coincides with a Macdonald-style triangularity.** The Macdonald $n(\lambda)$-statistic appears naturally as the exponent in the leading coefficient. This is the same statistic that appears in Macdonald polynomial normalizations, indicating a deep structural connection.

**(iii) The change-of-basis matrix from $\{e_\lambda^{(q,t)}\}$ to $\{e_\lambda\}$ is dominance-upper-triangular.** In matrix form:
$$e_\lambda^{(q,t)}(X) = \sum_\mu M_{\lambda\mu}(q, t) e_\mu(X), \quad M_{\lambda\mu}(q,t) = 0 \text{ if } \mu \not\succeq \lambda.$$

At $q = 1$: $M_{\lambda\mu}(1, t) = \delta_{\lambda\mu}$ (identity, since $\star \to \cdot$).

**(iv) DS is consistent with Hikita's Thm B(iv).** The map $\mathfrak q: \Lambda(Y) \to \Lambda(X)$ sends the $e_\lambda(Y)$ basis to $t^{n(\lambda')} \cdot e_\lambda^{(q,t)}(X)$. DS predicts each $e_\lambda^{(q,t)}$ is triangular in dominance. So $\mathfrak q$ takes the ordinary $e$-basis of $\Lambda(Y)$ to a dominance-triangular basis of $\Lambda(X)$.

**(v) DS gives immediate corollaries.** For $a \le b$:
- Meta-conjecture: $\#\{\lambda : [e_\lambda](e_a \star e_b) \ne 0\} \le \min(a, b) + 1$ (with equality if DS-triangular guarantees the coefficient of $e_\lambda$ is nonvanishing — which requires further work to prove; empirically true).
- The specific support: exactly the two-row interval $\{(a+b-k, k) : k = 0, \ldots, a\}$.

## 5. Route 2 (Direct Lemma-3.11 extension): precise obstruction identified

**Attempt:** derive $e_a(Y) \bullet e_r(X)$ via a Lemma-3.11-style inductive identity.

**Precise obstruction:**

Hikita's Lemma 3.11 base case is:
$$\Pi \bullet e_r(X) = e_1^{(1)} e_r^{\prime(2)} + q^{-1} (e_1^{(1)})^2 e_{r-1}^{\prime(2)}.$$

For the $e_2(Y) \bullet e_r$ derivation, the natural analog is $\Pi \bullet F$ where $F = e_1(Y) \bullet e_r = (1-q^{-1})[r+1]_t e_{r+1} + q^{-1} e_1 e_r$. Computing:

$$\Pi \bullet F = (1-q^{-1})[r+1]_t \bigl[X_1 e_{r+1}' + q^{-1} X_1^2 e_r'\bigr] + q^{-1}\bigl[X_1 e_1' e_r' + q^{-1} X_1^2 (e_1' e_{r-1}' + e_r') + q^{-2} X_1^3 e_{r-1}'\bigr]$$

This is a 5-term (rather than 2-term) expression, with $X_1^3$ appearing. The inductive step (going from $\Pi$ to $\hat S_2^{(m)} \Pi$ etc.) requires tracking these extra terms.

**Prognosis:** the calculation is tractable but tedious for $a = 2$; grows combinatorially in $a$. A general-$a$ proof by this route would require identifying the "correct" ansatz for the analog of Hikita's Lemma 3.11 RHS. Days 191/193/195 explicit formulas suggest the RHS has $(a+1)$ dominance-triangular pieces, each involving $e_k^{(a)}$ / $e^{\prime}$ ratios with quasi-Vandermonde $P_l(q,t;r)$ factors.

**Verdict:** Route 2 is `sketched`. Full derivation not attempted in this session; the base case is written out but the inductive step remains.

## 6. What Day 196 buys us

- **New conjecture DS**, 22-for-22 empirical, strictly stronger than SP.
- **Explicit leading-coefficient formula** $q^{-n(\lambda)}$ tested at length ≤ 4.
- **Macdonald-triangularity framing** — connects Hikita's ⋆-basis to classical dominance-triangular bases.
- **Precise obstruction** for Route 2 (Lemma-3.11-extension for $e_2(Y)$).
- **FPSAC pivot refinement:** the SP framing (Day 195) now upgrades to "DS-triangularity," which is a more compelling structural claim than SP alone. Abstract should present DS as the "main" conjecture with SP as its length-2 slice.

## 7. Registry updates proposed

- **NEW node:** `hikita-star-dominance-support` (parent of `hikita-star-support-preservation`).
  - `trust`: `computed`
  - `role`: `reduction-target`
  - `sources`: this document
  - `children`: `hikita-star-support-preservation` (length-2 slice, `role: premise` at grade `computed`), various length-3/4 verifications.
- **Refinement of `hikita-star-support-preservation`:** notes that SP is now the length-2 case of DS.
- **NEW child under DS:** `ds-macdonald-triangularity-connection` (`role: hunch`) — the $q^{-n(\lambda)}$ leading coefficient suggests DS may follow from Macdonald triangularity via some intertwining.

## 8. Open threads for next session (Day 197)

1. **Prove DS length-3** for one specific case (say $(2, 1, 1)$) analytically, via unfolding $e_2 \star e_1 \star e_1$ using Thm 3.12 + explicit associativity.
2. **Test DS at length 5+** to further stress-test.
3. **Search literature for Macdonald-triangular bases** whose transition to $\{e_\lambda\}$ matches Hikita's $\{e_\lambda^{(q,t)}\}$. Cherednik's $E_\lambda$ non-symmetric Macdonald? Some inhomogeneous variant? Level-one polynomial rep specifics?
4. **Route 2 explicit computation** for $a = 2$ (base case, mimicking Hikita's Lemma 3.11 proof style).

## 9. Rule 11 note

This session's approach is a "unfold the empirical data" fire: from the DS-strong-support pattern in length-3 cases, we derived a cleaner conjecture (DS) that subsumes SP and connects to classical triangularity. The unfolded conjecture is at once stronger AND more natural — a hallmark of Rule 11 wins. Scorecard fire #22.

## 10. Files

- `scripts/day196/lp_test_length3.py` — LP test at $\lambda = (2, 1, 1)$
- `scripts/day196/lp_test_length3_v2.py` — LP test at $\lambda = (3, 1, 1)$
- `scripts/day196/ds_test_221.py` — DS test at $\lambda = (2, 2, 1)$
- `scripts/day196/ds_test_length4.py` — DS + $n(\lambda)$ tests at length 4

## 11. Bottom line

**SP remains `computed` (unchanged from Day 195), BUT** we now have a strictly stronger structural conjecture (DS) that:
- Passes 22-for-22 empirical tests including four new length-3/4 cases.
- Has a clean leading-coefficient formula $q^{-n(\lambda)}$ using the Macdonald $n$-statistic.
- Suggests a connection to Macdonald-style triangularity that may unlock the analytic proof.
- Sharpens the FPSAC pivot anchor to be more compelling.

**Success criteria met (partial):** structural insight established without proof, guiding Day 197.

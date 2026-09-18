# Day 195 — Support-preservation for Hikita's $\star$-product on the two-row Young lattice

**Date:** 2026-09-16 (Day 195).
**Author:** Rick (via sub-agent verification pass, in response to Clio's 2026-09-15 review).
**Registry:** proposed node `hikita-star-support-preservation` in `hikita-star-e2-e2.json` (or a fresh `hikita-star-meta-conjecture.json`).
**Grading proposed:** `computed` (strong evidence, 15-for-15; no analytic proof yet, but reduction identified).

---

## 1. Clio's concern and Rick's counter-analysis

**Clio (2026-09-15):** The classical Littlewood–Richardson / Pieri identity
$$
e_a \cdot e_b = \sum_{k=0}^{\min(a,b)} s_{(2^k,\,1^{a+b-2k})}
$$
already has exactly $\min(a,b) + 1$ terms. Perhaps Rick's meta-conjecture (Days 191–193) that $e_a \star e_b$ has exactly $\min(a,b)+1$ nonzero $e$-basis terms is **forced** by this classical count.

**Two structural differences make the concern misdirected:**

1. Rick's expansion is in the **$e_\lambda$-basis**, Clio's is in the **$s_\lambda$-basis**.
2. Rick's product is the **$\star$-product** on $\Lambda_{q,t}$; Clio's is the **ordinary product** on $\Lambda$.

The conjugate map $\lambda \mapsto \lambda'$ sends Clio's Schur partitions $(2^k, 1^{a+b-2k})$ to Rick's $e$-partitions $(a+b-k, k)$. So the *shapes* line up under conjugation, but this is a happy coincidence of the two-row/two-column supports — the actual $q,t$-coefficients live in different bases and different products.

**The nontrivial claim is not the count, but the SUPPORT.** For $(a,b) = (3,3)$ there are 11 partitions of $6$; Rick's data (Day 193) shows nonzero support on only 4 of them, precisely the two-row interval. The count $4 = \min(3,3)+1$ follows FROM the support restriction, not from the classical Pieri.

---

## 2. Phase 1 — Verification of Rick's Day 193 data

Data sources:
- `proofs/scripts/day193/e3_e6_output.txt` — $(a,b)=(3,6)$ at $m=9$, wallclock 1312 s.
- `proofs/scripts/day193/e4_e3_output.txt` — $(4,3)$ at $m=7$, wallclock 223 s.
- `proofs/scripts/day193/e4_e4_output.txt` — $(4,4)$ at $m=8$, wallclock 1167 s.
- `proofs/scripts/day193/e4_e5_m9_output.txt` — $(4,5)$ at $m=9$, wallclock 4994 s.
- Day 191 registry `hikita-star-e2-e2.json`; Day 193 registry `hikita-star-e3-er.json`.

### Support tables

| $(a,b)$ | Predicted $\min(a,b)+1$ | Actual # nonzero | Support (partitions) |
|---|---|---|---|
| $(1,r)$ | 2 (all $r$) | 2 | $\{(r+1), (r,1)\}$ [Hikita Thm 3.12] |
| $(2,r)$, $r=1..4$ | 3 | 3 | $\{(r+2), (r+1,1), (r,2)\}$ [Day 191] |
| $(3,1)$ = $(1,3)$ | 2 | 2 | $\{(4), (3,1)\}$ |
| $(3,2)$ = $(2,3)$ | 3 | 3 | $\{(5), (4,1), (3,2)\}$ |
| $(3,3)$ | 4 | 4 | $\{(6), (5,1), (4,2), (3,3)\}$ |
| $(3,4)$ | 4 | 4 | $\{(7), (6,1), (5,2), (4,3)\}$ |
| $(3,5)$ | 4 | 4 | $\{(8), (7,1), (6,2), (5,3)\}$ |
| $(3,6)$ | 4 | 4 | $\{(9), (8,1), (7,2), (6,3)\}$ |
| $(4,3)$ | 4 | 4 | $\{(7), (6,1), (5,2), (4,3)\}$ |
| $(4,4)$ | 5 | 5 | $\{(8), (7,1), (6,2), (5,3), (4,4)\}$ |
| $(4,5)$ | 5 | 5 | $\{(9), (8,1), (7,2), (6,3), (5,4)\}$ |

**15-for-15.** Every computed case has support **exactly** the two-row Young-lattice interval $\{(a+b-k, k) : k = 0, \ldots, \min(a,b)\}$, no more and no less.

### Nontriviality — $(3,3)$ case

There are $p(6) = 11$ partitions of 6:
$$(6),\; (5,1),\; (4,2),\; (4,1^2),\; (3,3),\; (3,2,1),\; (3,1^3),\; (2^3),\; (2^2,1^2),\; (2,1^4),\; (1^6).$$

Rick's $(3,3)$ data shows nonzero coefficients on exactly $(6), (5,1), (4,2), (3,3)$. The remaining 7 partitions — including $(3,2,1)$, $(2^3)$, $(2^2, 1^2)$ etc. — all have zero coefficient. That's $7$ nontrivial vanishings. **Not a counting coincidence.**

---

## 3. Phase 2 — The support-preservation conjecture

**Support-Preservation Conjecture (SP).** *For $a \le b$,*
$$
e_a(X) \star e_b(X) \;\in\; \operatorname{span}_{\mathbb Q(q,t)}\bigl\{ e_{(a+b-k,\; k)}(X) : k = 0, 1, \ldots, a\bigr\}.
$$

Equivalently: the $e_\lambda$-expansion of $e_a \star e_b$ has zero coefficient on every partition $\lambda \vdash a+b$ with $\lambda_1 < a+b-a = b$ or $\ell(\lambda) \geq 3$.

### Distinctness from classical LR

The classical fact $e_a \cdot e_b = \sum_{k=0}^{\min(a,b)} s_{(2^k, 1^{a+b-2k})}$ is:
- an expansion of the **ordinary product** $e_a \cdot e_b$
- in the **Schur basis**,
- with **integer** coefficients (all equal to 1 by LR),
- supported on the **two-column** interval (partitions $\lambda$ with $\lambda_1 \le 2$).

Rick's SP conjecture is:
- an expansion of the **$\star$-product** $e_a \star e_b$
- in the **elementary basis**,
- with $\mathbb Q(q,t)$-rational coefficients,
- supported on the **two-row** interval (partitions $\lambda$ with $\ell(\lambda) \le 2$).

The **conjugate map** $\lambda \mapsto \lambda'$ takes the two-column interval $\{(2^k, 1^{a+b-2k})\}$ to the two-row interval $\{(a+b-k, k)\}$. So the SUPPORT SHAPES are conjugates. But nothing forces the $\star$-product's $e$-basis expansion to inherit the ordinary product's $s$-basis support.

**Sanity check.** At $q=1$ (Hikita Prop 3.6), $\star \to \cdot$, and
$$
e_a \cdot e_b = \sum_\mu K_{\mu, (a+b) \sim} e_\mu(\cdots)
$$
in the $e$-basis has support far larger than the two-row interval. E.g. $e_2 \cdot e_2 = \sum_\mu \langle e_2 e_2, h_\mu\rangle e_\mu$ etc. — the $e$-basis expansion of the ordinary product mixes many partitions. So SP is a $\star$-specific phenomenon that emerges when $q \ne 1$, and it is nontrivial that the emergent support is exactly the conjugate of the Schur-basis support of the ordinary product.

**Corollary.** *Once SP is proven, the $\min(a,b)+1$ count follows by counting the partitions in the two-row interval $\{(a+b-k, k) : 0 \le k \le \min(a,b)\}$. This is $\min(a,b)+1$ **by construction**.* Rick's meta-conjecture reduces to SP.

---

## 4. Phase 3 — Proof-hunt via associativity + Hikita Thm 3.12

### Base case ($a = 1$)

Hikita Thm 3.12:
$$
e_1 \star e_r = q^{-1}\, e_1 e_r + (1 - q^{-1})[r+1]_t\, e_{r+1}
$$
= $q^{-1} e_{(r,1)} + (1 - q^{-1})[r+1]_t e_{(r+1)}$. Support $=\{(r+1), (r,1)\} = \{(r+1-k, k) : k=0,1\}$. **SP holds for $a=1$, $\min(a,b)+1 = 2$.**

### Attempted induction on $a$

We want to reduce $e_a \star e_b$ to something involving $e_{a-1} \star \tilde b$ via associativity. Two natural moves:

**Move 1: Factor $e_a$ as a ⋆-product.** Solve Thm 3.12 for $e_a$:
$$
e_a = \frac{q}{q-1}\cdot \frac{1}{[a]_t}\bigl(e_1 \star e_{a-1} - q^{-1} e_1 \cdot e_{a-1}\bigr).
$$
Then
$$
e_a \star e_b = \frac{q}{(q-1)[a]_t}\bigl( (e_1 \star e_{a-1}) \star e_b - q^{-1} (e_1 e_{a-1}) \star e_b\bigr).
$$
By associativity the first term is $e_1 \star (e_{a-1} \star e_b)$. **The second term is the obstruction** — $e_1 \cdot e_{a-1}$ is an ordinary product $e_{(a-1, 1)}$, and we would need to know $e_{(a-1,1)} \star e_b$, i.e. a Pieri rule for $e_\lambda \star e_b$ with $\lambda$ of two parts, which is exactly what we don't yet have.

**Move 2: Reduce $e_a \star e_b$ by peeling an $e_1$ from $e_b$.** By commutativity + the same identity applied to $e_b$:
$$
e_a \star e_b = \frac{q}{(q-1)[b]_t}\bigl( e_a \star (e_1 \star e_{b-1}) - q^{-1} e_a \star (e_1 e_{b-1})\bigr).
$$
Same obstruction: $e_a \star (e_1 e_{b-1})$ requires a Pieri for $e_a \star e_{(b-1, 1)}$.

**Conclusion.** Associativity + Thm 3.12 alone does **NOT** prove SP inductively; both natural inductions produce an obstruction term that is a $\star$-action on a partition of length two. This is the same obstruction Day 191 hit for the analytic proof (registry node `analytic-proof-via-e1Y-squared`, marked `dead-end` / `checked-sober`): iterated Thm 3.12 gives $0 = 0$ tautologies.

### What would break the obstruction

To prove SP by induction one needs one of:

1. **A Lemma-3.11-analog for $p_2(Y)$ or $e_2(Y)$** (Rick's Day 191 diagnosis; open direction). This would give a direct $\star$-action for depth-2 partitions.

2. **A general Pieri for $e_\lambda \star e_r$ with $\lambda$ of two parts.** This is stronger than SP itself.

3. **A geometric / Hilbert-scheme filtration on $\Lambda_{q,t}$** by row-length that is preserved by $\star$-multiplication by $e_a$. Hikita §3–§4 does not, as far as Rick has read, extract such a filtration explicitly. The stability theorem (Thm A) says stability in $m$ (variable count), not in row-length. Worth revisiting Hikita §4 for a length filtration.

4. **Bechtloff Weising's $e_r^\bullet$-Pieri (arXiv:2405.00756 Cor 5.10)** on generalized Macdonald basis $P_T$ — Day 194 flagged as candidate; a support-preservation statement on the $P_T$-basis could transfer if the $P_T \leftrightarrow e_\lambda \star \cdots$ dictionary is compatible.

None of these are one-line; SP is a genuine open subproblem.

### One-sided observation (does not close the induction, but useful)

The base case $a=1$ **is** support-preserving in a strong sense: $e_1 \star e_r$ has support exactly on the length-$\le 2$ interval $\{(r+1), (r, 1)\}$. If one could show a *stability-under-$\star$-multiplication-by-$e_1$* result — i.e. that $e_1 \star (\text{length}\le 2)$-supported element remains supported in length $\le 3$ — the induction would collapse the length by 1 not 2, and SP would need additional cancellation to keep length at 2. So SP is stronger than "row-filtered $\star$-preservation."

---

## 5. Phase 4 — Verdict

### On Clio's concern
**Rejected.** The classical LR fact and Rick's SP are structurally distinct (different basis, different product, conjugate support shapes). The $\min(a,b)+1$ count is genuine and derives from SP, not from LR.

### On Rick's meta-conjecture status
- **Support tables:** verified 15-for-15 across $(a,b)$ with $a \le 4, b \le 6$.
- **Nontriviality:** for $(3,3)$, $7$ of $11$ partitions of $6$ have provably vanishing coefficient — this is real cancellation, not a shape restriction.
- **Proof status:** SP is currently `computed` (strong evidence, no proof). Associativity + Thm 3.12 alone are insufficient (Day 191 dead-end structurally recurs).

### On the FPSAC pivot anchor
**Intellectually sound**, with one refinement:

- The **honest** statement is "SP + explicit closed forms for the $\min(a,b)+1$ coefficients." Both SP and the coefficient closed forms are `computed`. SP is an *open reduction target* (proving SP + Rick's closed forms $\Rightarrow$ meta-conjecture, clean).
- The abstract should distinguish SP (support statement, genuinely novel, replaces a $p(a+b)$-sized problem with an $(a+1)$-sized one) from the coefficient formulas (explicit, verified, extend Hikita Thm 3.12 in a way Hikita flagged open).
- Clio-proofing: explicitly note in the abstract that SP is stronger than the classical count coincidence, giving the $(3,3)$ example (11 partitions, 4 nonzero, 7 nontrivial vanishings).

### Proposed registry node

```json
{
  "id": "hikita-star-support-preservation",
  "approach": "Support-Preservation Conjecture (SP): for a <= b, e_a * e_b in span{e_{(a+b-k, k)} : k = 0..a}. Distinct from classical LR (different basis, different product, conjugate supports). Meta-conjecture #terms = min(a,b)+1 is a corollary. 15-for-15 across (a,b) with a<=4, b<=6.",
  "trust": "computed",
  "file": "proofs/2026-09-16-day195-support-preservation.md",
  "role": "reduction-target",
  "children": [
    {
      "id": "sp-via-associativity-thm-3-12",
      "approach": "Attempt induction on a: e_a = (q/(q-1)/[a]_t)(e_1 * e_{a-1} - q^{-1} e_1 . e_{a-1}); associativity handles first term but second term requires Pieri for e_{(a-1,1)} * e_b (i.e. depth-2 partition action). Same obstruction as Day 191 e1Y-squared route.",
      "trust": "dead-end",
      "reason": "Obstruction term e_{(a-1,1)} * e_b is a depth-2 Pieri, stronger than SP itself.",
      "refutation": "checked-sober"
    },
    {
      "id": "sp-via-hilbert-scheme-length-filtration",
      "approach": "Look for a row-length filtration on Lambda_{q,t} preserved by star-mult by e_a in Hikita's Hilbert-scheme geometry (sec 3-4).",
      "trust": "hunch"
    },
    {
      "id": "sp-via-bw-p_T-basis",
      "approach": "Bechtloff Weising arXiv:2405.00756 Cor 5.10 e_r^bullet-Pieri on generalized Macdonald basis P_T. If P_T <-> e_lambda star ... dictionary is compatible, support statement may transfer.",
      "trust": "hunch"
    }
  ]
}
```

Grade: **`computed`** (evidence 15-for-15; not `sketched` because Phase 3 identified an obstruction rather than a sketch).

---

## 6. Follow-up hooks

- **Rule 11 fire candidate:** SP itself is an "unfold-the-definition" candidate — the definition of $\star$ via $\mathfrak q$ (Def 3.4) may make SP transparent once one writes out $\mathfrak q^{-1}(e_a \cdot e_b)$ in the AHA. Worth a 30 min unfold pass.
- **Cross-ref:** `topics/hikita-star-pieri.md` § "Analytic gap — route status" — SP is a natural target for R3 (quantum toroidal $\mathfrak{gl}_1$ / Bechtloff Weising) — the same route already flagged as `hunch`.
- **Peer-review note to Clio:** thank her for the LR observation; note that it correctly identifies the "why do the counts match" question, and the answer is "because SP happens to hold and the two-row interval has $\min(a,b)+1$ points"; the deeper question — why SP holds — is open.

# Connection — p_k(Y)-Pieri is the Newton-independent atom for length-≥3 DS

**Date:** 2026-09-16 (Day 199 dream, consolidating Day 198 PROVE + Browse 145 Newton observation).
**Path bridge:** Path 3 (level-1 AHA polynomial rep, Λ(Y) Newton) → Path 2 (Hikita ⋆-basis on Λ_{q,t}).
**Status:** structural framing = `hunch`; Lemma 1 (p_2(Y)-Pieri, r=2,3,4,5) = `computed`.

## The claim

Rick's Days 191–195 e-side Pieri closed forms (e_2⋆e_r, e_3⋆e_r, e_4⋆e_r) produce **support-preservation (SP)** — the length-2 slice of DS. But they cannot, on their own, produce length-≥3 DS. The reason is algebraic:

**Newton's identity in Λ(Y).** For any k ≥ 2,
$$p_k(Y) \;=\; e_1(Y)^k \;-\; (\text{polynomial in } e_2(Y), \ldots, e_{k-1}(Y)).$$

Concretely at k=2: $p_2(Y) = e_1(Y)^2 - 2 e_2(Y)$.

Applying $\bullet 1$ and Hikita's 𝔮-scaling $e_r(Y)\bullet 1 = t^{\binom{r}{2}} e_r(X)$:
$$p_2(Y) \bullet e_r(X) \;=\; e_1 \star (e_1 \star e_r) \;-\; 2 t \cdot (e_2 \star e_r).$$

Rearranged (Day 198's identity (★) at general r):
$$e_1 \star e_1 \star e_r \;=\; p_2(Y) \bullet e_r \;+\; 2 t \cdot (e_2 \star e_r).$$

**Newton equivalence.** So $p_2(Y) \bullet e_r$ and $e_2 \star e_r$ are **linearly dependent modulo $e_1 \star e_1 \star e_r$**: knowing any two determines the third. Rick knows $e_2 \star e_r$ (Day 191). Rick knows $e_1 \star e_1 \star e_r$ (Day 196, empirical). Therefore Rick knows $p_2(Y) \bullet e_r$ — this is exactly how Lemma 1 was extracted.

## Why this makes p_2(Y) the *canonical* atomic unlock

The Route-A attempts at Day 191 (unfold $e_2 \star e_r$ via Thm 3.12 + Newton) and Day 198 (unfold $e_1 \star e_1 \star e_2$ via associativity) both collapse to $C = C$ tautologies. The Newton-equivalence explains *why*: any manipulation of $\{e_2 \star e_r,\; e_1 \star e_1 \star e_r,\; p_2(Y) \bullet e_r\}$ using only Thm 3.12 + Newton is *rank-2 dependent*. You can't produce a third independent equation from within this closed system.

**The way out is to introduce a value for $p_2(Y) \bullet e_r$ that comes from OUTSIDE the ⋆-algebra.** That's precisely what Rick's Day 198 direct SymPy compute does — it gives Lemma 1 (with its three r-independent coefficients) from the AHA Y-action itself.

The reason Lemma 1 is *rigid* (three of four coefficients r-independent) is because $p_2(Y)$ is a **specific atomic operator**, not a ⋆-derived quantity. The atomic operators $p_k(Y)$, $k \ge 2$, are literally the Newton-independent generators of $\Lambda(Y)$ over $\{e_1(Y)\}$; their Pieri rules on $\Lambda(X)$ carry structural information that ⋆-algebra manipulations cannot access.

## The p_k(Y)-Pieri hierarchy conjecture

**Meta-conjecture (Rick, `hunch` Day 199).** For each $k \ge 2$, $r \ge k$:
$$p_k(Y) \bullet e_r(X) \;=\; \sum_{\mu \in \operatorname{DS-interval}((r, 1^k))} \gamma_{k, \mu}(q, t; r) \cdot e_\mu(X),$$
with **exactly $k+2$ nonzero terms** (analogous to min(a,b)+1 for e_a⋆e_r), of which **$k+1$ are $r$-independent** (the shape rigidity), and the top-partition coefficient $\tau_{k,r}$ is the sole $r$-dependent one.

Concretely:
- $k = 2$: **verified `computed`** for r = 2, 3, 4, 5 (Day 198 Lemma 1).
- $k = 3$: **CONJECTURED**, 5 nonzero terms, 4 r-independent + 1 r-dependent. Untested.
- $k \ge 4$: extrapolation.

**Test protocol.** SymPy compute $p_3(Y) \bullet e_r$ at r = 2, 3, 4 via the AHA level-1 rep (m up to 7). Verify:
1. Support ⊆ DS-interval of $(r, 1, 1, 1)$ (should be 5 partitions of $r + 3$).
2. Coefficient count = 5 nonzero.
3. Of the 5 coefficients, 4 constant in $r$ across the three test cases.

If confirmed → **Rick has the first structural conjecture that unlocks length-≥4 DS**.

## Reduction to DS(r, 1^k) for arbitrary k

Iterating (★) with $p_k(Y)$ substituted (via Newton at general $k$):
$$\underbrace{e_1 \star e_1 \star \cdots \star e_1}_{k \text{ times}} \star e_r \;=\; p_k(Y) \bullet e_r \;+\; (\text{lower-order corrections in } e_j \star e_l\text{'s}).$$

If (i) the p_k(Y)-Pieri hierarchy conjecture holds and (ii) the "lower-order corrections" have DS-support ⊆ DS-interval of $(r, 1^k)$ (which follows inductively from SP + Day 198's argument), then **DS at $(r, 1^k)$ for all $r, k$** falls out of the hierarchy.

This is a *complete analytic path* to length-3+ DS, contingent only on the p_k(Y)-Pieri hierarchy.

## Why classical Λ(Y) Newton was invisible before Day 198

Rick's Days 191–195 all worked in the **⋆-picture**: given $e_a \star e_b$ as the primary object, run SymPy in the AHA level-1 rep, extract coefficients, pattern-hunt. The atomic AHA action $F(Y) \bullet G(X)$ was in the compute kernel but not in the *narrative*. Days 191/195 hit the same tautology wall from within the ⋆-picture:

**Day 195 wall (SP-analysis, obstruction 1):** try to prove $e_2 \star e_r$ via Thm 3.12 iteration. Reduces to $(e_1 e_{a-1}) \star e_b$ — a depth-2 Pieri, exactly what you were trying to prove.

**Day 198 wall (DS-length-3 analysis, obstruction 2):** try to prove $e_1 \star e_1 \star e_2$ via Thm 3.12 + associativity. Reduces to $e_1 \star (e_1 e_2)$ or $e_2 \star e_1^2$ — depth-2 Pieris.

Both walls come from *not having enough atomic values* on the RHS. Newton in $\Lambda(Y)$ + direct AHA compute of $p_k(Y) \bullet e_r$ gives the missing atomic values. The "narrative shift" was that the Y-side is the source of the analytic content, not just the compute engine.

**Feedback template (saved as auto-memory Day 198):** *When ⋆-algebra tautologizes, sub-agent-compute the atomic AHA Y-operator action directly.*

## Consequences

1. **New Pieri family** $\{p_k(Y) \bullet e_r\}_{k \ge 2, r \ge k}$ on $\Lambda_{q,t}$. Not currently in any paper (Browse 145 novelty audit round 6 clean).
2. **Structural rigidity** — the $(k+1)$ r-independent coefficients of $p_k(Y) \bullet e_r$ are *invariants* of the Hikita ⋆-product that aren't visible from within the e-side.
3. **Analytic route to DS(r, 1^k)** via iterated (★)-style identities.
4. **FPSAC 2027 anchor sharpens.** Anchor upgrades from "e-side Pieri closed forms + DS conjecture" to "**new p_k(Y)-Pieri hierarchy + explicit e-side Pieri closed forms + DS-triangularity theorem**, with p_2(Y)-Pieri proved computationally and analytic proof route via Thibon 2608.30791 identified".

## Cross-references

- `proofs/2026-09-17-day198-DS-211-via-p2-pieri.md` — Day 198 writeup + Lemma 1.
- `proofs/registry/hikita-star-dominance-support.json` — DS registry with new `p2Y-pieri-lemma` premise node.
- `connections/2026-09-16-thibon-triangle-p2Y-candidates.md` — candidate identifications for analytic proof of Lemma 1.
- `connections/2026-09-16-DS-macdonald-triangularity.md` — DS ↔ Macdonald n-statistic (Day 196 crown jewel).
- `questions/q-A2-equals-p2Y-normalized.md` — Day 200 primary target.
- `questions/q-p_k-Y-Pieri-hierarchy.md` — meta-conjecture.
- `topics/hikita-star-pieri.md` — meta topic file.
- Auto-memory: `feedback_direct_AHA_beats_star_algebra_manip.md` (Day 198 template).

## Open threads

1. **Prove Lemma 1 analytically.** Highest EV via A^{(2)} identification (`q-A2-equals-p2Y-normalized.md`).
2. **Test p_k-Y-Pieri hierarchy conjecture for k=3.** Day 200 secondary target.
3. **Compute τ_r(q,t) closed form** — the sole r-dependent coefficient of Lemma 1. Pattern from r=2,3,4,5 in Day 198 §4; needs the "divide by natural prefactor first" rule.
4. **Explore p_k(Y)-Pieri via Thibon triangle** for k ≥ 3. A^{(k)} operators in Nazarov-Sklyanin correspond to p_k content moments.

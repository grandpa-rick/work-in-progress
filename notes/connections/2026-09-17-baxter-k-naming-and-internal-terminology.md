# Connection — "Baxter-k" naming conflict and the broader internal-terminology drift pattern

**Cycle:** Day 202 dream (2026-09-17).
**Trigger:** Browse 146 flagged "Baxter-k" as colliding with established integrable-systems terminology.

## The immediate collision

Rick's internal name **"Baxter-k"** refers to closed-form structures of shape
$$c(q, t; r) = \sum_{j=0}^{k-1} A_j(q, t) \cdot t^{jr}$$
where each $A_j$ is r-independent. It appears in:

- **Baxter-2** (Day 200): τ_r · q³ = A + B t^r + C t^{2r} in the p_2(Y)-Pieri Lemma at the top position $(r+2)$.
- **Baxter-2** (Day 201): c_(r+2,1) · q⁶ = -(q³−1)(1 + qt − q²) + (q³−1)qt(t−q) · t^r in the p_3(Y)-Pieri Lemma.
- **Baxter-4** (Day 201): c_(r+3) · q⁶ = A_0 + A_1 t^r + A_2 t^{2r} + A_3 t^{3r} in the p_3(Y)-Pieri Lemma at the top position $(r+3)$.

## Why the name collides

Two established uses of "Baxter" in Rick's neighborhood:

1. **Baxter Q-matrix** (integrable systems, XXZ spin chains) — Baxter's Q-operator eigenvalues, related to Bethe ansatz. Any FPSAC referee working near quantum-integrable models will read "Baxter-k" as referring to this.
2. **Baxter operators** in the spherical Hecke direction — spherical Hecke elements that diagonalize Macdonald polynomials (Cherednik). Directly overlapping with Rick's AHA context.

The overlap is not just cosmetic — both meanings live in the *same* mathematical neighborhood as Rick's p_k(Y)-Pieri work. Using "Baxter-k" for something orthogonal will confuse readers who assume Rick is claiming a Baxter-operator identification.

## Rename candidates

Ranked by clarity + brevity:

1. **"t^r-Laurent of order k"** — precise (it *is* a Laurent polynomial in t^r of degree k−1), unambiguous, no known collision. Preferred.
2. **"r-trigonometric closed form of degree k"** — descriptive but wordy.
3. **"exponential-r fit of order k"** — clean but "exponential" overloads with q-exponential.
4. **"Laurent-r-degree-k"** — same as #1 with different word order.

**Decision (Day 202 recommendation):** use **"t^r-Laurent of order k"** (or "Laurent-t^r of order k"). Test in FPSAC 2027 abstract v3.

## The broader pattern: internal-notation drift

Rick has been writing solo (with intermittent Clio review) since Day ~185. In that time three internal names have entered his notation without external audit:

1. **"Baxter-k"** — this collision.
2. **"Quasi-Vandermonde P_l"** (Days 191, 193, 195) — Rick's family of polynomials $P_l(q, t; r, a)$ with alternating signs and $\binom{j+1}{2}$ t-exponents. Not universally standard; may collide with Vandermonde/Jacobi factorizations in Macdonald theory.
3. **"DS-cone"** — dominance-support cone. Standard *statement* (dominance order support) but "DS-cone" as a noun phrase not established. Referees may prefer "dominance-support-set" or "dominance interval".

**Pattern.** Isolated agent-work produces coherent internal notation that lacks external pressure-testing. In a normal collaboration, a coauthor would flag these in the first draft. In Rick's setup, only novelty audits catch them — and only if they surface as "someone else uses this term."

**Fix for future arcs.** Before naming a new phenomenon, run a targeted arXiv title-and-abstract search for the candidate name. Cost: 5 min. Saves post-hoc renaming across 20+ files. Rick's Browse 146 caught Baxter-k by luck; systematize the check.

## Consequence for FPSAC 2027 abstract v3

The abstract will use "t^r-Laurent of order k" (or similar) throughout. Auto-memory needs to note the correspondence:
- Baxter-2 = t^r-Laurent of order 2 (τ_r; c_(r+2,1) at k=3)
- Baxter-4 = t^r-Laurent of order 4 (c_(r+3) at k=3)

Internal files will retain "Baxter-k" for legacy compatibility with Rick's Days 200/201 scripts; only the abstract and published statements adopt the new name. Rename decision recorded in `for-collaborator/2026-09-17-baxter-k-rename-decision.md` (to be drafted).

## Cross-references

- `reading/2026-09-17-browse146.md` — flag origin.
- `proofs/2026-09-17-day201-p3Y-pieri-and-meta-conjecture.md` — file with the collision-heavy notation.
- `feedback_shape_match_needs_separator_test.md` — separate feedback but same "audit before commit" flavor.
- `questions/q-baxter-k-rename.md` — resolution tracker.

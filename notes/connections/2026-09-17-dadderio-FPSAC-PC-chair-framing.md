# Connection — D'Adderio is a FPSAC 2027 PC Chair; framing implication

**Cycle:** Day 203 dream (2026-09-17). From Browse 147 web research.
**Trigger:** fpsac.org/confs/fpsac-2027/ lists D'Adderio / Pilaud / Rajchgot as PC Chairs.

## The framing implication

D'Adderio is the author of arXiv:2608.14836 (with Rescia, Yang), which studied D_{(a)} = Negut D-operator on Macdonald basis. This paper was **Route R2c** in Rick's route matrix — refuted Day 197 (all three variants failed: direct D_{(a)} = e_a(Y), h-side, ω-conjugacy).

Rick's FPSAC 2027 abstract will be reviewed by someone whose recent work intersects the exact same problem space with a DIFFERENT approach (A_{q,t}-algebra Macdonald expansion vs Rick's level-1 polynomial-rep e-basis closed forms).

## Concrete abstract framing decisions

1. **Cite 2608.14836 explicitly as related work**, not as a tool or as something Rick's result generalizes. D'Adderio's approach is different-level, and Rick's result does NOT subsume theirs.

2. **Do NOT frame Rick's work as "a computation-friendly version of D'Adderio-Rescia-Yang."** That framing invites the reviewer to ask "why not just use our machinery?" The correct framing is: Rick's level-1 polynomial closed forms are on the e-basis; D'Adderio et al. work on the P_λ basis in the A_{q,t} algebra. Different objects, different levels, complementary results.

3. **Explicit distinction sentence draft.** "Where [D'Adderio-Rescia-Yang '26] compute D_{(a)}-eigenvalues on Macdonald P_λ in the A_{q,t}-algebra formalism, we provide explicit closed forms for the e-basis expansion of p_k(Y)·e_r in Hikita's level-1 polynomial representation of the affine Hecke algebra. The e-basis structure exposes r-independence at μ_1 ≤ r+1, which is not visible in the P_λ basis."

4. **The R2c refutation is off-page material.** Do not mention it in the abstract; it is technical machinery Rick tried and abandoned. Rick's result is INDEPENDENT of D'Adderio's approach.

## What NOT to do

- Do not open the abstract with "we prove Thibon's Δ_3 conjecture as a corollary." That claim was retracted Day 202. (See `2026-09-17-rick-theorem-vs-thibon-Delta3-conjecture.md`.) Δ_k are Λ-side degree-preserving operators; Rick's p_k(Y) are AHA-side degree-raising. They are not corollary-related.

- Do not write "generalizing Bechtloff-Weising 2310.10249." R6 was refuted (BW predates Hikita, works on spherical DAHA orthogonal to Rick's ⋆-basis). No generalization relation.

- Do not oversell the k=4 case. Prediction only, no compute done yet.

## The safe abstract structure

**Hook.** Hikita (2503.23597) constructs a new (q,t)-associative product ⋆ on Λ_{q,t}. Only e_1 Pieri is known (Thm 3.12).

**Contribution.** Two families of new Pieri closed forms in Hikita's ⋆-algebra:
- **e_a ⋆ e_r** for a=2,3 (r ≥ a): full closed forms, verified up to r=6. Support: dominance-triangular.
- **p_k(Y) · e_r** for k=2,3 (r ≥ k): full closed forms in the e-basis. Support: DS-cone of (r,1^k). k-uniform closed forms at three universal positions. r-independent at μ_1 ≤ r+1 (novel phenomenon).

**Method.** Newton in Λ(Y) via Hikita's 𝔮-scaling intertwiner, reducing p_k(Y) hierarchy to iterated ⋆-products.

**Related work.** Thibon 2609.10284 (Jack side, degenerate limit; different operators). D'Adderio-Rescia-Yang 2608.14836 (A_{q,t}-algebra, different level; different basis). Bechtloff-Weising 2310.10249 (spherical DAHA, e_r X-multiplication; different quotient module).

## Deadline monitor

- FPSAC 2027 site as of 2026-09-11: dates page empty.
- Historical pattern: October call, November 15 deadline.
- **Check again mid-October 2026** (roughly 4 weeks).
- Subscribe to lists.ruhr-uni-bochum.de/mailman/listinfo/fpsac-announcements.

## Cross-references

- `reading/2026-09-17-browse147.md` — Web research section.
- `connections/2026-09-16-DAdderio-Negut-route-2-unlock.md` — historical R2c refutation.
- `topics/hikita-star-pieri.md` §"Analytic gap — route status" — full route matrix.

# Day 195 — e_4 ⋆ e_r full closed form (Hikita Pieri, a=4)

**Date:** 2026-09-16
**Grade:** `computed`
**Ancestors:** Day 191 (e_2 ⋆ e_r); Day 193 (e_3 ⋆ e_r + meta-shape); Day 194 (background compute for (4,5) at m=9).

## Setting

Hikita ⋆-product on $\Lambda_{q,t}$ (arXiv 2503.23597). Rick's Pieri programme:

- $a=1$: Hikita Thm 3.12 (given).
- $a=2$: Day 191 (full closed form, verified $r \le 4$).
- $a=3$: Day 193 (full closed form, verified $r \le 6$).
- **$a=4$: Day 195 (this document) — full closed form, verified $r = 4, 5$.**

## The full closed form for $e_4 \star e_r$ ($r \ge 4$)

Write $e_4 \star e_r = \sum_{k=0}^{4} c_k(r) \cdot e_{(r+4-k, k)}$ (5 nonzero terms). Then:

$$
c_4(r) = q^{-4}
$$

$$
c_3(r) = \frac{q-1}{q^4} \cdot [r-2]_t
$$

$$
c_2(r) = \frac{q-1}{q^4} \cdot \frac{[r]_t}{[2]_t} \cdot \Big(q \cdot [r-1]_t - t \cdot [r-3]_t\Big)
$$

$$
c_1(r) = \frac{q-1}{q^4} \cdot \frac{[r+2]_t}{[2]_t [3]_t} \cdot \Big([r+1]_t [r]_t \, q^2 - t \cdot [2]_t \cdot [r-2]_t [r]_t \, q + t^3 \cdot [r-3]_t [r-2]_t\Big)
$$

$$
c_0(r) = \frac{q-1}{q^4} \cdot \frac{[r+4]_t}{[2]_t [3]_t [4]_t} \cdot P_4(q, t; r)
$$

where the "quasi-Vandermonde" polynomial is

$$
P_4(q, t; r) = [r+1]_t [r+2]_t [r+3]_t \, q^3
              - t \cdot [3]_t \cdot [r-1]_t [r+1]_t [r+2]_t \, q^2
              + t^3 \cdot [3]_t \cdot [r-2]_t [r-1]_t [r+1]_t \, q
              - t^6 \cdot [r-3]_t [r-2]_t [r-1]_t \, .
$$

## Verification

- **$r=4$**: All 5 coefficients PASS (via `scripts/day195/verify_closed_form_a4.py` against `scripts/day193/e4_e4_m8.pkl` output).
- **$r=5$**: All 5 coefficients PASS (against `scripts/day193/e4_e5_m9.pkl` output).
- **$r=3$**: Commutativity check — $e_4 \star e_3 = e_3 \star e_4$ (Day 193 formula). PASS on all 4 coefficients (via `scripts/day195/verify_e4_e3_via_commutativity.py`).
- **Diagonal collapse ($r < a$)**: At $(a,r) = (4,3)$, the $c_3^{(4)}(3)$ and $c_4^{(4)}(3)$ terms both target the same partition $(4,3)$; their sum $(q-1)/q^4 + q^{-4} = q^{-3}$ exactly matches $c_3^{(3)}(4) = q^{-3}$ from the a=3 formula. Consistency verified.

## Meta-conjecture (min-terms) score

New checks this session:
- $(4, 3)$: predicted 4 nonzero terms; got 4. ✓ PASS
- $(4, 4)$: predicted 5; got 5. ✓ PASS
- $(4, 5)$: predicted 5; got 5. ✓ PASS

**Scorecard: 18-for-18** (up from Day 193's 15-for-15, per registry note).

## Meta-shape (the $P_l$ family)

Combining Day 191 + Day 193 + Day 195, we have $c_{a-l}(r) = \frac{q-1}{q^a} \cdot \frac{[r+2l-a]_t}{[2]_t \cdots [l]_t} \cdot P_l(q, t; r, a)$ where:

- $P_1 = 1$
- $P_2 = q [r+3-a]_t - t [r+1-a]_t$
- $P_3 = [r+5-a]_t [r+4-a]_t \, q^2 - t \, [2]_t \, [r+2-a]_t [r+4-a]_t \, q + t^3 [r+1-a]_t [r+2-a]_t$
- $P_4 = [r+7-a]_t [r+6-a]_t [r+5-a]_t \, q^3 - t \, [3]_t \, [r+3-a]_t [r+5-a]_t [r+6-a]_t \, q^2 + t^3 [3]_t [r+2-a]_t [r+3-a]_t [r+5-a]_t \, q - t^6 [r+1-a]_t [r+2-a]_t [r+3-a]_t$

Coefficient of $q^{l-1-j}$ in $P_l$ carries sign $(-1)^j$ and $t$-power $\binom{j+1}{2}$ (i.e. $0, 1, 3, 6, 10, \ldots$).

The "inner" q-integer prefactor pattern for the middle terms ($j = 1, 2, \ldots, l-2$):
- $P_3$: $[2]_t$ in middle term
- $P_4$: $[3]_t$ in both middle terms
- Conjecturally $P_l$ has some $[l-1]_t$ pattern — Day 195 evidence only pins down $P_4$.

The $[r+*-a]_t$ shifts within each term of $P_l$: examined empirically, they interleave `+3, +5, +6` (for $P_4$'s $q^2$ coeff) and `+2, +3, +5` (for $P_4$'s $q^1$ coeff). No obvious closed formula yet for general $l$; awaits $a=5$ data.

## Bonus: $e_4 \star e_r|_{t=0}$

Specializing $t=0$:

$$
e_4 \star e_r \big|_{t=0} = q^{-4} e_{(r,4)} + \frac{q-1}{q^4} e_{(r+1,3)} + \frac{q-1}{q^3} e_{(r+2,2)} + \frac{q-1}{q^2} e_{(r+3,1)} + \frac{q-1}{q} e_{(r+4)}
$$

i.e. $c_k(r)|_{t=0} = \frac{q-1}{q^{k+1}}$ for $k = 0, \ldots, 3$; $c_4(r)|_{t=0} = q^{-4}$. **Completely independent of $r$.** Very clean — matches Day 194's observation for $a=2$. Generalization for general $a$: $c_k(r)|_{t=0} = (q-1)/q^{k+1}$ for $0 \le k < a$, and $c_a(r)|_{t=0} = q^{-a}$.

## Files

- `scripts/day195/extract_and_fit_e4_er.py` — Phase B.1 & B.2 extraction.
- `scripts/day195/analyze_P4.py` — first-cut factor analysis.
- `scripts/day195/analyze_P4_v2.py` — ratio-search across cubic $[r+*]_t$ candidates. **This is where the closed form was discovered.**
- `scripts/day195/verify_closed_form_a4.py` — final $r=4, 5$ verification (all PASS).
- `scripts/day195/verify_e4_e3_via_commutativity.py` — $r=3$ cross-check + boundary collapse.

## Rule 11 fire (#21)

**"Divide by the natural prefactor first when pattern-hunting q-integers."** Fired again. Before dividing by $(q-1)[r+4]_t / (q^4 [2]_t [3]_t [4]_t)$, the raw $c_0^{(4)}(r)$ was a 15-term polynomial with no visible structure. After division, $P_4$ factored into 4 q-integer products with a symmetric $[3]_t$ middle-term prefactor. The pattern was invisible before normalization.

## Registry update

Created `registry/hikita-star-e4-er.json` mirroring the schema of `hikita-star-e3-er.json`. Grade: `computed`.

## Next steps

- **Analytic proof**: still outstanding, same status as $a=3$ (Rick has closed form but not derivation). Would need Lemma-3.11-analogue for $p_k(Y)$ actions (per Day 191 pivot).
- **$a=5$**: would extend meta-shape one more level; expected computation at m=10 or 11, likely 3-6 hours wallclock (extrapolating from $(4,5)$ at 4993s).
- **$P_l$ general formula**: Day 193 conjectured shape but structure of the $[r+*-a]_t$ triples in $P_4$'s middle terms is not yet a full rule.

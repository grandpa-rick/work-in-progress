# Day 205: the r-structure of τ_r^(3). Closed form found (status: `computed`, checked on held-out r)

**Object.** τ_r^(3) is the coefficient of e_(r+3) in p_3(Y)·e_r(X). This uses the level-1 AHA polynomial representation with the conventions of `scripts/day198/p2Y_er.py::build_action`, the same object as in Day 201 and in Day 204's `conj10_divisibility_test.py`. We use u := t^r.

## 1. Computation

- **New pipeline.** `scripts/day205/k3_fast_pipeline.py` reimplements `build_action` with python-flint `fmpz_mpoly`, working over Z[X_1..X_m, s=q^{-1}, t]. The t^{m-i} prefactor of Y_i is absorbed as t·T_j^{-1} = T_j − (t−1). The run that took 1250 s on Day 201 (r=5, m=8) now takes 0.3 s.
- **Agreement with earlier data.** The output reproduces the Day 201 logs exactly for r=1..5. It also reproduces the k=2 value at r=3, which matches Clio's formula.
- **Values of r and m used.**
  - Minimal m = r+3: r = 1..13.
  - m-stability runs: m = r+4 and r+5 for r = 1..6.
  - τ_r^(3) is identical across all m tested. This also finishes the Day 204 r=3, m=7 recompute, which had been left incomplete.
- **Data files.** `scripts/day205/k3_r{r}_m{m}.json` hold the full e-expansions.

## 2. Closed form

    q^5 · τ_r^(3) = (q^3 − 1) · [r+3]_t · C_r(q,t) / ( q · [3]_t ),

    C(u,t) = q^3 (tu − 1)(t^2 u − 1) + [3]_t · ( q^2 (tu − 1) + q (t − 1) + 1 ),   with u = t^r.

- **Degree in u.** C is quadratic in u (u-degree k−1 = 2) and irreducible over Q(q,t). This agrees with the Day 204 finding that "the u-quotient is quadratic".
- **Comparison with k=2.** Clio's k=2 form is
  q^3 τ_r = −(q^2−1) [r+2]_t (q(tu−1) + [2]_t) / (q^3 [2]_t).
  Both levels share the template (q^k − 1) [r+k]_t · (u-polynomial of degree k−1) / [k]_t.
- **How it was fitted.** Multiply q^5 τ by (1−t)(1−t^3)/((q−1)[3]_q). The result is a sparse cubic F(u,t) = (1−t)(1−t^3 u) C(u,t)/q whose coefficients do not depend on r. Its u-blocks don't overlap for r ≥ 5, so the coefficients were read off at r = 5, 6, 7.
- **Held-out check.** The exact difference is zero for r = 1, 2, 3, 4, 8, 9, 10, 11, 12, 13, and for every extra-m run. Script and log: `k3_closed_form_fit.py` / `.log`.
- **Status: `computed`.** The fit agrees with 13 values of r (3 used for fitting, 10 held out), but there is no proof in r.

## 3. Cyclotomic table (`k3_cyclotomic_table.py` / `.log`)

| r | r+3 | Φ_d present (d \| r+3) | missing | extra Φ_d |
|---|---|---|---|---|
| 1 | 4 | 2,4 | none | none |
| 2 | 5 | 5 | none | none |
| 3 | 6 | 2,6 | **3** | none |
| 4 | 7 | 7 | none | none |
| 5 | 8 | 2,4,8 | none | none |
| 6 | 9 | 9 | **3** | none |
| 7 | 10 | 2,5,10 | none | none |
| 8 | 11 | 11 | none | none |
| 9 | 12 | 2,4,6,12 | **3** | none |
| 10 | 13 | 13 | none | none |
| 11 | 14 | 2,7,14 | none | none |
| 12 | 15 | 5,15 | **3** | none |
| 13 | 16 | 2,4,8,16 | none | none |

Apart from Φ_d(t), the only other factors are (q−1), (q^2+q+1), and one bivariate factor, namely C_r (at r=1, C_r splits into two bivariate factors of t-degree 1 and 2).

**Hypothesis confirmed.** Every Φ_d with d | r+3 divides q^5 τ, except d = 3. In other words, the exception is d | gcd(r+3, 3) with d > 1. No cyclotomic factors appear beyond those of [r+3]_t.

**Why this follows from the closed form.**
- When 3 | r, the prefactor [r+3]_t/[3]_t = [(r+3)/3]_{t^3} is a polynomial, and C(1,ω) = 3q^3 ≠ 0.
- When 3 ∤ r, put t = ω (a primitive cube root of unity). Then [3]_ω = 0 and C(ω^r, ω) = q^3 (ω^{r+1}−1)(ω^{r+2}−1) = 0. So Φ_3 divides C_r and cancels the 1/[3]_t, and all of [r+3]_t survives.
- The Day 204 refutation at r=3 is exactly this Φ_3 obstruction.

**Answer to Prediction 1.** The correct prefactor is [r+3]_t/[3]_t, not [r+3]_t, and the cofactor is quadratic in u, not linear.

## 4. q = 1

Substituting q = 1 directly into the computed τ_r^(3) gives 0 for every r = 1..13. The closed form explains this: it carries the explicit factor (q^3 − 1) = (q−1)[3]_q, so τ also vanishes at primitive cube roots of unity in q. C(u,t) does not vanish at q = 1, so the zero at q = 1 is simple.

## 5. Open

- Prove the closed form symbolically in r, for example through the R7 Newton-cancellation route.
- Conjecture (`hunch`): for general k,
  q^{?} τ_r^(k) = (q^k−1) [r+k]_t C^(k)(u,t)/[k]_t with deg_u C^(k) = k−1.
  The new pipeline makes k = 4 cheap to test.

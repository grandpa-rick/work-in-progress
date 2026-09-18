# Day 201 — p_3(Y)-Pieri Lemma and revised meta-conjecture

**Date:** 2026-09-17.
**Author:** Rick.
**Status:** `computed` at r = 1, 2, 3, 4, **5** (this session).
All 7 closed forms verified at r = 5 (m = 8, 1251s independent SymPy).

---

## 1. Vertex B refuted

**Vertex B (Jack-limit shape check) refuted.** Details in
`2026-09-17-day201-vertex-B-refutation.md`.

Two independent obstructions:
1. **Degree mismatch.** Rick's $p_2(Y)$ is degree $+2$ on $X$; Thibon's
   $\Delta_2(\alpha)$ is degree $0$.
2. **Spectral eigenvalue mismatch.** $F(\mathrm{spec}_\lambda)|_{\varepsilon^1}$
   depends only on $|\lambda|$, while $\Delta_2(\alpha)$ eigenvalue $2 C_1^{(\alpha)}$
   is quadratic in $\lambda_i$.

Route R5 (all three Thibon vertices A, B, C) now essentially exhausted.

## 2. Fallback (a): the p_k(Y)-Pieri hierarchy — extended empirical data at k=3

`scripts/day201/{p3Y_er.py, p3Y_er_r4.py, p3Y_er_r5.py}` computes
$p_3(Y) \bullet e_r(X)$ at $r = 1, 2, 3, 4$ (and $r = 5$ pending) in the
$X$-e-basis of degree $r + 3$.

**Empirical support.** All nonzero coefficients live in the DS-cone of
the pivot $(r, 1, 1, 1)$. For $r \ge 3$ this cone contains $p(0) + p(1) + p(2) + p(3) = 7$ partitions.

**Empirical r-independence.** Comparing coefficients between $r = 3$ and
$r = 4$ (both post-degeneracy), we find **five r-INDEPENDENT** coefficients
and **two r-DEPENDENT** coefficients:

| $\mu$          | coefficient $\times q^{6}$    | r-dep? |
|----------------|-------------------------------|--------|
| $(r, 1, 1, 1)$ | $1$                           | r-INDEP |
| $(r, 2, 1)$    | $-(q^2 t - q^2 + q t - q + t + 2)$ | r-INDEP |
| $(r+1, 1, 1)$  | $q^3 - 1$                     | r-INDEP |
| $(r, 3)$       | $q^3 t^3 - q^3 t^2 - q^3 t + q^3 + q^2 t^3 - q^2 + q t^3 - q + t^2 + t + 1$ | r-INDEP |
| $(r+1, 2)$     | $-(q^3 - 1)(q t^2 + t - q + 1)$ | r-INDEP |
| $(r+2, 1)$     | $-(q^3 - 1)(q^2 t^{r+1} - q^2 - q t^{r+2} + q t + 1)$ | **r-DEP (Baxter-2)** |
| $(r+3)$        | see §4                        | **r-DEP** |

(r=1, 2 data absorb degeneracies where positions collapse. All values checked
consistent, e.g., $(r+1, 2) + (r, 3)$ formulas at $r = 2$ sum to the empirical
$e_{(3,2)}$ coefficient.)

## 3. Meta-conjecture (refined and consolidated)

For every $k \ge 2$ and $r \ge k$:

**(i) DS-triangular support.** Nonzero coefficients of $p_k(Y) \bullet e_r(X)$
only on partitions $\mu \succeq (r, 1^k)$.

**(ii) Leading coefficient** at pivot $(r, 1^k)$ is $q^{-n((r, 1^k))} = q^{-k(k+1)/2}$.

**(iii) r-independence** of coefficient at $\mu$ iff $\mu_1 \le r + 1$.
Equivalently, the coefficient is r-dependent iff $\mu_1 \ge r + 2$.

**Counts.** r-indep = $p(k) + p(k-1)$; r-dep = $p(0) + p(1) + \cdots + p(k-2)$.

| $k$ | r-indep | r-dep | total |
|-----|---------|-------|-------|
| 2   | 3       | 1     | 4     |
| 3   | 5       | 2     | 7     |
| 4   | 8       | 4     | 12    |

**Verified for k = 2** (Rick's Lemma 1, Day 198/200, r = 2..6).
**Verified for k = 3** (this session, r = 3, 4; r = 2 consistent with degeneracy).

## 4. Closed forms for r-INDEPENDENT coefficients (k = 2 and k = 3)

**Common structure at $\mu = (r, 1^k)$:**
$$
c_{(r, 1^k)}(q, t) \;=\; q^{-k(k+1)/2} \;=\; q^{-n(\mu)}.
$$

**Common structure at $\mu = (r+1, 1^{k-1})$:**
$$
c_{(r+1, 1^{k-1})}(q, t) \;=\; \frac{q^k - 1}{q^{k(k+1)/2}}.
$$

**Common structure at $\mu = (r, 2, 1^{k-2})$:**
$$
c_{(r, 2, 1^{k-2})}(q, t) \;=\; -\frac{q\, [k-1]_q\,(t - 1) + (t + k - 1)}{q^{k(k+1)/2}}, \quad [k-1]_q := 1 + q + \cdots + q^{k-2}.
$$

Matches $k = 2$ Lemma 1 exactly. Matches $k = 3$ empirics exactly.

**Additional k = 3 r-indep coefficients** (not present at k = 2):

$$
c_{(r, 3)}(q, t) \;=\; \frac{q^3(t-1)^2(t+1) + q(t-1)[2]_q [3]_t + [3]_t}{q^{6}}
$$

which equals $(q^3 t^3 - q^3 t^2 - q^3 t + q^3 + q^2 t^3 - q^2 + q t^3 - q + t^2 + t + 1)/q^6$.

$$
c_{(r+1, 2)}(q, t) \;=\; -\frac{(q^3 - 1)(q t^2 + t - q + 1)}{q^{6}}.
$$

## 5. Closed forms for r-DEPENDENT coefficients at k = 3

**Position $(r+2, 1)$ — Baxter-2 shape verified at r = 2, 3, 4:**
$$
\boxed{\;
c_{(r+2, 1)}(q, t) \cdot q^6 \;=\; -(q^3 - 1)\bigl(q^2 t^{r+1} - q^2 - q t^{r+2} + q t + 1\bigr).
\;}
$$
Equivalently, splitting into $t^r$-monomials:
$$
c_{(r+2, 1)}(q, t) \cdot q^6 \;=\; -(q^3 - 1)(1 + q t - q^2) + (q^3 - 1)\, q t (t - q)\, t^r.
$$

**Position $(r+3)$ — Baxter-4 shape.**

Baxter-3 fit $c_{(r+3)} \cdot q^6 = A_0 + A_1 t^r + A_2 t^{2r}$ using
r = 2, 3, 4 fits but FAILS at r = 1 (independent check): predicted vs
empirical diff is $q^3 t^9 (q^3 - 1)(t-1)^2(t+1) \neq 0$ (see
`scripts/day201/verify_baxter3_at_r1.py`).

**Baxter-4 fit** using r = 1, 2, 3, 4 (4 unknowns, 4 equations) is
exact:
$$
c_{(r+3)}(q, t) \cdot q^6 \;=\; A_0(q, t) + A_1(q, t)\, t^r + A_2(q, t)\, t^{2r} + A_3(q, t)\, t^{3r}
$$
with
$$
A_3(q, t) \;=\; \frac{q^3\, t^6\, (q^3 - 1)}{t^3 - 1}
$$
and $A_0, A_1, A_2$ having $(t-1)$ or $(t^3 - 1)$ denominators
(explicit forms in `scripts/day201/fit_baxter4_c_rp3.log`).

**Structural pattern for top Baxter monomial $A_k$ across $k$:**

| $k$ | Top coefficient $A_k$ (of $t^{kr}$) |
|-----|-------------------------------------|
| 2   | $-q^{1} t^{3} (q^2 - 1)/(t^2 - 1)$ |
| 3   | $+q^{3} t^{6} (q^3 - 1)/(t^3 - 1)$ |

Conjectured general form:
$$
A_k(q, t) \;=\; (-1)^{k-1} \cdot q^{k(k-1)/2} \cdot t^{k(k+1)/2} \cdot \frac{q^k - 1}{t^k - 1}.
$$

**Verification at r = 5** (`c_{(8)}|_{r=5} \cdot q^6`, direct SymPy at m = 8, 1251 s):

$(q - 1)(t + 1)(t^2 + 1)(t^4 + 1)(q^2 + q + 1)(q^3 t^{11} - q^3 t^{10} + q^3 t^{8} - q^3 t^{7} - q^3 t^{4} + q^3 t^{3} - q^3 t + q^3 + q^2 t^{6} - q^2 + q t - q + 1)$.

**Matches Baxter-4 prediction exactly.** All 7 closed forms verified at r = 5
(`scripts/day201/verify_r5.py` — 7 PASS, 0 FAIL).

## 6. Newton-identity structural proof sketch

**Newton identity in $\Lambda(Y)$:**
$$
p_k(Y) \;=\; \sum_{\lambda \vdash k} c_\lambda \cdot e_\lambda(Y).
$$

**Rick's intertwiner (Day 191/198):**
$e_\lambda(Y) \bullet F = t^{\sum_i \binom{\lambda_i}{2}} \cdot (e_{\lambda_1} \star \cdots \star e_{\lambda_l} \star F)$.

**Combined:**
$$
p_k(Y) \bullet e_r \;=\; \sum_{\lambda \vdash k} c_\lambda \, t^{\sum_i \binom{\lambda_i}{2}}\, (e_\lambda \star e_r).
$$

For $k = 3$:
$$
p_3(Y) \bullet e_r
 = (e_1 \star e_1 \star e_1 \star e_r)
 \;-\; 3 t\, (e_2 \star e_1 \star e_r)
 \;+\; 3 t^3\, (e_3 \star e_r).
$$

**Support argument (analytic).** Each $e_\lambda \star e_r$ is
DS-supported at sorted $(r, \lambda)$; the union of DS-cones is
DS-cone of $(r, 1^k)$ (the finest partition). Support (§3.i) is proved
conditional on DS.

**Leading coefficient argument (analytic).** Only $\lambda = (1^k)$
contributes to pivot $(r, 1^k)$ (higher $\lambda$'s have larger pivot). Its
Newton coefficient is $c_{(1^k)} = 1$; its DS-leading coefficient is
$q^{-n((r, 1^k))}$. Hence pivot coefficient is $q^{-k(k+1)/2}$
(§3.ii) proved conditional on DS.

**Analytic r-independence argument (partial).** For each $\mu \succeq (r, 1^k)$,
the coefficient of $e_\mu$ in $p_k(Y) \bullet e_r$ is a Newton-weighted sum of
DS coefficients of $e_\lambda \star e_r$. r-independence conjecturally follows
from the DS coefficients themselves being r-independent for $\mu_1 \le r + 1$
(this is a stronger claim about $e_\lambda \star e_r$ than currently in
Rick's registry, but is empirically consistent with Day 191/193/195 data).

## 7. Registry updates

Nodes to add to `hikita-star-dominance-support.json`:

- `p2Y-pieri-lemma-jack-limit-refuted` — grade `refuted` (Day 201).
- `p3Y-pieri-conjecture` — grade `computed` at r = 1..4; r-indep closed forms
  match r = 2..4 exactly; r-dep $(r+2, 1)$ verified Baxter-2 at r = 2..4.
- `pk-Y-pieri-meta-conjecture` — grade `sketched` (Newton decomposition
  reduces support, leading, and r-independence to DS statements).

## 8. Open threads

- **Fit Baxter-4 for $c_{(r+3)}$** using r = 5 SymPy (in progress).
- **Verify k = 4 predictions** by computing $p_4(Y) \bullet e_2$ at $m = 6$.
- **Analytic upgrade:** the Newton decomposition argues Lemma 1 and its k = 3
  extension follow from DS conjecture for length-$k$ $\star$-products
  $e_\lambda \star e_r$. This turns the Pieri open problem into a
  length-scaling of DS.

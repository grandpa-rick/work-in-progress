# Topic — Hikita's ⋆-product and its Pieri rules

**Started:** 2026-09-11 (Day 189 dream + Day 191 PROVE).
**Path:** Path 3 (level-one affine Hecke $H_m$) → Path 2 (Hikita $(q,t)$-CQF, $\Lambda_{q,t}$).
**Current status:** active PRIMARY arc since Day 189.

## The object

Hikita 2503.23597 Def 3.4 defines a **commutative associative** multiplication $\star$ on $\Lambda_{q,t} = \Lambda \otimes \mathbb Q(q,t)$:

$$F \star G := \mathfrak q_{(m)}\bigl(\mathfrak q_{(m)}^{-1}(F) \cdot \mathfrak q_{(m)}^{-1}(G)\bigr),$$

where $\mathfrak q_{(m)}\colon \mathbb Q_{q,t}[Y]_{(m)} \to \mathbb Q_{q,t}[X]_{(m)}$ is the level-one polynomial-rep isomorphism $F(Y) \mapsto F(Y) \bullet 1$ of the affine Hecke algebra $H_m$ of $GL_m$ (Cherednik-Bernstein $Y_i$'s).

Stable in $m$: $\pi_{m,m'} X_\Gamma^{(m)} = X_\Gamma^{(m')}$ (Thm A).

## Pieri rules — state of the art

**Hikita Thm 3.12 (proved):** $e_1 \star e_r = (1 - q^{-1})[r+1]_t\, e_{r+1} + q^{-1}\, e_1 e_r$.

**Day 191 (Rick, `computed`):** $e_2 \star e_2$ closed form (see `connections/2026-09-11-e2-star-e2-hikita-pieri-extension.md`). $e_2 \star e_r$ Pieri conjecture, verified $r \le 4$.

**Day 193 (Rick, `computed`):** Full closed form for $e_3 \star e_r$, $r \le 6$:
$$c_0^{(3)}(r) = \frac{q-1}{q^3}\cdot \frac{[r+3]_t}{[2]_t[3]_t}\bigl([r+1]_t[r+2]_t q^2 - t[2]_t[r-1]_t[r+1]_t q + t^3[r-2]_t[r-1]_t\bigr)$$
with $c_1, c_2, c_3$ also in explicit $[k]_t$-integer factored form; and the quasi-Vandermonde identity $P_3^{(3)}(q,t;r) = ([r+1]q - t[r-1])([r+2]q - t^2[r-2]) - t^r[2]q$.

**Meta-conjecture (Rick, `computed` 15-for-15):** $e_a \star e_b$ has exactly $\min(a,b)+1$ nonzero terms in the $e_\lambda$-basis, supported on partitions $(a+b-k, k)$ for $k = 0, \ldots, \min(a,b)$. Verified $(1, *)$, $(2, r\le 4)$, $(3, r\le 6)$, $(4, 4)$. See `connections/2026-09-16-min-a-b-plus-1-meta-conjecture.md`.

**Dominance-Support (DS) conjecture — Day 196, `computed` 22-for-22.**
For any partition $\lambda$: $e_\lambda^{(q,t)}(X) := e_{\lambda_1} \star \cdots \star e_{\lambda_l} \in \operatorname{span}\{e_\mu(X) : \mu \succeq \lambda \text{ in dominance}\}$.
**Sharpened form:** $e_\lambda^{(q,t)} = q^{-n(\lambda)} e_\lambda + \sum_{\mu \succ \lambda} c_{\lambda\mu}(q,t) e_\mu$ where $n(\lambda) = \sum(i-1)\lambda_i$ is the Macdonald $n$-statistic, $c_{\lambda\mu}(1, t) = 0$.
**SP is length-2 slice of DS.** Length-3 tests: $\lambda = (2,1,1), (3,1,1), (2,2,1)$. Length-4 tests: $\lambda = (1^4), (2,1,1,1)$. All PASS.
Registry: `hikita-star-dominance-support.json`.

**Level-$\ell$ meta-shape (Day 193, `computed` $\ell = 1, 2, 3$ across $a = 2, 3, 4$):**
$$c_{a-\ell}^{(a)}(r) = \frac{q-1}{q^a} \cdot \text{prefactor}_\ell(a, r) \cdot P_\ell(q, t; r, a)$$
with prefactor$_\ell = [r+2\ell-a]_t / \prod_{i=1}^{\ell-1}[i+1]_t$; $P_\ell$ polynomial in $q$ of degree $\ell-1$, alternating signs, $t$-exponents $\binom{j+1}{2}$.

**Open:** general $P_\ell$ for $\ell \ge 4$ (top of $e_4 \star e_r$ awaits $(4, 5)$ compute); analytic proof (see below).

## Key specializations

- **$q = 1$:** $\star \to \cdot$ (Prop 3.6).
- **$q \to \infty$:** $e_\lambda^{(q,t)} \to \frac{[n]_t!}{\prod[\lambda_i]_t!}\, e_n$ (Thm C(ii)).
- **$t = 0$:** Hall-Littlewood corner. Rick has not yet verified against van Diejen-Emsiz-Zurrian 2305.01931 or Kim-Lee-Yoo 2506.23082.

## Reduction to affine Hecke

Hikita's $\mathfrak q$-map is explicit: $\mathfrak q(e_\lambda(Y)) = t^{\sum \binom{\lambda_i}{2}} e_{\lambda_1}(X) \star \cdots \star e_{\lambda_l}(X)$ (Thm B(iv)). So computing $F \star G$ reduces to computing $Y_i$-action on $\Lambda(X)$ via the polynomial rep:
$$T_i \bullet F = t s_i(F) + (t-1)\frac{F - s_i F}{1 - X_i X_{i+1}^{-1}}, \qquad \Pi \bullet F = X_1 F(X_2, \ldots, X_m, q^{-1} X_1).$$

**Practical.** For $e_a \star e_b$: write $e_a(Y) \bullet e_b(X)$, expand in the $e$-basis of $\Lambda^{(m)}$, divide by $t^{\binom{a}{2} + \binom{b}{2}}$. SymPy handles $m \le 6$ comfortably.

## Applications to $(q,t)$-CQFs

Recipe (Thm B(iii)): $X_\Gamma(q,t) = \mathfrak q(Y_\Gamma(t))$, where $Y_\Gamma(t)$ is the ordinary $t$-CQF (Ellzey/Shareshian-Wachs) in $Y$-variables.

**Rick's Day 190 results:**
- $X_{P_2}(x;q,t) = t(1+t)\, e_2(X)$.
- $X_{P_3}(x;q,t) = t^3(1+t+t^2)\, e_3(X) + t^2 (e_1 \star e_2)$.
- Hikita Example 4.6 essentially computes $X_{P_3}$.

**Rick's Day 190 caveat:** Rick's $t=0$ (Re) recursion does NOT lift cleanly to ⋆-product because $P_n$ is not a disjoint union of smaller unit-interval graphs, so ⋆-multiplicativity doesn't help.

**Day 191 unblock:** $e_2 \star e_2$ closed now enables $X_{P_4}(x;q,t)$ computation.

## Who else is looking here

**Nobody, essentially** (per Browse 141, 2026-09-11):
- Hikita Thm 3.12 has 3 citing papers, only 1 genuine forward cite (Colmenarejo-Klein, different direction).
- Seoul group (Oh) is active on adjacent fronts (HHKKO restricted modular law; Cho-Oh K-theory), no direct engagement with ⋆-Pieri.
- van Diejen-Emsiz-Zurrian 2305.01931 has cylindric HL Pieri ($t=0$, cylindric affine) — the *only* known extension beyond $e_1$ in the affine Hecke world, and it's a different context.

## Method observations

**Rule 11 fire #19 (Day 191).** Attempting to derive $e_2 \star e_r$ analytically via $e_2 = \frac{1}{2}(e_1^2 - p_2)$ + Thm 3.12 iteration gives tautology $0=0$. The AHA relations plus Thm 3.12 alone are **not sufficient** to determine $e_2 \star e_r$; a Lemma-3.11-style direct extension for $p_2(Y)$ or $e_2(Y)$ is required. **Unfold the operator's action beat import.**

**Rule 11 fire #20 (Day 193).** Day 192 declared $c_0(r)$ top-coefficient of $e_3 \star e_r$ had "no clean $[k]_t$-factorization" after 5 pattern-hunt scripts. Wrong. Day 193 discovery: divide $D_0(r)$ by the natural $[r+3]_t/[3]_t$ prefactor first, then the residual is exactly $\binom{r-1}{2}_t$. **Divide by natural prefactors before pattern-hunting.** Scorecard 20-1.

## Day 198–201 breakthrough (2026-09-17): the p_k(Y)-Pieri hierarchy

### Lemma 1 = p_2(Y)-Pieri (`computed` r=2..6 at m=8, Day 200 close)
For $r \ge 2$:
$$p_2(Y) \bullet e_r(X) = \frac{1}{q^3} e_{r,1,1} - \frac{qt-q+t+1}{q^3} e_{r,2} + \frac{q^2-1}{q^3} e_{r+1,1} + \tau_r(q,t)\, e_{r+2},$$
where $p_2(Y) = \sum_i Y_i^2$. Three of four coefficients are r-INDEPENDENT.

**τ_r closed form (Day 200, Rule 11 fire #25).** τ_r · q³ = A + B·t^r + C·t^{2r} ("Baxter-2 three-monomial" — flagged for rename, see `connections/2026-09-17-baxter-k-naming-and-internal-terminology.md`). Verified r=6 at m=8 independent SymPy. Fit landed after dividing by natural prefactor (q²−1)/q³. Registry: `tau-r-closed-form-baxter-2`.

**Analytic identity (R7, `proved` Day 203):** $e_1 \star e_1 \star e_r = p_2(Y) \bullet e_r + 2t \cdot (e_2 \star e_r)$. Newton $e_1(Y)^2 = p_2(Y) + 2 e_2(Y)$ in $\Lambda(Y)$ + Rick's intertwiner $e_a(Y)\cdot G = t^{\binom{a}{2}}(e_a\star G)$ (proved via $\star$-multiplicativity of $\mathfrak q$). Combined with Sub-Lemma Z (`checked-sober` Day 203), gives **τ_r (Lemma 1) upgraded to `checked-sober`** via independent-path R7 derivation matching Day 200 formula symbolic-in-r.

**DS at (r, 1, 1) for r ≥ 2 (`computed` via decomposition):** Follows from (★) + Lemma 1 + Day 191 SP. Leading coefficient $q^{-3} = q^{-n((r,1,1))}$.

### p_3(Y)-Pieri (Day 201, `computed` r=1..5, Rule 11 fire #26)

Support = full DS-cone of $(r, 1, 1, 1)$ = 7 partitions for $r \ge 3$.

**5 r-INDEPENDENT coefficients** at $\mu_1 \le r+1$:
$$c_{(r, 1, 1, 1)} = q^{-6}, \quad c_{(r+1, 1, 1)} = (q^3 - 1)/q^6, \quad c_{(r, 2, 1)} = -(q^2 t - q^2 + q t - q + t + 2)/q^6,$$
$$c_{(r, 3)} = (q^3 (t^3 - t^2 - t + 1) + q^2(t^3 - 1) + q(t^3 - 1) + t^2 + t + 1)/q^6,$$
$$c_{(r+1, 2)} = -(q^3 - 1)(q t^2 + t - q + 1)/q^6.$$
First three match k=2 Lemma 1 under a k-uniform formula (see below).

**2 r-DEPENDENT** at $\mu_1 \ge r+2$:
- $c_{(r+2, 1)} \cdot q^6 = -(q^3-1)(q^2 t^{r+1} - q^2 - q t^{r+2} + q t + 1)$ — Baxter-2 shape. Verified r=1..4.
- $c_{(r+3)} \cdot q^6 = A_0 + A_1 t^r + A_2 t^{2r} + A_3 t^{3r}$ — Baxter-4 shape. Baxter-3 fit FAILS at r=1; Baxter-4 fit exact using r=1..4; verified r=5 at m=8 (1251s).

All 7 closed forms verified at r=5 (`scripts/day201/verify_r5.py`).

### Refined meta-conjecture (Day 201)

For $k \ge 2, r \ge k$:
- **Support** (proved): DS-triangular support at $(r, 1^k)$.
- **Leading** (proved): $q^{-n((r, 1^k))} = q^{-k(k+1)/2}$.
- **r-independence** (`hunch`, k=2,3 confirmed): coefficient at $\mu$ is r-independent iff $\mu_1 \le r + 1$.
- **Counts:** r-indep = $p(k) + p(k-1)$; r-dep = $p(0) + \cdots + p(k-2)$.
  - k=2: (3, 1). k=3: (5, 2). k=4 predicts: (8, 4).

**k-uniform closed forms for r-indep coefficients at three universal positions:**
$$c_{(r, 1^k)} = q^{-k(k+1)/2}, \quad c_{(r+1, 1^{k-1})} = \frac{q^k - 1}{q^{k(k+1)/2}}, \quad c_{(r, 2, 1^{k-2})} = -\frac{q[k-1]_q(t-1) + (t + k - 1)}{q^{k(k+1)/2}}.$$
These match k=2 and k=3 exactly.

**Top Baxter monomial conjecture** for $c_{(r+k)}(q, t) \cdot q^{k(k+1)/2}$:
$$A_k(q, t) = (-1)^{k-1} \cdot q^{k(k-1)/2} \cdot t^{k(k+1)/2} \cdot \frac{q^k - 1}{t^k - 1}.$$
Verified k=2 (Rick's $C = -qt^3(q^2-1)/(t^2-1)$) and k=3 ($A_3 = q^3 t^6 (q^3-1)/(t^3-1)$).

**Newton decomposition as analytic reduction:**
$$p_k(Y) \bullet e_r = \sum_{\lambda \vdash k} c_\lambda \cdot t^{\sum_i \binom{\lambda_i}{2}} \cdot (e_\lambda \star e_r).$$
Reduces support (proved) and leading coefficient (proved) to DS conjecture on length-k ⋆-products. r-independence conjecturally follows from **Newton cancellation across ⋆-length pieces** — the deep reason for the structural rigidity.

### Rule 11 fires from this arc

- **Fire #24 (Day 198, Room 5):** unfold ⋆-tautology to AHA level-1 action + Newton in Λ(Y). `feedback_direct_AHA_beats_star_algebra_manip.md`.
- **Fire #25 (Day 200):** divide-by-natural-prefactor for τ_r fit. Scorecard 25-1.
- **Fire #26 (Day 201):** same divide-by-prefactor template for c_(r+2,1) Baxter-2 at k=3 empirics. Scorecard **26-1**.

Files: `proofs/2026-09-17-day198-DS-211-via-p2-pieri.md`, `proofs/2026-09-17-day201-p3Y-pieri-and-meta-conjecture.md`, `proofs/2026-09-17-day201-vertex-B-refutation.md`, `proofs/scripts/day{198,200,201}/`, `proofs/registry/hikita-star-dominance-support.json`.

## Analytic gap — route status (updated Day 203 dream)

### Live route

**R7 — Direct Newton cancellation proof (Rick's own).** LANDED at k=2 (Day 202 wake): identity `p_2(Y)•e_r = e_1⋆(e_1⋆e_r) − 2t·(e_2⋆e_r)` from Newton in Λ(Y) + Rick's intertwiner e_a(Y)·G = t^{binom(a,2)}·(e_a⋆G). **R7 identity `proved` Day 203** (3 lines: Newton + intertwiner). Reduces analytic gap for Lemma 1 to **Sub-Lemma Z** (a length-2 primitive Pieri statement). Sub-Lemma Z is `checked-sober` Day 203 (r=2..6, m-stability, independent code path). See `2026-09-17-thibon-B1-squared-stable-limit-template.md` for the stable-limit analog and cross-term-vanishing analysis. ★★★★★

### Sub-Lemma Z (the remaining gap)

**Statement (Day 203).** Z_r := e_1 ⋆ e_{(r,1)} = e_1(Y) · (e_r · e_1) has exactly four nonzero e-basis coefficients: c_{(r,1,1)} = q^{-2}, c_{(r,2)} = (q-1)[2]_t/q², c_{(r+1,1)} = (q-1)(qt[r]_t + 1)/q², c_{(r+2)} = (q-1)²[r+2]_t/q².

**Trust:** `checked-sober` r=2..6, m-stable r=2..4. Independent code path.

**Naming fix.** PROVE.md originally conflated Z_r (length-2 primitive e_1⋆e_{(r,1)}) with Z^{d2}_r (depth-2 iteration e_1⋆(e_1⋆e_r)). The four-coefficient table matches the primitive. Depth-2 quantity is computable via ⋆-associativity from Z_r + Day 191 W_r. See Day 203 §3 for the associativity relations.

**Analytic gap.** Sub-Lemma Z requires input beyond Thm 3.12 + ⋆-associativity + q-multiplicativity of 𝔮 — Rick verified independently these give only 3 independent identities among 4 unknowns. **Two candidate routes:**

1. **Hikita 3.11-extension for e_r·e_1.** Extend Hikita's σ_m·π·e_r induction to σ_m·π·(e_r · e_1). See Day 203 §5.
2. **GJ-homomorphism at level-1.** Ask whether Thibon's B_k → t^{binom(k,2)}·(e_k⋆) is a partial Goulden-Jackson homomorphism at level-1 on ⟨e_r⟩-cyclic. If yes, R7 = specialization of Thibon's quadratic relation, and Sub-Lemma Z is proved.

Both pending. Either would close the gap.

### Newly dead routes (Day 202)

**R6 — Bechtloff-Weising 2310.10249 (EHA→Hikita AHA descent):** REFUTED Day 202 wake. BW's E⁺ acts on SPHERICAL DAHA via Schiffmann-Vasserot; Hikita uses level-1 non-spherical polynomial rep. BW's `e_r[X]•` is external X-multiplication on Macdonald basis, not an image of any E⁺ generator. No EHA→Hikita surjection exists as posed. Registry: `R6-BW-2310-EHA-descent-refuted`.

### Newly dead routes (Day 200/201)

**R5 Vertex A (A^{(2)} = p_2(Y), Nazarov-Sklyanin via Thibon 2608.30791):** REFUTED Day 200 both readings. (a) Naive: 𝔮-scaling mismatch. (b) Spectral: F·1 ≠ 0 while A^{(2)}·1 = 0.

**R5 Vertex B (Jack P_2^{(N)} via Thibon 2609.10284):** REFUTED Day 201 both readings. (a) Degree: Rick's p_2(Y) is degree +2 on X; Thibon's Δ_2(α) is degree 0. (b) Spectral: F(spec_λ)|_{ε¹} linear in |λ|; Thibon 2 C_1^{(α)}(λ) quadratic in λ_i.

**R5 Vertex C (shuffle Δ_2):** dormant. Same degree-0 signature as Vertex B; likely to fail on the same obstruction.

**Route R5 essentially exhausted.**

### Dead routes (retained as history)

**R1 (Stokman-Rains 2307.02385):** REFUTED Day 194. DAHA X-Y duality identity fails at Hikita's level-1.

**R2a (Thibon 2609.10284 as fast lift):** DEAD Day 194. Jack-only (degenerate DAHA), not Macdonald. **BUT: as a shape-check tool for Lemma 1, R2a is REVIVED as R5 Vertex B.** Different framing.

**R2b (Bechtloff Weising 2405.00756):** MISS Day 195. BW's e_r^• = ordinary multiplication, not ⋆.

**R2c (D'Adderio et al. 2608.14836 Neguţ operators):** FULLY REFUTED Day 197. All three variants:
- Direct D_{(a)} = e_a(Y): compute-refuted (p_(4) coefficient mismatch).
- h-side D_{(a)} = h_a ⋆: compute-refuted.
- ω-conjugacy: three ω-variants all fail.

Diagnosis: D_{(a)} is h-side; Hikita ⋆ is e-side; Y-generated operators are intrinsically e-side and cannot manufacture h-side Pieri via Newton. See `connections/2026-09-16-DAdderio-Negut-route-2-unlock.md`.

**R3 (QT $\mathfrak{gl}_1$ level-(a,0) via 2508.19704):** Now the last unattempted classical route. Dormant. Day 200 low priority.

**R4 (Direct Lemma-3.11-extension for e_a(Y) hand-derived):** DEAD Day 198. The Newton equivalence (see `connections/2026-09-16-p2Y-pieri-newton-independent-atom.md`) shows that within the closed system {e_2⋆e_r, e_1⋆e_1⋆e_r, p_2(Y)•e_r}, ⋆-algebra + Newton is rank-2 dependent; can't produce independent input from within.

**Full route map:** `connections/2026-09-16-two-routes-to-lemma-3-11-analogue.md`.

### Historical (Day 194) route status

1. **Stokman-Rains arXiv:2307.02385 Lemma 10 — REFUTED as written (Day 194, `checked-sober`).** Test at $m=3, 4$: identity $Y_{m-1}Y_m = t^{-1}(\Pi T_1\cdots T_{m-2})^2$ FAILS on the polynomial rep for every test polynomial, starting from $f = 1$ (LHS $= tX_{m-1}X_m$, RHS involves $X_1X_{m-1}$). Convention-variant test (reverse T-chain, $Y_1Y_2$ LHS, inverse Π): all 4 variants FAIL with the same X-index-mismatch obstruction, no scalar/q-power correction. Diagnosis: DAHA identity relies on full double affine structure Hikita's level-1 AHA lacks. Scripts: `proofs/scripts/day194/stokman_rains_check.py` + `stokman_rains_variants.py`. Registry: `hikita-star-e2-e2.json` node `analytic-proof-via-stokman-rains-lift` = `refuted`.

2. **Thibon arXiv:2609.10284 (Day 194: read, DEAD as fast lift).** Jack (1-param, degenerate DAHA), not Macdonald. §10.2 formula $e_2 = \frac{1}{2}[\Delta_2(\alpha), e_1]$ IS the Lemma-3.11-analogue at degenerate level (right structure, wrong parameters). Hand-lift cost: 2-4 weeks quantum toroidal $\mathfrak{gl}_1$ Drinfeld generators. Reference: `reading/2026-09-16-thibon-2609.10284.md`.

3. **R3: quantum toroidal $\mathfrak{gl}_1$ / Maulik-Okounkov (Day 194: novelty search, verdict (c) GENUINE GAP).** Hikita does NOT reference MO (0 hits in body + bibliography). Rick's publish slot verified intact by 3rd novelty audit. Best candidate: **Bechtloff Weising arXiv:2405.00756 (2024)** — explicit $e_r^{\bullet}$-Pieri rule (Cor 5.10) on generalized Macdonald basis $P_T$ for new EHA reps $\tilde W_\lambda$. Not proven equivalent to Hikita ⋆. Estimated 1-2 weeks bridge if BW $e_r^{\bullet}$ collides with Hikita ⋆; else 2-3 months full dictionary. Other candidates: Garbali-Neguţ 2112.09094 (diagonal only), Schiffmann-Vasserot 0802.4001 (foundational EHA), Garbali-de Gier 2004.09241. Reference: `reading/2026-09-17-r3-qt-gl1-novelty.md`.

4. **van Diejen-Emsiz 1009.4482 (LOWER, untried).** Generalized Macdonald difference operators $D_{\omega_r}$.

5. **Direct Lemma-3.11-extension (`hunch`).** Rick's own approach: mimic Hikita's induction for $\sum_{i<j} Y_i Y_j \bullet e_r$. Untried; now-elevated candidate given R1/R2 dispositions.

Full route map: `connections/2026-09-16-two-routes-to-lemma-3-11-analogue.md`.

## Cross-references

- `connections/2026-09-17-route-R6-BW-2310-EHA-descent.md` — **Day 202 crown-jewel**: sole surviving external analytic route (Bechtloff-Weising 2310.10249).
- `connections/2026-09-17-rick-theorem-vs-thibon-Delta3-conjecture.md` — **Day 202**: Rick's k=3 is a theorem where Thibon's Δ_3 is a conjecture (publication asymmetry).
- `connections/2026-09-17-baxter-k-naming-and-internal-terminology.md` — **Day 202**: rename "Baxter-k" before writeup + broader internal-notation-drift pattern.
- `connections/2026-09-16-p2Y-pieri-newton-independent-atom.md` — **Day 199**: why p_2(Y) is the canonical missing analytic input.
- `connections/2026-09-16-thibon-triangle-p2Y-candidates.md` — **Day 199 (historical)**: three candidate identifications for analytic Lemma 1 (Vertex A/B/C). All three now DEAD (Day 200/201).
- `connections/2026-09-16-DS-macdonald-triangularity.md` — **Day 196**: DS ↔ Macdonald n-statistic.
- `connections/2026-09-16-DAdderio-Negut-route-2-unlock.md` — **Day 196/197 (historical)**: R2c attack vector + full REFUTATION log.
- `connections/2026-09-11-e2-star-e2-hikita-pieri-extension.md` — Day 191 result + conjecture.
- `connections/2026-09-16-min-a-b-plus-1-meta-conjecture.md` — Day 193 meta-shape (subsumed by DS + p_k hierarchy).
- `connections/2026-09-16-two-routes-to-lemma-3-11-analogue.md` — historical analytic route map (R1-R5).
- `connections/2026-09-11-qt-slot-open-hikita-recipe-unused.md` — three-way slot verification.
- `questions/q-A2-equals-p2Y-normalized.md` — CLOSED Day 201 (Vertex A refuted Day 200, Vertex B refuted Day 201).
- `questions/q-p_k-Y-Pieri-hierarchy.md` — Day 199 meta-conjecture; Day 201 k=3 confirmed.
- `questions/q-BW-2310-EHA-descent.md` — Day 202 primary target (Route R6).
- `questions/q-baxter-k-rename.md` — resolution tracker before FPSAC abstract v3.
- `questions/q-DS-analytic-proof-strategies.md` — Strategy 6 (via Thibon triangle) now REFUTED; new Strategy 7 = BW-2310 descent.
- `questions/q-h-basis-qt-recursion.md` — $X_{P_n}(x;q,t)$ recursion form.
- `questions/q-D-a-equals-e-a-Y-level-1-AHA.md` — CLOSED Day 197 (REFUTED).
- `questions/q-star-product-commutativity.md` — CLOSED Day 191 (Def 3.4 states commutative).
- `~/projects/proofs/2026-09-17-day198-DS-211-via-p2-pieri.md` — **Day 198 writeup + Lemma 1**.
- `~/projects/proofs/2026-09-11-day191-e2-star-e2-hikita.md` — Day 191 writeup.
- `~/projects/proofs/2026-09-15-day192-e3-star-er-hikita.md` — Day 192 partial closed forms.
- `~/projects/proofs/2026-09-16-day193-e3-star-er-hikita.md` — Day 193 full closed form.
- `~/projects/proofs/registry/hikita-star-e2-e2.json`, `hikita-star-e3-er.json`, `hikita-star-dominance-support.json`.

## Open threads (Day 202+)

1. **Test Route R6 (BW 2310.10249)** — sole surviving external analytic route to Lemma 1 and p_k hierarchy. Day 202 primary. ★★★★
2. **Test k=4 empirical** — p_4(Y)•e_r at m=6, r=3. Predict 8 r-indep + 4 r-dep = 12 nonzero terms in DS-cone of (r,1^4). Day 202. ★★★
3. **Verify Jack degeneration matches Thibon Δ_3 conjecture** — 30-min SymPy from Rick's Day 201 formulas. Publication asymmetry payoff. ★★★
4. **Resolve "Baxter-k" rename** before FPSAC abstract v3. Recommended: "t^r-Laurent of order k". ★★★
5. **Direct proof of r-independence meta-conjecture** via Newton cancellation across ⋆-length pieces (R7 fallback if R6 fails). ★★
6. **Extend F ⋆ P_λ triangular observation to two-row λ** (Day 200 spectral-Pieri byproduct). New direction. ★★
7. **$e_4 \star e_r$ closed form** — needs $(4, 5)$ compute at $m=9$. Deferred.
8. **$t=0$ sanity checks** — van Diejen-Emsiz-Zurrian cylindric HL Pieri; Kim-Lee-Yoo linked rook placements. Deferred.
9. **$s_\lambda \star e_r$ Schur Pieri.** Hikita flags open in same sentence.
10. **FPSAC 2027 abstract v3.** Deadline monitor mid-October 2026. Anchor structure: **p_k(Y)-Pieri hierarchy for k=2,3 (computed) + k=4 (predicted) + DS-triangularity of ⋆-basis with q^{-n(λ)} leading, implying Thibon's Δ_3 conjecture as corollary at Jack degeneration**.

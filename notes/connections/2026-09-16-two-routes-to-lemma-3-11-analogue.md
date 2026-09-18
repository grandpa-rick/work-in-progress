# Connection — Two analytic routes to the Lemma-3.11-analogue

**Date:** 2026-09-16 (Day 193 dream / Browse 143).
**Path bridge:** Path 3 (DAHA / degenerate AHA machinery) → Path 2 (Hikita level-1 AHA).
**Status:** both routes `hunch`; priority order set.

## The analytic gap

Rick has closed forms for $e_2 \star e_r$ (Day 191, `computed` $r \le 4$) and $e_3 \star e_r$ (Day 193, `computed` $r \le 6$). Analytic proof fails via the natural Newton's-identity route: writing $e_2(Y) = \tfrac{1}{2}(e_1(Y)^2 - p_2(Y))$ and using Hikita's Thm 3.12 ($e_1 \star e_r$ Pieri) iteratively hits **tautology $0 = 0$**. Reason: Thm 3.12 + AHA relations alone do *not* determine $p_2(Y) \bullet e_r(X)$ — this action is new data.

The missing ingredient is a **Lemma-3.11-analogue** in Hikita's level-1 AHA: an explicit formula for either $e_k(Y) \bullet e_r(X)$ or $p_k(Y) \bullet e_r(X)$ for $k \ge 2$.

## Route 1: Stokman-Rains Lemma 10 (PRIORITY 1, 30 min)

**Source:** Stokman-Rains arXiv:2307.02385 Lemma 10, bisymmetric Macdonald Pieri.

**Result in DAHA.** For the double affine Hecke algebra of $GL_N$ with generators $\omega, T_i, Y_i$:
$$e_r(Y_1, \ldots, Y_N) = \frac{1}{[N-r]_t![r]_t!}\, S^t_N\, Y_{N-r+1} \cdots Y_N$$
for ALL $r \ge 1$, where $S^t_N$ is the $t$-symmetrizer. Proof uses only the elementary DAHA commutation
$$Y_{N-r+1}\cdots Y_N = t^{-r(r-1)/2}\bigl(\omega\, T_1 \cdots T_{N-r}\bigr)^r.$$

**Lift attempt to Hikita.** Hikita's level-1 AHA of $GL_m$ has $\Pi$ and $T_i$ with $\Pi T_i = T_{i-1}\Pi$, structurally analogous to DAHA's $\omega T_i = T_{i+1}\omega$ (up to direction). Key check:
$$Y_{m-1} Y_m \;\stackrel{?}{=}\; t^{-1}\bigl(\Pi\, T_1 \cdots T_{m-2}\bigr)^2 \quad\text{in Hikita's level-1 AHA.}$$

**If YES:** immediate consequence $e_2(Y_1, \ldots, Y_m) = \frac{1}{[m-2]_t![2]_t!} S^t_m Y_{m-1} Y_m$, and $e_2(Y) \bullet e_r(X)$ closed form falls out with essentially no new work. This gives an **unconditional proof of Day 191's $e_2 \star e_r$ conjecture**.

**If NO:** learn what commutation obstruction breaks it; likely a level-1-specific correction term.

**Cost.** 30 min SymPy at $m = 3, 4$. This is the FIRST thing Day 194 wake should try.

## Route 2: Thibon $\Delta_2(\alpha)$ content operator (BACKUP, 1 hr read + 30 min derivation)

**Source:** Thibon arXiv:2609.10284 "Jack Content Operators and the Deformed $W_{1+\infty}$ Algebra" (posted 2026-09-09, one week old at Day 193).

**Result in degenerate DAHA.** Defines $\Delta_2(\alpha)$, the degree-2 power-sum content operator in the spherical degenerate DAHA ($t \to 1$ limit of Hikita's level-1 AHA). Key relation:
$$\psi_3 = 3\Delta_2(\alpha) + 2(\alpha - 1) E,$$
where $\psi_3$ is a mode of the affine Yangian of $\mathfrak{gl}_1$.

**Lift attempt to Hikita.** If $\Delta_2(\alpha)$ extends to the $(q, t)$ level-1 AHA — i.e., if $p_2(Y) = \sum Y_i^2$ has a closed form in terms of $E, \psi_3$-analogs in Hikita's setting — then Newton's identity
$$e_2(Y) = \frac{1}{2}\bigl(e_1(Y)^2 - p_2(Y)\bigr) = \frac{1}{2}\bigl(\Pi^2 - p_2(Y)\bigr)$$
gives the analytic proof.

**Cost.** 1 hr read §§2-3 of Thibon + 30 min derivation attempt. Priority 2 for Day 194.

## Route 3: van Diejen-Emsiz $D_{\omega_r}$ (LOWER PRIORITY)

**Source:** van Diejen-Emsiz arXiv:1009.4482, "Hecke algebraic Pieri" — all-$e_r$ Pieri in symmetric Macdonald via generalized Macdonald difference operators.

**Lift attempt.** Construct a level-1 $D_{\omega_2}$ analog in Hikita's algebra. Structural, requires significant translation. Not attempted; save for after Routes 1 & 2 exhausted.

## Why these matter (path bridge)

Both routes are **Path 3 → Path 2 lifts**: extract an identity or operator in a well-studied (D)AHA setting, translate to Hikita's level-1 AHA, apply via the $\mathfrak q$-map to close the analytic proof. This is the same seed-bridge structure that produced the Day 191, 192, 193 closed forms via compute — but now the *proof* rather than the *formula*.

**If Route 1 works:** the FPSAC 2027 abstract upgrades from "computed Pieri conjecture" to "**proved** Pieri Theorem" for $e_2 \star e_r$. This is the analytic capstone.

## Route status matrix (end of Day 196 + Browse 144)

| Route | Cost | Priority | Grade |
|-------|------|----------|-------|
| R1: Stokman-Rains commutation lift | 30 min | — | **REFUTED** (Day 194, `checked-sober`; 4 variants FAIL at $m=3,4$) |
| R2a: Thibon $\Delta_2(\alpha)$ lift | 2-4 weeks | ★ | **DEAD as fast lift** (Day 194; Jack-only, needs QT gl_1 lift by hand) |
| R2b: Bechtloff-Weising 2405.00756 collision | 30 min collision test | — | **MISS** (Day 195; BW's $e_r^\bullet$ = ordinary multiplication, not ⋆) |
| **R2c: D'Adderio et al. 2608.14836 Neguţ formula** | **30 min SymPy check** | **★★★★** | **`hunch` — Day 197 primary target** |
| R3: QT $\mathfrak{gl}_1$ level-(a,0) template (2508.19704) | 30 min read | ★★★ | `hunch` (Day 197 secondary) |
| R4: Row-length filtration in Hikita §3–4 | 3+ hr | ★ | `hunch` (structural, deferred) |
| R5: van Diejen-Emsiz $D_{\omega_r}$ | 3+ hr | ★ | `hunch` (deferred) |
| R6: Direct Lemma-3.11 extension (Rick's own) | est 3-6 hr | ★★ | **base case sketched** (Day 196: 5-term expression with $X_1^3$; inductive step remains) |

## Day 196 / Browse 144 disposition: **R2c (D'Adderio) is the primary attack vector**

D'Adderio-Interdonato-Iraci-Pagaria **arXiv:2608.14836** (Aug 2026) gives explicit linear-time formula for Neguţ operators $D_\gamma$ in $A_{q,t}$:
$$D_\gamma F = d_- (-y_1)^{\gamma_1-1} \hat z_1 (-y_1)^{\gamma_2} \hat z_1 \cdots (-y_1)^{\gamma_l} \hat z_1 d_+ F$$
with $D_{(m)} \cdot F = e_m \cdot F$ at $q = 1$. **If $D_{(a)} = e_a(Y)$ in Hikita's level-1 polynomial rep**, then D'Adderio's formula IS the Lemma-3.11-extension for all $a \ge 2$, and Rick's entire Pieri program falls out analytically.

- Base evidence: $D_{(1)}$ matches Thm 3.12 at both $q=1$ and general $q$ up to normalization ingredients.
- A_{q,t} contains Hikita's ⋆-algebra structure (via Griffin-Mellit et al. 2504.06936).
- 30-min SymPy check at $m=3$, $a=2$ decides.

Full analysis: `connections/2026-09-16-DAdderio-Negut-route-2-unlock.md`.

## Day 194 disposition of Route 1 (`checked-sober` REFUTATION)

Test: $Y_{m-1}Y_m \stackrel{?}{=} t^{-1}(\Pi T_1 \cdots T_{m-2})^2$ on the polynomial rep.

- $m=3$, $f=1$: LHS $= tX_2X_3$; RHS $= X_1 X_3$. Difference $X_3(tX_2 - X_1)$.
- $m=4$, $f=1$: LHS $= tX_3X_4$; RHS $= t^2 X_1 X_3$. Difference $-X_3 t(tX_1 - X_4)$.

Convention-variant tests all FAIL:
- Reverse T-chain ($T_{m-2}\cdots T_1\Pi$): identical X-monomial support on LHS, scalar rescale on RHS; still $X_{m-1}X_m$ vs "wrong indices".
- $Y_1 Y_2$ on LHS with forward T-chain: LHS $\sim X_1 X_2 t^{?}$, RHS $\sim X_1 X_3$; index-mismatch.
- $Y_1 Y_2$ + reverse T-chain, and $\Pi$-on-right variants: identical obstructions.

**Structural read.** The obstruction is X-index-mismatch (LHS involves $X_{m-1}X_m$ or $X_1X_2$, RHS involves shifted-index products). No $q$-power or $t$-power correction fixes an $X_2/X_1$ ratio. DAHA identity relies on the full double-affine structure ($Y_i$ invertible, X-Y duality) that Hikita's level-1 AHA lacks — the $q$-central element in level-1 breaks the DAHA symmetry that Stokman-Rains uses.

**Scripts.** `proofs/scripts/day194/stokman_rains_check.py`, `stokman_rains_variants.py`.

**Consequence.** Route 1 slot is now closed. Analytic-proof effort concentrates on Route 2 (Thibon) and the fallback: Rick's own direct Lemma-3.11-extension.

## Cross-references

- `topics/hikita-star-pieri.md` — meta topic file.
- `connections/2026-09-16-min-a-b-plus-1-meta-conjecture.md` — parallel Day 193 crown-jewel connection.
- `questions/q-stokman-rains-lift.md` — priority-1 SymPy check.
- `questions/q-p2-Y-content-operator-lift.md` — priority-2 Thibon read.
- `reading/2026-09-16.md` — Browse 143 with full paper summaries.
- `proofs/registry/hikita-star-e3-er.json` — `hikita-star-e3-er-analytic-proof` node currently `hunch`.

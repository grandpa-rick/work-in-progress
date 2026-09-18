# Day 195 — BW ($e_r^\bullet$ Pieri) vs Rick's ($e_2 \star e_r$) collision test

Date: 2026-09-16.
Source: Bechtloff Weising, arXiv:2405.00756v2, SIGMA 22 (2026), 076.
Reading notes: `/home/agent/projects/reading/2026-09-16-bechtloff-weising-2405.00756.md`.

## Setup

Rick's Day 191 formula for the Hikita star product $\star$ on $\Lambda_{q,t}$ (verified $r \le 4$):
$$e_2 \star e_r = \frac{1}{q^2}\, e_2 e_r + \frac{1-q^{-1}}{q}\,[r]_t\, e_1 e_{r+1} + (1-q^{-1})\,\frac{[r+2]_t}{[2]_t}\!\left([r+1]_t - \frac{t\,[r-1]_t}{q}\right)\! e_{r+2}.$$

Meta-conjecture: $e_a \star e_b$ has exactly $\min(a,b) + 1$ non-zero $e$-basis terms, supported on partitions $(a+b-k, k)$, $k = 0, \ldots, \min(a,b)$. (15-for-15 through Day 193.)

## BW's structure

- Algebra: positive EHA $\mathcal{E}^+$ (Burban–Schiffmann).
- Rep: $\widetilde{W}_\lambda$ for $\lambda \in \mathbb{Y}$. Case $\lambda = \emptyset$: $\widetilde{W}_\emptyset \cong \Lambda_{q,t}$ with the standard action.
- Basis: $\mathfrak{P}_T$ for $T \in \Omega(\lambda)$. For $\lambda = \emptyset$, $\Omega(\emptyset) \leftrightarrow \mathbb{Y}$ and $\mathfrak{P}_\mu = P_\mu(x; q^{-1}, t)$ up to nonzero scalar (Prop 4.19).
- Multiplication: $e_r^\bullet := \lim_n e_r(X_1, \ldots, X_n)$ — this is **ordinary $e_r$-multiplication** in $\Lambda_{q,t}$ when $\lambda = \emptyset$.
- Pieri rule (Cor 5.10 = Thm 5.7 limit):
  $$e_r^\bullet \mathfrak{P}_T = \sum_S \mathfrak{d}_{S,T}^{(r)} \mathfrak{P}_S,$$
  where $S$ is obtained from $T$ by increasing labels in $r$ distinct boxes by $+1$ while preserving RSSYT (strict-column) structure. In the $\lambda=\emptyset$ case this reduces to $\nu / \mu$ being a horizontal $r$-strip — the classical Macdonald Pieri rule.

## Collision test

### Test A — Do the products coincide?

**BW's $e_r^\bullet$** on $\widetilde{W}_\emptyset = \Lambda_{q,t}$ = ordinary multiplication by $e_r$.

**Rick's $\star$** (Hikita 2503.23597) is a *different* associative product $\star$ on $\Lambda_{q,t}$, coming from the quantum toroidal action / Hilbert-scheme geometry. It is **not** ordinary multiplication.

In particular: $e_2 \cdot e_r$ (ordinary product) expanded in $e$-basis is $e_2 e_r$ — a single term (it IS a monomial in $e$'s). Rick's $e_2 \star e_r$ has THREE terms. So the products are not equal.

**A = FAIL** (products differ structurally). No coincidence of Hikita $\star$ with BW $e_r^\bullet$.

### Test B — Does BW's coefficient shape match Rick's?

Rick's coefficients (Day 191) live in $\mathbb{Q}(q, t)$ with prefactors like $\frac{1}{q^2}, \frac{1-q^{-1}}{q}, (1-q^{-1})$, and $t$-integers $[r]_t, [r+2]_t/[2]_t$, etc.

BW's $d_{S,T}^{(r)}$ (Thm 5.7) are big products over inversion sets $\mathrm{Inv}(\tau)$ and $\mathrm{Inv}(\Psi^r(\tau))$ summed over $\tau \in \mathrm{PSYT}_{\ge 0}(\lambda; T)$. Specialized to $\lambda = \emptyset$ these must reduce to the standard Macdonald Pieri coefficients (Macdonald book, Ch VI, (6.24)) — ratios like $\prod_{s \in R_{\lambda/\mu}} \frac{b_\lambda(s)}{b_\mu(s)}$ where $b_\mu(s) = (1 - q^{a_\mu(s)} t^{l_\mu(s)+1})/(1 - q^{a_\mu(s)+1} t^{l_\mu(s)})$.

These are **not the shape** of Rick's coefficients. Rick's live on a different basis (elementary $e$-basis, not Macdonald $P$-basis) and have a much simpler $(q, t)$-rational structure with the $\min(a,b)+1$-term collapse.

**B = FAIL** (coefficient shape mismatch).

### Test C — Is there a specialization/isomorphism?

Could there be a hidden intertwiner between $(\Lambda_{q,t}, \star)$ (Hikita) and $(\Lambda_{q,t}, \cdot)$ (BW)? This is asking whether Hikita's star-product algebra structure is isomorphic to the ordinary commutative algebra $(\Lambda_{q,t}, \cdot)$ as $\mathbb{Q}(q,t)$-algebras.

They cannot be isomorphic in any natural way that preserves the $e_a$ generators, because:
- $(\Lambda_{q,t}, \cdot)$ has $e_a \cdot e_b = e_b \cdot e_a$ trivially, and $e_a \cdot e_b$ has a **single** $e$-monomial support.
- $(\Lambda_{q,t}, \star)$: Rick has verified $e_a \star e_b$ has $\min(a,b)+1$ $e$-basis terms.

A change of basis / rescaling would either preserve the $\min(a,b)+1$ count or destroy it — no way to interpolate. So the two products encode different associative-algebra data on $\Lambda_{q,t}$.

**C = FAIL** (no plausible specialization).

## Verdict: MISS (clear)

BW's paper does NOT give a route to Rick's meta-conjecture. Specifically:
- **Wrong product.** BW's $e_r^\bullet$ is ordinary $e_r$-multiplication on $\Lambda_{q,t}$ (in the standard rep). Hikita's $\star$ is a different associative product — it's what encodes the Hikita quantum-toroidal action.
- **Wrong basis.** BW works with the Macdonald $P_\mu$ basis; Rick works with the elementary $e_\lambda$ basis. The Day 191 collapse ($\min(a,b)+1$ terms) is an $e$-basis phenomenon.
- **Wrong shape of coefficients.** BW's $d_{S,T}^{(r)}$ are Macdonald arm/leg products; Rick's coefficients are simple $(q, t)$-rationals.

## What BW gives instead (may still be useful)

1. **A different $\mathcal{E}^+$-module family.** For $\lambda \ne \emptyset$, $\widetilde{W}_\lambda$ is a *new* rep of $\mathcal{E}^+$ with distinguished basis $\mathfrak{P}_T$. Hikita's $\Lambda_{q,t}$ with $\star$-product might correspond (under some correspondence) to a *specific* $\widetilde{W}_\lambda$ — but this would need to be checked and is not obvious.
2. **Combinatorial framework for stability limits.** BW's PSYT/RSSYT/$\Omega(\lambda)$ combinatorics + Dunkl–Luque vv Macdonald polynomials give a controlled way to take $n \to \infty$ limits of DAHA reps. If a Hikita-style $\star$-product admits a similar $n \to \infty$ presentation, BW's Prop 3.14 stability tools may be adaptable.
3. **Non-vanishing lemmas (§5.3).** The technique of proving $\mathfrak{d}_{T', T}^{(1)} \ne 0$ using $q \to \infty, t \to 0$ limits (Thm 5.16 proof) is a nice example of how to prove non-vanishing of coefficients in a Pieri-type rule. Could be a template for proving Rick's coefficients are nonzero.
4. **A concrete distinct construction.** Confirms that "positive EHA action on $\Lambda_{q,t}$" is *one* structure, and Hikita's $\star$ is a *different* algebra structure. So Rick's result is genuinely orthogonal to the Dunkl–Luque / BW / Schiffmann–Vasserot line.

## Novelty confirmation

BW's paper (SIGMA 22 (2026), submitted Jul 2025) does **not** address:
- Hikita's $\star$-product on $\Lambda_{q,t}$.
- $e_a \star e_b$ for $a \ge 2$.
- Any $\min(a,b)+1$-term collapse.

So Rick's Day 191 formula and meta-conjecture remain novel and are not subsumed by BW. This is Rick's third failed-collision confirmation (after Stokman–Rains Day 194 R1 and Thibon Day 194 R2).

## Next-step recommendation

**Recommendation:** Drop the BW-bridge route. The next best analytic candidates already flagged in memory are (i) direct quantum-toroidal / spherical-DAHA lifts of Hikita level-1 with attention to the X/Y-duality break template (feedback_daha_x_y_duality_breaks_at_level_1.md), and (ii) test whether the Carlsson–Mellit $\mathbb{A}_{q,t}$ algebra (BW's reference [5–6]) — which mediates between DAHA and $\Lambda_{q,t}$ in the shuffle theorem proof — might give a hidden intertwiner between Hikita's $\star$ and ordinary multiplication after a nontrivial change of basis. If neither pans out, retreat to `/expository` on the exact definition of Hikita's $\star$-product and hunt for a direct combinatorial identity.

## Verdict summary line

**MISS.** BW gives ordinary $e_r$-Pieri (horizontal-strip sum) on Macdonald basis in $\Lambda_{q,t}$. Hikita's $\star$ is a fundamentally different associative product; Rick's Day 191 formula is on the $e$-basis with $\min(a,b)+1$-term collapse. No bridge. BW's stability & non-vanishing techniques may serve as templates but do not touch the substance.

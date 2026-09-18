# Day 201 — Vertex B (Jack-limit shape check) REFUTED

**Date:** 2026-09-17.
**Author:** Rick.

## Problem

PROVE.md Day 201 primary target: verify that in the Jack limit
$q = t^{\alpha}$, $t \to 1$ (equivalently Sekiguchi's $\alpha = \theta$
degeneration), Rick's operator $p_2(Y)$ from Hikita's level-1 AHA
representation matches Thibon 2609.10284's degenerate-DAHA operator
$\Delta_2(\alpha)$ (§7, eq 136) acting on Jack polynomials, in the
stable $N \to \infty$ limit.

Target identity was
$$
\lim_{t \to 1,\, q = t^{\theta}} p_2(Y) \bullet e_r(X)
    \; \stackrel{?}{=} \; \Delta_2(\alpha = \theta) \bullet e_r^{\text{Jack}}(X).
$$

## Diagnosis

**Two independent structural obstructions, each fatal.**

### Obstruction 1 — degree mismatch

Rick's $p_2(Y)$ (verified directly, e.g. m=3 SymPy in
`scripts/day201/vertex_B_refutation.py`) satisfies
$Y_i \cdot 1 = X_i$, so each $Y_i$ raises the $X$-degree by $1$ and
$p_2(Y) = \sum Y_i^2$ raises degree by $2$.

Thibon's $\Delta_2(\alpha)$ (eq 136):
$$
\Delta_2(\alpha) = \alpha \sum_{i,j \ge 1} p_i p_j D_{i+j}
                 + \sum_{i,j \ge 1} p_{i+j} D_i D_j
                 + (\alpha - 1) \sum_{k \ge 1} (k-1) p_k D_k
$$
is **degree $0$** (each summand has matching $p$-degree and $D$-degree).

Hence the two operators live on different graded pieces:
$p_2(Y) \bullet e_r$ has degree $r+2$, while $\Delta_2(\alpha) \bullet e_r$
has degree $r$. No direct identification is possible.

### Obstruction 2 — spectral eigenvalue mismatch

This obstruction applies to the *reformulated* Vertex B via the
Hikita $\star$-product: interpret $p_2(Y)$ via the symmetric function
$F := (p_2(Y) \cdot 1)(X)$ and ask whether $F \star (-)$ (which is
degree-preserving in the $\star$-product) has Jack eigenvalue matching
$\Delta_2(\alpha)$.

Direct SymPy at $m=3$:
$$
F(X) \;=\; \frac{p_2(X) - (t-1)(q+1)\, e_2(X)}{q}.
$$

The $\star$-eigenvalue of $F \star (-)$ on $P_{\lambda}$ is $F$ evaluated
at the Macdonald spectrum
$\mathrm{spec}_{\lambda,i} = q^{\lambda_i} t^{-(i-1)}$.

Substituting $q = t^{\alpha}$, $t = 1+\varepsilon$, we get
$\mathrm{spec}_{\lambda,i} = 1 + \varepsilon(\alpha \lambda_i - (i-1)) + O(\varepsilon^2)$.

The $\varepsilon^0$ term of $F(\mathrm{spec}_\lambda)$ is $m$ (independent of $\lambda$).
The $\varepsilon^1$ term, computed at $m=5$:

| $\lambda$ | $F(\mathrm{spec}_\lambda) \vert_{\varepsilon^1}$ | Thibon $2 C_1^{(\alpha)}(\lambda)$ |
|-----------|--------|--------|
| ()        | $-5\alpha - 40$ | $0$ |
| (1)       | $-3\alpha - 40$ | $0$ |
| (2)       | $-\alpha - 40$  | $2\alpha$ |
| (1,1)     | $-\alpha - 40$  | $-2$      |
| (3)       | $\alpha - 40$   | $6\alpha$ |
| (2,1)     | $\alpha - 40$   | $2\alpha - 2$ |
| (1,1,1)   | $\alpha - 40$   | $-6$      |
| (2,2)     | $3\alpha - 40$  | $4\alpha - 4$ |
| (3,1)     | $3\alpha - 40$  | $6\alpha - 2$ |

**Key observation.** $F(\mathrm{spec}_\lambda)|_{\varepsilon^1}$ depends
only on $|\lambda|$. Explicitly:
$$
F(\mathrm{spec}_\lambda)\bigl|_{\varepsilon^1}
  \;=\; \alpha(2|\lambda| - m) - 2m + \text{const}(m)
  \;=\; 2\alpha|\lambda| + C_m
$$
where $C_m$ is a constant depending only on $m$. So $F$'s $\varepsilon^1$
information about $\lambda$ is one-dimensional (only $|\lambda|$).

Thibon's $2 C_1^{(\alpha)}(\lambda) = \alpha(\sum \lambda_i^2 - |\lambda|)
- 2 n(\lambda)$ is *quadratic* in the parts $\lambda_i$, distinguishing
partitions of the same size (compare $(2)$ vs $(1,1)$: same $|\lambda|=2$
but different eigenvalues $2\alpha$ vs $-2$).

**No constant rescaling of $F$ or its spectrum can convert a linear-in-$\lambda_i$
quantity into a quadratic-in-$\lambda_i$ one.** Vertex B fails at
$\varepsilon^1$ across every shape of the same size.

## Verification (SymPy)

`scripts/day201/vertex_B_refutation.py` computes both sides at $m=5$ for
$\lambda \in \{\emptyset, (1), (2), (1,1), (3), (2,1), (1,1,1), (2,2), (3,1)\}$.
Log preserved at `.log` sibling.

## Consequence

Vertex B — the sole surviving R5 analytic route for Lemma 1 —
is dead. Rick's $p_2(Y)$ operator in Hikita's level-1 rep does not
degenerate to Thibon's $\Delta_2(\alpha)$ under any (Hall-adjoint /
spectral / naive) reading. The two obstructions (degree, spectral)
are independent — either alone suffices.

**Route R5 (all three Thibon vertices A, B, C) essentially exhausted**:
- Vertex A (A$^{(2)}$, Nazarov–Sklyanin): refuted Day 200.
- Vertex B (Jack, Thibon 2609.10284): refuted Day 201 (this note).
- Vertex C (shuffle $\Delta_2$): dormant. No a priori reason to expect it fares
  better since it has the same degree-0 signature.

Lemma 1 (`p2Y-pieri-lemma`) remains at grade `computed` (r=2..6 at m=8).
An analytic upgrade requires a **different** approach — most likely a
direct AHA argument (Rule 11 Room 5 style, but generalized).

## Registry

Add node `p2Y-pieri-lemma-jack-limit-refuted` under `p2Y-pieri-lemma`
in `hikita-star-dominance-support.json`, grade `refuted`.

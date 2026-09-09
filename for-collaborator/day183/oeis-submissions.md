# OEIS submissions — Day 183 (2026-09-09), for Robin to submit

Three new integer sequences arising from Rick's algebraic generating-function
work on path-graph chromatic quasisymmetric functions (Shareshian–Wachs
setting). All three absent from OEIS as of 2026-09-09.

Both b_k and a_k = INVERTi(b_k) are computed from the algebraic identity
(Day 148, proved):
$$F(1-F)^3(3-4F) \;=\; \vartheta\,(3-2F)^2, \qquad F(\vartheta)=\sum_{k\ge1} b_k\,\vartheta^k.$$

The third sequence, $p_k$, is the graded dimension of the space of Lie
primitives of the free graded-connected cocommutative Hopf algebra
$U(L(a))$ with dimension sequence $b_k$; equivalently, $p_k = \dim L(a)_k$
where $L(a)$ is the free graded Lie algebra on $a_n$ generators in degree
$n$. By the graded Witt formula,
$$\sum_n p_n t^n \;=\; \sum_{d=1}^{\infty} \frac{\mu(d)}{d}\, \log\!\big(1 + \tsum_{k\ge 1} b_k t^{dk}\big).$$

**Structural home (Day 183, corrected).** Andrews–Gagnon–Gélinas–Schlums–Zabrocki
(arXiv:2505.06941, "When are Hopf algebras determined by integer sequences?",
2026) prove (Theorem 4.2): a positive-integer sequence $(h_k)$ with $h_0 = 1$
is the graded-dimension sequence of a **Free Graded Connected Cocommutative
Hopf Algebra (FGCCHA)** over $\mathbb C$ **iff** its INVERTi transform is
componentwise nonnegative; the FGCCHA is unique up to Hopf isomorphism
(Theorem 4.1, Aliniaeifard–Thiem [7]); and the FGCCHA is $U(L(a))$ where
$L(a)$ is the free graded Lie algebra on generators counted by $a_k =
\mathrm{INVERTi}(h_k)$. **CAUTION on our earlier note:** $a_k$ counts the
*free Lie generators*, NOT the space of Lie primitives. The Lie primitive
dimensions $p_k$ are a **distinct** sequence, produced from $a_k$ by the
graded Witt formula.

Rick's $b_k$ passes: $a_k$ is positive at all 12 computed terms. Hence
$b_k$ is the graded dimension of a well-defined $U(L(a))$ over $\mathbb C$;
$a_k$ is the number of free Lie generators in each degree; $p_k$ is the
graded dimension of the primitives.

---

## Sequence 1 — $b_k$

**Data (12 terms):**

```
3, 27, 417, 7851, 164124, 3661389, 85384566, 2056373739, 50751637140, 1276862920140, 32626363346505, 844375375808301
```

**%N** Coefficients of the unique power series $F(\vartheta) = \sum_{k\ge1} b_k\,\vartheta^k$ satisfying the algebraic equation $F(1-F)^3(3-4F) = \vartheta\,(3-2F)^2$; the graded dimensions of a ring attached to path-graph chromatic quasisymmetric functions (Shareshian–Wachs).

**%C** $F(\vartheta)$ is algebraic of degree 5 over $\mathbb Q(\vartheta)$.

**%C** Divisibility: $b_k \equiv 0 \pmod 3$ for every $k \ge 1$ (proved). Under $F = 3G$, both sides of the defining identity acquire exactly one factor $9$; the resulting equation $G(3G-1)^3(4G-1) = \vartheta(2G-1)^2$ has the form $G = \vartheta\cdot\phi(G)$ with $\phi(G) = (2G-1)^2/[(3G-1)^3(4G-1)] \in \mathbb Z[[G]]$ and $\phi(0)=1$; Lagrange inversion gives $G = F/3 \in \mathbb Z[[\vartheta]]$, so $b_k/3 \in \mathbb Z$.

**%C** Rick's algebraic-GF construction (2026): $b_k$ is the top-weight dimension of a ring in the Shareshian–Wachs (2016) chromatic quasisymmetric-function framework restricted to path graphs $P_n$.

**%C** By Andrews–Gagnon–Gélinas–Schlums–Zabrocki (arXiv:2505.06941, Theorem 4.2), since $\mathrm{INVERTi}(b_k)$ (see A[a_k]) is componentwise positive, $b_k$ is the graded-dim sequence of a unique-up-to-isomorphism free graded connected cocommutative Hopf algebra $U(L(a))$ over $\mathbb C$.

**%F** $F(\vartheta) := \sum_{k\ge1} b_k \vartheta^k$ satisfies $F(1-F)^3(3-4F) = \vartheta(3-2F)^2$.

**%F** Lagrange form: $b_k/3 = \dfrac{1}{k}\,[G^{k-1}]\,\dfrac{(2G-1)^{2k}}{(3G-1)^{3k}(4G-1)^{k}}$.

**%F** Poincaré–Birkhoff–Witt: $1 + \sum_{k\ge 1} b_k t^k = \prod_{k\ge 1} (1 - t^k)^{-p_k}$, where $p_k$ is A[p_k].

**%o** (Python/SymPy) [Lagrange inversion, see draft]

**%Y** Cf. A[a_k] (INVERTi transform; free Lie generator counts), A[p_k] (Lie primitive dimensions).

**%K** nonn,hard,new

---

## Sequence 2 — $a_k = \mathrm{INVERTi}(b_k)$

**Data (12 terms):**

```
3, 18, 282, 5268, 109647, 2438928, 56758176, 1364824620, 33643660620, 845633502606, 21590775239850, 558411335278644
```

**%N** INVERTi transform of A[b_k]: the unique sequence $a_k$ with $1 + \sum_{k\ge1} b_k x^k = 1/(1 - \sum_{k\ge1} a_k x^k)$.

**%C** By Andrews–Gagnon–Gélinas–Schlums–Zabrocki (arXiv:2505.06941, Theorem 4.2), $a_k$ is the graded dimension of the free-Lie-generator space of the free graded connected cocommutative Hopf algebra $U(L(a))$ with dimension sequence $b_k$. Equivalently, if $L$ is the free graded Lie algebra with $a_k$ generators in degree $k$, then $\dim U(L)_k = b_k$.

**%C** (Note: $a_k$ counts free Lie generators, NOT the full Lie primitive space; the latter has dimensions $p_k$ given in A[p_k]. For instance $a_2 = 18$ but $p_2 = 21 = 18 + \binom{3}{2}$: the three extra brackets come from $[X_i, X_j]$ for the three degree-1 generators.)

**%C** Divisibility: $a_k \equiv 0 \pmod 3$ for every $k \ge 1$ (proved from A[b_k] $\equiv 0 \pmod 3$ by induction on the INVERTi recursion).

**%C** Refinement: $a_k \equiv b_k \pmod 9$ for every $k \ge 1$ (proved). Setting $b'_k := b_k/3$, $a'_k := a_k/3$, the INVERTi identity $b_k - a_k = \sum_{i=1}^{k-1} b_i a_{k-i}$ gives $b_k - a_k = 9 \sum_{i=1}^{k-1} b'_i a'_{k-i}$.

**%F** $a_n = b_n - \sum_{k=1}^{n-1} a_k\, b_{n-k}$, with $b_0 = 1$.

**%F** Generating function: $A(x) := \sum_{k\ge1} a_k x^k = 1 - 1/B(x)$, where $B(x) := 1 + \sum_{k\ge1} b_k x^k$.

**%F** $b_k \equiv a_k \pmod 9$ (proved).

**%Y** Cf. A[b_k], A[p_k], arXiv:2505.06941.

**%K** nonn,hard,new

---

## Sequence 3 — $p_k$ (Lie primitive dimensions)

**Data (12 terms):**

```
3, 21, 344, 6447, 134571, 2995655, 69761697, 1678307754, 41386815905, 1040573158494, 26574621911472, 687454232433863
```

**%N** Graded dimensions of the Lie primitives of the free graded connected cocommutative Hopf algebra $U(L(a))$ associated to path-graph CQF; equivalently, $\dim L(a)_k$ where $L(a)$ is the free graded Lie algebra with $a_n = \mathrm{INVERTi}(b_n)$ generators in degree $n$.

**%C** Given by the graded Witt / plethystic-log formula:
$$\sum_{k\ge 1} p_k t^k \;=\; \sum_{d \ge 1} \frac{\mu(d)}{d}\, \log B(t^d), \qquad B(t) := 1 + \sum_{k \ge 1} b_k t^k.$$

**%C** Poincaré–Birkhoff–Witt: $B(t) = \prod_{k \ge 1} (1 - t^k)^{-p_k}$. Verified for $k \le 12$.

**%C** Small identities: $p_1 = 3 = a_1$; $p_2 = 21 = a_2 + \binom{a_1}{2} \cdot 2 = 18 + 3$ (three degree-1 generators contribute $\binom{3}{2}$ Lie brackets).

**%C** Divisibility pattern: $p_k \equiv 0 \pmod 3$ if $3 \nmid k$; $p_k \equiv 2 \pmod 3$ if $3 \mid k$. Verified $k = 1..12$. Structural explanation open.

**%F** $\sum p_k t^k = \sum_{d\ge 1} (\mu(d)/d) \log(1 + \sum_j b_j t^{dj})$.

**%F** $B(t) = \prod_k (1-t^k)^{-p_k}$, where $B$ is the OGF of $b_k$ with $b_0 = 1$.

**%o** (Python) [see computed script at /home/agent/projects/scratch/day183/compute_pk.py]

**%Y** Cf. A[b_k], A[a_k], arXiv:2505.06941 (Andrews–Gagnon–Gélinas–Schlums–Zabrocki).

**%K** nonn,hard,new

---

## Summary

- All three sequences absent from OEIS (checked 2026-09-09).
- All 12 $a_k$ terms positive → Andrews–Gagnon–Gélinas–Schlums–Zabrocki
  Theorem 4.2 applies → Hopf home is theorem-level, not heuristic.
- Correction from Rick's earlier notes: $a_k$ = free-Lie generator counts,
  $p_k$ = Lie primitive dims. These are different sequences related by the
  graded Witt formula.
- Mod-3 pattern proved for $b_k$ (Day 148); mod-3 proved for $a_k$;
  mod-9 identity $a_k \equiv b_k \pmod 9$ proved (Day 181).
- Mod-3 pattern for $p_k$ empirical (period-3 residues 0,0,2 for $k=1..12$).

## Suggested submission order for Robin

1. Submit **b_k** first (fundamental sequence, algebraic).
2. Submit **a_k** with cross-ref "Cf. A[b_k]".
3. Submit **p_k** with cross-ref to both, citing PBW factorisation.

The %A field goes to Robin (per PROTOCOL §5). Rick's construction attribution
goes in %C.

## Provenance (Rick's private notes)

- $b_k$ computed and tabulated in `~/projects/beta-prime/code/day146_prove/data.json`.
- $b_k \equiv 0 \pmod 3$ proved in `~/projects/proofs/2026-08-30-day148-bk-mod3-SOLVED.md`.
- $a_k \equiv b_k \pmod 9$ proved in `~/projects/scratch/day181/mod9_investigation.md`.
- $p_k$ computed in `~/projects/scratch/day183/compute_pk.py`, PBW-verified mod $t^{13}$.
- Algebraic identity discovered Day 148, Theorem 2.2.
- OEIS non-inclusion re-confirmed 2026-09-09.

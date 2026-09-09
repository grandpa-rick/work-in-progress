# OEIS submissions — Day 181

Two new integer sequences arising from Rick's algebraic generating-function
work on path-graph chromatic quasisymmetric functions (Shareshian–Wachs
setting). Neither sequence appears in OEIS as of 2026-09-08 (Browse 135).

Both are computed from the algebraic identity (Day 148, proved):
$$F(1-F)^3(3-4F) \;=\; \vartheta\,(3-2F)^2, \qquad F(\vartheta)=\sum_{k\ge1} b_k\,\vartheta^k.$$

12 terms of $b_k$ are on hand from `~/projects/beta-prime/code/day146_prove/data.json`.
From those, 12 terms of $a_k = \mathrm{INVERTi}(b_k)$ are computed below.

---

## Sequence 1 — $b_k$

**Data (12 terms):**

```
3, 27, 417, 7851, 164124, 3661389, 85384566, 2056373739, 50751637140, 1276862920140, 32626363346505, 844375375808301
```

**%N** Coefficients of the unique power series $F(\vartheta) = \sum_{k\ge1} b_k\,\vartheta^k$ satisfying the algebraic equation $F(1-F)^3(3-4F) = \vartheta\,(3-2F)^2$; equivalently, the graded dimension sequence of a ring attached to path-graph chromatic quasisymmetric functions (Shareshian–Wachs).

**%C** $F(\vartheta)$ is algebraic of degree 5 over $\mathbb Q(\vartheta)$.

**%C** Divisibility: $b_k \equiv 0 \pmod 3$ for every $k \ge 1$ (proved). Under the substitution $F = 3G$, both sides of the defining identity acquire exactly one factor $9$; the resulting equation $G(3G-1)^3(4G-1) = \vartheta(2G-1)^2$ has the form $G = \vartheta\cdot\phi(G)$ with $\phi(G) = (2G-1)^2/[(3G-1)^3(4G-1)] \in \mathbb Z[[G]]$ and $\phi(0)=1$, so by Lagrange inversion $G = F/3 \in \mathbb Z[[\vartheta]]$. Hence $b_k/3 = [\vartheta^k]G \in \mathbb Z$.

**%C** Rick's algebraic-GF construction (2026): $b_k$ is the top-weight dimension of a ring in the Shareshian–Wachs (2016) chromatic quasisymmetric-function framework restricted to path graphs $P_n$. Related to Hikita's independence-of-$q$ theorem for path-graph $e$-coefficients.

**%C** Empirically not a linear-recurrence sequence; the algebraic curve has degree 5.

**%F** $F(\vartheta) := \sum_{k\ge1} b_k \vartheta^k$ satisfies $F(1-F)^3(3-4F) = \vartheta(3-2F)^2$.

**%F** Lagrange form: $b_k/3 = \dfrac{1}{k}\,[G^{k-1}]\,\dfrac{(2G-1)^{2k}}{(3G-1)^{3k}(4G-1)^{k}}$.

**%o** (Python/SymPy)
```python
from sympy import symbols, series, solve, Poly, Rational, Symbol
from sympy import Function, Symbol, expand, series, O
# Compute via Lagrange inversion for G = F/3
# G = theta * phi(G), phi(G) = (2G-1)^2 / [(3G-1)^3 (4G-1)]
# Compute [theta^k] G via power series
from sympy import symbols, series, expand, O
th = symbols('th')
G = 0
for _ in range(15):
    # (2G-1)^2 / ((3G-1)^3 (4G-1))
    num = expand((2*G - 1)**2)
    den = expand((3*G - 1)**3 * (4*G - 1))
    # series division mod th^N
    N = 15
    from sympy import series as ser
    phi = ser(num/den, th, 0, N).removeO()
    G_new = expand(th * phi + O(th**N)).removeO()
    if G_new == G:
        break
    G = G_new
F = 3 * G
print([F.coeff(th, k) for k in range(1, 13)])
```

**%o** (Sage)
```sage
R.<th> = PowerSeriesRing(QQ, default_prec=15)
# Iterate the Lagrange fixed point for G = F/3
G = R(0)
for _ in range(20):
    phi = (2*G - 1)^2 / ((3*G - 1)^3 * (4*G - 1))
    G = (th * phi).add_bigoh(15)
F = 3 * G
print([F[k] for k in range(1, 13)])
```

**%Y** Cf. A(a_k) [companion INVERTi transform sequence, submitted alongside].

**%K** nonn,hard,new

**%A** Robin Langer (langer.robin(AT)gmail.com), Sep XX 2026

---

## Sequence 2 — $a_k = \mathrm{INVERTi}(b_k)$

**Data (12 terms):**

```
3, 18, 282, 5268, 109647, 2438928, 56758176, 1364824620, 33643660620, 845633502606, 21590775239850, 558411335278644
```

**%N** INVERTi transform of A(b_k): the unique sequence $a_k$ with $1 + \sum_{k\ge1} b_k x^k = 1/(1 - \sum_{k\ge1} a_k x^k)$, where $b_k$ is A(b_k).

**%C** Conjectural interpretation: predicted number of degree-$k$ primitive generators of a graded connected free noncommutative-cocommutative Hopf algebra whose graded-dimension sequence is $b_k$. If such a Hopf structure exists on the ring underlying $b_k$, then $a_k \ge 0$ and the multiplicities of primitive generators in degree $k$ are $a_k$ (Andrews–Gagnon–Gélinas–Schlums–Zabrocki, arXiv:2505.06941, "Zabrocki test").

**%C** All 12 computed terms are positive integers, consistent with the graded-free NC-cocommutative Hopf structure hypothesis.

**%C** Divisibility: empirically $a_k \equiv 0 \pmod 3$ for every $k$ in the range computed ($1 \le k \le 12$). Structural explanation open; likely inherited from $b_k \equiv 0 \pmod 3$ via the INVERTi recursion, but not established here.

**%C** Refined data: $a_k \bmod 9 = 3,0,3,3,0,0,0,0,3,3,0,3$ for $k=1..12$ (no obvious period).

**%F** $a_n = b_n - \sum_{k=1}^{n-1} a_k\, b_{n-k}$, with $b_0 = 1$.

**%F** Generating function: $A(x) := \sum_{k\ge1} a_k x^k = 1 - 1/B(x)$, where $B(x) := 1 + \sum_{k\ge1} b_k x^k$.

**%o** (Python)
```python
b = [1, 3, 27, 417, 7851, 164124, 3661389, 85384566, 2056373739,
     50751637140, 1276862920140, 32626363346505, 844375375808301]
a = [0]
for n in range(1, len(b)):
    a.append(b[n] - sum(a[k] * b[n-k] for k in range(1, n)))
print(a[1:])
```

**%o** (Sage)
```sage
R.<x> = PowerSeriesRing(QQ, default_prec=13)
b_coeffs = [1, 3, 27, 417, 7851, 164124, 3661389, 85384566, 2056373739,
            50751637140, 1276862920140, 32626363346505, 844375375808301]
B = R(b_coeffs)
A = 1 - 1/B
print([A[k] for k in range(1, 13)])
```

**%Y** Cf. A(b_k) [inverse INVERT gives back b_k].

**%Y** Cf. arXiv:2505.06941 (Andrews–Gagnon–Gélinas–Schlums–Zabrocki) for the INVERTi test for free NC-cocommutative Hopf algebras.

**%K** nonn,hard,new

**%A** Robin Langer (langer.robin(AT)gmail.com), Sep XX 2026

---

## Summary of numerical verification

- Given check values: $a_1..a_5 = 3, 18, 282, 5268, 109647$ — all match.
- Extended to 12 terms using $b_1..b_{12}$ (already on disk from Day 146).
- Mod-3 pattern for $a_k$ holds through $k=12$ (no break).
- Mod-9 pattern for $a_k$: $\{3,0,3,3,0,0,0,0,3,3,0,3\}$; no obvious period.
- Mod-9 pattern for $b_k$: $\{3,0,3,3,0,0,0,0,3,3,0,3\}$ — identical to $a_k$
  through $k=12$. Suggests $a_k \equiv b_k \pmod 9$ (worth checking further).

## Provenance

- $b_k$ computed and tabulated in `~/projects/beta-prime/code/day146_prove/data.json`
  (12 terms, from symbolic $\Psi_b$ recursion via `core.py`).
- $b_k \equiv 0 \pmod 3$ proved in `~/projects/proofs/2026-08-30-day148-bk-mod3-SOLVED.md`.
- Algebraic identity discovered Day 148, Theorem 2.2.
- OEIS non-inclusion checked in Browse 135, `~/projects/memory/reading/`
  (2026-09-08).

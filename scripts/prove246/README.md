# prove246 — scripts behind §8 of `rick-research/proofs/2026-10-18-q488-d-infinity.md` (Theorem D∞)

Engine: `fast.py` (exact Laurent series in ε over `Fraction`s, with w_b = 1/(x(1−ε)), so the residue on
H_ab = {w_a w_b = 1} is the ε^{−1} coefficient; random rational points, two seeds 3 and 11).
Subsets of [N] are increasing tuples. `geq(M, I)` is **top-aligned dominance**: |M| ≥ |I| and the top |I|
elements of M dominate I entrywise.

## The census domain and "allowed"

**Domain.** All (λ, μ, ν; a, b) with λ, μ, ν ⊆ [N], (λ, μ) unordered (λ before μ in the `subsets(N)` order),
ν ≥ λ and ν ≥ μ (top-aligned dominance), jump |ν| − max(|λ|, |μ|) ≥ 2, and 1 ≤ a < b ≤ N.

**Definition (one sentence).** A case (λ, μ, ν; a, b) of the domain is *allowed* iff
E(λ, μ, ν; ab) := {M″ ⊆ [N]∖{a, b} : M″ ∪ {a, b} ≤ ν, M″ ≥ λ, M″ ≥ μ} is non-empty
(all comparisons top-aligned dominance); this is the right-hand side of Theorem D∞.

As coded (identical in `census.py` and `wcensus.py`):

```python
def allowed(l,m,nu,a,b):
    rest=[c for c in range(1,N+1) if c not in (a,b)]
    for k in range(len(rest)+1):
        for Mpp in itertools.combinations(rest,k):
            M=tuple(sorted(Mpp+(a,b)))
            if geq(nu,M) and geq(Mpp,l) and geq(Mpp,m): return True
    return False
```

Domain sizes (triples × pairs = cases; allowed / not allowed):

| N | triples | pairs | cases | allowed | not allowed |
|---|---|---|---|---|---|
| 3 | 13 | 3 | 39 | 25 | **14** |
| 4 | 127 | 6 | 762 | **408** | **354** |
| 5 | 1136 | 10 | 11360 | **5578** | 5782 |

## Which script prints which §8 number

| §8 number | script | output |
|---|---|---|
| N=3: zeros = non-allowed, **14/39** | `python3 census.py 3` (~1 s) | `census_N3.txt`: pair (1,3) 5 + pair (2,3) 9 = 14 cases `(False, True, True)`; all 25 allowed are `(True, False, False)` |
| N=4: H12 127/127 nonzero at ∞, single q=0 zero; other pairs **354** non-allowed, Z0 = Z∞ = {E = ∅}; **408/408** allowed nonzero | `python3 census.py 4` (~30 s) | `census_N4.txt`: 25+70+64+88+107 = 354; allowed 127+102+57+63+39+20 = 408 |
| Weighted census, N=4, 408 allowed | e.g. `python3 wcensus.py 4 -10,-8,-6,-4,-2,0,2,4` (any comma list of c) | prints `wt <0|1> c <c> allowed 408 zeros ...` |
| N=5, **5578** allowed; F-system (wt 1) 0 zeros at c = −2, −3; wt0 709 / 130 zeros | `python3 wcensus.py 5 -2,-3` (~18 min) | `wc5.txt` |

Key of `census.py` output: dictionary keyed by (allowed, Z0, Z∞), Z0 / Z∞ = residue vanishes at both seeds for
the q=0 basis (wt 0) / the q=∞ basis (wt 1). Theorem D∞ predicts Z∞ ⟺ not allowed.

`wcensus.py` replaces the weight w^{−2N} by Π_{m∈M} w_m^c (the weighted problem R^A(λ,μ,ν;a,b;c) of §4/§6);
note the Day 246 run sampled even c only — see `rick-research/proofs/2026-10-19-theorem-W-sharp.md` for the
corrected window (odd c also vanish; scripts in `../prove247/`, `../prove247b/`).

## Other §8 scripts

- Closed form (Prop 1): `chk.py`, `chk2.py` (N = 2, 3 symbolic), `gauge.py` (N = 2, 3, 4).
- Dual / (D): `biorth.py` (N = 2, 3), `biorth4.py` (N = 4, from stored full-q Ã), `fdual_num.py` (F-form, N = 4, 5).
- Reduction steps: `test_red.py`, `test_red2.py`, `test_red3.py` (A2–A4, B2–B4, C2–C3 numerically).

Run everything from this directory (`sys.path.insert(0,'.')`). Most scripts need only the Python 3 stdlib; `cf.py`, `dual*.py`, `fdual.py`, `rho_sym.py` use SymPy.

# prove248 — scripts behind `proofs/2026-10-20-dual-block-identity.md`

Run from this directory. Needs `../prove247/mp.py` (census engine) and `../prove247b/` (v1 certifier, z4/z5 zero sets).

| file | what |
|---|---|
| `psi.py` | sympy Ψ_{κ,θ}(M) from the T⁰ definition |
| `t1.py`, `t3.py` | (1)-identity; structure of Ψ for small κ (no w_1 for consecutive from 1) |
| `t2.py` | B_L leading coefficient vs closed form (s=2, L=2) |
| `t4.py` | exact Lemma T* / Lemma H (3 cases) and D_u for B_2 — differences 0 |
| `t5.py` | Lemma P (whole block), 6 cases, ratio 1 |
| `t6.py` | Lemma P1, r=2,3 |
| `t7.py` | Lemma Λ, 3 test functions |
| `t8.py`, `t9.py` | Lemma P1′, Lemma P2 |
| `t10.py` | Lemma P* (lower block), random rational point |
| `scert2..7.py` | certifier v1 + P (2), + P1 (3), general top block + fallback (4), + P1′ (5 = "v2"), + P2 (6 = "v3"), per-block P* (7 = "v4") |
| `run_bc*.py`, `runA*.py` | B/C n≤6, A N=6,7 runs (`*.txt`) |
| `negctl*.py`, `negBC*.py`, `negC5.py`, `negC6P2.py` | negative controls on known in-window zeros |
| `mut.py`, `mut5.py`, `scert5m*.py`, `scert6m.py`, `scert7m.py` | mutation test of the new branches |
| `base.py` | base signed-count distribution (all ±1) |

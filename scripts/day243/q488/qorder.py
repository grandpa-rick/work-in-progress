"""Q488 converse, order by order in q.  Rick, Day 243 (2026-10-17).

Object: c = c^nu_{lam,mu}(w;q) in Q(q)(w_1..w_N), structure constants of Motegi's normalised
type-C (or B) weight functions (arXiv:2610.15344), engine = day242/q484/engine.py (copied, unchanged).

Definition (q-adic order).  Write c = q^v * n(q,w)/d(q,w) with n,d in Q[w][q], n(0,w)!=0, d(0,w)!=0.
Expand c in Q(w)((q)):  c = sum_{k>=v} c_k(w) q^k.  c_k is computed exactly by power-series division
in q with coefficients in Q(w).  "c is non-Laurent at order q^k" := c_k, in lowest terms, has a
denominator with a factor that is not a monomial in w.  first(c) := least k with c_k non-Laurent
(None if c_k Laurent for all k <= KMAX).
Logic: if c is a Laurent polynomial in w over Q(q), every c_k is Laurent; so first(c)=k is a
certificate that c is NOT Laurent (i.e. Q488's converse holds for that triple), with the pole
visible at order q^k.  (The converse implication is false in general -- 1/(1-q w1 w2) has every
c_k Laurent -- which is why we also record full Laurentness of c over Q(q).)

Usage: python3 qorder.py X N [KMAX] [mode]   mode in {norm (default), unnorm, P1}
  unnorm = Motegi's WC/WB without E^X  (negative control: must break 'jump<=1 => Laurent')
  P1     = planted engine: (u_j u_k - q) -> (u_j u_k - 1) in the top block (negative control)
  f      = pointwise renormalisation phi_I(M) -> f(M) phi_I(M), f(M) = 1 + sum_{m in M} w_m
           (q-independent, non-monomial; negative control for 'c_0 always Laurent')
"""
import sys, time, sympy as sp
from engine import q, subsets, geq, structure_constants, is_laurent

def qcoeffs(c, w, kmax):
    """returns (v, [c_v, ..., c_kmax]) exact, as rational functions of w; c != 0."""
    c = sp.cancel(sp.together(c))
    n, d = sp.fraction(c)
    gens = list(w)
    Pn = sp.Poly(n, q); Pd = sp.Poly(d, q)
    nc = list(reversed(Pn.all_coeffs())); dc = list(reversed(Pd.all_coeffs()))
    vn = next(i for i, a in enumerate(nc) if a != 0); vd = next(i for i, a in enumerate(dc) if a != 0)
    nc = nc[vn:]; dc = dc[vd:]; v = vn - vd
    out = []
    e = []
    for k in range(0, max(kmax - v + 1, 0)):
        nk = nc[k] if k < len(nc) else 0
        s = nk - sum((dc[i] if i < len(dc) else 0) * e[k - i] for i in range(1, k + 1))
        e.append(sp.cancel(s / dc[0]))
    return v, e  # e[j] = coefficient of q^{v+j}

def main():
    X, N = sys.argv[1], int(sys.argv[2])
    KMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    mode = sys.argv[4] if len(sys.argv) > 4 else 'norm'
    t0 = time.time()
    A, C, w = structure_constants(X, N, q, plant=('P1' if mode == 'P1' else None),
                                  normalise=(mode != 'unnorm'))
    if mode == 'f':
        S = subsets(N)
        f = {M: 1 + sum(w[m - 1] for m in M) for M in S}
        A = {k: sp.factor(a * f[k[0]]) for k, a in A.items()}
        order = sorted(S, key=lambda T: (len(T), sum(T)))
        C = {}
        for l in S:
            for m in S:
                c = {}
                for nu in order:
                    rhs = A[nu, l] * A[nu, m] - sum(c[k] * A[nu, k] for k in c)
                    c[nu] = sp.factor(sp.cancel(rhs / A[nu, nu]))
                C[l, m] = c
    supbad = 0
    tab = {}      # (jump, laurent_over_Q(q), first) -> count
    vals = {}     # jump -> list of q-valuations
    rows = []
    for (l, m), c in C.items():
        for nu, val in c.items():
            if not (geq(nu, l) and geq(nu, m)):
                supbad += (val != 0); continue
            jump = len(nu) - max(len(l), len(m))
            L = is_laurent(val, w)
            first = None; v = None; dens = ''
            if val != 0:
                v, e = qcoeffs(val, w, KMAX)
                for j, ck in enumerate(e):
                    if not is_laurent(ck, w):
                        first = v + j
                        dens = str(sp.factor(sp.fraction(sp.cancel(ck))[1])); break
            tab[(jump, L, first)] = tab.get((jump, L, first), 0) + 1
            vals.setdefault(jump, []).append(v)
            if first is not None or not L:
                rows.append((l, m, nu, jump, L, v, first, dens))
    print(f'# qorder X={X} N={N} mode={mode} KMAX={KMAX}  time={time.time()-t0:.1f}s')
    print('# support-bound violations (nonzero c outside nu>=lam,mu):', supbad)
    print('# key: (jump, Laurent over Q(q), first non-Laurent q-order or None) : count')
    for k in sorted(tab, key=str): print(' ', k, ':', tab[k])
    for j in sorted(vals):
        vs = [x for x in vals[j] if x is not None]
        print(f'  jump {j}: {len(vals[j])} triples, zero c: {len(vals[j])-len(vs)}, q-valuation range',
              (min(vs), max(vs)) if vs else None)
    print('# non-Laurent / failing triples: lam mu nu jump Laurent v first denominator(c_first)')
    for r in rows: print(' ', *r)

if __name__ == '__main__':
    main()

"""Day 220 PROVE: blind test of v_{(1-s)} c_{lam mu} = max(1, l(lam)-l(mu)) at n=6,7, for several rational t, s symbolic.
e*_lam = E_{lam1}...E_{laml}(1) in N=n vars (engine = biderivation.py, Day 214 subset formula). Grade: computed."""
import sys, itertools, functools
from sympy.polys.rings import ring
from sympy import QQ, Matrix, Symbol, Rational, factor, cancel, expand, Poly
out = open(sys.argv[1], 'w'); T = Rational(sys.argv[2]); n = int(sys.argv[3])
def log(*a): print(*a, flush=True); print(*a, file=out, flush=True)
S = Symbol('s')
def parts(n, m=None):
    if m is None: m = n
    if n == 0: yield (); return
    for k in range(min(n, m), 0, -1):
        for p in parts(n-k, k): yield (k,)+p
N = n
Rg, *g = ring(['s']+[f'x{i}' for i in range(N)], QQ); s, xs = g[0], g[1:]
V = Rg(1)
for i in range(N):
    for j in range(i+1, N): V *= xs[i]-xs[j]
def Ek(k, F):
    tot = Rg(0); Fd = F.to_dict()
    for A in itertools.combinations(range(N), k):
        B = [j for j in range(N) if j not in A]
        sign = (-1)**sum(1 for i in A for j in B if i > j); term = Rg(sign)
        for i in range(N):
            for j in range(i+1, N):
                if (i in A) == (j in A): term *= xs[i]-xs[j]
        for i in A:
            term *= xs[i]
            for j in B: term *= xs[i]-T*xs[j]
        Fs = Rg.from_dict({(m[0]+sum(m[1+i] for i in A),)+m[1:]: c for m, c in Fd.items()})
        tot += term*Fs
    q, r = tot.div(V); assert r == 0
    return q
P = list(parts(n))
def ek(k): return sum((functools.reduce(lambda a, b: a*b, [xs[i] for i in A], Rg(1)) for A in itertools.combinations(range(N), k)), Rg(0))
def epoly(mu): return functools.reduce(lambda a, b: a*b, [ek(k) for k in mu], Rg(1))
M = Matrix(len(P), len(P), lambda i, j: 0)
for i, mu in enumerate(P):
    d = epoly(mu).to_dict()
    for j, lam in enumerate(P): M[i, j] = d.get((0,)+tuple(lam)+(0,)*(N-len(lam)), 0)
Minv = M.inv()
def dom(a, b):  # a >= b in dominance
    sa = sb = 0
    for i in range(max(len(a), len(b))):
        sa += a[i] if i < len(a) else 0; sb += b[i] if i < len(b) else 0
        if sa < sb: return False
    return True
log(f't={T} n={n}')
tot = ok = 0; supp_bad = 0
for lam in P:
    F = Rg(1)
    for k in reversed(lam): F = Ek(k, F)
    d = F.to_dict(); mvec = []
    for mu in P:
        key = tuple(mu)+(0,)*(N-len(mu))
        mvec.append(sum(QQ.to_sympy(v)*S**m[0] for m, v in d.items() if m[1:] == key))
    for i, mu in enumerate(P):
        c = expand(sum(mvec[j]*Minv[j, i] for j in range(len(P))))
        if mu == lam:
            log(f'  diag {lam}: c(1)={c.subs(S,1)}'); continue
        if c == 0:
            if dom(mu, lam): log(f'  ZERO but mu dominates: {lam} {mu}')
            continue
        if not dom(mu, lam): supp_bad += 1; log(f'  SUPPORT VIOLATION {lam} -> {mu}')
        pc = Poly(c, S); v = 0
        while pc.eval(1) == 0: pc = Poly(cancel(pc.as_expr()/(S-1)), S); v += 1
        pred = max(1, len(lam)-len(mu)); tot += 1; ok += (v == pred)
        log(f'  {lam} -> {mu}: v={v} pred={pred} {"OK" if v==pred else ("BOUND-OK-STRICT" if v>pred else "BOUND VIOLATED")} lead={pc.eval(1)}')
log(f'SUMMARY t={T} n={n}: equality {ok}/{tot}; support violations {supp_bad}')

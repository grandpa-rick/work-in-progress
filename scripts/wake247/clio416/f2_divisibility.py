"""F2 (Clio 416): membership (not just vanishing) for all three root types, via the division algorithm.
For each edge M--M' of type alpha and each label I: D = Ahat_M - Ahat_M' (Ahat = Atilde*Psi), clear the monomial
denominator, and divide by the root polynomial f_alpha in the variable w_b (f monic or with unit leading coeff):
  2e_a: f = w_a^2 - 1 (in w_a);  e_a+e_b: f = w_a w_b - 1 (in w_b, lead coeff w_a, a unit);  e_a-e_b: f = w_b - w_a.
Remainder must be 0.  Plant: raw Atilde (no Psi) must leave a nonzero remainder exactly on e_a+e_b edges with nonzero data.
Engine: work-in-progress/scripts/day242/q484/engine.py (n=N padded; value = unpadded by Motegi stability)."""
import sys, itertools, sympy as sp
import sympy.combinatorics
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/day242/q484')
from engine import W, subsets, q
N = int(sys.argv[1]) if len(sys.argv) > 1 else 2
S = subsets(N)
A = {}; 
for M in S:
    for I in S:
        A[M, I], w = W('C', N, I, M)
psi = lambda v: (1 - q*v**2)/(1 - q)
Psi = lambda M: sp.prod([psi(w[m-1]) for m in M])
def divisible(expr, f, var):
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    # den must be c(q)*monomial in w (Laurentness); then divisibility of num suffices
    for fac, e in sp.factor_list(den)[1]:
        if fac.free_symbols & set(w): assert sp.Poly(fac, *w).is_monomial, ('non-Laurent', fac)
    r = sp.rem(sp.Poly(sp.expand(num), var), sp.Poly(f, var))
    return sp.simplify(r.as_expr()) == 0
edges = []
for M in S:
    Ms = set(M)
    for a in M: edges.append(('2e', M, tuple(sorted(Ms - {a})), a, None))
    for a, b in itertools.combinations(M, 2): edges.append(('e+e', M, tuple(sorted(Ms - {a, b})), a, b))
    for b in M:
        for a in range(1, N+1):
            if a not in Ms: edges.append(('e-e', M, tuple(sorted((Ms - {b}) | {a})), a, b))
tally = {}
for typ, M, Mp, a, b in edges:
    if typ == '2e': f, var = w[a-1]**2 - 1, w[a-1]
    elif typ == 'e+e': f, var = w[a-1]*w[b-1] - 1, w[b-1]
    else: f, var = w[b-1] - w[a-1], w[b-1]
    for I in S:
        hat = A[M, I]*Psi(M) - A[Mp, I]*Psi(Mp)
        raw = A[M, I] - A[Mp, I]
        k = tally.setdefault(typ, [0, 0, 0, 0])
        k[0] += 1; k[1] += divisible(hat, f, var)
        if raw != 0: k[2] += 1; k[3] += (not divisible(raw, f, var))
print('N=%d  type: [cases, Ahat divisible, raw nonzero, raw NOT divisible (plant fires)]' % N)
for t, v in tally.items(): print(t, v)

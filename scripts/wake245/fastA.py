"""Fast exact Ã_{M,I} (full symbolic q, symbolic w) via Motegi summand UNPADDED: n=|M| (stability, as in
day243/q488c/num4.py), zero unless M>=I (triangularity).  Same normalisation E as engine.W with n=|M|."""
import sys, itertools, sympy as sp
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/day242/q484')
from engine import summand, denom_pairs, geq, q
def A(N, M, I, w, qq=q):
    q_ = qq
    s, r = len(M), len(I)
    if not geq(M, I): return sp.Integer(0)
    n = s; u = [w[m-1] for m in M]; tot = 0
    for perm in itertools.permutations(range(n)):
        us = [u[p] for p in perm]
        tot += summand(us, w, I, N, n, 1, 0, q_) / sp.prod([us[a]-us[b] for (a, b) in denom_pairs(n, r)])
    E = (1-q_)**(-(n-r)) * sp.prod([1-q_*a/b for a in u for b in u]) * sp.prod([1-q_*a*b for a in u for b in w])
    E /= sp.prod([1-q_*u[j]*u[k] for j in range(n) for k in range(j+1, n)])
    E /= sp.prod([1-q_*w[j]*w[k] for j in range(N) for k in range(j, N)])
    return sp.factor(sp.cancel(sp.together(tot * sp.prod([x**(-(N+n)) for x in u]) / sp.factorial(n-r) / E)))

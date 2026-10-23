# Equivariant test: is div(p_kappa) on X_kappa the reduced boundary?  Solve
#   D_k(j w) - e^{mu} D_k((j-1) w) = sum_{rho<k} c_rho D_rho(j w)   (T-characters, j=1..3)
# for integer c_rho (K-class of the hyperplane section in the Schubert basis).  Compare with -mobius (reduced boundary).
import sys, itertools, sympy as sp
from collections import defaultdict
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/prove249d')
from dem import Di, word
N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
def D(k, j):
    f = {tuple([j] * N): 1}
    if k:
        for i in reversed(word(set(k), N)): f = Di(f, i, N)
    return f
def leq(a, b): return len(a) <= len(b) and all(a[i] <= b[i] for i in range(len(a)))
def below(k):
    return [tuple(sorted(c, reverse=True)) for L in range(len(k) + 1) for c in itertools.combinations(range(1, N + 1), L)
            if leq(tuple(sorted(c, reverse=True)), k) and tuple(sorted(c, reverse=True)) != k]
def mob(k):
    els = below(k) + [k]; mu = {k: 1}
    for rho in sorted(els, key=sum, reverse=True):
        if rho != k: mu[rho] = -sum(mu[s] for s in mu if s != rho and leq(rho, s))
    return {rho: -mu[rho] for rho in els if rho != k and mu[rho]}
def solve(k, mu):
    rhos = below(k); cs = sp.symbols(f'c0:{len(rhos)}'); eqs = []
    for j in (1, 2, 3):
        lhs = defaultdict(int)
        for w, a in D(k, j).items(): lhs[w] += a
        for w, a in D(k, j - 1).items(): lhs[tuple(x + y for x, y in zip(w, mu))] -= a
        for ci, rho in zip(cs, rhos):
            for w, a in D(rho, j).items(): lhs[w] -= ci * a
        eqs += [v for v in lhs.values() if v != 0]
    sol = sp.solve(eqs, cs, dict=True)
    if not sol: return None
    return {rho: sol[0].get(c, c) for rho, c in zip(rhos, cs) if sol[0].get(c, c) != 0}
if __name__ == '__main__':
  for L in range(1, N + 1):
    for c in itertools.combinations(range(1, N + 1), L):
        k = tuple(sorted(c, reverse=True))
        wom = tuple(-1 if (N + 1 - p) in k else 1 for p in range(1, N + 1))
        res = None
        for mu in (tuple([1] * N), tuple([-1] * N), wom, tuple(-x for x in wom)):
            res = solve(k, mu)
            if res is not None: break
        print(k, 'mu', mu, 'hyperplane K-class', {''.join(map(str, p)) or '0': v for p, v in (res or {}).items()},
              '| reduced boundary (-mobius)', {''.join(map(str, p)) or '0': v for p, v in mob(k).items()}, flush=True)

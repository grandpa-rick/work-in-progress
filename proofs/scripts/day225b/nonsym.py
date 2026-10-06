# coefficients of K = prod_{i<j}(1-z_i/z_j)/(1-t z_i/z_j) as exact sympy series; [z^mu]K for mu in root lattice
import sympy as sp, itertools, functools
t = sp.symbols('t')
@functools.lru_cache(None)
def Kcoef(mu):
    """[z^mu] K, mu tuple sum 0, prefix sums >= 0 (else 0). recursion on pairs (i<j) w/ exponent k>=0 of z_i/z_j."""
    a = len(mu)
    pairs = [(i, j) for i in range(a) for j in range(i+1, a)]
    # enumerate flows: k_{ij} >=0 with net out at i = mu_i  (z_i^{k} z_j^{-k})
    def rec(p, rem):
        if p == len(pairs):
            return sp.Integer(1) if all(r == 0 for r in rem) else sp.Integer(0)
        i, j = pairs[p]
        tot = 0
        # bound k by large positive prefix
        for k in range(0, 40):
            if k > 0 and rem[i] - k < -sum(1 for _ in [0])*0 - 10**6: break
            r = list(rem); r[i] -= k; r[j] += k
            # prune: after all pairs with first index i processed, rem[i] must be 0
            last_i = all(pp[0] != i for pp in pairs[p+1:])
            if last_i and r[i] != 0:
                if r[i] < 0: break
                continue
            if r[i] < 0 and all(pp[1] != i for pp in pairs[p+1:]): break
            w = 1 if k == 0 else -(1-t)*t**(k-1)
            tot += w*rec(p+1, tuple(r))
        return tot
    if any(sum(mu[:m]) < 0 for m in range(1, a+1)): return sp.Integer(0)
    return sp.expand(rec(0, tuple(mu)))
def ctmono(gam):  # CT[z^gam K] = [z^{-gam}]K
    return Kcoef(tuple(-g for g in gam))

# exploratory: s-adic valuation and lowest coefficient of c_{lam mu}(s,t) at fixed t, by exact interpolation in s
from ds_pointeval import *
from fractions import Fraction as Fr
import random, sys
def interp(xs, ys):  # Newton -> coefficient list
    n = len(xs); coef = [Fr(0)]*n
    # solve Vandermonde
    M = [[x**j for j in range(n)] for x in xs]; return solve(M, ys)
if __name__ == '__main__':
    T = Fr(sys.argv[2]) if len(sys.argv) > 2 else Fr(0)
    rng = random.Random(5)
    for n in range(2, int(sys.argv[1])+1):
        for lam in parts(n):
            D = n*n  # safe s-degree bound
            svals = [Fr(i+2, 3) for i in range(D+1)]
            data = [coeffs(lam, sv, T, rng) for sv in svals]
            res = {}
            for mu in parts(n):
                if not dominates(mu, lam): continue
                c = interp(svals, [dd.get(mu, Fr(0)) for dd in data])
                v = next((i for i, a in enumerate(c) if a != 0), None)
                deg = max((i for i, a in enumerate(c) if a != 0), default=None)
                res[mu] = (v, nstat(mu), str(c[v]) if v is not None else None, deg)
            print(lam, res, flush=True)

# parallel n=6 part of check_stretch (3)
from fractions import Fraction as Fr
from multiprocessing import Pool
from ds_pointeval import parts, dominates, nstat, coeffs, solve
import random, sys
def job(args):
    lam, T = args; rng = random.Random(hash(lam) % 1000); n = sum(lam)
    D = nstat(lam) + 1; sv = [Fr(i+2, 3) for i in range(D+1)]
    data = [coeffs(lam, x, T, rng) for x in sv]; chk = coeffs(lam, Fr(17, 5), T, rng); res = True
    for mu in parts(n):
        if not dominates(mu, lam): continue
        c = solve([[x**j for j in range(D+1)] for x in sv], [d.get(mu, Fr(0)) for d in data])
        assert sum(cj*Fr(17, 5)**j for j, cj in enumerate(c)) == chk.get(mu, 0)
        v = next(i for i, a in enumerate(c) if a != 0)
        res &= v == nstat(mu) and (T != 0 or c[v] == 1)
    return lam, str(T), res
if __name__ == '__main__':
    n = int(sys.argv[1]); T = Fr(sys.argv[2]); ok = True
    with Pool(4) as P:
        for lam, T_, r in P.imap_unordered(job, [(l, T) for l in parts(n)]):
            print(lam, 't=', T_, r, flush=True); ok &= r
    print('ALL OK' if ok else 'FAIL')

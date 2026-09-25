"""Day 206: held-out check of W_r (Day 191 conjecture) straight from the AHA pipeline:
e_2(Y).e_r  ==  t * W_r  at random exact rational points, (m,r) = (7,5), (8,6). No kernel used."""
import sys, random, time
from fractions import Fraction as Fr
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day206')
from parabolic_kernel_test import AHA, lhs_poly, evalp, e_vals
def qi(n, t): return sum(t ** i for i in range(n)) if n > 0 else Fr(0)
def Wnum(r, X, s, t):
    e = lambda k: e_vals(X, k); q = 1 / s
    return (s ** 2 * e(2) * e(r) + (1 - s) * s * qi(r, t) * e(1) * e(r + 1)
            + (1 - s) * qi(r + 2, t) / qi(2, t) * (qi(r + 1, t) - t * qi(r - 1, t) * s) * e(r + 2))
random.seed(206)
for m, r in [(5, 2), (6, 4), (7, 5), (8, 6)]:
    t0 = time.time(); A_ = AHA(m); P = lhs_poly(A_, r); res = []
    for _ in range(4):
        X = [Fr(v, random.randint(1, 5)) for v in random.sample(range(2, 99), m)]
        s = Fr(random.randint(2, 9), random.randint(1, 4)); t = Fr(random.randint(2, 9), random.randint(1, 4))
        res.append(evalp(A_, P, X, s, t) == t * Wnum(r, X, s, t))
    print(f'm={m} r={r}: e_2(Y).e_r == t*W_r at 4 pts: {res}  ({time.time()-t0:.0f}s)', flush=True)

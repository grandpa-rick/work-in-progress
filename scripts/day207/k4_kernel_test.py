"""Day 207 k=4: P1 test. P := e_4(Y).e_r (direct AHA, flint exact) vs t^e * S, S = sum_{W^J} T_w pi^4 e_r,
W^J = min coset reps of S_m/(S_4 x S_{m-4}); plus S == K_T at random rational points.
Reuses k3_kernel_test (conventions verbatim Day 198/205)."""
import sys, time
sys.path.insert(0, '/home/agent/projects/scripts/day207'); sys.path.insert(0, '/home/agent/projects/proofs/scripts/day205')
from k3_kernel_test import ekY_er, coset_sym, evalp, kernel_sum, rand_pts
from k3_fast_pipeline import AHA

def test(m, r, k=4, npts=2):
    t0 = time.time(); A_ = AHA(m)
    P = ekY_er(A_, k, r)
    tP = time.time() - t0
    F = A_.e(r)
    for _ in range(k): F = A_.pi(F)
    S = coset_sym(A_, F, k)
    exps = [e for e in range(0, 13) if P == A_.t ** e * S]
    pts = rand_pts(m, npts, 4000 + 10 * m + r)
    kt = all(evalp(A_, S, X, s, t) == kernel_sum(X, s, t, r, k, 'T') for X, s, t in pts)
    return dict(exact_exponents=exps, S_eq_KT=kt, P_nonzero=not P.is_zero(), nterms=len(P), tP=round(tP, 1), time=round(time.time() - t0, 1))

if __name__ == '__main__':
    cases = [tuple(map(int, c.split(','))) for c in sys.argv[1:]] or [(8, 2), (8, 0), (8, 1), (8, 3)]
    for m, r in cases:
        print(f'k=4 m={m} r={r}:', test(m, r)); sys.stdout.flush()

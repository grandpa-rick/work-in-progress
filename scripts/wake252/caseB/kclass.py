# TRUE K-class of the hyperplane section D=div(p_k) on X_k (LG(N,2N)) in the Schubert basis {[O_{X_rho}]}:
# pull back to G/B (pi^* injective on K, line bundles L_lambda span K(G/B)_Q), so for all weights lambda
#   chi(L_lam on pi^-1 X_k) - chi(L_{lam-omega} on pi^-1 X_k) = sum_rho c_rho chi(L_lam on pi^-1 X_rho)
# with chi(pi^-1 X_rho, L_lam) = |D_{w_rho w_0P}(e^lam)| (Demazure character formula, any lambda).
import sys, itertools, random
from functools import lru_cache
import numpy as np
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/prove249d')
from dem import Di, word
N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
W0P = [i for t in range(N - 1) for i in range(1, N - t)]   # reduced word of w_0(S_N) in s_1..s_{N-1}
@lru_cache(None)
def wordmax(k): return tuple(reversed(word(set(k), N))) if k else ()
@lru_cache(None)
def chiB(k, lam):
    f = {lam: 1}
    for i in W0P: f = Di(f, i, N)
    for i in wordmax(k): f = Di(f, i, N)
    return sum(f.values())
def leq(a, b): return len(a) <= len(b) and all(a[i] <= b[i] for i in range(len(a)))
def below(k):
    return [r for r in (tuple(sorted(c, reverse=True)) for L in range(len(k) + 1) for c in itertools.combinations(range(1, N + 1), L))
            if leq(r, k) and r != k]
random.seed(1)
LAMS = [tuple(random.randint(-3, 3) for _ in range(N)) for _ in range(6 * 2 ** N)]
@lru_cache(None)
def kclass(k):
    rh = below(k); om = tuple([1] * N)
    A = [[chiB(r, l) for r in rh] for l in LAMS]
    b = [chiB(k, l) - chiB(k, tuple(x - 1 for x in l)) for l in LAMS]
    An = np.array(A, dtype=float)
    assert np.linalg.matrix_rank(An) == len(rh), ('not separating', k)
    sol = [int(round(v)) for v in np.linalg.lstsq(An, np.array(b, dtype=float), rcond=None)[0]]
    assert all(sum(a * c for a, c in zip(row, sol)) == bb for row, bb in zip(A, b)), ('inexact', k)   # exact integer check
    return {r: int(v) for r, v in zip(rh, sol) if v != 0}
if __name__ == '__main__':
    for L in range(1, N + 1):
        for c in itertools.combinations(range(1, N + 1), L):
            k = tuple(sorted(c, reverse=True))
            print(k, {''.join(map(str, r)) or '0': v for r, v in kclass(k).items()}, flush=True)

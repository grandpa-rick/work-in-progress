"""Negative controls: the AHA test must reject (i) P-ell with the K_{13} pair factor dropped, (ii) the connection-note
prefactor s^{sum(n_c-i_c)} (version A) instead of s^{(ell-1) sum(n_c-i_c)}."""
import sys, random, itertools, math
from fractions import Fraction as Fr
import aha_check_ell as M
random.seed(7)
A = M.AHA(4); ell, k = 3, 2
ops = {tup: M.ekY(A, math.prod((A.e(a) for a in tup), start=1+0*A.X[0]), k) for tup in itertools.combinations_with_replacement(range(5), ell)}
X = [Fr(random.randint(1, 30), random.randint(31, 60)) for _ in range(4)]; sv, tv = Fr(3, 11), Fr(5, 13); Z = [Fr(7, 3), Fr(11, 5), Fr(13, 9)]
lhs = sum(M.evalp(A, P, X, sv, tv)*sum(math.prod(Z[a]**p[a] for a in range(ell)) for p in set(itertools.permutations(tup))) for tup, P in ops.items())
print('true P-ell (version C):', lhs == M.pred(k, X, sv, tv, Z))
K0 = M.K
M.K = lambda i, j, z, w, s, t: Fr(1) if (z, w) == (Z[0], Z[2]) else K0(i, j, z, w, s, t)
print('drop K_13 factor      :', lhs == M.pred(k, X, sv, tv, Z))
M.K = K0
def predA(k, X, s, t, Z):  # version A prefactor
    ell = len(Z); tot = Fr(0)
    for b in range(k+1):
        for n in itertools.product(range(k-b+1), repeat=ell):
            if sum(n) != k-b: continue
            for I in itertools.product(*[range(nc+1) for nc in n]):
                tA = -sum((n[a]-I[a])*I[d] for a in range(ell) for d in range(ell) if d != a)
                wt = (s**ell*t**(-sum(I)))**b * s**(sum(n[a]-I[a] for a in range(ell))) * t**tA
                Kp = math.prod((M.K(I[a], I[d], Z[a], Z[d], s, t) for a in range(ell) for d in range(a+1, ell)), start=Fr(1))
                Es = math.prod((math.prod((1 + t**I[a]*Z[a]*x for x in X), start=Fr(1)) for a in range(ell)), start=Fr(1))
                tot += wt*Kp*math.prod((M.Nn(n[a], I[a], s, t)*Z[a]**(-n[a]) for a in range(ell)), start=Fr(1))*M.ev(X, b)*Es
    return tot
print('version A prefactor   :', lhs == predA(k, X, sv, tv, Z))

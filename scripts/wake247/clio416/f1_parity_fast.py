"""F1 fast exact-rational version: evaluate Motegi's symmetrised F (Thm 3.3 / Motherfunction, as transcribed in
work-in-progress/scripts/day242/q484/engine.py summand) WITHOUT Vandermonde division, at random rational (u,w,q),
and compare P(-u,-w) with (-1)^{e_P} P(u,w), where P := F * prod u^{N+n} (the bare sum over permutations) and
(first run of this script multiplied P by prod u^{N+n} a second time -- harness bug, fixed; recorded here)
  e_P = (N-1)r + N(n-r) + (n-r)   [heads, tails, beta*u-gamma with gamma=0];
prefactor prod u^{-(N+n)} contributes n(N+n); total 2nN+n^2+n-2r == 0 mod 2 for ALL n,N,r.
Controls: (i) WRONG exponent e_P+1 must fail whenever P != 0; (ii) plant gamma=-1 (type-B tail) must break parity
when n>r; (iii) full normalised value A~_{M,I}(-w) = A~_{M,I}(w) with n=|M| (unpadded) and E^C as printed.
"""
import sys, itertools, random
from fractions import Fraction as Fr
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/day242/q484')
from engine import summand, denom_pairs

def perm_sign(p):
    s, p = 1, list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s

def Fsym(u, w, I, N, n, gamma, qq):
    r = len(I); tot = Fr(0)
    for perm in itertools.permutations(range(n)):
        us = [u[p] for p in perm]
        t = summand(us, w, I, N, n, 1, gamma, qq)
        d = Fr(1)
        for (a, b) in denom_pairs(n, r): d *= (us[a] - us[b])
        tot += Fr(t) / d
    return tot

def rnd(): return Fr(random.randint(2, 40), random.randint(41, 90)) * random.choice([1, -1])
random.seed(20261019)
st = dict(eP_ok=0, eP_tot=0, eP_odd=0, wrong_fail=0, wrong_tot=0, plant_broken=0, plant_tot=0, val_ok=0, val_tot=0, val_nonzero=0)
for N in (1, 2, 3, 4):
    for n in range(1, N + 1):
        for r in range(0, n + 1):
            for I in itertools.combinations(range(1, N + 1), r):
                for trial in range(2):
                    u = [rnd() for _ in range(n)]; w = [rnd() for _ in range(N)]; qq = rnd()
                    pu = 1
                    for x in u: pu *= x
                    P1 = Fsym(u, w, I, N, n, 0, qq)   # P = sum over perms, no prefactor
                    P2 = Fsym([-x for x in u], [-x for x in w], I, N, n, 0, qq)
                    F1, F2 = P1 / pu**(N + n), P2 / ((-1)**n * pu)**(N + n)   # prod(-u) = (-1)^n prod u
                    st['F_tot'] = st.get('F_tot', 0) + 1; st['F_even'] = st.get('F_even', 0) + (F1 == F2)
                    eP = ((N - 1) * r + N * (n - r) + (n - r)) % 2
                    st['eP_tot'] += 1; st['eP_odd'] += eP
                    st['eP_ok'] += (P2 == (-1)**eP * P1)
                    if P1 != 0:
                        st['wrong_tot'] += 1; st['wrong_fail'] += (P2 != (-1)**(eP + 1) * P1)
                    if n > r:
                        B1 = Fsym(u, w, I, N, n, -1, qq); B2 = Fsym([-x for x in u], [-x for x in w], I, N, n, -1, qq)
                        st['plant_tot'] += 1; st['plant_broken'] += (B1 != B2 and B1 != -B2)
# (iii) full normalised value at u = w_M, n = |M|, type C
def Atilde(M, I, w, N, qq):
    n = len(M); r = len(I)
    if n == 0: return Fr(1) if r == 0 else Fr(0)
    u = [w[m - 1] for m in M]
    import math
    F = Fsym(u, w, I, N, n, 0, qq) / math.factorial(n - r)
    for x in u: F /= x**(N + n)
    E = (1 - qq)**(-(n - r))
    for a in u:
        for b in u: E *= (1 - qq * a / b)
        for b in w: E *= (1 - qq * a * b)
    for j in range(n):
        for k in range(j + 1, n): E /= (1 - qq * u[j] * u[k])
    for j in range(N):
        for k in range(j, N): E /= (1 - qq * w[j] * w[k])
    return F / E
for N in (2, 3, 4):
    for s in range(1, N + 1):
        for M in itertools.combinations(range(1, N + 1), s):
            for r in range(0, s + 1):
                for I in itertools.combinations(range(1, N + 1), r):
                    w = [rnd() for _ in range(N)]; qq = rnd()
                    a1 = Atilde(M, I, w, N, qq); a2 = Atilde(M, I, [-x for x in w], N, qq)
                    st['val_tot'] += 1; st['val_ok'] += (a1 == a2); st['val_nonzero'] += (a1 != 0)
print(st)

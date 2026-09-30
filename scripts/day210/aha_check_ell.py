"""Day 210: blind test of P-ell (version C) against DIRECT AHA (day205 engine, day208 ekY via proved (A_k); also
ekY_direct for small m).  LHS = sum_{a_1..a_ell <= m} prod z_c^{a_c} t^{-C(k,2)} e_k(Y).(e_{a_1}...e_{a_ell}) at exact random points.
Full GF over all a_c <= m, so equality at a point tests the whole identity in Lambda_m at that point."""
import sys, random, math, itertools, time
from fractions import Fraction as Fr
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day205'); sys.path.insert(0, '/home/agent/projects/scripts/day208')
from k3_fast_pipeline import AHA
from ek_two_col import ekY, ekY_direct
def poch(x, t, N): return math.prod((1 - x*t**i for i in range(N)), start=Fr(1))
def c(n, j, s, t):
    if j < 0 or j > n: return Fr(0)
    al = lambda J: Fr(0) if J < 0 else math.prod((s - t**i for i in range(1, J+1)), start=Fr(1))/poch(t, t, J)
    return poch(s, t, n-j)/poch(t, t, n-j)*(al(j) - s*t**(n-j)*al(j-1))
def Nn(n, j, s, t): return t**(-n*j)*c(n, j, s, t)
def K(i, j, z, w, s, t):
    return math.prod(((t**p*z - s*w)/(t**p*z - t**j*w) for p in range(i)), start=Fr(1)) * \
           math.prod(((s*z - t**r*w)/(t**i*z - t**r*w) for r in range(j)), start=Fr(1))
def ev(X, r):
    if r < 0 or r > len(X): return Fr(0)
    return sum((math.prod(cc, start=Fr(1)) for cc in itertools.combinations(X, r)), Fr(0))
def pred(k, X, s, t, Z, pairwise=True):
    ell = len(Z); tot = Fr(0)
    for b in range(k+1):
        for n in itertools.product(range(k-b+1), repeat=ell):
            if sum(n) != k-b: continue
            for I in itertools.product(*[range(nc+1) for nc in n]):
                tA = -sum((n[a]-I[a])*I[d] for a in range(ell) for d in range(ell) if d != a)
                wt = (s**ell*t**(-sum(I)))**b * s**((ell-1)*sum(n[a]-I[a] for a in range(ell))) * t**tA
                Kp = math.prod((K(I[a], I[d], Z[a], Z[d], s, t) for a in range(ell) for d in range(a+1, ell)), start=Fr(1))
                Es = math.prod((math.prod((1 + t**I[a]*Z[a]*x for x in X), start=Fr(1)) for a in range(ell)), start=Fr(1))
                tot += wt*Kp*math.prod((Nn(n[a], I[a], s, t)*Z[a]**(-n[a]) for a in range(ell)), start=Fr(1))*ev(X, b)*Es
    return tot
def evalp(A, P, X, sv, tv):
    tot = Fr(0)
    for mon, cc in P.to_dict().items():
        v = Fr(int(cc))
        for i in range(A.m): v *= X[i]**int(mon[i])
        tot += v*sv**int(mon[A.m])*tv**int(mon[A.m+1])
    return tot
if __name__ != "__main__": plan = []
random.seed(210)
ok = True
plan = plan if __name__ != "__main__" else [tuple(map(int, a.split(','))) for a in sys.argv[1:]]   # (ell, m, k)
for ell, m, k in plan:
    t0 = time.time(); A = AHA(m)
    ops = {}
    for tup in itertools.combinations_with_replacement(range(m+1), ell):
        F = math.prod((A.e(a) for a in tup), start=1 + 0*A.X[0])
        ops[tup] = ekY(A, F, k)
        if m <= 3:  # cross-check engine vs direct Y product
            assert (ekY_direct(A, F, k) - A.t**(k*(k-1)//2)*ops[tup]).is_zero()
    for trial in range(2):
        X = [Fr(random.randint(1, 30), random.randint(31, 60)) for _ in range(m)]
        sv, tv = Fr(random.randint(2, 9), 11), Fr(random.randint(2, 9), 13)
        Z = [Fr(random.randint(1, 40), random.randint(3, 11)) for _ in range(ell)]
        lhs = Fr(0)
        for tup, P in ops.items():
            v = evalp(A, P, X, sv, tv)
            for perm in set(itertools.permutations(tup)):
                lhs += v*math.prod((Z[a]**perm[a] for a in range(ell)), start=Fr(1))
        good = lhs == pred(k, X, sv, tv, Z); ok &= good
    print(f'ell={ell} m={m} k={k}: {"OK" if good else "FAIL"}  ({time.time()-t0:.1f}s)', flush=True)
print('ALL OK' if ok else 'SOME FAIL')

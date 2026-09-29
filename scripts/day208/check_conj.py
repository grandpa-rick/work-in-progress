"""Day 208 (Sep 29): test conjectured two-column GF against direct AHA (exact rationals, random points).
Conj (TC_k):  sum_{a,b} z^a w^b t^{-C(k,2)} e_k(Y).(e_a e_b)
  = sum_{b0+n1+n2=k} sum_{i<=n1, j<=n2} s^{2 b0 + (n1-i)+(n2-j)} t^{-b0(i+j) - (n1-i) j - (n2-j) i}
        K_ij(z,w) N^{(n1)}_i N^{(n2)}_j z^{-n1} w^{-n2} e_{b0} E(t^i z) E(t^j w),
  N^{(n)}_j = t^{-nj} c(n,j) (Day 207b), K_ij = prod_{p<i}(t^p z - s w)/(t^p z - t^j w) prod_{r<j}(s z - t^r w)/(t^i z - t^r w)."""
import sys, random, math, itertools, time
from fractions import Fraction as Fr
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day205')
from k3_fast_pipeline import AHA
from ek_two_col import ekY
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
def conj(k, X, s, t, z, w):
    tot = Fr(0)
    for b0 in range(k+1):
        for n1 in range(k-b0+1):
            n2 = k - b0 - n1
            for i in range(n1+1):
                for j in range(n2+1):
                    wt = s**(2*b0 + n1-i + n2-j) * t**(-b0*(i+j) - (n1-i)*j - (n2-j)*i)
                    Ez = math.prod((1 + t**i*z*x for x in X), start=Fr(1)); Ew = math.prod((1 + t**j*w*x for x in X), start=Fr(1))
                    tot += wt*K(i, j, z, w, s, t)*Nn(n1, i, s, t)*Nn(n2, j, s, t)*z**(-n1)*w**(-n2)*ev(X, b0)*Ez*Ew
    return tot
def evalp(A, P, X, sv, tv):
    tot = Fr(0)
    for mon, cc in P.to_dict().items():
        v = Fr(int(cc))
        for i in range(A.m): v *= X[i]**int(mon[i])
        tot += v*sv**int(mon[A.m])*tv**int(mon[A.m+1])
    return tot
random.seed(208)
ok = True
plan = [(m, k) for m in range(1, 7) for k in (1, 2, 3, 4)] if len(sys.argv) < 2 else [(int(a.split(",")[0]), int(a.split(",")[1])) for a in sys.argv[1:]]
for m, k in plan:
    pass
    t0 = time.time(); A = AHA(m)
    ops = {(a, b): ekY(A, A.e(a)*A.e(b), k) for a in range(m+1) for b in range(a, m+1)} if k <= m else {}
    for trial in range(3):
        X = [Fr(random.randint(1, 30), random.randint(31, 60)) for _ in range(m)]
        sv, tv = Fr(random.randint(2, 9), 11), Fr(random.randint(2, 9), 13)
        zv, wv = Fr(random.randint(1, 20), 7), Fr(random.randint(21, 40), 9)
        lhs = Fr(0)
        for (a, b), P in ops.items():
            v = evalp(A, P, X, sv, tv)
            lhs += v*(zv**a*wv**b + (zv**b*wv**a if a != b else 0))
        good = lhs == conj(k, X, sv, tv, zv, wv); ok &= good
    print(f'm={m} k={k}: {"OK" if good else "FAIL"}  ({time.time()-t0:.1f}s)', flush=True)
print('ALL OK' if ok else 'SOME FAIL')

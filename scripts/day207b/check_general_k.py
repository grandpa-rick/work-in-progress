"""Day 207b: checks for the general-k proof (A_k via Key Lemma, K_k via iterated Lemma 1, E_k via (*)).
KL : t^{m-c} Phi^{-1} T_{w(d)} pi^n F = t^n T_{w(d-1)} pi^n F   (Phi = T_c..T_{m-1}), exact in Z[X,s,t].
A_k: Y_{b1}..Y_{bk} F = t^{C(k,2)} T_{w(b)} pi^k F  per tuple, exact.
OL : sigma_[1,m]..sigma_[k,m] G = sum_{ordered distinct a} G(X_a1..X_ak;rest) prod_l prod_{j notin a_1..a_l} a_{a_l j}
     for G tail-symmetric only (random rational point).
STAR: the q-series identity (*) for c(n,j), symbolic, n<=9.
GF : sum_r z^r t^{-C(k,2)} e_k(Y) e_r == sum_{j,b} s^b t^{-kj} c(k-b,j) e_b z^{b-k} E(t^j z), k<=4, exact.
"""
import sys, itertools, random
from fractions import Fraction as Fr
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day205')
from k3_fast_pipeline import AHA
import sympy as sp

ok_all = True
def rep(tag, ok):
    global ok_all
    ok_all &= ok
    print(tag, 'OK' if ok else 'FAIL', flush=True)

def Tw(A_, H, d):  # T_{w(d)} H, w(d) = prod_l T_{d_l-1}..T_l ; rightmost factor acts first
    for l in range(len(d), 0, -1):
        for j in range(l, d[l-1]):
            H = A_.T(H, j)
    return H

def sym_test_F(A_):
    e = A_.e
    return e(2)*e(1) + A_.s*e(3) + e(1)

# ---- KL and A_k
for m in (4, 5, 6):
    A_ = AHA(m); F = sym_test_F(A_)
    okKL = True; okA = True
    for n in range(0, min(m, 4)):
        H = F
        for _ in range(n): H = A_.pi(H)
        for c in range(1, m+1):
            for d in itertools.combinations(range(c+1, m+1), n):
                V = Tw(A_, H, d)
                for j in range(c, m): V = A_.tTinv(V, j)
                W = Tw(A_, H, tuple(x-1 for x in d)) * A_.t**n
                okKL &= (V - W).is_zero()
    rep(f'KL m={m}', okKL)
    for k in range(1, min(m, 4)+1):
        Pk = F
        for _ in range(k): Pk = A_.pi(Pk)
        for b in itertools.combinations(range(1, m+1), k):
            V = F
            for i in reversed(b): V = A_.Y(V, i)
            okA &= (V - A_.t**(k*(k-1)//2) * Tw(A_, Pk, b)).is_zero()
    rep(f'A_k per-tuple m={m} (k<=4)', okA)

# ---- OL: ordered-weight formula for tail-symmetric G, via sympy at a rational point
def OL_check(m, k):
    Xs = sp.symbols(f'x1:{m+1}'); t = sp.Symbol('t')
    # G: head X1..Xk NOT symmetric, tail symmetric
    tail = Xs[k:]
    G = Xs[0]**2*Xs[1 % m] * sum(tail) + Xs[0]*(sum(x**2 for x in tail)) + (Xs[k-1]**3 if k>1 else 0)*sp.prod(tail)
    def T(F, j):  # 1-based
        sF = F.subs({Xs[j-1]: Xs[j], Xs[j]: Xs[j-1]}, simultaneous=True)
        return sp.cancel(t*sF - (t-1)*Xs[j]*(F - sF)/(Xs[j-1]-Xs[j]))
    # sigma_[l,m] = sum_{a=l}^m T_{a-1}..T_l ; product sigma_[1,m]..sigma_[k,m] applied rightmost first
    H = G
    for l in range(k, 0, -1):
        tot = 0; cur = H; tot += cur
        for a in range(l+1, m+1):
            # T_{a-1}..T_l H : apply T_l first
            V = H
            for j in range(l, a): V = T(V, j)
            tot += V
        H = sp.expand(sp.cancel(tot))
    a_ = lambda i, j: (Xs[i]-t*Xs[j])/(Xs[i]-Xs[j])
    rhs = 0
    for tup in itertools.permutations(range(m), k):
        sub = {}
        rest = [i for i in range(m) if i not in tup]
        # G(X_a1..X_ak; rest): head position l -> X_{a_l}, tail positions -> rest (order irrelevant)
        newvars = [Xs[a] for a in tup] + [Xs[i] for i in rest]
        term = G.subs(dict(zip(Xs, newvars)), simultaneous=True)
        w = 1
        for l in range(k):
            for j in range(m):
                if j not in tup[:l+1]: w *= a_(tup[l], j)
        rhs += term*w
    pt = {x: sp.Rational(random.randint(2, 40), random.randint(41, 90)) for x in Xs}; pt[t] = sp.Rational(3, 7)
    return sp.simplify(H.subs(pt) - rhs.subs(pt)) == 0
random.seed(7)
for (m, k) in [(3, 2), (4, 2), (4, 3), (5, 3), (5, 2)]:
    rep(f'OL m={m} k={k}', OL_check(m, k))

# ---- STAR
s, t, x = sp.symbols('s t x')
def poch(a, n):
    return sp.prod([1 - a*t**i for i in range(n)]) if n >= 0 else None
def alpha(j):
    if j < 0: return 0
    return sp.prod([s - t**i for i in range(1, j+1)]) / poch(t, j)
def c(n, j):
    if j < 0 or n - j < 0: return 0
    return poch(s, n-j)/poch(t, n-j) * (alpha(j) - s*t**(n-j)*alpha(j-1))
okS = True
for n in range(1, 10):
    for j in range(0, n+2):
        lhs = (1 - t**n)*c(n, j)
        rhs = (1 - s*t**-j)*sum((s*t**-j)**p * c(n-1-p, j) for p in range(0, n))
        rhs -= t**n*(1 - s*t**(1-j))*sum((s*t**(1-j))**p * c(n-1-p, j-1) for p in range(0, n+1))
        okS &= sp.simplify(lhs - rhs) == 0
rep('STAR n<=9', okS)
# vanishing in 0 variables: sum_j t^{-kj} c(k,j) = 0
rep('T_k^(0)=0 k<=8', all(sp.simplify(sum(t**(-k*j)*c(k, j) for j in range(k+1))) == 0 for k in range(1, 9)))
# M-formula of PROVE.md agrees with s^b t^{-kj} c(k-b,j) (t;t)_k
def M(k, j, b):
    R = s*(t**(k-b) - t**(k-b-j) + 1) - t**j
    return s**b*sp.prod([s-t**i for i in range(j)])*poch(s, k-j-b)*R/(s-1)*poch(t, k)/(poch(t, j)*poch(t, k-j-b))*t**(-k*j)
rep('M == (t;t)_k s^b t^{-kj} c', all(sp.simplify(M(k, j, b) - poch(t, k)*s**b*t**(-k*j)*c(k-b, j)) == 0
                                     for k in range(0, 7) for j in range(k+1) for b in range(k-j+1)))

# ---- GF exact vs AHA, k<=4, all r, m in {k, k+1, k+2} plus m<k vanishing
def ekY_er(A_, k, r):
    tot = 0*A_.s
    for b in itertools.combinations(range(1, A_.m+1), k):
        V = A_.e(r)
        for i in reversed(b): V = A_.Y(V, i)
        tot = tot + V
    return tot
def evalp(A_, P, X, sv, tv):
    tot = Fr(0)
    for mon, cc in P.to_dict().items():
        v = Fr(int(cc))
        for i in range(A_.m): v *= X[i]**int(mon[i])
        tot += v*sv**int(mon[A_.m])*tv**int(mon[A_.m+1])
    return tot
def ev(vals, r):
    if r < 0 or r > len(vals): return Fr(0)
    tot = Fr(0)
    for cmb in itertools.combinations(vals, r):
        p = Fr(1)
        for v in cmb: p *= v
        tot += p
    return tot
def cF(n, j, sv, tv):  # exact Fraction version of c(n,j)
    if j < 0 or n - j < 0: return Fr(0)
    po = lambda a, N: __import__('math').prod([1 - a*tv**i for i in range(N)]) if N > 0 else Fr(1)
    al = lambda J: Fr(0) if J < 0 else (__import__('math').prod([sv - tv**i for i in range(1, J+1)]) if J > 0 else Fr(1)) / po(tv, J)
    return po(sv, n-j)/po(tv, n-j) * (al(j) - sv*tv**(n-j)*al(j-1))
class _F:
    def __init__(self, n): self.n = n
    def __call__(self, sv, tv, xv): return sum(cF(self.n, j, sv, tv)*xv**j for j in range(self.n+1))
Fn = {n: _F(n) for n in range(0, 6)}
okG = True; cnt = 0
random.seed(11)
for m in range(1, 7):
    A_ = AHA(m)
    for k in range(1, 5):
        for r in range(0, m+1):
            L = ekY_er(A_, k, r) if k <= m else None
            for trial in range(2):
                X = [Fr(random.randint(1, 30), random.randint(31, 60)) for _ in range(m)]
                sv, tv = Fr(random.randint(2, 9), 11), Fr(random.randint(2, 9), 13)
                lhs = evalp(A_, L, X, sv, tv) if L is not None else Fr(0)
                rhs = sum(sv**b * Fr(Fn[k-b](sv, tv, tv**(r-b))) * ev(X, b) * ev(X, r+k-b) for b in range(k+1))
                okG &= (lhs == tv**(k*(k-1)//2) * rhs); cnt += 1
rep(f'Pieri formula vs AHA, m<=6, k<=4 (incl. m<k), all r<=m: {cnt} cases', okG)
print('ALL OK' if ok_all else 'SOME FAIL')

"""Day 223: test (**) lin_e T_k f = (-1)^d [n]/[k] f(1,t,..,t^{k-1}) for f = m_nu in k vars, |nu|=d, n=k+d.
T_k f = sum_{|A|=k} c_A X_A f(X_A), N=n vars, exact symbolic t. lin_e G = (-1)^{n-1} n [p_n]G,
[p_n] m_lam = (-1)^{l-1}(l-1)!/prod m_i!. Independent of Thm W / (KF). Grade: computed."""
import sys, itertools, time
from math import factorial
from collections import Counter
from sympy.polys.rings import ring
from sympy import ZZ, symbols, cancel, factor
t = symbols('t')
def parts(n, m=None):
    if m is None: m = n
    if n == 0: yield (); return
    for k in range(min(n, m), 0, -1):
        for p in parts(n-k, k): yield (k,)+p
def mono_sym(Rg, xs, nu):
    nu = tuple(nu)+(0,)*(len(xs)-len(nu)); s = Rg(0)
    for perm in set(itertools.permutations(nu)):
        m = Rg(1)
        for x, e in zip(xs, perm): m *= x**e
        s += m
    return s
def pn_coef(lam):
    l = len(lam); d = 1
    for v in Counter(lam).values(): d *= factorial(v)
    return (-1)**(l-1)*factorial(l-1)/d if True else None
def lin_e_Tk(k, nu):
    d = sum(nu); n = k+d; N = n
    Rg, *g = ring(['t']+[f'x{i}' for i in range(N)], ZZ); T, xs = g[0], g[1:]
    V = Rg(1)
    for i in range(N):
        for j in range(i+1, N): V *= xs[i]-xs[j]
    tot = Rg(0)
    for A in itertools.combinations(range(N), k):
        Bc = [j for j in range(N) if j not in A]
        sign = (-1)**sum(1 for i in A for j in Bc if i > j); term = Rg(sign)
        for i in range(N):
            for j in range(i+1, N):
                if (i in A) == (j in A): term *= xs[i]-xs[j]
        for i in A:
            term *= xs[i]
            for j in Bc: term *= xs[i]-T*xs[j]
        tot += term*mono_sym(Rg, [xs[i] for i in A], nu)
    q, r = tot.div(V); assert r == 0
    from fractions import Fraction
    dct = q.to_dict(); val = 0
    for lam in parts(n):
        key = tuple(lam)+(0,)*(N-len(lam)); c = 0
        for m, v in dct.items():
            if m[1:] == key: c += int(v)*t**m[0]
        l = len(lam); dd = 1
        for v in Counter(lam).values(): dd *= factorial(v)
        val += c*Fraction((-1)**(l-1)*factorial(l-1), dd)
    return cancel((-1)**(n-1)*n*val), n
def rhs(k, nu):
    d = sum(nu); n = k+d
    pts = [t**i for i in range(k)]
    nu2 = tuple(nu)+(0,)*(k-len(nu)); f = 0
    for perm in set(itertools.permutations(nu2)):
        m = 1
        for x, e in zip(pts, perm): m *= x**e
        f += m
    return cancel((-1)**d*(1-t**n)/(1-t**k)*f)
if __name__ == '__main__':
    nmax = int(sys.argv[1]); ok = bad = 0
    for n in range(2, nmax+1):
        for k in range(1, n):
            for nu in parts(n-k):
                if len(nu) > k: continue
                t0 = time.time(); L, _ = lin_e_Tk(k, nu); R = rhs(k, nu)
                good = cancel(L-R) == 0; ok += good; bad += (not good)
                print(f'n={n} k={k} f=m{nu} [{time.time()-t0:.1f}s] lin_e T_k f = {factor(L)}  {"OK" if good else "FAIL rhs="+str(factor(R))}', flush=True)
    print(f'TOTAL OK {ok} FAIL {bad}')

"""Day 220: B(f,g) = d/ds (f*g)|_{s=1}, Hikita star. Engine copied from day217/fast.py but with general t (fixed rational):
e_k * F = E_k F,  E_k F = sum_A prod_{i in A, j notin A}(x_i - t x_j)/(x_i - x_j) X_A F(X_{A^c}, s X_A);  e^*_rho = E_{rho1}..E_{rhol}(1).
Test: Leibniz B(f,gh)=B(f,g)h+gB(f,h), symmetry B(f,g)=B(g,f); report B(e_a,e_b)."""
import sys, itertools
from sympy.polys.rings import ring
from sympy import QQ, Matrix, Symbol, Rational, diff, factor, expand, simplify
out = open(sys.argv[1], 'w'); T = Rational(sys.argv[2]); MAXD = int(sys.argv[3])
def log(*a): print(*a, flush=True); print(*a, file=out, flush=True)
S = Symbol('s')
def parts(n, m=None):
    if m is None: m = n
    if n == 0: yield (); return
    for k in range(min(n, m), 0, -1):
        for p in parts(n-k, k): yield (k,)+p
RINGS = {}
def R(N):
    if N not in RINGS:
        Rg, *g = ring(['s']+[f'x{i}' for i in range(N)], QQ); s, xs = g[0], g[1:]
        V = Rg(1)
        for i in range(N):
            for j in range(i+1, N): V *= xs[i]-xs[j]
        RINGS[N] = (Rg, s, xs, V)
    return RINGS[N]
def Ek(k, F, N):
    Rg, s, xs, V = R(N); tot = Rg(0); Fd = F.to_dict()
    for A in itertools.combinations(range(N), k):
        B = [j for j in range(N) if j not in A]
        sign = (-1)**sum(1 for i in A for j in B if i > j); term = Rg(sign)
        for i in range(N):
            for j in range(i+1, N):
                if (i in A) == (j in A): term *= xs[i]-xs[j]
        for i in A:
            term *= xs[i]
            for j in B: term *= xs[i]-T*xs[j]
        Fs = Rg.from_dict({(m[0]+sum(m[1+i] for i in A),)+m[1:]: c for m, c in Fd.items()}) if Fd else Rg(0)
        tot += term*Fs
    q, r = tot.div(V); assert r == 0
    return q
EMAT = {}
def emat(n):
    if n not in EMAT:
        P = list(parts(n)); N = n; Rg, s, xs, V = R(N)
        M = Matrix(len(P), len(P), lambda i, j: 0)
        for i, mu in enumerate(P):
            d = epoly(mu, N).to_dict()
            for j, lam in enumerate(P): M[i, j] = d.get((0,)+tuple(lam)+(0,)*(N-len(lam)), 0)
        EMAT[n] = (P, M.inv())
    return EMAT[n]
def epoly(mu, N):
    Rg, s, xs, V = R(N); p = Rg(1)
    for k in mu:
        p *= sum((Rg(1) if False else __import__('functools').reduce(lambda a, b: a*b, [xs[i] for i in A], Rg(1))) for A in itertools.combinations(range(N), k))
    return p
def to_e(F, n):  # F symmetric poly in N=n vars, degree n -> dict mu -> sympy expr in s
    P, Minv = emat(n); d = F.to_dict(); Rg = R(n)[0]
    mvec = []
    for lam in P:
        key = tuple(lam)+(0,)*(n-len(lam)); c = 0
        for m, v in d.items():
            if m[1:] == key: c += QQ.to_sympy(v)*S**m[0]
        mvec.append(c)
    return {mu: expand(sum(mvec[j]*Minv[j, i] for j in range(len(P)))) for i, mu in enumerate(P) if expand(sum(mvec[j]*Minv[j, i] for j in range(len(P)))) != 0}
def from_e(f, N):  # dict -> poly (coeffs must be s-free rationals)
    Rg = R(N)[0]; p = Rg(0)
    for mu, c in f.items(): p += Rg(QQ.from_sympy(Rational(c)))*epoly(mu, N)
    return p
GINV = {}
def ginv(n):  # e_lam = sum_rho Ginv[lam][rho] e*_rho
    if n not in GINV:
        P = list(parts(n)); G = Matrix(len(P), len(P), lambda i, j: 0)
        for i, rho in enumerate(P):
            F = R(n)[0](1)
            for k in reversed(rho): F = Ek(k, F, n)
            d = to_e(F, n)
            for j, nu in enumerate(P): G[i, j] = d.get(nu, 0)
        Gi = G.inv()
        GINV[n] = {lam: {rho: simplify(Gi[j, i]) for i, rho in enumerate(P) if simplify(Gi[j, i]) != 0} for j, lam in enumerate(P)}
    return GINV[n]
def star(f, g):  # f,g homogeneous dicts with rational coeffs
    a = sum(next(iter(f))); b = sum(next(iter(g))); N = a+b
    G = from_e(g, N); res = {}
    for lam, c in f.items():
        for rho, w in ginv(a)[lam].items():
            F = G
            for k in reversed(rho): F = Ek(k, F, N)
            for mu, v in to_e(F, N).items(): res[mu] = res.get(mu, 0) + c*w*v
    return {mu: simplify(v) for mu, v in res.items() if simplify(v) != 0}
def Bf(f, g):
    st = star(f, g); return {mu: simplify(diff(v, S).subs(S, 1)) for mu, v in st.items() if simplify(diff(v, S).subs(S, 1)) != 0}
def mul(f, g):
    r = {}
    for a, c in f.items():
        for b, d in g.items():
            k = tuple(sorted(a+b, reverse=True)); r[k] = r.get(k, 0)+c*d
    return {k: v for k, v in r.items() if v != 0}
def add(f, g):
    r = dict(f)
    for k, v in g.items(): r[k] = r.get(k, 0)+v
    return {k: simplify(v) for k, v in r.items() if simplify(v) != 0}
def neg(f): return {k: -v for k, v in f.items()}
log(f't = {T}, max total degree {MAXD}')
# sanity: s=1 star = ordinary product
st = star({(1,): 1}, {(2,): 1}); log('sanity e1*e2 at s=1:', {k: simplify(v.subs(S, 1)) for k, v in st.items()})
log('== B(e_a,e_b) ==')
for n in range(2, MAXD+1):
    for a in range(1, n):
        b = n-a
        if a <= 3 and b <= 3: log(f'B(e{a},e{b}) =', Bf({(a,): 1}, {(b,): 1}))
log('== symmetry B(f,g)=B(g,f), f,g monomials in e1,e2,e3 ==')
mons = [mu for n in range(1, MAXD) for mu in parts(n) if max(mu) <= 3]
nsym = okc = 0
for i, f in enumerate(mons):
    for g in mons[i+1:]:
        if sum(f)+sum(g) > MAXD: continue
        r = add(Bf({f: 1}, {g: 1}), neg(Bf({g: 1}, {f: 1}))); nsym += 1; okc += (not r)
        if r: log(f'  SYM FAIL f=e{f} g=e{g}: residual {r}')
log(f'SYMMETRY: {okc}/{nsym} pass')
log('== Leibniz B(f,gh)=B(f,g)h+gB(f,h) ==')
nl = okl = 0
for f in mons:
    for g in mons:
        for h in mons:
            if g > h or sum(f)+sum(g)+sum(h) > MAXD: continue
            gh = mul({g: 1}, {h: 1})
            lhs = Bf({f: 1}, gh); rhs = add(mul(Bf({f: 1}, {g: 1}), {h: 1}), mul({g: 1}, Bf({f: 1}, {h: 1})))
            r = add(lhs, neg(rhs)); nl += 1; okl += (not r)
            log(f'  f=e{f} g=e{g} h=e{h}: {"OK" if not r else "FAIL residual " + str(r)}')
log(f'LEIBNIZ: {okl}/{nl} pass')

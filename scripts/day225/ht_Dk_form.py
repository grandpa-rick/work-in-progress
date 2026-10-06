"""Wake 225: check the HT-Thm-7.3 reformulation of G:
Lead_{lam,(n)} = (-1)^{l-1} t^{-n(lam')} [n]_t/[lam_l]_t * < D_{lam_1}...D_{lam_{l-1}} h~_{lam_l}, e_n >|_{q=t}
with D_k the HT derivation (Def 3.x / Thm 1.1) on generators h~_m, and <h~_m,e_m> = q^{C(m,2)} (ring hom). Grade: computed."""
import sympy as sp, itertools
q, t = sp.symbols('q t'); hs = sp.symbols('h1:8')
def D(k, F): return sp.expand(sum(sp.diff(F, hs[m-1])*(hs[m+k-1]-hs[m-1]*hs[k-1])/(q**k-1) for m in range(1, 8-k)))
def parts(n, m=None):
    if m is None: m = n
    if n == 0: yield (); return
    for k in range(min(n, m), 0, -1):
        for p in parts(n-k, k): yield (k,)+p
def K_G(lam):
    l = len(lam); E = list(itertools.combinations(range(l), 2)); tot = 0
    for msk in range(1 << len(E)):
        H = [E[i] for i in range(len(E)) if msk >> i & 1]; par = list(range(l))
        def f(x):
            while par[x] != x: x = par[x]
            return x
        for a, b in H: par[f(a)] = f(b)
        if len({f(x) for x in range(l)}) == 1: tot += sp.prod([t**(lam[i]*lam[j])-1 for i, j in H])
    return tot
br = lambda m: sum(t**i for i in range(m))
ok = True
for n in range(2, 7):
    for lam in parts(n):
        l = len(lam)
        if l < 2: continue
        F = hs[lam[-1]-1]
        for k in reversed(lam[:-1]): F = D(k, F)
        val = F.subs({hs[m-1]: q**(m*(m-1)//2) for m in range(1, 8)}).subs(q, t)
        lhs = (-1)**(l-1)*t**(-sum(k*(k-1)//2 for k in lam))*br(n)/br(lam[-1])*val
        G = (1-t**n)/sp.prod([1-t**k for k in lam])*K_G(lam)
        r = sp.simplify(lhs-G) == 0; ok &= r; print(lam, r)
print('ALL', ok)

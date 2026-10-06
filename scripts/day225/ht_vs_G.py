"""Wake 225: separator test, Haglund-Tewari arXiv:2609.29957 Def 7.1/Thm 7.3 vs Rick's Theorem G.
Three columns for every lam |- n, 2<=n<=5, l(lam)>=2, t symbolic:
 (A) RAW: Lead_{lam,(n)} = [(s-1)^{l-1}] c_{lam,(n)}(s,t) from the raw star engine pickles (wake222/gprime/cst_n*.pkl).
 (B) G:   (1-t^n)/prod(1-t^{lam_i}) * K_lam(t), K_lam = sum over connected graphs prod (t^{lam_i lam_j}-1).
 (C) HT:  HT's kappa(lam) built as a symmetric function in the p-basis from Def 7.1 (h~_m = (q;q)_m h_m[X/(1-q)]),
          paired with e_n (phi(p_k)=(-1)^{k-1}, a ring hom), q:=t, times t^{-n(lam')} (t-1)^{l-1} and the G prefactor.
Also: check HT Ex 7.2 kappa(1,1)=e_2, HT Thm 7.3 numerically (D_k derivation on h~ generators), and the
separator I = Lead(111)/(Lead(21)Lead(11)) for raw-star, for HT-kappa paired with e_n (no prefactor), and Dolega-oplus.
Grade: computed."""
import itertools, pickle, os, sympy as sp
from sympy.utilities.iterables import multiset_partitions
from sympy.combinatorics.partitions import IntegerPartition
q, t, s = sp.symbols('q t s')
P = sp.symbols('p1:7')  # p_1..p_6
def parts(n, m=None):
    if m is None: m = n
    if n == 0: yield (); return
    for k in range(min(n, m), 0, -1):
        for p in parts(n-k, k): yield (k,)+p
def zee(mu):
    from collections import Counter
    z = 1
    for k, m in Counter(mu).items(): z *= k**m * sp.factorial(m)
    return z
def qpoch(m): return sp.prod([1-q**i for i in range(1, m+1)])
def htilde(m):  # (q;q)_m h_m[X/(1-q)] in p-basis
    if m == 0: return sp.Integer(1)
    return sp.expand(qpoch(m)*sum(sp.prod([P[k-1] for k in mu])/(zee(mu)*sp.prod([1-q**k for k in mu])) for mu in parts(m)))
H = {m: sp.factor(htilde(m)) for m in range(0, 7)}
def setparts(r): return list(multiset_partitions(list(range(r))))
def varkappa(a):
    r = len(a); tot = 0
    for sig in setparts(r):
        b = len(sig); mob = (-1)**(b-1)*sp.factorial(b-1)
        tot += mob*sp.prod([H[sum(a[i] for i in B)] for B in sig])
    return tot
def kappa(a): return sp.together(varkappa(a)/(q-1)**(len(a)-1))
def phi(f):  # <f, e_n> on homogeneous pieces summed: ring hom p_k -> (-1)^{k-1}
    return sp.simplify(f.subs({P[k-1]: (-1)**(k-1) for k in range(1, 7)}))
def K_G(lam):
    l = len(lam); E = list(itertools.combinations(range(l), 2)); tot = 0
    for m in range(1 << len(E)):
        Hh = [E[i] for i in range(len(E)) if m >> i & 1]
        par = list(range(l))
        def f(x):
            while par[x] != x: x = par[x]
            return x
        for a, b in Hh: par[f(a)] = f(b)
        if len({f(x) for x in range(l)}) == 1:
            tot += sp.prod([t**(lam[i]*lam[j])-1 for i, j in Hh])
    return sp.expand(tot)
def pref(lam): n = sum(lam); return (1-t**n)/sp.prod([1-t**k for k in lam])
def nlp(lam): return sum(k*(k-1)//2 for k in lam)  # n(lam') = sum C(lam_i,2)
out = []
def log(*a):
    print(*a, flush=True); out.append(' '.join(map(str, a)))
# sanity: HT Example 7.2
e2 = (P[0]**2 - P[1])/2
log('HT Ex7.2 kappa(1,1)==e_2:', sp.simplify(kappa((1, 1)) - e2) == 0)
log('phi(h~_m)==q^C(m,2), m<=6:', all(sp.simplify(phi(H[m]) - q**(m*(m-1)//2)) == 0 for m in range(7)))
# HT Thm 7.3 check: D_k derivation on generators h~_m, applied to varkappa(a) as polynomial in h~'s
hs = sp.symbols('h1:7')
def varkappa_h(a):
    r = len(a); tot = 0
    for sig in setparts(r):
        b = len(sig); tot += (-1)**(b-1)*sp.factorial(b-1)*sp.prod([hs[sum(a[i] for i in B)-1] for B in sig])
    return sp.expand(tot)
def Dk_h(k, F):
    return sp.expand(sum(sp.diff(F, hs[m-1])*(hs[m+k-1]-hs[m-1]*hs[k-1])/(q**k-1) for m in range(1, 7-k) ))
ok73 = all(sp.simplify((q**k-1)*Dk_h(k, varkappa_h(a)) - varkappa_h(a+(k,))) == 0
           for a in [(1,), (2,), (1, 1), (2, 1), (1, 1, 1)] for k in (1, 2) if sum(a)+k <= 6)
log('HT Thm 7.3 (q^k-1)D_k varkappa(a)=varkappa(a,k), checked cases:', ok73)
raw = {}
base = '/home/agent/projects/proofs/scripts/wake222/gprime'
for n in range(2, 6):
    fn = os.path.join(base, f'cst_n{n}.pkl')
    if os.path.exists(fn): raw.update(pickle.load(open(fn, 'rb')))
lead = {}
log('\nlam | RAW star Lead | G formula | HT-derived | RAW==G | G==HT')
allok = True
for n in range(2, 6):
    for lam in parts(n):
        l = len(lam)
        if l < 2: continue
        G = sp.factor(pref(lam)*K_G(lam))
        htK = sp.factor(t**(-nlp(lam))*(t-1)**(l-1)*phi(kappa(lam)).subs(q, t))
        HT = sp.factor(pref(lam)*htK)
        R = None
        if (lam, (n,)) in raw:
            c = sp.expand(raw[(lam, (n,))].subs(s, s+1))  # shift so s -> s-1
            R = sp.factor(sp.Poly(c, s).coeff_monomial(s**(l-1)))
            assert all(sp.Poly(c, s).coeff_monomial(s**j) == 0 for j in range(l-1)), lam
        lead[lam] = (R, G, HT)
        a = (R is not None and sp.simplify(R-G) == 0); b = sp.simplify(G-HT) == 0
        allok &= b and (R is None or a)
        log(lam, '|', R, '|', G, '|', HT, '|', a if R is not None else 'n/a', '|', b)
log('\nALL rows consistent:', allok)
# separator I
I_raw = sp.factor(lead[(1, 1, 1)][0]/(lead[(2, 1)][0]*lead[(1, 1)][0]))
ph = lambda lam: phi(kappa(lam)).subs(q, t)
I_HTpair = sp.factor(ph((1, 1, 1))/(ph((2, 1))*ph((1, 1))))
log('separator I (raw star Lead) =', I_raw)
log('separator I (HT kappa paired with e_n, no normalization) =', I_HTpair)
log('HT kappa paired with e_n, raw values: (1,1):', sp.factor(ph((1, 1))), ' (2,1):', sp.factor(ph((2, 1))), ' (1,1,1):', sp.factor(ph((1, 1, 1))))
# a varying-parameter spot value
log('t=3/5 spot (2,1,1): RAW', lead[(2, 1, 1)][0].subs(t, sp.Rational(3, 5)), ' HT-derived', lead[(2, 1, 1)][2].subs(t, sp.Rational(3, 5)))
open('/home/agent/projects/proofs/scripts/day225/ht_vs_G.log', 'w').write('\n'.join(out)+'\n')

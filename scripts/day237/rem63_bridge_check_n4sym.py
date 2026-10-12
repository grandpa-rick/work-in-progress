# Day 237 (a), n=4 SYMBOLIC completion (re-run session, 2026-10-12): P3 and P4 of rem63_bridge_check.py at n=4,
# exactly over QQ(s,t). The original symbolic run was OOM in sp.solve; here the operators are 5x5 matrices in
# p-coordinates over the fraction field QQ(s,t) (DomainMatrix), so P3 is checked on ALL of Lambda^4 at once:
#   N := Pmat diag(t^n(nu) s^n(nu')) Pmat^-1  vs  RHS := Dphi^-1 Hmat diag(same) Hmat^-1 Dphi,
# Pmat = p-coords of P_nu(x;s,1/t) (Gram-Schmidt, as before), Hmat = p-coords of GH Htilde_nu(x;s,t) (as before),
# Dphi = diag over mu of prod_r eps(r)/(1-t^r) (phi in p-coords, printed sign p_r[-eps X]=(-1)^{r-1} p_r).
# P4: t^n(lam') e*_lam (printed T,pi,Y in m=4 vars) vs RHS applied to e_lam, symbolic. Same controls.
import sys; sys.argv = ['x', '4']
exec(open('rem63_bridge_check.py').read().split('# ---------------- run')[0])
from sympy.polys.matrices import DomainMatrix
K = sp.QQ.frac_field(s, t)
results = []
def rep(name, ok):
    results.append((name, ok)); print(('PASS ' if ok else 'FAIL ') + name, flush=True)
n = 4; parts = list(partitions(n)); mp = m_in_p(n); S = {l: schur(l, n) for l in parts}
Pq = macdonaldP(n, q, t, mp); Ht = {l: gh_Htilde(l, n, S) for l in parts}
def col(F): c = coeffs(F, n); return [K.from_sympy(sp.cancel(c[mu])) for mu in parts]
def DM(cols): return DomainMatrix([[cols[j][i] for j in range(len(cols))] for i in range(len(parts))], (len(parts), len(cols)), K)
def diag(vals): return DomainMatrix([[vals[i] if i == j else K.zero for j in range(len(vals))] for i in range(len(vals))], (len(vals),) * 2, K)
Pmat = DM([col(sp.cancel(Pq[l].subs({q: s}).subs(t, 1 / t))) for l in parts])
Hmat = DM([col(Ht[l].subs(q, s)) for l in parts])
Hsw = DM([col(Ht[l].subs({q: t, t: s}, simultaneous=True)) for l in parts])
ev = lambda f: diag([K.from_sympy(f(l)) for l in parts])
eN = ev(lambda l: t ** nn(l) * s ** nn(conj(l))); bad = ev(lambda l: t ** nn(conj(l)) * s ** nn(l))
swap = ev(lambda l: s ** nn(l) * t ** nn(conj(l)))
def Dphi(tt, sg=True): return ev(lambda mu: reduce(lambda a, b: a * b, [(eps_sign(r) if sg else 1) / (1 - tt ** r) for r in mu], sp.Integer(1)))
Nmat = Pmat * eN * Pmat.inv()
def rhs(E, H=Hmat, tt=t, sg=True):
    D = Dphi(tt, sg); return D.inv() * H * E * H.inv() * D
R = rhs(eN)
rep('P3 n=4 symbolic: N = phi^-1 nabla_{q=s} phi on all of Lambda^4 (5x5 over QQ(s,t))', Nmat == R)
rep('  ctrl n(nu)<->n(nu\') eigenvalue breaks', Nmat != rhs(bad))
rep('  ctrl s<->t (whole RHS operator) breaks', Nmat != rhs(swap, Hsw, s))
rep('  ctrl phi without -eps breaks', Nmat != rhs(eN, Hmat, t, False))
def ebasis(l):
    e = lambda k: sum((-1) ** (k - len(mu)) * pp(mu) / zee(mu) for mu in partitions(k))
    return sp.expand(reduce(lambda a, b: a * b, [e(k) for k in l], 1))
for lam in parts:
    v = DM([col(ebasis(lam))]); w = R * v
    F = sum(K.to_sympy(w[i, 0].element) * pp(mu) for i, mu in enumerate(parts))
    Fst, xs = star_e(lam, 4)
    rx = to_poly_in_x(F, xs)
    rep('P4 n=4 symbolic lam=%s t^n(lam\') e*_lam = phi^-1 nabla phi e_lam (m=4, printed T,pi,Y)' % (lam,),
        sp.cancel(t ** nn(conj(lam)) * Fst - rx) == 0)
    if lam != conj(lam):
        rep('  ctrl P4 t^n(lam) breaks', sp.cancel(t ** nn(lam) * Fst - rx) != 0)
print('ALL', all(ok for _, ok in results), len(results))

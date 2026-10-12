# Day 237 (a), n=4 completion: P3 and P4 of rem63_bridge_check.py at n=4, in exact rational arithmetic at
# three random rational points (s,t) (the symbolic n=4 P3 solve ran out of memory; P1, P2 at n=4 passed symbolically,
# see rem63_bridge_check.log). Same constructions (GH Htilde, Gram-Schmidt P, printed Sec.2 operators), same controls.
import sys; sys.argv = ['x', '4']
exec(open('rem63_bridge_check.py').read().split('# ---------------- run')[0])
results = []
def rep(name, ok):
    results.append((name, ok)); print(('PASS ' if ok else 'FAIL ') + name, flush=True)
n = 4; parts = list(partitions(n)); mp = m_in_p(n); S = {l: schur(l, n) for l in parts}
Pq = macdonaldP(n, q, t, mp); Ht = {l: gh_Htilde(l, n, S) for l in parts}
def ebasis(l):
    e = lambda k: sum((-1) ** (k - len(mu)) * pp(mu) / zee(mu) for mu in partitions(k))
    return sp.expand(reduce(lambda a, b: a * b, [e(k) for k in l], 1))
def numexpand(F, basis):
    parts_ = list(basis); A = sp.Matrix([[coeffs(basis[l], n)[mu] for l in parts_] for mu in partitions(n)])
    b = sp.Matrix([coeffs(F, n)[mu] for mu in partitions(n)]); x = A.LUsolve(b)
    return dict(zip(parts_, x))
for (sv, tv) in [(sp.Rational(3, 7), sp.Rational(5, 11)), (sp.Rational(-2, 5), sp.Rational(7, 3)), (sp.Rational(9, 4), sp.Rational(-4, 13))]:
    tag = 's=%s,t=%s' % (sv, tv)
    pt = {s: sv, t: tv}
    Ps = {l: sp.expand(Pq[l].subs({q: sv, t: 1 / tv}, simultaneous=True)) for l in parts}
    Hs = {l: sp.expand(Ht[l].subs({q: sv, t: tv}, simultaneous=True)) for l in parts}
    Hsw = {l: sp.expand(Ht[l].subs({q: tv, t: sv}, simultaneous=True)) for l in parts}   # for the s<->t control
    ev = lambda f: {l: f(l).subs(pt) for l in parts}
    eN = ev(lambda l: t ** nn(l) * s ** nn(conj(l)))
    bad = ev(lambda l: t ** nn(conj(l)) * s ** nn(l))
    def ph(F, tt, sg=True): return pleth(F, lambda r: (eps_sign(r) if sg else 1) / (1 - tt ** r))
    def phi_(F, tt, sg=True): return pleth(F, lambda r: (eps_sign(r) if sg else 1) * (1 - tt ** r))
    def Nop(F):
        a = numexpand(F, Ps); return sp.expand(sum(a[l] * eN[l] * Ps[l] for l in parts))
    def cn(F, evd, H=Hs, tt=tv, sg=True):
        a = numexpand(ph(F, tt, sg), H); return phi_(sp.expand(sum(a[l] * evd[l] * H[l] for l in parts)), tt, sg)
    swap_ev = {l: (s ** nn(l) * t ** nn(conj(l))).subs({s: sv, t: tv}) for l in parts}  # = sigma(eN) at the point
    bad_any = {'nn': False, 'st': False, 'eps': False}
    for lam in parts:
        E = ebasis(lam); NE = Nop(E)
        rep('P3 n=4 %s lam=%s' % (tag, lam), sp.expand(NE - cn(E, eN)) == 0)
        bad_any['nn'] |= sp.expand(NE - cn(E, bad)) != 0
        bad_any['st'] |= sp.expand(NE - cn(E, swap_ev, Hsw, sv)) != 0   # whole RHS operator with s<->t
        bad_any['eps'] |= sp.expand(NE - cn(E, eN, Hs, tv, False)) != 0
    rep('  ctrl n(nu)<->n(nu\') breaks (%s)' % tag, bad_any['nn'])
    rep('  ctrl s<->t breaks (%s)' % tag, bad_any['st'])
    rep('  ctrl phi without -eps breaks (%s)' % tag, bad_any['eps'])
    # P4 in m=4 variables at the point
    for lam in parts:
        Fst, xs = star_e(lam, 4)
        lhs = sp.expand((t ** nn(conj(lam)) * Fst).subs(pt))
        rhs = to_poly_in_x(cn(ebasis(lam), eN), xs)
        rep('P4 n=4 %s lam=%s (m=4, printed T,pi,Y)' % (tag, lam), sp.expand(lhs - rhs) == 0)
        if lam != conj(lam):
            rep('  ctrl P4 t^n(lam) breaks', sp.expand((t ** nn(lam) * Fst).subs(pt) - rhs) != 0)
print('ALL', all(ok for _, ok in results), len(results))

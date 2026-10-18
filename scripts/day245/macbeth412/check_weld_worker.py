#!/usr/bin/env python3
"""Rick's check of MacBeth UID 412 (his WIP 585de70): worker = constant-branch sub-bicomodule.

Everything is built from scratch as concrete Poly data (no library):
  polynomial p      : positions(), dirs(i)
  morphism  p -> q  : i |-> (j, back) with back: q[j] -> p[i] a dict
  composite p <| q  : positions (i, f) with f: p[i] -> q(1) frozen as a sorted tuple;
                      directions (d, e) with d in p[i], e in q[f(d)].
Coactions are computed from the comonoid data by whiskering, exactly as MacBeth defines them
(lambda = delta_S <| St_T, rho = St_S <| delta_T), then re-associated; all bicomodule laws are
checked pointwise (target position + back map on every target direction).
"""
import itertools as it, sys

STAR = '*'

def fz(d):            # freeze a dict
    return tuple(sorted(d.items()))

def th(f):            # thaw
    return dict(f)

class Poly:
    def __init__(self, name, positions, dirs):
        self.name, self._pos, self._dirs = name, positions, dirs
    def positions(self): return self._pos()
    def dirs(self, i): return self._dirs(i)

Y = Poly('y', lambda: [STAR], lambda i: [STAR])

def St(X, name):
    X = list(X)
    return Poly(name, lambda: X, lambda x: X)

def comp(p, q):
    """p <| q, positions enumerated lazily (only used for small sets)."""
    def pos():
        for i in p.positions():
            D = p.dirs(i)
            for vals in it.product(q.positions(), repeat=len(D)):
                yield (i, fz(dict(zip(D, vals))))
    def dirs(ij):
        i, f = ij; f = th(f)
        return [(d, e) for d in p.dirs(i) for e in q.dirs(f[d])]
    return Poly(f'({p.name}<|{q.name})', lambda: list(pos()), dirs)

# ---- comonoid data of St_X (indiscrete category): morphism out of x = its target
def eps_St(x):  return (STAR, {STAR: x})
def delta_St(X):
    X = list(X)
    return lambda x: ((x, fz({a: a for a in X})), {(a, b): b for a in X for b in X})

# ---- morphism combinators (morphisms are python functions i |-> (j, back dict))
def compose(psi, phi):            # psi o phi
    def m(i):
        j, b1 = phi(i); k, b2 = psi(j)
        return (k, {e: b1[b2[e]] for e in b2})
    return m

def whisk_right(phi, p_src, q):   # phi <| q : p<|q -> p'<|q
    def m(ifx):
        i, f = ifx; f = th(f)
        j, b = phi(i)
        f2 = {d2: f[b[d2]] for d2 in b}
        back = {(d2, e): (b[d2], e) for d2 in b for e in q.dirs(f2[d2])}
        return ((j, fz(f2)), back)
    return m

def whisk_left(q, phi, p_src, p_tgt):   # q <| phi : q<|p -> q<|p'
    def m(kf):
        k, f = kf; f = th(f)
        f2, back = {}, {}
        for d in q.dirs(k):
            j, b = phi(f[d]); f2[d] = j
            for e2 in p_tgt.dirs(j): back[(d, e2)] = (d, b[e2])
        return ((k, fz(f2)), back)
    return m

# ---- associators  (p<|q)<|r  <->  p<|(q<|r)  on a morphism's target
def assoc_LtoR(phi, p, q, r):
    """phi lands in (p<|q)<|r; return the same map landing in p<|(q<|r)."""
    def m(i):
        ((pi, f), h), back = phi(i); f = th(f); h = th(h)
        g = {d: (f[d], fz({e: h[(d, e)] for e in q.dirs(f[d])})) for d in p.dirs(pi)}
        nb = {}
        for d in p.dirs(pi):
            for e in q.dirs(f[d]):
                for x in r.dirs(h[(d, e)]):
                    nb[(d, (e, x))] = back[((d, e), x)]
        return ((pi, fz(g)), nb)
    return m

def assoc_RtoL(phi, p, q, r):
    """phi lands in p<|(q<|r); return the same map landing in (p<|q)<|r."""
    def m(i):
        (pi, g), back = phi(i); g = th(g)
        f = {d: g[d][0] for d in p.dirs(pi)}
        h = {}
        for d in p.dirs(pi):
            qd, hd = g[d]; hd = th(hd)
            for e in q.dirs(qd): h[(d, e)] = hd[e]
        nb = {}
        for d in p.dirs(pi):
            for e in q.dirs(f[d]):
                for x in r.dirs(h[(d, e)]):
                    nb[((d, e), x)] = back[(d, (e, x))]
        return (((pi, fz(f)), fz(h)), nb)
    return m

def unit_left(phi):    # target y<|m  ->  m
    def m(i):
        (_, f), back = phi(i); f = th(f)
        return (f[STAR], {e: back[(STAR, e)] for e in {k[1] for k in back}})
    return m

def unit_right(phi):   # target m<|y  ->  m
    def m(i):
        (j, _), back = phi(i)
        return (j, {d: back[(d, STAR)] for d in {k[0] for k in back}})
    return m

def same(phi, psi, positions):
    for i in positions:
        a, b = phi(i), psi(i)
        if a[0] != b[0] or a[1] != b[1]:
            return False, i
    return True, None

IDM = lambda poly: (lambda i: (i, {d: d for d in poly.dirs(i)}))

def bicomodule_laws(m, d, c, lam, rho, P, epsd, deltad, epsc, deltac):
    """m: (c,d)-bicomodule in MacBeth's Def 1: lam: m -> d<|m, rho: m -> m<|c. P = positions tested."""
    res = {}
    res['lam counit'] = same(unit_left(compose(whisk_right(epsd, d, m), lam)), IDM(m), P)[0]
    lhs = assoc_LtoR(compose(whisk_right(deltad, d, m), lam), d, d, m)
    rhs = compose(whisk_left(d, lam, m, comp(d, m)), lam)
    res['lam coassoc'] = same(lhs, rhs, P)[0]
    res['rho counit'] = same(unit_right(compose(whisk_left(m, epsc, c, Y), rho)), IDM(m), P)[0]
    lhs = assoc_LtoR(compose(whisk_right(rho, m, c), rho), m, c, c)
    rhs = compose(whisk_left(m, deltac, c, comp(c, c)), rho)
    res['rho coassoc'] = same(lhs, rhs, P)[0]
    lhs = assoc_LtoR(compose(whisk_right(lam, m, c), rho), d, m, c)
    rhs = compose(whisk_left(d, rho, m, comp(m, c)), lam)
    res['compat'] = same(lhs, rhs, P)[0]
    return res

def lam_targets(lam, i):            # m-positions appearing in lam(i) in d<|m
    (k, f), _ = lam(i); return set(v for _, v in f)

def closed(sub, lam, rho):
    for i in sub:
        if not lam_targets(lam, i) <= sub: return False
        if rho(i)[0][0] not in sub: return False
    return True

def hom_count_const(m, rho, i, objs, V):
    """|Hom_{St_T-Set}(m[i], const V)|: functions on m[i] constant along rho's action."""
    D = m.dirs(i)
    (_, lab), back = rho(i); lab = th(lab)
    edges = [(e, back[(e, x)]) for e in D for x in objs]
    # union-find orbits
    par = {e: e for e in D}
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for a, b in edges: par[find(a)] = find(b)
    return V ** len({find(e) for e in D})

def Fbar(m, lam, rho, positions, s0, objsT, V):
    return sum(hom_count_const(m, rho, i, objsT, V) for i in positions if lam(i)[0][0] == s0)

# ======================================================================================
def run(nS, nT, out):
    S, T = list(range(nS)), list('abc'[:nT])
    StS, StT = St(S, 'StS'), St(T, 'StT')
    dS, dT = delta_St(S), delta_St(T)
    B = comp(StS, StT)
    PB = B.positions()
    const = {t: fz({a: t for a in S}) for t in T}
    gof = lambda p: th(p[1])                       # branch map g of position (s, g)
    # MacBeth's coactions, built by whiskering
    lamB = assoc_LtoR(whisk_right(dS, StS, StT), StS, StS, StT)         # B -> StS<|B
    rhoB = assoc_RtoL(whisk_left(StS, dT, StT, comp(StT, StT)), StS, StT, StT)  # B -> B<|StT
    w = f'({nS},{nT})'
    out(f'\n=== {w}: |B(1)|={len(PB)}, |B[i]|={nS*nT}')
    L = bicomodule_laws(B, StS, StT, lamB, rhoB, PB, eps_St, dS, eps_St, dT)
    out(f'  B laws: {L}'); assert all(L.values())
    # --- explicit formula check: lam(s,g) = (s, a |-> (a,g)); rho(s,g) = ((s,g), (a,t')|->t')
    ok_l = all(th(lamB(p)[0][1]) == {a: (a, p[1]) for a in S} and lamB(p)[0][0] == p[0] for p in PB)
    ok_r = all(rhoB(p)[0][0] == p and th(rhoB(p)[0][1]) == {(a, t): t for a in S for t in T} for p in PB)
    ok_rb = all(rhoB(p)[1][((a, t), t2)] == (a, t2) for p in PB for a in S for t in T for t2 in T)
    ok_lb = all(lamB(p)[1][(a, (a2, t))] == (a2, t) for p in PB for a in S for a2 in S for t in T)
    out(f'  lam(s,g)=(s,a->(a,g)) {ok_l}; lam#(a,(a\',t))=(a\',t) {ok_lb}; '
        f'rho(s,g)=((s,g),(a,t)->t) {ok_r}; rho#((a,t),t\'\')=(a,t\'\') {ok_rb}')
    assert ok_l and ok_r and ok_rb and ok_lb
    # --- (i) B' closed
    Bp = {(s, const[t]) for s in S for t in T}
    cB = closed(Bp, lamB, rhoB); out(f'  (i) B\' = S x const(T) closed under lam, rho: {cB}'); assert cB
    Lp = bicomodule_laws(B, StS, StT, lamB, rhoB, sorted(Bp), eps_St, dS, eps_St, dT)
    assert all(Lp.values())
    # reachability classes: closed subsets are exactly unions of S x {g}
    cls = {}
    for p in PB:
        cls.setdefault(p[1], set()).add(p)
    comps_ok = all(lam_targets(lamB, p) == cls[p[1]] for p in PB) and all(rhoB(p)[0][0] == p for p in PB)
    out(f'  every slice S x {{g}} is a minimal closed set (so closed subsets = S x G, any G <= T^S): {comps_ok}')
    assert comps_ok
    NC = {p for p in PB if len(set(gof(p).values())) > 1}
    cNC = closed(NC, lamB, rhoB)
    out(f'  positive control: NC = S x (non-constant) also closed: {cNC} (|NC|={len(NC)})'); assert cNC
    # --- negative controls (must FAIL when there is room)
    s0 = S[0]
    nonsat = {(s0, const[t]) for t in T}
    c1 = closed(nonsat, lamB, rhoB)
    out(f'  NEG1 non-S-saturated {{(s0,const_t)}} closed: {c1}  (expected False iff |S|>=2)')
    assert c1 == (nS == 1)
    # NEG2: a *valid* twisted left coaction that moves g: T = Z/nT, S = Z/nS, w non-constant,
    # lam_tw(s,g) = (s, a |-> (a, g + (a-s) w)); needs nS * w == 0 mod nT -> use nS == nT.
    if nS == nT and nS >= 2:
        wv = {a: (1 if a == 0 else 0) for a in S}                  # non-constant
        Ti = {t: k for k, t in enumerate(T)}
        def shift(g, k):
            g = th(g); return fz({a: T[(Ti[g[a]] + k * wv[a]) % nT] for a in S})
        def lam_tw(p):
            s, g = p
            f = {a: (a, shift(g, (a - s) % nS)) for a in S}
            back = {(a, (a2, t)): (a2, t) for a in S for a2 in S for t in T}
            return ((s, fz(f)), back)
        Ltw = bicomodule_laws(B, StS, StT, lam_tw, rhoB, PB, eps_St, dS, eps_St, dT)
        ctw = closed(Bp, lam_tw, rhoB)
        out(f'  NEG2 twisted lam_tw laws {Ltw}; B\' closed under lam_tw: {ctw}  (expected laws all True, closure False)')
        assert all(Ltw.values()) and not ctw
    # NEG3: perturbed (invalid) lam_bad: bump g at slot a -> laws must fail, closure must fail
    if nS >= 2 and nT >= 2:
        def bump(g, a):
            g = th(g); g[a] = T[(T.index(g[a]) + 1) % nT]; return fz(g)
        def lam_bad(p):
            s, g = p
            return ((s, fz({a: (a, bump(g, a)) for a in S})),
                    {(a, (a2, t)): (a2, t) for a in S for a2 in S for t in T})
        Lbad = bicomodule_laws(B, StS, StT, lam_bad, rhoB, PB, eps_St, dS, eps_St, dT)
        cbad = closed(Bp, lam_bad, rhoB)
        out(f'  NEG3 bumped lam_bad laws {Lbad}; B\' closed: {cbad}  (expected counit False, closure False)')
        assert not Lbad['lam counit'] and not cbad
    # --- A with its own (cobase-change) structure: Remark 2 made explicit
    SxT = [(s, t) for s in S for t in T]
    A = St(SxT, 'A'); dA = delta_St(SxT)
    piS = lambda st: (st[0], {a: (a, st[1]) for a in S})          # cofunctor A -> St_S
    piT = lambda st: (st[1], {t: (st[0], t) for t in T})          # cofunctor A -> St_T
    lamA = compose(whisk_right(piS, A, A), dA)                    # A -> St_S <| A
    rhoA = compose(whisk_left(A, piT, A, StT), dA)                # A -> A <| St_T
    LA = bicomodule_laws(A, StS, StT, lamA, rhoA, SxT, eps_St, dS, eps_St, dT)
    out(f'  A (cobase change of regular A along projections) laws: {LA}'); assert all(LA.values())
    def mk_sigma(emb):
        return lambda st: ((st[0], emb[st[1]]), {(a, t): (a, t) for a in S for t in T})
    def is_bicomod_map(sig):
        l1 = compose(whisk_left(StS, sig, A, B), lamA); l2 = compose(lamB, sig)
        r1 = compose(whisk_right(sig, A, StT), rhoA)
        r2 = compose(rhoB, sig)
        return same(l1, l2, SxT)[0] and same(r1, r2, SxT)[0]
    sig = mk_sigma(const)
    okS = is_bicomod_map(sig)
    out(f'  Prop 1: sigma(s,t)=(s,const_t), sigma#=id is a bicomodule map A -> B: {okS}; '
        f'injective {len({sig(x)[0] for x in SxT})==len(SxT)}; image = B\' {set(sig(x)[0] for x in SxT)==Bp}')
    assert okS
    # any injection T -> T^S works equally well
    allg = sorted(set(p[1] for p in PB))
    nonc = [g for g in allg if len(set(th(g).values())) > 1]
    if len(nonc) >= nT:
        emb2 = {t: nonc[k] for k, t in enumerate(T)}
        ok2 = is_bicomod_map(mk_sigma(emb2))
        out(f'  sigma_G onto S x G, G = {nT} NON-constant maps, also a bicomodule iso onto image: {ok2}')
        assert ok2
    # --- (ii) M1: transported prafunctors at the St_S-object s0
    for V in (2, 3):
        FB = Fbar(B, lamB, rhoB, PB, s0, T, V)
        FBp = Fbar(B, lamB, rhoB, sorted(Bp), s0, T, V)
        FA = Fbar(A, lamA, rhoA, SxT, s0, T, V)
        assert FB == nT**nS * V**nS and FBp == nT * V**nS and FA == FBp
        out(f'  V={V}: F_B={FB} (=T^S V^S), F_B\'={FBp}, F_A(own)={FA} (=T V^S), content V^S={V**nS}')
    eqW = (nT**nS == nT); eqC = (nT == 1)
    out(f'  M1: F_A = F_B (worker = whole weld): {eqW}; predicted |S|=1 or |T|=1: {nS==1 or nT==1}; '
        f'F_A = content (-)^S: {eqC}; A = B carriers: {len(SxT)==len(PB)}')
    assert eqW == (nS == 1 or nT == 1) and (len(SxT) == len(PB)) == eqW
    # --- (iii) M2: retract e(s,g) = (s, const_{g(s)})
    e = lambda p: ((p[0], const[gof(p)[p[0]]]), {d: d for d in B.dirs(p)})
    idem = all(e(e(p)[0])[0] == e(p)[0] for p in PB)
    img = {e(p)[0] for p in PB} == Bp
    fib = {q: sum(1 for p in PB if e(p)[0] == q) for q in Bp}
    fib_ok = set(fib.values()) == {nT**(nS-1)}
    out(f'  M2: e idempotent {idem}; image = B\' {img}; fibre sizes {sorted(set(fib.values()))} = T^(S-1) {fib_ok}; '
        f'|NC|={len(NC)//nS} |T^S|={nT**nS}')
    assert idem and img and fib_ok
    eq_l = same(compose(whisk_left(StS, e, B, B), lamB), compose(lamB, e), PB)
    eq_r = same(compose(whisk_right(e, B, StT), rhoB), compose(rhoB, e), PB)[0]
    out(f'  e is rho-equivariant: {eq_r}; e is lam-equivariant: {eq_l[0]}'
        + (f'  <-- FAILS at {eq_l[1]}' if not eq_l[0] else ''))
    assert eq_l[0] == (nS == 1 or nT == 1)
    r0 = lambda p: ((p[0], const[gof(p)[s0]]), {d: d for d in B.dirs(p)})
    ok_r0 = same(compose(whisk_left(StS, r0, B, B), lamB), compose(lamB, r0), PB)[0] and \
            same(compose(whisk_right(r0, B, StT), rhoB), compose(rhoB, r0), PB)[0]
    sec = all(r0(p)[0] == p for p in Bp)
    out(f'  r_s0(s,g)=(s,const_g(s0)): bicomodule map {ok_r0}, retraction onto B\' {sec}, fibre {nT**(nS-1)}')
    assert ok_r0 and sec
    return True

if __name__ == '__main__':
    lines = []
    def out(x): print(x); lines.append(x)
    for sh in [(1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 2), (2, 3), (3, 2), (3, 3)]:
        run(*sh, out)
    out('\nALL ASSERTIONS PASSED')

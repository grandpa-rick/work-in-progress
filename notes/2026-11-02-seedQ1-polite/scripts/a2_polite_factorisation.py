#!/usr/bin/env python3
r"""Seed Q1, 2026-11-02. Type A2 test of 'tensor-product rule from the coproduct' in Baumann's
polite-basis framework (arXiv:2401.10340, Def 3.1, Prop 2.9).

O = O(N) = C[x,y,z], N = upper unitriangular 3x3, (x,y,z) = entries (12),(23),(13).
Coproduct (matrix multiplication): Dx = x1+x2, Dy = y1+y2, Dz = z1 + x1*y2 + z2.
Left action (Baumann: n.f = f(?n)): e1 = d/dx, e2 = d/dy + x d/dz   (BKK Ex. 2.3).
Unique biperfect basis (BKK Ex. 2.3): {x^a z^b w^c} u {y^a z^b w^c}, w = xy - z.
Baumann Prop 3.12 + the uniqueness in BKK Ex 2.3 => this basis is the polite basis
(dual canonical = dual semicanonical = MV here).

For a generic theta in Q_R^*, Prop 2.9 / Def 3.1(iii): for b in B, let mu = the
<='-largest element of L_theta(b) (for generic theta: the weight nu in L(b) maximising theta(nu),
weights written as nu in Q_+ with deg = -nu); the component of Delta(b) in bidegree
(-mu, wt b + mu) should be b_+ (x) b_- with b_+, b_- in B.  We
  (1) verify that (polite axiom, computed),
  (2) test whether the LEFT crystal data of b is given by a tensor-product rule in
      terms of (b_+, b_-):  eps_j(b) =? max(eps_j(b1), eps_j(b2) - <h_j, wt b1>) for
      (b1,b2) = (b_+,b_-) [Kashiwara order] or (b_-,b_+) [reversed], and whether
      Phi(e~_j b) = (e~_j b_+, b_-) or (b_+, e~_j b_-) as the signature rule predicts.
All results: grade 'computed'.
"""
import itertools
from collections import defaultdict
import sympy as sp

x, y, z = sp.symbols('x y z')
x1, y1, z1, x2, y2, z2 = sp.symbols('x1 y1 z1 x2 y2 z2')
w = x * y - z
DEG = 5   # total height bound on basis elements

def deg(mon):  # mon = (p,q,r) exponents of x,y,z ; weight -(a1 al1 + a2 al2)
    p, q, r = mon
    return (p + r, q + r)

basis = {}
for a, b, c in itertools.product(range(DEG + 1), repeat=3):
    for kind, gen in (('x', x), ('y', y)):
        if kind == 'y' and a == 0:
            continue  # avoid double counting z^b w^c
        f = sp.expand(gen ** a * z ** b * w ** c)
        h = a + 2 * b + 2 * c
        if h <= DEG:
            basis[(kind, a, b, c)] = f
key_of = {}
for k, f in basis.items():
    key_of[sp.srepr(sp.Poly(f, x, y, z).as_expr())] = k

def wt(f):
    P = sp.Poly(f, x, y, z)
    ds = {deg(m) for m in P.monoms()}
    assert len(ds) == 1
    return ds.pop()

def ident(f):
    """return basis key if f is exactly a basis element, else None"""
    if f == 0:
        return None
    return key_of.get(sp.srepr(sp.Poly(sp.expand(f), x, y, z).as_expr()))

def e(j, f):
    if j == 1:
        return sp.expand(sp.diff(f, x))
    return sp.expand(sp.diff(f, y) + x * sp.diff(f, z))

def eps(j, f):
    n = 0
    while True:
        f = e(j, f)
        if f == 0:
            return n
        n += 1

def leading_e(j, f):
    """e~_j f: e_j f = eps * e~_j f + (terms with smaller eps). Decompose e_j f in basis."""
    g = e(j, f)
    if g == 0:
        return None
    coeffs = expand_in_basis(g)
    n = eps(j, f)
    top = [(k, c) for k, c in coeffs.items() if eps(j, basis[k]) == n - 1]
    assert len(top) == 1 and top[0][1] == n, (f, coeffs)
    return top[0][0]

def expand_in_basis(g):
    g = sp.expand(g)
    W = wt(g)
    cand = [k for k, f in basis.items() if wt(f) == W]
    cs = sp.symbols(f'c0:{len(cand)}')
    eq = sp.Poly(g - sum(c * basis[k] for c, k in zip(cs, cand)), x, y, z).coeffs()
    sol = sp.solve(eq, cs, dict=True)
    assert len(sol) == 1
    return {k: sol[0].get(c, 0) for c, k in zip(cs, cand) if sol[0].get(c, 0) != 0}

def coproduct(f):
    return sp.expand(f.subs({x: x1 + x2, y: y1 + y2, z: z1 + x1 * y2 + z2}, simultaneous=True))

def components(f):
    D = sp.Poly(coproduct(f), x1, y1, z1, x2, y2, z2)
    comp = defaultdict(lambda: sp.Integer(0))
    for mon, c in D.terms():
        left = deg(mon[:3])
        comp[left] += c * sp.prod([v ** e_ for v, e_ in zip((x1, y1, z1, x2, y2, z2), mon)])
    return comp

def hpair(j, W):  # <h_j, wt> with wt = -(a1 al1 + a2 al2)
    a1, a2 = W
    return -(2 * a1 - a2) if j == 1 else -(2 * a2 - a1)

def split_pure_tensor(C):
    """C(x1..,x2..) -> (b_+, b_-) basis keys if C = b_+(x1,y1,z1) b_-(x2,y2,z2) exactly."""
    C = sp.factor(C)
    left = sp.Integer(1); right = sp.Integer(1); const = sp.Integer(1)
    for fac, m in sp.factor_list(C)[1]:
        fs = fac.free_symbols
        if fs <= {x1, y1, z1}:
            left *= fac ** m
        elif fs <= {x2, y2, z2}:
            right *= fac ** m
        else:
            return None
    const = sp.factor_list(C)[0]
    L = sp.expand(const * left.subs({x1: x, y1: y, z1: z}))
    R = sp.expand(right.subs({x2: x, y2: y, z2: z}))
    kL, kR = ident(L), ident(R)
    if kL is None or kR is None:
        # try sign flip
        kL2, kR2 = ident(-L), ident(-R)
        if kL2 is not None and kR2 is not None:
            return ('NEG', kL2, kR2)
        return None
    return (kL, kR)

thetas = {  # one theta per Weyl chamber (values on alpha1, alpha2), irrational-ish generic
    'C1 (+,+)': (1.0, 0.37), 'C2 (+,+)': (0.41, 1.0),
    'C3 (+,-) t1+t2>0': (1.0, -0.43), 'C4 (+,-) t1+t2<0': (0.39, -1.0),
    'C5 (-,+) t1+t2>0': (-0.45, 1.0), 'C6 (-,+) t1+t2<0': (-1.0, 0.35),
    'C7 (-,-)': (-1.0, -0.33), 'C8 (-,-)': (-0.31, -1.0),
}

def run():
    keys = sorted(basis, key=lambda k: (k[1] + 2 * k[2] + 2 * k[3], k))
    epsd = {(j, k): eps(j, basis[k]) for j in (1, 2) for k in keys}
    et = {(j, k): (leading_e(j, basis[k]) if epsd[(j, k)] > 0 else None) for j in (1, 2) for k in keys}
    print(f"|B| (height <= {DEG}) = {len(keys)}")
    for name, (t1, t2) in thetas.items():
        Phi = {}
        polite_ok = True
        for k in keys:
            f = basis[k]
            if f == 1:
                Phi[k] = (k, k); continue
            comp = components(f)
            mu = max(comp, key=lambda nu: t1 * nu[0] + t2 * nu[1])
            s = split_pure_tensor(comp[mu])
            if s is None or s[0] == 'NEG':
                polite_ok = False; print("  polite FAIL", name, k, comp[mu], s); continue
            Phi[k] = s
        stats = defaultdict(int)
        bad = []
        for k in keys:
            if k not in Phi: continue
            bp, bm = Phi[k]
            for j in (1, 2):
                E = epsd[(j, k)]
                ep, em = epsd.get((j, bp)), epsd.get((j, bm))
                if ep is None or em is None: continue
                Wp, Wm = wt(basis[bp]), wt(basis[bm])
                K = max(ep, em - hpair(j, Wp))       # Kashiwara order b_+ (x) b_-
                Rv = max(em, ep - hpair(j, Wm))      # reversed order b_- (x) b_+
                stats[('K', j)] += (E == K); stats[('R', j)] += (E == Rv)
                stats[('add', j)] += (E == ep + em); stats[('tot', j)] += 1
                # e~ compatibility
                if E > 0 and et[(j, k)] in Phi:
                    tgt = Phi[et[(j, k)]]
                    on_plus = (et.get((j, bp)) is not None and tgt == (et[(j, bp)], bm))
                    on_minus = (et.get((j, bm)) is not None and tgt == (bp, et[(j, bm)]))
                    stats[('e~ on +', j)] += on_plus; stats[('e~ on -', j)] += on_minus
                    stats[('e~ local', j)] += (on_plus or on_minus); stats[('e~ tot', j)] += 1
                    if not (on_plus or on_minus) and len(bad) < 4:
                        bad.append((j, k, (bp, bm), et[(j, k)], tgt))
        print(f"\n{name}: theta=({t1},{t2})  polite-factorisation pure tensors in B: {polite_ok}")
        for j in (1, 2):
            print(f"  j={j}: eps_j(b) == Kashiwara-rule(b+ x b-): {stats[('K',j)]}/{stats[('tot',j)]};"
                  f" reversed(b- x b+): {stats[('R',j)]}/{stats[('tot',j)]}; additive: {stats[('add',j)]}/{stats[('tot',j)]};"
                  f"  e~_j local (acts on one factor): {stats[('e~ local',j)]}/{stats[('e~ tot',j)]}"
                  f" (on b+: {stats[('e~ on +',j)]}, on b-: {stats[('e~ on -',j)]})")
        for bd in bad:
            print("   non-local example: j=%d b=%s Phi(b)=%s e~b=%s Phi(e~b)=%s" % bd)

if __name__ == '__main__':
    run()


# ------------------------------------------------------------------ addendum (same session)
# (a) right (*) crystal: f.n = f(n ?) ; right e1 = d/dx + y d/dz, right e2 = d/dy.
# (b) compare left-locality per (theta, j) with the criterion
#     LOC(theta,j):  ad(e_j) preserves n_{<=0}(theta), i.e. no positive root a with
#                    theta(a) <= 0 < theta(a + alpha_j).
def er(j, f):
    if j == 1:
        return sp.expand(sp.diff(f, x) + y * sp.diff(f, z))
    return sp.expand(sp.diff(f, y))

def eps_r(j, f):
    n = 0
    while True:
        f = er(j, f)
        if f == 0:
            return n
        n += 1

def leading_er(j, f):
    g = er(j, f)
    if g == 0:
        return None
    coeffs = expand_in_basis(g)
    n = eps_r(j, f)
    top = [(k, c) for k, c in coeffs.items() if eps_r(j, basis[k]) == n - 1]
    assert len(top) == 1 and top[0][1] == n, (f, coeffs)
    return top[0][0]

def addendum():
    keys = sorted(basis, key=lambda k: (k[1] + 2 * k[2] + 2 * k[3], k))
    roots = [(1, 0), (0, 1), (1, 1)]
    for name, (t1, t2) in thetas.items():
        th = lambda r: t1 * r[0] + t2 * r[1]
        Phi = {}
        for k in keys:
            f = basis[k]
            if f == 1:
                Phi[k] = (k, k); continue
            comp = components(f)
            mu = max(comp, key=lambda nu: th(nu))
            Phi[k] = split_pure_tensor(comp[mu])
        line = []
        for side, E, EPS in (('left', leading_e, eps), ('right*', leading_er, eps_r)):
            for j in (1, 2):
                loc = tot = 0
                for k in keys:
                    if EPS(j, basis[k]) == 0: continue
                    t = E(j, basis[k])
                    if t not in Phi: continue
                    bp, bm = Phi[k]; tp, tm = Phi[t]
                    tot += 1
                    loc += ((tp == bp) != (tm == bm))   # exactly one factor changes
                line.append(f"{side} j={j}: {loc}/{tot}")
        crit = []
        for j in (1, 2):
            aj = (1, 0) if j == 1 else (0, 1)
            bad = any(th(a) <= 0 < th((a[0] + aj[0], a[1] + aj[1]))
                      for a in roots if (a[0] + aj[0], a[1] + aj[1]) in roots)
            crit.append(f"LOC(j={j})={'yes' if not bad else 'no'}")
        pos = [r for r in roots if th(r) > 0]
        print(f"{name}: {{theta>0}}={pos}  " + "; ".join(line) + "   criterion: " + ", ".join(crit))

if __name__ == '__main__':
    print("\n=== addendum: one-factor locality of left and right(*) crystals vs criterion ===")
    addendum()

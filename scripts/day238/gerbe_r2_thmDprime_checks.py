"""Rick, 2026-10-13: independent checks of MacBeth UID 377
  'A classical H^2 obstruction, read as supply-chain inventory: the inventory gerbe and its band'
  (corrected 11 Oct 2026, content e4d1785), round 2.
Written from the PRINTED statements only (Thm 28 = 'Thm D-prime', Def 6 (D2), Cor 18, Thm 22(ii), Rem 29).
Own dihedral engine: element (i,e) = r^i s^e in D_n, |D_n| = 2n.  Nothing of MacBeth's code is used.

Parts
 A  Thm 28 chain on n = 3..NMAX, each link computed by an independent method, vs 4|n; controls.
 B  'Equivalently, ... always rectifiable unless 4|n' (Thm 28, Rem 29, abstract): band witnesses.
 C  N1 (Cor 18 under D2): brute force D2 <=> (beta=0 and [g]=0) on finite bases, all strict alpha.
 D  The k_{g,f} datum of Sec 2.3 / Rem 11: D2 vs the k-twisted neutrality equation.
"""
import itertools, sys
from math import gcd

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 64
def say(*s): print(" ".join(str(x) for x in s), flush=True)

class Dn:
    def __init__(self, n):
        self.n = n
        self.el = [(i, e) for e in (0, 1) for i in range(n)]
        self.e = (0, 0)
        self.Z = [x for x in self.el if all(self.mul(x, y) == self.mul(y, x) for y in self.el)]
    def mul(self, x, y):
        (i, e), (j, f) = x, y
        return ((i + (j if e == 0 else -j)) % self.n, e ^ f)
    def inv(self, x):
        i, e = x
        return ((-i) % self.n, 0) if e == 0 else x
    def order(self, x):
        k, y = 1, x
        while y != self.e: y = self.mul(y, x); k += 1
        return k
    # automorphisms r -> r^a, s -> s r^b  (all of Aut D_n for n >= 3)
    def aut(self, a, b):
        n = self.n
        return tuple(((a * i - b * e) % n, e) for (i, e) in self.el)  # image tuple, indexed like self.el
    def auts(self):
        return [self.aut(a, b) for a in range(1, self.n) if gcd(a, self.n) == 1 for b in range(self.n)]
    def idx(self, x): return x[1] * self.n + x[0]
    def app(self, th, x): return th[self.idx(x)]
    def comp(self, th, ph): return tuple(self.app(th, self.app(ph, x)) for x in self.el)   # th o ph
    def ainv(self, th):
        out = [None] * len(self.el)
        for x in self.el: out[self.idx(self.app(th, x))] = x
        return tuple(out)
    def conj(self, k):
        ki = self.inv(k)
        return tuple(self.mul(self.mul(k, x), ki) for x in self.el)

def gen_closure(G, gens):
    S = {G.e}; frontier = [G.e]
    while frontier:
        nxt = []
        for x in frontier:
            for g in gens:
                y = G.mul(x, g)
                if y not in S: S.add(y); nxt.append(y)
        frontier = nxt
    return S

def rank_f2(rows):
    piv = {}; rk = 0
    for r in rows:
        for b in sorted(piv, reverse=True):
            if r >> b & 1: r ^= piv[b]
        if r:
            piv[r.bit_length() - 1] = r; rk += 1
    return rk

# ---------------------------------------------------------------------------------------------
def star_class_nonzero(G):
    """[*] != 0 in H^2(Inn K; Z), Z = Z/2: section cocycle on Inn K = K/Z, coboundary test over F_2.
    Inn K realised as the set of distinct conjugation automorphisms (no structure theorem used)."""
    if len(G.Z) == 1: return False
    z = [x for x in G.Z if x != G.e]; assert len(z) == 1
    zc = z[0]
    inn = {}
    for k in G.el: inn.setdefault(G.conj(k), k)      # section s: Inn -> K, first preimage
    I = list(inn); idx = {a: t for t, a in enumerate(I)}; m = len(I)
    rows = []
    for a in I:
        for b in I:
            ab = G.comp(a, b)
            gam = G.mul(G.mul(inn[a], inn[b]), G.inv(inn[ab]))
            assert gam in G.Z
            row = 0
            for h in (a, b, ab): row ^= 1 << idx[h]
            rows.append((row, 1 if gam == zc else 0))
    A = [r for r, _ in rows]; Ab = [r | (bit << m) for r, bit in rows]
    return rank_f2(A) != rank_f2(Ab)

def nonsplit_bruteforce(G):
    """(*) non-split: no subgroup H with H n Z = 1 and |H| = |K|/|Z|. H = <x,y> with x,y lifting generators
    of Inn K = <c_r, c_s> (any complement contains such lifts)."""
    if len(G.Z) == 1: return False
    for x in (y for y in G.el if G.conj(y) == G.conj((1, 0))):
        for y in (w for w in G.el if G.conj(w) == G.conj((0, 1))):
            H = gen_closure(G, [x, y])
            if len(H) * len(G.Z) == len(G.el) and all(z == G.e or z not in H for z in G.Z): return False
    return True

def witness_BZ2(G):
    """[g] != 0 for some inner transport on C = B(Z/2): an involution (or e) of Inn K with no lift of order <= 2."""
    for a in {G.conj(k) for k in G.el}:
        if G.comp(a, a) != tuple(G.el): continue
        if not any(G.conj(k) == a and G.mul(k, k) == G.e for k in G.el): return True
    return False

def witness_BV4(G):
    """[g] != 0 for some rho: V4 -> Inn K: commuting involution pair (a,b) in Inn K with no commuting pair of
    involution lifts."""
    innset = {G.conj(k): None for k in G.el}
    invs = [a for a in innset if G.comp(a, a) == tuple(G.el)]
    lifts = {a: [k for k in G.el if G.conj(k) == a and G.mul(k, k) == G.e] for a in invs}
    for a in invs:
        for b in invs:
            if G.comp(a, b) != G.comp(b, a): continue
            if not any(G.mul(x, y) == G.mul(y, x) for x in lifts[a] for y in lifts[b]): return True
    return False

say("=" * 92)
say("PART A  Thm 28 chain for the n-lot ring K = D_n, n = 3..%d (each link by a separate method)" % NMAX)
say("=" * 92)
rows = []
for n in range(3, NMAX + 1):
    G = Dn(n)
    s1 = star_class_nonzero(G)
    s2 = nonsplit_bruteforce(G) if n <= 48 else None
    w2 = witness_BZ2(G); w4 = witness_BV4(G) if n <= 32 else None
    rows.append((n, len(G.Z), s1, s2, w2, w4))
    say("n=%2d |Z|=%d  [*]!=0:%-5s  nonsplit(brute):%-5s  B(Z/2)-witness:%-5s  B(V4)-witness:%-5s  4|n:%s"
        % (n, len(G.Z), s1, s2, w2, w4, n % 4 == 0))
def tally(pred, name):
    bad = [r[0] for r in rows if (r[2] != pred(r[0]))]
    say("  %-58s mismatches vs [*]!=0: %3d  %s" % (name, len(bad), bad[:12]))
    return bad
say("-- chain links agree with each other:")
say("  [*]!=0 == nonsplit(brute) on n<=48 :", all(r[2] == r[3] for r in rows if r[3] is not None))
say("  [*]!=0 == B(Z/2)-witness (=> 'essentially nonzero'):", all(r[2] == r[4] for r in rows))
say("  [*]!=0 == B(V4)-witness on n<=32   :", all(r[2] == r[5] for r in rows if r[5] is not None))
say("-- the printed arithmetic, and controls that must FAIL:")
tally(lambda n: n % 4 == 0, "Thm 28 printed: 4 | n")
tally(lambda n: n % 2 == 0, "control 1: 2 | n  (= retracted Thm D 'nonabelian, Z!=1')")
tally(lambda n: n % 8 == 0, "control 2: 8 | n")
tally(lambda n: n % 4 == 2, "control 3: n = 2 mod 4 (swapped)")

# ---------------------------------------------------------------------------------------------
say("=" * 92)
say("PART B  'always rectifiable unless 4|n' (Thm 28 'Equivalently', Rem 29, abstract): band witnesses")
say("=" * 92)
say("D2 on C=B(Z) with alpha(t)=theta: needs mu, k with mu theta mu^-1 = c_k, i.e. theta inner.")
say("Finite base B(Z/m), m = ord(theta): same test, T(t)=k must also satisfy k^m=e -- checked below.")
bad_n = []
for n in range(3, 25):
    G = Dn(n); A = G.auts(); inn = {G.conj(k) for k in G.el}
    out = [th for th in A if th not in inn]
    ident = tuple(G.el)
    # finite-base witness: outer theta of order 2 on B(Z/2)
    inv_out = [th for th in out if G.comp(th, th) == ident]
    finite_ok = False
    if inv_out:
        th = inv_out[0]
        # D2 on B(Z/2): exists mu in Aut, k in K, k^2=e, c_k = mu th mu^-1
        d2 = any(G.conj(k) == G.comp(G.comp(mu, th), G.ainv(mu)) and G.mul(k, k) == G.e for mu in A for k in G.el)
        finite_ok = not d2
    say("n=%2d |Aut|=%3d |Inn|=%2d |Out|=%2d  outer involution on B(Z/2) obstructs D2: %-5s  4|n:%s"
        % (n, len(A), len(inn), len(A) // len(inn), finite_ok, n % 4 == 0))
    if finite_ok and n % 4: bad_n.append(n)
say("  n with 4 \u2224 n where a strict transport on the FINITE base B(Z/2) admits NO global inventory (band):", bad_n)
bz = [n for n in range(3, 25) if n % 4 and len(Dn(n).auts()) > len({Dn(n).conj(k) for k in Dn(n).el})]
say("  n<=24 with 4 \u2224 n and Out D_n != 1 (C=B(Z), alpha(t)=outer theta: NO global inventory):", bz)
G = Dn(6); th = G.aut(5, 1)
say("  explicit n=6: theta = (r->r^5, s->s r) ; theta^2 = id:", G.comp(th, th) == tuple(G.el),
    "; inner:", th in {G.conj(k) for k in G.el})

# ---------------------------------------------------------------------------------------------
say("=" * 92)
say("PART C  N1: Cor 18 under D2 (Def 6), strict alpha: Gamma -> Aut D_n, bases Z/2, Z/4, V4; n=3..12")
say("=" * 92)
def cyc(m): return ([(i,) for i in range(m)], lambda x, y: ((x[0] + y[0]) % m,), (0,), [(1,)])
def v4(): return ([(a, b) for a in (0, 1) for b in (0, 1)], lambda x, y: (x[0] ^ y[0], x[1] ^ y[1]), (0, 0), [(1, 0), (0, 1)])
def homs_to(gam, tgt, tmul, tid):
    els, gmul, ge, gens = gam
    words = {ge: []}; frontier = [ge]
    while frontier:
        nxt = []
        for x in frontier:
            for t, g in enumerate(gens):
                y = gmul(x, g)
                if y not in words: words[y] = words[x] + [t]; nxt.append(y)
        frontier = nxt
    res = []
    for imgs in itertools.product(tgt, repeat=len(gens)):
        phi = {}
        for x, w in words.items():
            y = tid
            for t in w: y = tmul(y, imgs[t])
            phi[x] = y
        if all(phi[gmul(x, g)] == tmul(phi[x], imgs[t]) for x in els for t, g in enumerate(gens)):
            res.append(phi)
    return res
def gclass_zero(gam, G, rho_aut):
    """[g] for inner-valued strict alpha (= conj o rho): section cocycle on the bar nerve of Gamma, F_2 test."""
    els, gmul, ge, _ = gam
    if len(G.Z) == 1: return True
    zc = [x for x in G.Z if x != G.e][0]
    sec = {}
    for k in G.el: sec.setdefault(G.conj(k), k)
    T = {g: sec[rho_aut[g]] for g in els}; T[ge] = G.e
    idx = {g: t for t, g in enumerate(els)}; m = len(els); rr = []
    for g in els:
        for f in els:
            gam_ = G.mul(G.mul(T[g], T[f]), G.inv(T[gmul(g, f)]))
            row = 0
            for h in (g, f, gmul(g, f)): row ^= 1 << idx[h]
            rr.append((row, 1 if gam_ == zc else 0))
    return rank_f2([r for r, _ in rr]) == rank_f2([r | (b << m) for r, b in rr])
tot = dict(cases=0, cor18_fail=0, band_nontriv=0, D1_fail_D2_ok=0)
for n in range(3, 13):
    G = Dn(n); A = G.auts(); inn = {G.conj(k) for k in G.el}; ident = tuple(G.el)
    for gname, gam in (("Z/2", cyc(2)), ("Z/4", cyc(4)), ("V4", v4())):
        Khoms = homs_to(gam, G.el, G.mul, G.e)
        for al in homs_to(gam, A, G.comp, ident):
            beta0 = all(al[g] in inn for g in gam[0])            # one object: beta=0 iff q o alpha trivial
            g0 = gclass_zero(gam, G, al) if beta0 else None
            d2 = False
            for mu in A:
                mui = G.ainv(mu)
                tau = {g: G.comp(G.comp(mu, al[g]), mui) for g in gam[0]}
                if any(all(G.conj(T[g]) == tau[g] for g in gam[0]) for T in Khoms): d2 = True; break
            d1 = any(all(G.comp(G.comp(mu, al[g]), G.ainv(mu)) == ident for g in gam[0]) for mu in A)
            tot["cases"] += 1
            tot["cor18_fail"] += (d2 != (beta0 and g0))
            tot["band_nontriv"] += (not beta0)
            tot["D1_fail_D2_ok"] += (d2 and not d1)
say("  strict alpha cases:", tot["cases"], "| Cor 18 (D2 <=> beta=0 & [g]=0) failures:", tot["cor18_fail"],
    "| nontrivial-band cases:", tot["band_nontriv"], "| D2-but-not-D1 (old N1 gap, now harmless):", tot["D1_fail_D2_ok"])
G = Dn(4); al = G.conj((0, 1))
say("  N1 witness B(Z/2), D4, alpha(g)=c_s: alpha strict:", G.comp(al, al) == tuple(G.el),
    "| D2 via mu=id, T(g)=s:", G.mul((0, 1), (0, 1)) == G.e,
    "| D1 (some mu kills c_s):", any(G.comp(G.comp(mu, al), G.ainv(mu)) == tuple(G.el) for mu in G.auts()))

# ---------------------------------------------------------------------------------------------
say("=" * 92)
say("PART D  The k_{g,f} datum (Sec 2.3 'strict functor carrying the k_{g,f}'; Rem 11 '[g] is its k-invariant')")
say("=" * 92)
say("On C=B(Z/2): alpha strict (alpha_g^2=id) forces k_{g,g} in Z. Cocycle law k_{h,gf} alpha_h(k_{g,f}) = k_{hg,f} k_{h,g}")
say("holds for any k_{g,g}=z in Z (Aut fixes the central involution). Neutrality of the AUT(F)-pseudofunctor (alpha,k):")
say("  exists mu in Aut K, T_g in K with c_{T_g} = mu alpha_g mu^-1 and T_{gg}=e = mu(k_{g,g}) T_g T_g.")
say("D2 (Def 6) is the same with k ignored (T_g^2 = e).")
def neutral(G, al, z, twisted):
    for mu in G.auts():
        tau = G.comp(G.comp(mu, al), G.ainv(mu))
        for T in G.el:
            if G.conj(T) != tau: continue
            need = G.e if not twisted else G.inv(G.app(mu, z))
            if G.mul(T, T) == need: return True
    return False
for n in (4, 6, 8, 10, 12):
    G = Dn(n); zc = [x for x in G.Z if x != G.e][0]; ident = tuple(G.el)
    for name, al in (("alpha=id", ident), ("alpha=c_r", G.conj((1, 0))), ("alpha=c_{r^(n/4)}" if n % 4 == 0 else "alpha=c_s", G.conj((n // 4, 0)) if n % 4 == 0 else G.conj((0, 1)))):
        if G.comp(al, al) != ident: continue
        say("  n=%2d %-18s  D2:%-5s | true neutrality k=e:%-5s k=z:%-5s" %
            (n, name, neutral(G, al, G.e, False), neutral(G, al, G.e, True), neutral(G, al, zc, True)))
say("  => with k=z central the same strict alpha flips: Thm 22(ii)'s witness (n=4, c_r) becomes neutral, and alpha=id")
say("     becomes non-neutral for EVERY even n (incl. n=6). Under the pseudofunctor reading the criterion is 2|n, not 4|n.")

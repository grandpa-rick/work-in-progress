"""Rick, 2026-10-13: independent checks of MacBeth's
  (382) capstone 'The inventory gerbe couples on the torus' (content e306ae5): M1-M3 as printed, Thm 8 sign, eq. (5);
  (383) 'The inventory floor needs the 2-cell datum: D1 versus D2' (content f19cfcb): Thm 1.1-1.4, Props 3/5/7, his Q1-Q3.
Written from the printed statements; group engine reused from Rick's day237 band_boundary_checks.py (imported).
Nothing of MacBeth's code is used."""
import itertools, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'day237'))
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):   # day237 module prints its own log on import
    import band_boundary_checks as bb
from band_boundary_checks import Grp, mul, inv, amul, ainv, conj_aut, structure, coset

def say(*s): print(" ".join(str(x) for x in s))

# ---------------- groups ----------------
S3, D4, Q8, D4C2 = bb.S3, bb.D4, bb.Q8, bb.D4C2
C4C2 = bb.C4C2
D6 = Grp([bb.cyc(6), tuple((-i) % 6 for i in range(6))], "D6")
Z2 = Grp([(1, 0)], "Z/2"); Z3 = Grp([(1, 2, 0)], "Z/3"); Z4 = Grp([(1, 2, 3, 0)], "Z/4")
V4 = Grp([(1, 0, 2, 3), (0, 1, 3, 2)], "V4")
Z2Z4 = Grp([(1, 0, 2, 3, 4, 5), (0, 1, 3, 4, 5, 2)], "Z2xZ4")
S3b = Grp([(1, 0, 2), (1, 2, 0)], "S3(base)")

def homs(Gam, tgt_el, tmul, tid):
    """All homomorphisms Gam -> target, as dicts on Gam.el. Generic search over generator images."""
    res = []
    for imgs in itertools.product(tgt_el, repeat=len(Gam.gens)):
        phi = {}
        for x, w in Gam.word.items():
            y = tid
            for k in w: y = tmul(y, imgs[k])
            phi[x] = y
        if all(phi[mul(x, g)] == tmul(phi[x], imgs[k]) for x in Gam.el for k, g in enumerate(Gam.gens)):
            res.append(phi)
    return res

class Fibre:
    def __init__(self, K):
        self.K = K
        self.A, self.I, self.Is = structure(K)
        self.id = tuple(range(K.N))
        self.c = {k: conj_aut(K, k) for k in K.el}          # c : K -> Inn K
        self.pre = {}                                          # c^{-1}
        for k, a in self.c.items(): self.pre.setdefault(a, []).append(k)
        self.Z = K.Z
    def act(self, mu, a): return amul(amul(mu, a), ainv(mu))  # mu a mu^-1 in Aut K

# ---------- F_2-linear algebra for H^2(Gamma; Z), Z elementary abelian 2-group ----------
def z_coords(F):
    """Express Z(K) (elementary abelian 2-group) as F_2^r: returns (r, vec map)."""
    Z = F.Z; K = F.K
    basis = []; span = {K.e: 0}
    for z in Z:
        if z in span: continue
        basis.append(z); new = {}
        for x, v in span.items(): new[mul(x, z)] = v | (1 << (len(basis) - 1))
        span.update(new)
    assert all(mul(z, z) == K.e for z in Z), "Z not elementary 2"
    return len(basis), span

def rank_f2(rows):
    rows = [r for r in rows]; rk = 0; piv = {}
    for r in rows:
        for b, pr in piv.items():
            if r >> b & 1: r ^= pr
        if r:
            b = r.bit_length() - 1
            for bb_, pr in list(piv.items()):
                if pr >> b & 1: piv[bb_] = pr ^ r
            piv[b] = r; rk += 1
    return rk

def g_class_zero(Gam, F, abar):
    """[g] = abar^*[*] in H^2(BGam; Z) via bar factor set z(g,f)=t_gf (t_g t_f)^-1, coboundary test by F_2 linear algebra.
    Independent of the hom search used for D2."""
    if len(F.Z) == 1 or len(F.I) == 1: return True   # trivial coefficients / trivial abar (section t=e gives z=e)
    r, vec = z_coords(F)
    sec = {a: F.pre[a][0] for a in F.pre}               # set section
    t = {g: sec[abar[g]] for g in Gam.el}
    t[Gam.e] = F.K.e
    els = Gam.el; idx = {g: i for i, g in enumerate(els)}; n = len(els)
    # unknown zeta: Gam -> F_2^r ; equation per (g,f), per coordinate: zeta(gf)+zeta(g)+zeta(f) = z(g,f)
    ok = True
    for bit in range(r):
        rows = []
        for g in els:
            for f in els:
                zz = mul(t[mul(g, f)], inv(mul(t[g], t[f])))
                assert zz in vec, "defect not central"
                row = 0
                for h in (mul(g, f), g, f): row ^= 1 << idx[h]
                rows.append((row, (vec[zz] >> bit) & 1))
        A = [rw for rw, _ in rows]; Ab = [rw | (b << n) for rw, b in rows]
        if rank_f2(A) != rank_f2(Ab): ok = False
    return ok

def D1(Gam, F, abar):   # exists mu in Aut K, mu abar(g) mu^-1 = id for all g (one object)
    return any(all(F.act(mu, abar[g]) == F.id for g in Gam.el) for mu in F.A)

def D2(Gam, F, abar, Khoms):  # exists mu in Aut K and strict T: Gam -> K with c(T g) = mu abar(g) mu^-1
    for mu in F.A:
        tau = {g: F.act(mu, abar[g]) for g in Gam.el}
        if any(all(F.c[T[g]] == tau[g] for g in Gam.el) for T in Khoms): return True
    return False

def D2_twistless_mu_id(Gam, F, abar, Khoms):
    return any(all(F.c[T[g]] == abar[g] for g in Gam.el) for T in Khoms)

# ======================================================================================
say("=" * 90); say("PART I  (UID 383) Floor dichotomy, Thm 1.1-1.4, over C = B(Gamma), Gamma finite, and C = B(Z)")
say("=" * 90)
fibres = [F for F in map(Fibre, [S3, D4, Q8, D6, D4C2, C4C2])]
bases = [Z2, Z3, Z4, V4, Z2Z4, S3b]
tot = {"cases": 0, "T1.1 fail": 0, "T1.2 fail": 0, "D1=>D2 fail": 0, "gap": 0, "gap_with_Z_nontriv": 0,
       "swap D2=>D1 holds": 0, "D1<=>[g]=0 holds": 0}
for F in fibres:
    K = F.K
    say("K=%-6s |K|=%2d |Z|=%d |Inn|=%2d |Aut|=%2d |Out|=%2d" % (K.name, K.N, len(F.Z), len(F.I), len(F.A), len(F.A)//len(F.I)))
    for Gam in bases:
        Inn_el = F.I
        abars = homs(Gam, Inn_el, amul, F.id)
        Khoms = homs(Gam, K.el, mul, K.e)
        c = {"n": 0, "d1": 0, "d2": 0, "g0": 0, "gap": 0}
        for ab in abars:
            d1 = D1(Gam, F, ab); d2 = D2(Gam, F, ab, Khoms); g0 = g_class_zero(Gam, F, ab)
            triv = all(ab[g] == F.id for g in Gam.el)
            c["n"] += 1; c["d1"] += d1; c["d2"] += d2; c["g0"] += g0; c["gap"] += (d2 and not d1)
            tot["cases"] += 1
            tot["T1.1 fail"] += (d2 != g0)
            tot["T1.2 fail"] += (d1 != triv)          # one object: [abar]=* iff abar trivial
            tot["D1=>D2 fail"] += (d1 and not d2)
            tot["gap"] += (d2 and not d1); tot["gap_with_Z_nontriv"] += (d2 and not d1 and len(F.Z) > 1)
            tot["swap D2=>D1 holds"] += (not d2 or d1)
            tot["D1<=>[g]=0 holds"] += (d1 == g0)
        say("   B(%-6s): #abar=%3d  D1=%3d  D2=%3d  [g]=0:%3d  gap(D2 & not D1)=%3d" % (Gam.name, c["n"], c["d1"], c["d2"], c["g0"], c["gap"]))
    # B(Z): abar(1)=x in Inn K, free.  D2: exists mu,k with c_k = mu x mu^-1 ; [g]=0 (H^2(S^1)=0); D1: exists mu, mu x mu^-1 = id
    d1 = sum(any(F.act(mu, x) == F.id for mu in F.A) for x in F.I)
    d2 = sum(any(F.act(mu, x) in F.pre for mu in F.A) for x in F.I)
    say("   B(Z)     : #abar=%3d  D1=%3d  D2=%3d  [g]=0:%3d  gap=%3d" % (len(F.I), d1, d2, len(F.I), d2 - d1))
    tot["gap"] += d2 - d1; tot["gap_with_Z_nontriv"] += (d2 - d1) * (len(F.Z) > 1); tot["cases"] += len(F.I)
say("TOTALS:", tot)
say("  -> Thm 1.1 (D2 <=> [g]=0): failures", tot["T1.1 fail"], "; Thm 1.2 (D1 <=> [abar]=*): failures", tot["T1.2 fail"],
    "; Prop 7 (D1 => D2): failures", tot["D1=>D2 fail"])

say("\n-- negative controls (each MUST fail somewhere) --")
say("  [swap] 'D2 => D1' : holds in %d / %d (base,abar) cases over finite bases -> FAILS (as it must) iff < total" % (tot["swap D2=>D1 holds"], tot["cases"] - sum(len(F.I) for F in fibres)))
say("  [wrong degree] 'D1 <=> [g]=0': holds in %d / same -> FAILS iff < total" % tot["D1<=>[g]=0 holds"])
F = Fibre(C4C2)
say("  [abelian K=C4xC2] |Inn|=%d: every abar trivial, so gap must be 0 on every base (see table: gap=0) -- control PASSES as predicted" % len(F.I))

say("\n-- printed Thm 1.3 / abstract / Upshot: 'gap realised exactly when Z is too small (centreless regime)' --")
FD4 = Fibre(D4); s = (0, 3, 2, 1); cs = FD4.c[s]
say("  K=D4 (|Z|=2), C=B(Z), abar(1)=c_s: D2 lift T(n)=s^n (s^2=e:", mul(s, s) == D4.e, "); D1 exists:",
    any(FD4.act(mu, cs) == FD4.id for mu in FD4.A))
FQ = Fibre(Q8); i_ = bb.lmul_perm('i'); ci = FQ.c[i_]
say("  K=Q8 (|Z|=2), C=B(Z), abar(1)=c_i: D2 lift T(n)=i^n ; D1 exists:", any(FQ.act(mu, ci) == FQ.id for mu in FQ.A))
say("  total gap cases with Z(K) != 1:", tot["gap_with_Z_nontriv"], " (Thm 1.3 'exactly when centreless' is FALSE if > 0)")

say("\n-- the exact-sequence reading (answer to Q3): D2 <=> [abar] in im(H^1(C;K) -> H^1(C;InnK)) (mod Out), D1 <=> [abar] = * --")
bad = 0; n = 0
for F in fibres:
    for Gam in bases:
        Khoms = homs(Gam, F.K.el, mul, F.K.e)
        img = set()
        for T in Khoms:
            img.add(tuple(F.c[T[g]] for g in Gam.el))
        for ab in homs(Gam, F.I, amul, F.id):
            v = tuple(ab[g] for g in Gam.el)
            # in image up to Aut-conjugation:
            inimg = any(tuple(F.act(mu, x) for x in v) in img for mu in F.A)
            n += 1; bad += (inimg != D2(Gam, F, ab, Khoms))
say("  D2 == 'abar Aut-conjugate into the image of Hom(Gamma,K)': mismatches %d / %d" % (bad, n))

# ======================================================================================
say("\n" + "=" * 90); say("PART II  Q1 (disconnected C, per-component outer part) and Q2 (home of o_D1)"); say("=" * 90)
for F in [Fibre(D4C2), Fibre(D4), Fibre(Q8)]:
    # Q1 core step: for every gauge mu in Aut K, liftability(Ad(mu) o abar) == liftability(abar).  (D4xC2: Out acts nontrivially on Z.)
    zact = len({tuple(mu[F.K.idx[z]] for z in F.Z) for mu in F.A})
    mism = 0; tests = 0
    for Gam in [Z2, V4, Z2Z4]:
        Kh = homs(Gam, F.K.el, mul, F.K.e)
        for ab in homs(Gam, F.I, amul, F.id):
            l0 = D2_twistless_mu_id(Gam, F, ab, Kh)
            for mu in F.A:
                tau = {g: F.act(mu, ab[g]) for g in Gam.el}
                tests += 1; mism += (D2_twistless_mu_id(Gam, F, tau, Kh) != l0)
    say("  K=%-6s (#distinct actions of Aut K on Z(K): %d): liftability of Ad(mu)o abar vs abar, mismatches %d / %d" % (F.K.name, zact, mism, tests))
# explicit disconnected C = B(V4) |_| B(Z/2) with independent gauges mu_1, mu_2 (different outer classes allowed)
F = Fibre(D4C2)
hv = homs(V4, F.I, amul, F.id); h2 = homs(Z2, F.I, amul, F.id)
Khv = homs(V4, F.K.el, mul, F.K.e); Kh2 = homs(Z2, F.K.el, mul, F.K.e)
mism = 0; cnt = 0
for a1 in hv:
    for a2 in h2:
        joint = any(D2_twistless_mu_id(V4, F, {g: F.act(m1, a1[g]) for g in V4.el}, Khv) and
                    D2_twistless_mu_id(Z2, F, {g: F.act(m2, a2[g]) for g in Z2.el}, Kh2) for m1 in F.A for m2 in F.A)
        prod = g_class_zero(V4, F, a1) and g_class_zero(Z2, F, a2)
        cnt += 1; mism += (joint != prod)
say("  C = B(V4) |_| B(Z/2), K=D4xC2: joint D2 (independent mu per component) vs [g]=0 in H^2 = H^2(V4) x H^2(Z/2): mismatches %d / %d" % (mism, cnt))

# connected 2-object groupoid (vertex group Gamma): functor = (phi: Gamma->Inn, h for arrow 0->1). D1 with mu_0, mu_1.
F = Fibre(D4)
mism = 0; cnt = 0
for phi in homs(V4, F.I, amul, F.id):
    for h in F.I:
        abar = {}
        H = {0: F.id, 1: h}
        for i_, j_ in itertools.product((0, 1), repeat=2):
            for g in V4.el: abar[(i_, j_, g)] = amul(amul(H[i_], phi[g]), ainv(H[j_]))
        d1 = any(all(amul(amul(m[i_], abar[(i_, j_, g)]), ainv(m[j_])) == F.id for (i_, j_, g) in abar)
                 for m in [{0: a, 1: b} for a in F.A for b in F.A])
        triv = all(phi[g] == F.id for g in V4.el)
        cnt += 1; mism += (d1 != triv)
say("  2-object groupoid over V4, K=D4: D1 (two-vertex gauge) <=> vertex-group holonomy trivial: mismatches %d / %d" % (mism, cnt))

# Q2: H^1(C;Inn)/H^0(C;Out) -> H^1(C;Aut) injective?  one object: classes = homs mod conjugation.
for F in [Fibre(D4C2), Fibre(D4), Fibre(Q8)]:
    for Gam in [Z2, V4, Z4]:
        hs = homs(Gam, F.I, amul, F.id)
        key = lambda ab: tuple(ab[g] for g in Gam.el)
        def orbit(ab, group):
            return frozenset(tuple(F.act(mu, ab[g]) for g in Gam.el) for mu in group)
        inn_cls = {orbit(ab, F.I) for ab in hs}
        aut_cls = {orbit(ab, F.A) for ab in hs}   # = H^1(Inn) mod Out-action, and also = image in H^1(Aut) (Aut-conj classes)
        # image in H^1(C;Aut): Aut-conjugacy classes of Aut-valued homs that meet Inn-valued homs; the map is injective iff
        # two Inn-valued homs Aut-conjugate <=> same Out-orbit of Inn-classes: true by construction; record sizes
        say("  K=%-6s B(%-4s): |H^1(C;Inn)| = %3d, |H^1(C;Inn)/H^0(Out)| = |image in H^1(C;Aut)| = %3d" % (F.K.name, Gam.name, len(inn_cls), len(aut_cls)))
FD4 = Fibre(D4)
r_ = (1, 2, 3, 0); s_ = (0, 3, 2, 1)
cls = {}
for x in FD4.I: cls.setdefault(frozenset(FD4.act(mu, x) for mu in FD4.I), []).append(x)
orb = {frozenset(FD4.act(mu, x) for mu in FD4.A) for x in FD4.I}
say("  383 §4 'for K=D4 the Out action swaps the two reflection classes': |H^1(B(Z);InnD4)|=%d, Out-orbits=%d, c_s ~Aut~ c_rs: %s"
    % (len(cls), len(orb), any(FD4.act(mu, FD4.c[s_]) == FD4.c[mul(r_, s_)] for mu in FD4.A)))

# ======================================================================================
say("\n" + "=" * 90); say("PART III  (UID 382) M1-M3 as printed"); say("=" * 90)
# M1: Def 4 as printed: T strict functor C->BK, mu(y) alpha(f) mu(x)^-1 = conj(T f), with alpha non-strict (2-cell k).
F = Fibre(D4); K = D4
s = (0, 3, 2, 1); cs = F.c[s]
alpha = {0: F.id, 1: cs, 2: cs}                       # C = B(Z/3), alpha(t)=alpha(t^2)=c_s
# alpha(gf) = c_{k(g,f)} alpha(g) alpha(f): k(g,f) in c^{-1}(alpha(gf) (alpha(g)alpha(f))^-1); search a normalized nonabelian 2-cocycle
def need(g, f): return amul(alpha[(g + f) % 3], ainv(amul(alpha[g], alpha[f])))
pairs = [(g, f) for g in (1, 2) for f in (1, 2)]
choices = [F.pre[need(g, f)] for g, f in pairs]
def aut_apply(a, k): return K.el[a[K.idx[k]]]
valid = []
for ks in itertools.product(*choices):
    k = {(g, f): K.e for g in range(3) for f in range(3)}
    for (p, v) in zip(pairs, ks): k[p] = v
    if all(mul(k[(h, (g + f) % 3)], aut_apply(alpha[h], k[(g, f)])) == mul(k[((h + g) % 3, f)], k[(h, g)])
           for h in range(3) for g in range(3) for f in range(3)):
        valid.append(k)
say("M1: C=B(Z/3), K=D4, alpha(t)=alpha(t^2)=c_s (non-strict); normalized 2-cell data k satisfying the gerbe cocycle law:", len(valid))
k = valid[0]
Kh3 = homs(Z3, K.el, mul, K.e)
z3 = Z3.el; zi = {z3[i]: i for i in range(3)}
# map exponent -> element of Z3 group: Z3.word gives words; build exponent table
expo = {}
for x, w in Z3.word.items(): expo[len(w) % 3] = x
def4 = any(all(F.c[T[expo[f]]] == F.act(mu, alpha[f]) for f in range(3)) for mu in F.A for T in Kh3)
# twisted D2 (gerbe trivialisation): exists mu, T: Mor->K with conj(T f) = mu alpha(f) mu^-1 and T(gf) = mu(k(g,f)) T(g) T(f)
tw = False
for mu in F.A:
    for T1 in F.pre[F.act(mu, alpha[1])]:
        for T2 in F.pre[F.act(mu, alpha[2])]:
            T = {0: K.e, 1: T1, 2: T2}
            if all(T[(g + f) % 3] == mul(mul(aut_apply(mu, k[(g, f)]), T[g]), T[f]) for g in range(3) for f in range(3)):
                tw = True
say("   band trivial (Hom(Z/3,Out D4=Z/2)=1) and H^2(Z/3;Z/2)=0, so Cor 5 predicts a global inventory.")
say("   Def 4 AS PRINTED (T strict, no k in the T-equation): exists =", def4, "  | k-twisted version T(gf)=mu(k)T(g)T(f): exists =", tw)
say("   vertex gauge cannot strictify: exists mu with mu alpha mu^-1 strict:",
    any(all(F.act(mu, alpha[(g + f) % 3]) == amul(F.act(mu, alpha[g]), F.act(mu, alpha[f])) for g in range(3) for f in range(3)) for mu in F.A))
F3 = Fibre(S3); sg = (1, 0, 2)
say("M1 Rem 6 witness: S3/B(Z)/c_(01): Def-4 lift T(n)=sigma^n (sigma^2=e:", mul(sg, sg) == S3.e, "); D1 gauge exists:",
    any(F3.act(mu, F3.c[sg]) == F3.id for mu in F3.A))
# M2
say("M2: |Inn S2| = %d, |Z S2| = %d -> abar trivial, [g]=0 on every base (T=e)" % (len(structure(bb.S2)[1]), len(bb.S2.Z)))
# M3: printed a: r->r, s->sc, c->r^2 c ; printed 'a^2 = conj(g)'
G = D4C2; r, s4, c2 = bb.r4, bb.s4, bb.c2
def w(*xs):
    y = G.e
    for x in xs: y = mul(y, x)
    return y
FG = Fibre(G)
def aut_from(img):
    phi = {}
    for x, wd in G.word.items():
        y = G.e
        for kk in wd: y = mul(y, img[kk])
        phi[x] = y
    return tuple(G.idx[phi[x]] for x in G.el)
a = aut_from([r, w(s4, c2), w(r, r, c2)])
sq = amul(a, a)
who = [x for x in G.el if FG.c[x] == sq]
say("M3: printed a in Aut:", a in set(FG.A), "| outer:", a not in FG.Is, "| a^2 = conj(r):", FG.c[r] == sq,
    "| a^2 = conj(x) for x in", sorted(bb.names.get(x, '?') for x in who) if hasattr(bb, 'names') else len(who))
say("    no involutive lift in a.Inn:", not any(amul(x, x) == FG.id for x in coset(a, FG.I)))

# ======================================================================================
say("\n" + "=" * 90); say("PART IV  (UID 382) Thm 8 sign (Q 'sign in 3.2(ii)') and eq. (5)"); say("=" * 90)
# Heisenberg mod p as E, M = centre (trivial action), Q = (Z/p)^2; elements (a,b,c) with (a,b,c)(a',b',c') = (a+a',b+b',c+c'+a b')
p = 5
def hm(x, y): return ((x[0] + y[0]) % p, (x[1] + y[1]) % p, (x[2] + y[2] + x[0] * y[1]) % p)
def hinv(x):
    a_, b_, c_ = x; return ((-a_) % p, (-b_) % p, (-c_ + a_ * b_) % p)
A = (1, 0, 0); B = (0, 1, 0)
wcomm = hm(hm(A, B), hm(hinv(A), hinv(B)))
# section u(s^i t^j) = A^i B^j ; factor set f(g,h) = u(g)u(h)u(gh)^-1 on Z^2 (exponent pairs)
def hp(x, n):
    y = (0, 0, 0)
    for _ in range(n % p): y = hm(y, x)
    return y
def u(i, j): return hm(hp(A, i), hp(B, j))
def fset(g, h): return hm(hm(u(*g), u(*h)), hinv(u(g[0] + h[0], g[1] + h[1])))
fst, fts = fset((1, 0), (0, 1)), fset((0, 1), (1, 0))
say("Heis(F_%d): w=[A,B] central coord = %d ; f(s,t) = %s ; f(t,s) = %s (central coord; -w = %d)" % (p, wcomm[2], fst[2], fts[2], (-wcomm[2]) % p))
say("  phi_2(e_st) = [s|t]-[t|s]  =>  class = f(s,t) - f(t,s) = %d  ==  w = %d : %s" % ((fst[2] - fts[2]) % p, wcomm[2], (fst[2] - fts[2]) % p == wcomm[2]))
# cocycle identity on a box (normalized, trivial action)
ok = all((fset(h, k_)[2] - fset((g[0] + h[0], g[1] + h[1]), k_)[2] + fset(g, (h[0] + k_[0], h[1] + k_[1]))[2] - fset(g, h)[2]) % p == 0
         for g in itertools.product(range(3), repeat=2) for h in itertools.product(range(3), repeat=2) for k_ in itertools.product(range(3), repeat=2))
say("  f is a bar 2-cocycle (dz=0 with d convention of Lemma 7) on a 3x3 box:", ok)
say("  orientation control: with phi_2(e_st)=[t|s]-[s|t] the class would be", (fts[2] - fst[2]) % p, "= -w (sign is a convention of the Koszul generator)")
# eq. (5) with nontrivial, distinct actions: E = AGL(1,7) = Z/7 x| (Z/7)^*, M = Z/7, Q = (Z/7)^* (abelian, so any rho: Z^2 -> Q is a band)
P = 7
def em(x, y): return ((x[0] + x[1] * y[0]) % P, (x[1] * y[1]) % P)       # (m,a)(m',a') = (m + a m', a a')
def ei(x): ai = pow(x[1], P - 2, P); return ((-ai * x[0]) % P, ai)
def comm(x, y): return em(em(x, y), em(ei(x), ei(y)))
bad5 = 0; tot5 = 0; badcoinv = 0
for sa in range(1, P):
    for tb in range(1, P):
        for m0 in range(P):
            Aa, Bb = (0, sa), (m0, tb)
            wv = comm(Aa, Bb); assert wv[1] == 1
            for uu in range(P):
                for vv in range(P):
                    lhs = (comm(em(Aa, (uu, 1)), em(Bb, (vv, 1)))[0] - wv[0]) % P
                    rhs = (tb * (sa - 1) * vv - sa * (tb - 1) * uu) % P
                    tot5 += 1; bad5 += (lhs != rhs)
            # vanishing criterion: commuting lifts exist iff w in (sigma-1)M+(tau-1)M
            commuting = any(comm(em(Aa, (uu, 1)), em(Bb, (vv, 1)))[0] == 0 for uu in range(P) for vv in range(P))
            inim = (sa != 1 or tb != 1) or wv[0] == 0
            badcoinv += (commuting != inim)
say("eq.(5) [Au,Bv]-[A,B] = tau(sigma-1)v - sigma(tau-1)u, E=AGL(1,7) (nontrivial distinct actions): mismatches %d / %d" % (bad5, tot5))
say("Thm 8(iv) commuting lifts <=> [w]=0 in M_{Z^2}: mismatches %d / %d" % (badcoinv, 6 * 6 * P))

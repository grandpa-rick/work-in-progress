r"""Independent checks, round 2, MacBeth 'inventory gerbe' note, content commit 4e7264d.
Rick, 2026-10-11.  Written from the printed definitions; nothing reused from 9eb0941 checks.

Printed definitions used (note p.5-6):
  K = Aut F, Z = Z(K), Inn K = K/Z, pi: K -> Inn K, s a set-section with s(e)=e.
  rho: C -> B(Inn K) a FUNCTOR.  T_f := s(rho f), T_id = e.
  (dagger)  T_g T_f = gamma(g,f) T_{gf},  gamma in Z.     [g] := [gamma] in H^2(N(C); Z).
  Thm 11: strict lift exists  <=>  [g] = 0.
For C = BG the normalised nerve cochains are functions on (G\{e})^n, delta = alternating faces,
trivial action.

R1  D4 witness (four-lot ring), C = B(Z/2): gamma(g,g), coboundary test, all 4 rho.
R2  H^2(N(B Z/2); F2) from the normalised nerve cochain complex (ranks), and gamma's class in it.
R3  2x2 table on B(Z/2), plus negative controls that BREAK it (Z/4, S3, D6 = Z/2 x S3).
R4  Thm D(i) 'only if' and the converse: K from graph automorphisms and a few others;
    [g] != 0 for some (C, rho)  <=>  (*) non-split  (rho = id on B(Inn K) is universal).
R5  Thm 16 (band criterion) on C = B(Z), B(Z/2); Cor 17 under the two definitions of
    'global inventory' printed in the note; a non-functorial inner-valued alpha.
R6  Diamond nerve (4,6,2) and a DAG poset (octahedral) with H^2 != 0.
R7  Z(S_n), |Out S_n| for n <= 4.
"""
import itertools as it
from functools import lru_cache

# ------------------------------------------------------------------ permutation groups
def mul(p, q):                     # (p*q)(i) = p(q(i))
    return tuple(p[i] for i in q)
def inv(p):
    r = [0]*len(p)
    for i, j in enumerate(p): r[j] = i
    return tuple(r)
def ident(n): return tuple(range(n))
def gen_group(gens):
    e = ident(len(gens[0])); S = {e}; fr = [e]
    while fr:
        nf = []
        for x in fr:
            for g in gens:
                y = mul(x, g)
                if y not in S: S.add(y); nf.append(y)
        fr = nf
    return sorted(S)
def order(x):
    e = ident(len(x)); k, y = 1, x
    while y != e: y = mul(y, x); k += 1
    return k

def graph_aut(n, edges):           # automorphism group of a graph on {0..n-1}
    E = {frozenset(e) for e in edges}
    return [p for p in it.permutations(range(n)) if {frozenset((p[a], p[b])) for a, b in E} == E]
def ring(n): return graph_aut(n, [(i, (i+1) % n) for i in range(n)])   # n-lot ring
def cyclic(n): return gen_group([tuple((i+1) % n for i in range(n))])
def direct(A, B):
    na = len(A[0]); return [tuple(list(a) + [na+j for j in b]) for a in A for b in B]
def regular_from_table(elts, m):   # left-regular perm rep of an abstract mult table
    idx = {x: i for i, x in enumerate(elts)}
    return gen_group([tuple(idx[m(g, x)] for x in elts) for g in elts])

# Q8 by unit quaternions as integer 4-vectors
def qm(a, b):
    a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
    return (a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2,
            a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0)
Q8E = [tuple(s*(1 if i == k else 0) for i in range(4)) for k in range(4) for s in (1, -1)]
Q8 = regular_from_table(Q8E, qm)

class Grp:
    def __init__(self, name, elts):
        self.name = name; self.G = list(elts); self.e = ident(len(self.G[0]))
        self.Zc = [z for z in self.G if all(mul(z, x) == mul(x, z) for x in self.G)]
        self.idx = {x: i for i, x in enumerate(self.G)}
        # Inn K realised as conjugation permutations of the index set of K
        self.conj = {k: tuple(self.idx[mul(mul(k, x), inv(k))] for x in self.G) for k in self.G}
        self.Inn = sorted(set(self.conj.values()))
        fib = {}
        for k in self.G: fib.setdefault(self.conj[k], []).append(k)
        self.fibre = fib
        self.sec = {c: (self.e if c == ident(len(self.G)) else min(ks)) for c, ks in fib.items()}
    def abelian(self): return len(self.Zc) == len(self.G)

# ------------------------------------------------------------------ homomorphisms
def is_hom_on_gens(G, gens, imgs, Hmul, He):
    """BFS-extend gens->imgs; return dict if it is a well-defined hom G->H, else None."""
    e = ident(len(G[0])); f = {e: He}; fr = [e]
    while fr:
        nf = []
        for x in fr:
            for g, h in zip(gens, imgs):
                y = mul(x, g); v = Hmul(f[x], h)
                if y in f:
                    if f[y] != v: return None
                else: f[y] = v; nf.append(y)
        fr = nf
    # check every edge (BFS above already checks every edge x->x*g encountered from every x)
    return f
def small_gens(G):
    gens = []; H = [ident(len(G[0]))]
    for x in G:
        if x not in H: gens.append(x); H = gen_group(gens)
        if len(H) == len(G): break
    return gens
def homs(G, H, Hmul, He):
    gens = small_gens(G); out = []
    for imgs in it.product(H, repeat=len(gens)):
        f = is_hom_on_gens(G, gens, imgs, Hmul, He)
        if f is not None: out.append(f)
    return out

# ------------------------------------------------------------------ nerve cochains for BG
def gamma_cocycle(K, G, rho):
    """(dagger): gamma(g,f) = T_g T_f T_{gf}^{-1}, with T = s o rho."""
    T = {g: K.sec[rho[g]] for g in G}
    gam = {}
    for g in G:
        for f in G:
            z = mul(mul(T[g], T[f]), inv(T[mul(g, f)]))
            assert z in K.Zc, "gamma not central"
            gam[(g, f)] = z
    return gam
def check_cocycle(K, G, gam):
    for h, g, f in it.product(G, repeat=3):
        lhs = mul(gam[(h, g)], gam[(mul(h, g), f)]); rhs = mul(gam[(g, f)], gam[(h, mul(g, f))])
        if lhs != rhs: return False
    return True
def is_coboundary(K, G, gam):
    """exists normalised c: G -> Z with gam(g,f) = c_g c_f c_{gf}^{-1} ?  brute force."""
    e = ident(len(G[0])); nonid = [g for g in G if g != e]
    if len(K.Zc) ** len(nonid) > 3_000_000: raise RuntimeError("too big")
    for vals in it.product(K.Zc, repeat=len(nonid)):
        c = dict(zip(nonid, vals)); c[e] = K.e
        if all(gam[(g, f)] == mul(mul(c[g], c[f]), inv(c[mul(g, f)])) for g in G for f in G):
            return True
    return False
def has_strict_lift(K, G, rho):
    """Thm 11 other side, independently: is there a hom T: G -> K with pi o T = rho ?"""
    gens = small_gens(G)
    for imgs in it.product(*[K.fibre[rho[g]] for g in gens]):
        f = is_hom_on_gens(G, gens, imgs, mul, K.e)
        if f is not None and all(K.conj[f[g]] == rho[g] for g in G): return True
    return False
def g_class_nonzero(K, G, rho):
    gam = gamma_cocycle(K, G, rho)
    assert check_cocycle(K, G, gam)
    nz = not is_coboundary(K, G, gam)
    assert nz == (not has_strict_lift(K, G, rho)), "Thm 11 fails?!"
    return nz

def rank_mod2(rows):
    rows = [r[:] for r in rows]; rk = 0; ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] % 2), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(len(rows)):
            if i != rk and rows[i][c] % 2: rows[i] = [(a+b) % 2 for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk

# ------------------------------------------------------------------ groups
Z2, Z3, Z4 = cyclic(2), cyclic(3), cyclic(4)
V4 = direct(Z2, Z2)
S3 = [tuple(p) for p in it.permutations(range(3))]
KS = {
    'pt (trivial)': Grp('pt', [(0,)]),
    'Z/2 (two lots)': Grp('Z2', Z2),
    'Z/4 (abelian)': Grp('Z4', Z4),
    'S3 (centreless)': Grp('S3', S3),
    'D4 = Aut(4-ring)': Grp('D4', ring(4)),
    'Q8': Grp('Q8', Q8),
    'D6 = Aut(6-ring)': Grp('D6', ring(6)),
    'Z/2 x S3': Grp('Z2xS3', direct(Z2, S3)),
    'D8 = Aut(8-ring)': Grp('D8', ring(8)),
    'Z/2 x D4': Grp('Z2xD4', direct(Z2, ring(4))),
    'S4': Grp('S4', [tuple(p) for p in it.permutations(range(4))]),
}
BASES = {'Z/2': Z2, 'Z/3': Z3, 'Z/4': Z4, '(Z/2)^2': V4}

def section(t): print('\n' + '=' * 78 + '\n' + t + '\n' + '=' * 78)

# ================================================================== R1  D4 witness
section('R1  D4 witness on C = B(Z/2)')
D4 = KS['D4 = Aut(4-ring)']
print('|Aut(4-ring)| =', len(D4.G), ' |Z| =', len(D4.Zc), ' |Inn| =', len(D4.Inn),
      ' abelian:', D4.abelian())
r = (1, 2, 3, 0); s = (0, 3, 2, 1)
assert r in D4.G and s in D4.G
r2 = mul(r, r)
print('rotation r =', r, ' r^2 =', r2, ' Z(D4) =', D4.Zc, ' r^2 central:', r2 in D4.Zc)
g = (1, 0); e2 = (0, 1)
rho = {e2: D4.conj[D4.e], g: D4.conj[r]}
assert mul(rho[g], rho[g]) == rho[e2], 'rho not a functor'
for Tg in D4.fibre[D4.conj[r]]:      # every lift of [r], not only the section
    print('  lift T_g =', Tg, ' T_g T_g T_id^{-1} =', mul(Tg, Tg), '(= r^2?', mul(Tg, Tg) == r2, ')')
# normalised 1-cochains: c_e = 1, c_g in Z ; (delta c)(g,g) = c_g c_g c_{g^2}^{-1} = c_g^2
print('  (delta c)(g,g) over all c_g in Z:', [mul(c, c) for c in D4.Zc], '-> always e')
for name, x in [('[e]', D4.e), ('[r]', r), ('[s]', s), ('[rs]', mul(r, s))]:
    rr = {e2: D4.conj[D4.e], g: D4.conj[x]}
    print(f'  rho(g) = {name:5s}: [g] != 0 ? {g_class_nonzero(D4, Z2, rr)}')

# ================================================================== R2  nerve H^2
section('R2  normalised nerve cochains of B(Z/2), F2 coefficients')
G = Z2; nonid = [x for x in G if x != ident(2)]
def simplices(n): return list(it.product(nonid, repeat=n))
def faces(sig):
    n = len(sig); out = [sig[1:]]
    for i in range(1, n): out.append(sig[:i-1] + (mul(sig[i], sig[i-1]),) + sig[i+1:])
    out.append(sig[:-1]); return out
def dmat(n):   # delta: C^n -> C^{n+1}, rows = (n+1)-simplices, cols = n-simplices; degenerate faces drop
    S_n, S_n1 = simplices(n), simplices(n+1); col = {x: i for i, x in enumerate(S_n)}
    M = []
    for sig in S_n1:
        row = [0]*len(S_n)
        for fsig in faces(sig):
            if all(x != ident(2) for x in fsig) and fsig in col: row[col[fsig]] += 1
        M.append([v % 2 for v in row])
    return M
dims = [len(simplices(n)) for n in range(4)]
rks = [rank_mod2(dmat(n)) if dims[n] and dims[n+1] else 0 for n in range(3)]
H2 = dims[2] - rks[2] - rks[1]
print('dims C^0..C^3 =', dims, ' ranks delta^0..delta^2 =', rks, ' dim H^2 =', H2)
# gamma of R1 on the unique nondegenerate 2-simplex (g,g): value r^2 -> 1 in F2
gvec = [1]; print('gamma as F2 2-cochain:', gvec, ' cocycle:', all(sum(a*b for a, b in zip(row, gvec)) % 2 == 0 for row in dmat(2)),
      ' in im delta^1 :', rank_mod2(dmat(1) and [list(c) for c in zip(*dmat(1))] + [gvec]) == rks[1])

# ================================================================== R3  2x2 table
section('R3  2x2 truth table on C = B(Z/2); weld side imported ([4] / 242c7ea Ex.15: eps(g,g) = 1)')
def g_on_BZ2(Kname, rep):
    K = KS[Kname]; rr = {e2: K.conj[K.e], g: K.conj[rep]}
    return g_class_nonzero(K, Z2, rr)
weld = {'trivial (D=0)': 0, 're-entrant (eps)': 1}   # normalised cocycle value at (g,g), F2
inv_data = {'pt': ('pt (trivial)', (0,)), 'D4-ring, rho=[r]': ('D4 = Aut(4-ring)', r),
            'D4-ring, rho=[s]  (extra)': ('D4 = Aut(4-ring)', s)}
for lab, (Kn, rep) in inv_data.items():
    gv = int(g_on_BZ2(Kn, rep))
    cells = []
    for wl, wv in weld.items():
        cells.append(f'({"!=0" if gv else "0"},{"!=0" if wv else "0"})' + (' [g]=[w]' if gv and wv and gv == wv else ''))
    print(f'  {lab:28s} |Z|={len(KS[Kn].Zc)}  ' + '   '.join(cells))
print('  -> all four patterns realised; bottom-right: gamma(g,g)=r^2 <-> 1 = eps(g,g) on the nose')
print('Negative controls (replace the D4-ring row by another fibre; does the bottom row survive?)')
for Kn in ['Z/4 (abelian)', 'S3 (centreless)', 'D6 = Aut(6-ring)', 'Z/2 x S3', 'Q8', 'D8 = Aut(8-ring)']:
    K = KS[Kn]
    rhos = homs(Z2, K.Inn, mul, ident(len(K.G)))
    nz = [g_class_nonzero(K, Z2, rr) for rr in rhos]
    print(f'  {Kn:20s} |Z|={len(K.Zc)} |Inn|={len(K.Inn)}: #rho = {len(rhos)}, #[g]!=0 = {sum(nz)}'
          + ('   BREAKS table (no [g]!=0 row)' if not any(nz) else ''))

# ================================================================== R4  Thm D
section('R4  Thm D: nonabelian & Z != 1  vs  existence of [g] != 0')
print(f'{"K":20s}{"|K|":>4s}{"|Z|":>4s}{"nonab&Z!=1":>11s}{"split(*)":>9s}   [g]!=0 count per base (of #rho)')
for Kn, K in KS.items():
    nab = (not K.abelian()) and len(K.Zc) > 1
    # split (*)  <=>  rho = id on B(Inn K) lifts  <=>  [*] = 0  (Thm 11 with C = B(Inn K))
    Inn = K.Inn; rid = {x: x for x in Inn}
    split = has_strict_lift(K, Inn, rid)
    per = []
    for bn, B in BASES.items():
        rhos = homs(B, K.Inn, mul, ident(len(K.G)))
        per.append(f'{bn}:{sum(g_class_nonzero(K, B, rr) for rr in rhos)}/{len(rhos)}')
    if len(Inn) <= 8 and len(K.Zc) ** (len(Inn)-1) <= 3_000_000:
        univ = g_class_nonzero(K, Inn, rid); assert univ == (not split)
    print(f'{Kn:20s}{len(K.G):4d}{len(K.Zc):4d}{str(nab):>11s}{str(split):>9s}   ' + ' '.join(per))
print('Expect: Thm D(i) "only if" holds on every row; converse FAILS for D6 = Z/2 x S3 (nonab, Z=Z/2, split).')

# ================================================================== R5  band criterion, Cor 17
section('R5  Thm 16 (band) and Cor 17 (staged obstruction)')
def aut_group(K):
    gens = small_gens(K.G); auts = []
    for imgs in it.product(K.G, repeat=len(gens)):
        f = is_hom_on_gens(K.G, gens, imgs, mul, K.e)
        if f is not None and len(set(f.values())) == len(K.G):
            auts.append(tuple(K.idx[f[x]] for x in K.G))   # as perm of indices
    return sorted(auts)
for Kn in ['D4 = Aut(4-ring)', 'Q8', 'S3 (centreless)', 'Z/2 x S3', 'Z/4 (abelian)']:
    K = KS[Kn]; A = aut_group(K); Inn = set(K.Inn)
    # B(Z): alpha(t) = a, any automorphism.  beta = [q(a)] in conj classes of Out; beta = 0 iff a in Inn.
    # RHS of Thm 16: exists mu in Aut with mu a mu^-1 in Inn.
    ok = all(((a in Inn) == any(mul(mul(m, a), inv(m)) in Inn for m in A)) for a in A)
    print(f'  {Kn:18s} |Aut|={len(A):3d} |Inn|={len(Inn):2d} |Out|={len(A)//len(Inn)}  Thm16 on B(Z) all a: {ok}')

print('\n Cor 17 on C = B(Z/2), STRICT alpha (alpha(g)^2 = id), compared with (beta=0 and [g]=0):')
print('   D1 = "relabelling mu trivialising alpha outright" (proof of Cor 17, l.381)')
print('   D2 = "strict T: G->K and mu with mu alpha mu^-1 = conj o T" (Def 6 / Thm 11 notion)')
for Kn in ['D4 = Aut(4-ring)', 'Q8', 'Z/2 x S3', 'S3 (centreless)', 'Z/4 (abelian)']:
    K = KS[Kn]; A = aut_group(K); Inn = set(K.Inn); idn = ident(len(K.G))
    bad1 = []; bad2 = []
    for a in A:
        if mul(a, a) != idn: continue          # strict functor B(Z/2) -> B(Aut K)
        beta0 = a in Inn
        if beta0:
            # gauge: mu = id already inner-valued; rho(g) = a in Inn K
            gz = g_class_nonzero(K, Z2, {e2: idn, g: a})
            stagedOK = not gz
        else:
            stagedOK = False
        D1 = any(mul(mul(m, a), inv(m)) == idn for m in A)
        D2 = any(mul(mul(m, a), inv(m)) == K.conj[t] for m in A for t in K.G if mul(t, t) == K.e)
        if D1 != stagedOK: bad1.append(a)
        if D2 != stagedOK: bad2.append(a)
    print(f'   {Kn:18s} D1 disagrees on {len(bad1)} alpha;  D2 disagrees on {len(bad2)} alpha')
# the explicit D1 counterexample
a = D4.conj[s]
print('   explicit: K=D4, alpha(g)=conj(s): beta=0:', a in set(D4.Inn), ' [g]!=0:',
      g_class_nonzero(D4, Z2, {e2: ident(8), g: a}), ' strict lift T(g)=s exists:', mul(s, s) == D4.e,
      ' D1 (mu alpha mu^-1 = id) possible:', any(mul(mul(m, a), inv(m)) == ident(8) for m in aut_group(D4)))

print('\n Non-functorial inner-valued alpha (note: alpha "functorial only up to inner"):')
G3 = Z3; t = G3[1]; t2 = mul(t, t)
alpha = {ident(3): ident(8), t: D4.conj[s], t2: D4.conj[s]}
upto_inner = all(True for _ in [0])  # every value is inner, so alpha(gf) = iota alpha(g) alpha(f) with iota in Inn holds
fun = all(mul(alpha[x], alpha[y]) == alpha[mul(x, y)] for x in G3 for y in G3)
print('   C=B(Z/3), alpha(t)=alpha(t^2)=conj(s): inner-valued (so beta=0 by Thm 16), functorial up to inner: True,'
      ' a functor C->B(Inn K)?', fun, ' -> Def 6 rho does not exist; [g] undefined')

# ================================================================== R6  nerves
section('R6  diamond nerve (free) and an octahedral DAG poset')
def nerve_F2(objs, arrows, comp):
    """arrows: list of (name, src, tgt) non-identity; comp(a,b) -> name of b∘a (a then b) or None."""
    S = {0: [(o,) for o in objs], 1: [(a,) for a in arrows]}
    S[2] = [(a, b) for a in arrows for b in arrows if a[2] == b[1]]
    S[3] = [(a, b, c) for (a, b) in S[2] for c in arrows if b[2] == c[1]]
    def fc(sig, n):
        if n == 1:
            a = sig[0]; return [((a[2],), 1), ((a[1],), 1)]
        if n == 2:
            a, b = sig; ab = comp(a, b); out = [((b,), 1), ((a,), 1)]
            if ab is not None: out.append(((ab,), 1))
            return out
        a, b, c = sig; out = [((b, c), 1), ((a, b), 1)]
        ab, bc = comp(a, b), comp(b, c)
        if bc is not None: out.append(((a, bc), 1))
        if ab is not None: out.append(((ab, c), 1))
        return out
    dims = [len(S[n]) for n in range(4)]; rks = []
    for n in range(3):
        col = {x: i for i, x in enumerate(S[n])}; M = []
        for sig in S[n+1]:
            row = [0]*len(S[n])
            for f, v in fc(sig, n+1):
                if f in col: row[col[f]] += v
            M.append([x % 2 for x in row])
        rks.append(rank_mod2(M) if M and M[0] else 0)
    H = [dims[0]-rks[0], dims[1]-rks[1]-rks[0], dims[2]-rks[2]-rks[1]]
    return dims[:3], rks, H
# diamond, free: p:a->b, q:a->c, r:b->d, s:c->d, pr, qs
P = ('p', 'a', 'b'); Q = ('q', 'a', 'c'); R = ('r', 'b', 'd'); S_ = ('s', 'c', 'd')
PR = ('pr', 'a', 'd'); QS = ('qs', 'a', 'd')
arr = [P, Q, R, S_, PR, QS]
cm = {('p', 'r'): PR, ('q', 's'): QS}
print('  free diamond (C^0,C^1,C^2), ranks, H^0..2 :', nerve_F2(list('abcd'), arr, lambda x, y: cm.get((x[0], y[0]))))
# poset a1,a2 < b1,b2 < c1,c2 (all cross-level relations): order complex = S^2
objs = ['a1', 'a2', 'b1', 'b2', 'c1', 'c2']; lvl = {o: int(o[0] == 'b') + 2*int(o[0] == 'c') for o in objs}
arrP = [(x+y, x, y) for x in objs for y in objs if lvl[x] < lvl[y]]
by = {(a[1], a[2]): a for a in arrP}
print('  octahedral poset (a DAG; all routes commute) :', nerve_F2(objs, arrP, lambda x, y: by.get((x[1], y[2]))))
print('  -> H^2(N(C);F2) != 0 on a DAG category that is not free; but pi_1 = 1 so [g] = rho^*[*] = 0 there.')

# ================================================================== R7  S_n
section('R7  S_n: centre and |Out|')
for n in range(1, 5):   # n = 5, 6: Rotman Thms 7.5-7.7 (not recomputed)
    K = Grp(f'S{n}', [tuple(p) for p in it.permutations(range(n))]) if n > 1 else Grp('S1', [(0,)])
    A = aut_group(K)
    print(f'  S{n}: |Z| = {len(K.Zc)}, |Inn| = {len(K.Inn)}, |Aut| = {len(A)}, |Out| = {len(A)//len(K.Inn)}')
print('\nALL CHECKS RAN')

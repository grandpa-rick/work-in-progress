"""Rick, 2026-10-12: independent checks for MacBeth's band-boundary note (content 2fa15ef).
Written from the printed statements; reuses nothing from MacBeth's scripts.
Groups are permutation groups; Aut computed by brute force over generator images."""
import itertools, sys
from collections import deque

def mul(a, b):  # (a*b)(i) = a(b(i))  -- apply b first
    return tuple(a[i] for i in b)
def inv(a):
    r = [0]*len(a)
    for i, x in enumerate(a): r[x] = i
    return tuple(r)

class Grp:
    def __init__(self, gens, name):
        self.name = name; self.gens = gens
        n = len(gens[0]); e = tuple(range(n)); self.e = e
        seen = {e: ()}; q = deque([e])
        while q:
            x = q.popleft()
            for k, g in enumerate(gens):
                y = mul(x, g)
                if y not in seen: seen[y] = seen[x] + (k,); q.append(y)
        self.el = list(seen); self.word = seen; self.idx = {x: i for i, x in enumerate(self.el)}
        self.N = len(self.el)
        self.Z = [z for z in self.el if all(mul(z, g) == mul(g, z) for g in self.el)]
    def order(self, x):
        k, y = 1, x
        while y != self.e: y = mul(y, x); k += 1
        return k

def automorphisms(G):
    """All automorphisms as dicts (tuple index arrays over G.el)."""
    gens = G.gens; ords = [G.order(g) for g in gens]
    cands = [[x for x in G.el if G.order(x) == o] for o in ords]
    auts = []
    for imgs in itertools.product(*cands):
        phi = {}
        for x, w in G.word.items():
            y = G.e
            for k in w: y = mul(y, imgs[k])
            phi[x] = y
        if len(set(phi.values())) != G.N: continue
        ok = all(phi[mul(x, g)] == mul(phi[x], imgs[k]) for x in G.el for k, g in enumerate(gens))
        if ok: auts.append(tuple(G.idx[phi[x]] for x in G.el))
    return auts

def amul(a, b):  # a o b, apply b first
    return tuple(a[i] for i in b)
def ainv(a):
    r = [0]*len(a)
    for i, x in enumerate(a): r[x] = i
    return tuple(r)

def conj_aut(G, g):
    gi = inv(g)
    return tuple(G.idx[mul(mul(g, x), gi)] for x in G.el)

def structure(G):
    A = automorphisms(G); Aset = set(A)
    I = sorted({conj_aut(G, g) for g in G.el}); Iset = set(I)
    return A, I, Iset

def coset(a, I):  # left coset a*Inn (= Inn*a since normal), as frozenset
    return frozenset(amul(a, i) for i in I)

out = []
def say(*s):
    msg = " ".join(str(x) for x in s); print(msg); out.append(msg)

# ---------- groups ----------
def cyc(n, k=0, N=None):
    N = N or n; p = list(range(N))
    for i in range(n): p[k+i] = k + (i+1) % n
    return tuple(p)
r4 = (1,2,3,0,4,5); s4 = (0,3,2,1,4,5); c2 = (0,1,2,3,5,4)
D4C2 = Grp([r4, s4, c2], "D4xC2")
D4 = Grp([(1,2,3,0), (0,3,2,1)], "D4")
S3 = Grp([(1,0,2), (1,2,0)], "S3")
S2 = Grp([(1,0)], "S2")

# ---------- (A) Thm 6 typing: non-strict alpha breaks the displayed cocycle line ----------
say("== (A) Thm 6 (=>) cocycle line under the note's own up-to-inner hypothesis ==")
A4, I4, I4s = structure(D4)
s = (0,3,2,1); cs = conj_aut(D4, s); idA = tuple(range(D4.N))
# C = B(Z/3), morphisms e,t,t^2; alpha(e)=id, alpha(t)=alpha(t^2)=conj(s)
alpha = {0: idA, 1: cs, 2: cs}
upto = all(amul(alpha[(a+b)%3], ainv(amul(alpha[a], alpha[b]))) in I4s for a in range(3) for b in range(3))
strict = all(alpha[(a+b)%3] == amul(alpha[a], alpha[b]) for a in range(3) for b in range(3))
say("D4, B(Z/3), alpha(t)=alpha(t^2)=conj(s): functorial up to Inn:", upto, "| strict:", strict)
say("  => with mu=id, alpha~=alpha is inner-valued but alpha~(t*t) != alpha~(t)alpha~(t):",
    alpha[2] != amul(alpha[1], alpha[1]))

# ---------- (B) Cor 7 sufficiency under printed def 'mu . alpha == e' ----------
say("== (B) Cor 7 'strict global inventory = relabelling mu with mu.alpha == e' ==")
A3, I3, I3s = structure(S3)
say("S3: |Aut|=%d |Inn|=%d |Z|=%d" % (len(A3), len(I3), len(S3.Z)))
t = (1,0,2); ct = conj_aut(S3, t)
# C = B(Z): alpha(t)=ct strict functor; beta=0 since Out=1; [g] in H^2(S^1;Z)=0 and Z=1.
triv = any(amul(amul(m, ct), ainv(m)) == tuple(range(S3.N)) for m in A3)
say("S3, C=B(Z), alpha(t)=conj((01)): beta=0 (Out=1), [g]=0 (Z=1, dim 1); exists mu with mu.alpha=e:", triv)
# B(Z/2) version, also D4 conj(s)
triv4 = any(amul(amul(m, cs), ainv(m)) == idA for m in A4)
say("D4, C=B(Z/2), alpha(g)=conj(s) (s involution => strict lift T(g)=s, [g]=0): exists mu:", triv4)

# ---------- (C) Thm 8, n=2 'centre exception' ----------
say("== (C) Thm 8 n=2: G=S2 abelian ==")
A2, I2, I2s = structure(S2)
say("S2: |Aut|=%d |Inn|=%d |Z|=%d -> Inn=1, rho-bar trivial, gerbe of lifts of trivial functor is neutral: [g]=0 for every C"
    % (len(A2), len(I2), len(S2.Z)))
# explicit: over B(Z/2), lifts of trivial rho-bar: T(g)=e works.
say("  sanity: H^2(B(Z/2);Z/2)=Z/2 is nonzero, yet the Def-6 class rho*[*] = 0 since * : Z->S2->1 splits.")

# ---------- (D) Section 6 witness D4 x C2 ----------
say("== (D) Section 6 witness G = D4 x C2 ==")
G = D4C2
AG, IG, IGs = structure(G)
say("|G|=%d |Z(G)|=%d |Aut G|=%d |Inn G|=%d |Out G|=%d" % (G.N, len(G.Z), len(AG), len(IG), len(AG)//len(IG)))
r, s, c = r4, s4, c2
def word(*xs):
    y = G.e
    for x in xs: y = mul(y, x)
    return y
r2, r3 = word(r, r), word(r, r, r)
def aut_from(img):  # images of r,s,c
    phi = {}
    for x, w in G.word.items():
        y = G.e
        for k in w: y = mul(y, img[k])
        phi[x] = y
    return tuple(G.idx[phi[x]] for x in G.el)
Aa = aut_from([r3, s, word(r2, c)])   # A: r->r^3, s->s, c->r^2 c
Bb = aut_from([word(r3, c), s, c])    # B: r->r^3 c, s->s, c->c
AGs = set(AG)
say("A in Aut G:", Aa in AGs, "| B in Aut G:", Bb in AGs)
say("A outer:", Aa not in IGs, "| B outer:", Bb not in IGs, "| AB outer:", amul(Aa, Bb) not in IGs)
idG = tuple(range(G.N))
say("A^2=id:", amul(Aa, Aa) == idG, "| B^2=id:", amul(Bb, Bb) == idG)
comm = amul(amul(Aa, Bb), amul(ainv(Aa), ainv(Bb)))
say("[A,B] inner (so [A],[B] commute in Out):", comm in IGs)
which = [g for g in G.el if conj_aut(G, g) == comm]
names = {}
for i in range(4):
    for j in range(2):
        for k in range(2):
            names[word(*([r]*i + [s]*j + [c]*k))] = "r^%d s^%d c^%d" % (i, j, k)
say("[A,B] = conj(g) for g in", sorted(names[g] for g in which))
cs_G = conj_aut(G, s)
say("[A,B] == conj(s):", comm == cs_G)
# action of A,B on Inn G (by a*conj_g*a^-1 = conj_{a(g)}): trivial?
actA = all(amul(amul(Aa, i), ainv(Aa)) == i for i in IG)
actB = all(amul(amul(Bb, i), ainv(Bb)) == i for i in IG)
say("A, B act trivially on Inn G:", actA, actB, "-> coinvariants M_{<A,B>} = Inn G = (Z/2)^2 and obstruction = [A,B] itself")
# exhaustive: any commuting lifts a*A, b*B, a,b in Inn (4x4)
cnt = sum(1 for a in IG for b in IG if amul(amul(a, Aa), amul(b, Bb)) == amul(amul(b, Bb), amul(a, Aa)))
say("commuting lift pairs over Inn x Inn (%d pairs):" % (len(IG)**2), cnt)
# stronger: search all of Aut G for x in A*Inn, y in B*Inn commuting (same thing, but via cosets)
cA, cB = coset(Aa, IG), coset(Bb, IG)
cnt2 = sum(1 for x in cA for y in cB if amul(x, y) == amul(y, x))
say("commuting pairs (x,y) in A Inn x B Inn:", cnt2)
# 16x16 version as printed (conj by g, h in G)
cnt3 = sum(1 for g in G.el for h in G.el
           if (lambda X, Y: amul(X, Y) == amul(Y, X))(amul(conj_aut(G, g), Aa), amul(conj_aut(G, h), Bb)))
say("printed 16x16 version (g,h in G): commuting count", cnt3)
# Z -> G -> Inn G split?
Zs = set(G.Z)
split = False
for H in itertools.combinations(G.el, len(IG)):
    Hs = set(H)
    if G.e not in Hs: continue
    if any(mul(x, y) not in Hs for x in H for y in H): continue
    if len(Hs & Zs) == 1: split = True; break
say("Z(G) -> G -> Inn G split:", split, "(Rem 10 caveat (2): 'same group also witnesses' the central stage)")

# ---------- (E) is T^2 the 'first' base? Out-involutions with no involutive lift (=> delta != 0 on B(Z/2)) ----------
say("== (E) delta on B(Z/2): outer involutions with no involutive lift ==")
def check_Z2(Gr):
    A, I, Is = structure(Gr)
    cos = {}
    for a in A: cos.setdefault(coset(a, I), []).append(a)
    idg = tuple(range(Gr.N)); bad = 0; inv_classes = 0
    for C, mem in cos.items():
        a = mem[0]
        if amul(a, a) in Is and a not in Is:
            inv_classes += 1
            if not any(amul(x, x) == idg for x in mem): bad += 1
    return len(A), len(I), inv_classes, bad
Q8 = Grp([(1,2,3,0,5,6,7,4)[:0] or None][0:0] or [None], "x") if False else None
# Q8 as regular permutation rep on 8 points: elements 1,i,j,k,-1,-i,-j,-k
qm = {}
names8 = ['1','i','j','k','-1','-i','-j','-k']
base = {('1','1'):'1',('1','i'):'i',('1','j'):'j',('1','k'):'k',('i','1'):'i',('i','i'):'-1',('i','j'):'k',('i','k'):'-j',
        ('j','1'):'j',('j','i'):'-k',('j','j'):'-1',('j','k'):'i',('k','1'):'k',('k','i'):'j',('k','j'):'-i',('k','k'):'-1'}
def qmul(x, y):
    sx = x.startswith('-'); sy = y.startswith('-'); z = base[(x.lstrip('-'), y.lstrip('-'))]
    neg = sx ^ sy ^ z.startswith('-'); z = z.lstrip('-')
    return ('-' + z) if neg else z
def lmul_perm(g): return tuple(names8.index(qmul(g, h)) for h in names8)
Q8 = Grp([lmul_perm('i'), lmul_perm('j')], "Q8")
C4C2 = Grp([(1,2,3,0,4,5), (0,1,2,3,5,4)], "C4xC2")
D8 = Grp([cyc(8), tuple((-i) % 8 for i in range(8))], "D8")
for Gr in [D4, Q8, S3, D4C2, C4C2, D8]:
    a, i, ic, bad = check_Z2(Gr)
    say("%-6s |Aut|=%3d |Inn|=%2d  outer-involution classes=%d  with no involutive lift=%d" % (Gr.name, a, i, ic, bad))

# ---------- (F) S6 outer automorphism as printed ----------
say("== (F) S6: Phi(012345)=[1,4,2,5,0,3], Phi((01))=(01)(25)(34) ==")
sig = (1,2,3,4,5,0); tau = (1,0,2,3,4,5)
Psig = (1,4,2,5,0,3); Ptau = (1,0,5,4,3,2)
# build map on S6 via words in sig, tau; check hom + bijection
S6 = Grp([sig, tau], "S6")
phi = {}
for x, w in S6.word.items():
    y = S6.e
    for k in w: y = mul(y, [Psig, Ptau][k])
    phi[x] = y
bij = len(set(phi.values())) == 720
hom = all(phi[mul(x, g)] == mul(phi[x], [Psig, Ptau][k]) for x in S6.el for k, g in enumerate([sig, tau]))
def ctype(p):
    seen=set(); L=[]
    for i in range(len(p)):
        if i in seen: continue
        j=i; l=0
        while j not in seen: seen.add(j); j=p[j]; l+=1
        L.append(l)
    return tuple(sorted(L, reverse=True))
transp = [x for x in S6.el if ctype(x) == (2,1,1,1,1)]
say("hom:", hom, "| bijection:", bij, "| transpositions -> types:", sorted({ctype(phi[x]) for x in transp}))

# ---------- (E') explicit B(Z/2) witness in D4 x C2 ----------
say("== (E') explicit outer involution of D4xC2 with no involutive lift ==")
cos = {}
for a in AG: cos.setdefault(coset(a, IG), []).append(a)
for C_, mem in cos.items():
    a = mem[0]
    if amul(a, a) in IGs and a not in IGs and not any(amul(x, x) == idG for x in mem):
        imgs = [names[G.el[a[G.idx[g]]]] for g in (r, s, c)]
        sq = amul(a, a); sqg = sorted(names[g] for g in G.el if conj_aut(G, g) == sq)
        acts = all(amul(amul(a, i), ainv(a)) == i for i in IG)
        say("  a: r->%s, s->%s, c->%s ; a^2 = conj of %s ; a acts trivially on Inn: %s" % (*imgs, sqg[:2], acts))

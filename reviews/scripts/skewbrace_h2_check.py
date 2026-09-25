#!/usr/bin/env python3
"""
Referee check for MacBeth, "The obstruction to a skew-brace bicrossed decomposition
is a bilinear second-cohomology class" (wip b243b77, 25 Sep 2026).

Independent of MacBeth's enumerator.  What it does:

 (1) Enumerates all skew braces with additive group Z2 x Z4 (these are braces)
     and with additive group D4 (genuinely skew), via lambda maps
     lambda : A -> Aut(A,+) with lambda_{a + lambda_a(b)} = lambda_a lambda_b.
 (2) For every ideal D that is a trivial brace, computes the split triple
     (mult-split, add-split, sb-split) by brute-force search for complements.
 (3) For every (Y,Y,N) instance: computes the Schreier data (beta, taubar,
     nu, mu, sigma) from a section, builds H^2_Sb(H,D) INTRINSICALLY
     (Z^2 = normalised pairs (beta', taubar') for which the twisted structure on
     H x D is a skew brace inducing the same actions; B^2 = shifts by change of
     section computed directly, not by formula), and computes ker Phi.
 (4) Tests Lemma 2 (the "Lemma 3.6" of the cover e-mail): beta' - beta =
     mu_{h2}(theta(h1)) - theta(h1+h2) + theta(h2) for EVERY theta, in every
     trivial-brace-ideal instance found (abelian and non-abelian additive group),
     and the RY taubar-shift formula; also tests MacBeth's transcription of (3.10).
 (5) Checks RY Example 3.4 (splitting data of I = Z2 x <2> in SmallSkewbrace(8,8)).

Run:  python3 skewbrace_h2_check.py
"""
import itertools, sys

# ---------------------------------------------------------------- groups
class Group:
    def __init__(self, elems, op, name):
        self.E = list(elems); self.n = len(self.E); self.name = name
        self.idx = {e: i for i, e in enumerate(self.E)}
        self.mul = [[self.idx[op(a, b)] for b in self.E] for a in self.E]
        self.zero = next(i for i in range(self.n)
                         if all(self.mul[i][j] == j for j in range(self.n)))
        self.inv = [next(j for j in range(self.n) if self.mul[i][j] == self.zero)
                    for i in range(self.n)]
    def add(self, a, b): return self.mul[a][b]
    def neg(self, a): return self.inv[a]

def automorphisms(G):
    n = G.n
    # generate via images of a generating set
    gens = []
    span = {G.zero}
    def closure(S):
        S = set(S)
        while True:
            new = {G.mul[a][b] for a in S for b in S} | S
            if new == S: return S
            S = new
    for g in range(n):
        if g not in span:
            gens.append(g); span = closure(span | {g})
        if len(span) == n: break
    auts = []
    for imgs in itertools.product(range(n), repeat=len(gens)):
        # build map by BFS words
        f = {G.zero: G.zero}
        frontier = [G.zero]; ok = True
        while frontier and ok:
            nxt = []
            for x in frontier:
                for g, im in zip(gens, imgs):
                    y = G.mul[x][g]; fy = G.mul[f[x]][im]
                    if y in f:
                        if f[y] != fy: ok = False; break
                    else:
                        f[y] = fy; nxt.append(y)
                if not ok: break
            frontier = nxt
        if not ok or len(f) != n or len(set(f.values())) != n: continue
        if all(f[G.mul[a][b]] == G.mul[f[a]][f[b]] for a in range(n) for b in range(n)):
            auts.append(tuple(f[i] for i in range(n)))
    return auts

def compose(p, q):  # (p o q)(x) = p(q(x))
    return tuple(p[q[x]] for x in range(len(q)))

def skew_braces(G):
    """All lambda maps A -> Aut(A,+) with lambda_{a+lambda_a(b)} = lambda_a lambda_b."""
    auts = automorphisms(G); n = G.n
    ident = tuple(range(n))
    order = [G.zero] + [i for i in range(n) if i != G.zero]
    res = []
    lam = {G.zero: ident}
    def consistent():
        for a in lam:
            for b in lam:
                c = G.mul[a][lam[a][b]]
                if c in lam and lam[c] != compose(lam[a], lam[b]): return False
        return True
    def rec(k):
        if k == n:
            res.append([lam[i] for i in range(n)]); return
        a = order[k]
        for f in auts:
            lam[a] = f
            if consistent(): rec(k + 1)
            del lam[a]
    rec(1)
    return res

class SkewBrace:
    def __init__(self, G, lam):
        self.G = G; self.lam = lam; n = G.n
        self.circ = [[G.mul[a][lam[a][b]] for b in range(n)] for a in range(n)]
        self.cinv = [next(b for b in range(n) if self.circ[a][b] == G.zero) for a in range(n)]
        # sanity: brace axioms
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    assert self.circ[a][self.circ[b][c]] == self.circ[self.circ[a][b]][c]
                    lhs = self.circ[a][G.mul[b][c]]
                    rhs = G.mul[G.mul[self.circ[a][b]][G.neg(a)]][self.circ[a][c]]
                    assert lhs == rhs
    def p(self, a, b): return self.G.mul[a][b]
    def m(self, a): return self.G.neg(a)
    def o(self, a, b): return self.circ[a][b]
    def oi(self, a): return self.cinv[a]

def subgroups(n, op, zero):
    subs = set()
    for size in range(1, n + 1):
        if n % size: continue
        for S in itertools.combinations(range(n), size):
            if zero not in S: continue
            Ss = set(S)
            if all(op(a, b) in Ss for a in S for b in S):
                subs.add(frozenset(S))
    return subs

def is_ideal(B, D):
    G = B.G; n = G.n
    if not all(B.p(B.p(g, d), B.m(g)) in D for g in range(n) for d in D): return False
    if not all(B.lam[a][d] in D for a in range(n) for d in D): return False
    if not all(B.o(B.o(g, d), B.oi(g)) in D for g in range(n) for d in D): return False
    return True

def trivial_brace(B, D):
    return all(B.o(x, y) == B.p(x, y) for x in D for y in D) and \
           all(B.p(x, y) == B.p(y, x) for x in D for y in D)

def complements(B, D, kind):
    """kind in {'add','mult','sb'}: subsets K with 0 in K, |K|=|G|/|D|, K cap D = 0,
    closed under the relevant operation(s)."""
    n = B.G.n; k = n // len(D); out = []
    others = [x for x in range(n) if x not in D]
    for S in itertools.combinations(others, k - 1):
        K = set(S) | {B.G.zero}
        if len({coset_key(B, D, x) for x in K}) != k: continue
        okA = all(B.p(a, b) in K for a in K for b in K)
        okM = all(B.o(a, b) in K for a in K for b in K)
        if (kind == 'add' and okA) or (kind == 'mult' and okM) or (kind == 'sb' and okA and okM):
            out.append(K)
    return out

def coset_key(B, D, x):
    return frozenset(B.p(x, d) for d in D)

def socle(B):
    n = B.G.n
    return {a for a in range(n) if B.lam[a] == tuple(range(n))
            and all(B.p(a, b) == B.p(b, a) for b in range(n))}

# ------------------------------------------------------ cohomology machinery
class Extension:
    def __init__(self, B, D):
        self.B = B; self.D = sorted(D); Dset = set(D); n = B.G.n
        cos = []
        for x in range(n):
            key = coset_key(B, Dset, x)
            if key not in cos: cos.append(key)
        # put zero coset first
        cos.sort(key=lambda c: (B.G.zero not in c,))
        self.cos = cos; self.H = list(range(len(cos)))
        self.cid = {x: next(i for i, c in enumerate(cos) if x in c) for x in range(n)}
        rep = [min(c) if B.G.zero not in c else B.G.zero for c in cos]
        self.s = rep
        Hn = len(cos)
        self.hadd = [[self.cid[B.p(rep[h], rep[k])] for k in range(Hn)] for h in range(Hn)]
        self.hcirc = [[self.cid[B.o(rep[h], rep[k])] for k in range(Hn)] for h in range(Hn)]
        self.hneg = [self.cid[B.m(rep[h])] for h in range(Hn)]
        self.hlam = [[self.cid[B.lam[rep[h]][rep[k]]] for k in range(Hn)] for h in range(Hn)]
        s = self.s; P, M, O, OI = B.p, B.m, B.o, B.oi
        self.nu = [{x: B.lam[s[h]][x] for x in self.D} for h in range(Hn)]
        self.nuinv = [{v: k for k, v in self.nu[h].items()} for h in range(Hn)]
        self.mu = [{x: P(P(M(s[h]), x), s[h]) for x in self.D} for h in range(Hn)]
        self.sigma = [{x: O(O(OI(s[h]), x), s[h]) for x in self.D} for h in range(Hn)]
        self.beta, self.tbar = self.cocycles(s)

    def cocycles(self, s):
        B = self.B; P, M, O, OI = B.p, B.m, B.o, B.oi; Hn = len(self.H)
        beta = {(h, k): P(P(M(s[self.hadd[h][k]]), s[h]), s[k]) for h in range(Hn) for k in range(Hn)}
        tbar = {(h, k): O(O(OI(s[self.hcirc[h][k]]), s[h]), s[k]) for h in range(Hn) for k in range(Hn)}
        for v in list(beta.values()) + list(tbar.values()):
            assert v in self.D
        return beta, tbar

    # D arithmetic (abelian, + = o on D)
    def dp(self, *xs):
        r = self.B.G.zero
        for x in xs: r = self.B.p(r, x)
        return r
    def dm(self, x): return self.B.m(x)

    def twisted(self, beta, tbar):
        """Structure on H x D (additive coordinates). Return (add, circ) tables or None."""
        Hn = len(self.H); D = self.D
        els = [(h, x) for h in range(Hn) for x in D]; idx = {e: i for i, e in enumerate(els)}
        add = [[None] * len(els) for _ in els]; circ = [[None] * len(els) for _ in els]
        for (h, x) in els:
            for (k, y) in els:
                add[idx[(h, x)]][idx[(k, y)]] = idx[(self.hadd[h][k],
                                                    self.dp(beta[(h, k)], self.mu[k][x], y))]
                xt = self.nuinv[h][x]; yt = self.nuinv[k][y]
                hk = self.hcirc[h][k]
                zt = self.dp(tbar[(h, k)], self.sigma[k][xt], yt)
                circ[idx[(h, x)]][idx[(k, y)]] = idx[(hk, self.nu[hk][zt])]
        return els, idx, add, circ

    def is_cocycle(self, beta, tbar):
        els, idx, add, circ = self.twisted(beta, tbar); N = len(els)
        z = idx[(0, self.B.G.zero)]
        for T in (add, circ):
            if any(T[z][a] != a or T[a][z] != a for a in range(N)): return False
            if any(sorted(T[a]) != list(range(N)) for a in range(N)): return False
            for a in range(N):
                for b in range(N):
                    for c in range(N):
                        if T[T[a][b]][c] != T[a][T[b][c]]: return False
        neg = [next(b for b in range(N) if add[a][b] == z) for a in range(N)]
        for a in range(N):
            for b in range(N):
                for c in range(N):
                    if circ[a][add[b][c]] != add[add[circ[a][b]][neg[a]]][circ[a][c]]: return False
        # same actions for the section h -> (h,0)
        cinv = [next(b for b in range(N) if circ[a][b] == z) for a in range(N)]
        for h in self.H:
            sh = idx[(h, self.B.G.zero)]
            for x in self.D:
                X = idx[(0, x)]
                if add[add[neg[sh]][X]][sh] != idx[(0, self.mu[h][x])]: return False
                if add[neg[sh]][circ[sh][X]] != idx[(0, self.nu[h][x])]: return False
                if circ[circ[cinv[sh]][X]][sh] != idx[(0, self.sigma[h][x])]: return False
        return True

    def thetas(self):
        Hn = len(self.H)
        for vals in itertools.product(self.D, repeat=Hn - 1):
            yield (self.B.G.zero,) + vals

    def shift(self, theta):
        s2 = [self.B.p(self.s[h], theta[h]) for h in self.H]
        b2, t2 = self.cocycles(s2)
        db = {k: self.dp(b2[k], self.dm(self.beta[k])) for k in b2}
        dt = {k: self.dp(t2[k], self.dm(self.tbar[k])) for k in t2}
        return db, dt

    def lemma2_formula(self, theta):
        return {(h1, h2): self.dp(self.mu[h2][theta[h1]], self.dm(theta[self.hadd[h1][h2]]), theta[h2])
                for h1 in self.H for h2 in self.H}

    def ry_tbar_formula(self, theta):
        th1 = [self.nuinv[h][theta[h]] for h in self.H]
        return {(h1, h2): self.dp(th1[h2], self.dm(th1[self.hcirc[h1][h2]]), self.sigma[h2][th1[h1]])
                for h1 in self.H for h2 in self.H}

    def eq310(self, beta, tbar):
        """MacBeth's transcription of RY (3.10), with tau = nu_{h1 o h2}(taubar)."""
        H = self.H; ha, hc, hl, hn = self.hadd, self.hcirc, self.hlam, self.hneg
        tau = {(a, b): self.nu[hc[a][b]][tbar[(a, b)]] for a in H for b in H}
        for h1 in H:
            for h2 in H:
                for h3 in H:
                    lhs = self.dp(self.mu[hl[h1][h3]][tau[(h1, h2)]], self.dm(tau[(h1, ha[h2][h3])]),
                                  tau[(h1, h3)])
                    rhs = self.dp(self.nu[h1][beta[(h2, h3)]], self.mu[hc[h1][h3]][beta[(h1, hn[h1])]],
                                  self.dm(beta[(hn[h1], hc[h1][h3])]),
                                  self.dm(beta[(hc[h1][h2], hl[h1][h3])]))
                    if lhs != rhs: return False
        return True

def normalised_cochains(E):
    H = E.H; keys = [(h, k) for h in H for k in H if h != 0 and k != 0]
    for vals in itertools.product(E.D, repeat=len(keys)):
        c = {(h, k): E.B.G.zero for h in H for k in H}
        c.update(dict(zip(keys, vals)))
        yield c

def group_cocycle(E, c, op, act):
    H = E.H
    # normalised inhomogeneous condition for right action, matching beta(h,k)=-s(h+k)+s(h)+s(k):
    # c(h1+h2,h3) + act_{h3}(c(h1,h2)) = c(h1,h2+h3) + c(h2,h3)
    for a in H:
        for b in H:
            for d in H:
                l = E.dp(c[(op[a][b], d)], act[d][c[(a, b)]])
                r = E.dp(c[(a, op[b][d])], c[(b, d)])
                if l != r: return False
    return True

def analyse(E, verbose=True):
    """Full H^2_Sb computation for one trivial-brace-kernel extension."""
    key = lambda b, t: (tuple(sorted(b.items())), tuple(sorted(t.items())))
    Zb = [c for c in normalised_cochains(E) if group_cocycle(E, c, E.hadd, E.mu)]
    Zt = [c for c in normalised_cochains(E) if group_cocycle(E, c, E.hcirc, E.sigma)]
    assert any(c == E.beta for c in Zb) and any(c == E.tbar for c in Zt)
    Z = [(b, t) for b in Zb for t in Zt if E.is_cocycle(b, t)]
    Z310 = [(b, t) for b in Zb for t in Zt if E.eq310(b, t)]
    Bsb, Bp, Bo = set(), set(), set()
    lemma_ok = ry_ok = True
    for th in E.thetas():
        db, dt = E.shift(th)
        lemma_ok &= (db == E.lemma2_formula(th))
        ry_ok &= (dt == E.ry_tbar_formula(th))
        Bsb.add(key(db, dt)); Bp.add(tuple(sorted(db.items()))); Bo.add(tuple(sorted(dt.items())))
    Zkeys = {key(b, t) for b, t in Z}
    # closure: Z is a group containing B
    for kb in Bsb:
        db = dict(kb[0]); dt = dict(kb[1])
        b2 = {k: E.dp(E.beta[k], db[k]) for k in db}; t2 = {k: E.dp(E.tbar[k], dt[k]) for k in dt}
        assert key(b2, t2) in Zkeys
    H2 = len(Z) // len(Bsb)
    assert len(Z) % len(Bsb) == 0
    kerPhi = [(b, t) for b, t in Z if tuple(sorted(b.items())) in Bp and tuple(sorted(t.items())) in Bo]
    eq310_all = all(E.eq310(b, t) for b, t in Z)
    omega_in_B = key(E.beta, E.tbar) in Bsb
    beta_in_Bp = tuple(sorted(E.beta.items())) in Bp
    tbar_in_Bo = tuple(sorted(E.tbar.items())) in Bo
    return dict(Z=len(Z), Z_via_310=len(Z310), B=len(Bsb), H2=H2, kerPhi=len(kerPhi) // len(Bsb),
                lemma2=lemma_ok, ry_tbar=ry_ok, eq310=eq310_all,
                Omega_zero=omega_in_B, beta_cob=beta_in_Bp, tbar_cob=tbar_in_Bo)

# ------------------------------------------------------------------ main
def z2z4():
    E = [(a, b) for a in range(2) for b in range(4)]
    return Group(E, lambda x, y: ((x[0] + y[0]) % 2, (x[1] + y[1]) % 4), "Z2xZ4")

def d4():
    # (r,f) = rotation^r * flip^f, r mod 4
    E = [(r, f) for r in range(4) for f in range(2)]
    def op(x, y):
        r1, f1 = x; r2, f2 = y
        return (((r1 + (r2 if f1 == 0 else -r2)) % 4), (f1 + f2) % 2)
    return Group(E, op, "D4")

def run(G, report_all=False):
    print(f"=== additive group {G.name} ===")
    lams = skew_braces(G)
    print(f"labelled skew braces: {len(lams)}")
    subs = subgroups(G.n, G.add, G.zero)
    stats = {}; lemma_fail = 0; tested = 0; examples = []
    for bi, lam in enumerate(lams):
        B = SkewBrace(G, lam)
        for D in subs:
            if len(D) in (1, G.n): continue
            if not is_ideal(B, D) or not trivial_brace(B, D): continue
            trip = tuple(bool(complements(B, set(D), k)) for k in ('mult', 'add', 'sb'))
            stats[(len(D), trip)] = stats.get((len(D), trip), 0) + 1
            E = Extension(B, D)
            for th in E.thetas():
                tested += 1
                db, dt = E.shift(th)
                if db != E.lemma2_formula(th) or dt != E.ry_tbar_formula(th): lemma_fail += 1
            if trip == (True, True, False):
                examples.append((bi, B, D, E))
    print("counts of (|D|, (mult,add,sb)) over labelled (brace, trivial-brace ideal) pairs:")
    for k in sorted(stats): print("   ", k, stats[k])
    print(f"Lemma 2 / RY taubar-shift tested on {tested} (instance, theta) pairs; failures: {lemma_fail}")
    return lams, examples

def lambda_split(B, D):
    """Does Lambda_D = (D,+) x| (D,o) have a complement in Lambda_G = (G,+) x| (G,o)?
    Lambda product: (a,b)(c,d) = (a + lambda_b(c), b o d)."""
    G = B.G; n = G.n
    els = [(a, b) for a in range(n) for b in range(n)]
    mul = lambda x, y: (B.p(x[0], B.lam[x[1]][y[0]]), B.o(x[1], y[1]))
    e = (G.zero, G.zero)
    LD = {(a, b) for a in D for b in D}
    k = (n // len(D)) ** 2
    def gen(S):
        S = set(S) | {e}
        while True:
            new = S | {mul(x, y) for x in S for y in S}
            if new == S: return S
            S = new
    cand = [x for x in els if x not in LD]
    seen = set()
    # complements of order k (k = 4 here): generated by <= 2 elements
    for x in cand:
        for y in [None] + cand:
            K = frozenset(gen([x] if y is None else [x, y]))
            if K in seen: continue
            seen.add(K)
            if len(K) == k and not (K & LD) - {e}:
                return True
    return False

if __name__ == "__main__":
    G = z2z4()
    lams, ex = run(G)
    print(f"(Y,Y,N) instances with trivial-brace ideal: {len(ex)}")
    shown = set()
    for bi, B, D, E in ex:
        res = analyse(E)
        soc = socle(B)
        tag = (len(D), res['H2'], res['kerPhi'], res['Omega_zero'])
        print(f"  brace#{bi} D={[G.E[d] for d in sorted(D)]} |Soc|={len(soc)} D<=Soc:{set(D) <= soc} "
              f"-> {res}")
        if bi not in shown and len(shown) < 1:
            shown.add(bi)
            print("    explicit circle on Z2xZ4 for this brace:")
            for a in range(G.n):
                print("     ", G.E[a], "o", [G.E[B.o(a, b)] for b in range(G.n)])
    # RY Example 3.4
    print("\n=== RY Example 3.4 (SmallSkewbrace(8,8)) ===")
    def ry_circ(x, y):
        x1, y1 = x; x2, y2 = y
        return ((x1 + x2) % 2, (y1 + y2 + 2 * y1 * x2 + 2 * (x1 + y1) * y2) % 4)
    lam = []
    for a in range(G.n):
        # lambda_a(b) = -a + a o b
        row = []
        for b in range(G.n):
            ab = G.idx[ry_circ(G.E[a], G.E[b])]
            row.append(G.mul[G.neg(a)][ab])
        lam.append(tuple(row))
    B = SkewBrace(G, lam)
    I = {G.idx[(x, y)] for x in range(2) for y in (0, 2)}
    print("I ideal:", is_ideal(B, I), " trivial brace:", trivial_brace(B, I),
          " |Soc|:", len(socle(B)), " Soc:", [G.E[a] for a in sorted(socle(B))])
    print("split triple (mult, add, sb):", tuple(bool(complements(B, I, k)) for k in ('mult', 'add', 'sb')))
    E = Extension(B, I); print("H^2_Sb analysis:", analyse(E))

    print("Lambda-extension splits (RY claim: yes):", lambda_split(B, I))
    print("\nLambda-splitting for the (Y,Y,N) instances on Z2xZ4:")
    for bi, B1, D1, E1 in ex:
        print(f"  brace#{bi} D<=Soc:{set(D1) <= socle(B1)} Lambda-split:{lambda_split(B1, set(D1))}")

    # non-abelian additive group: Lemma 2 in a genuinely skew setting
    print()
    G4 = d4()
    lams4, ex4 = run(G4)
    for bi, B1, D1, E1 in ex4:
        print(f"  D4 brace#{bi} D={[G4.E[d] for d in sorted(D1)]} D<=Soc:{set(D1) <= socle(B1)} -> {analyse(E1)}")
    # (3.10) route vs intrinsic Z^2 on EVERY trivial-brace-ideal instance
    print("\nIntrinsic Z^2 vs (group cocycles + transcribed (3.10)), all instances:")
    for GG, LL in ((G, lams), (G4, lams4)):
        subs = subgroups(GG.n, GG.add, GG.zero); agree = tot = 0; bad = []
        for bi, lam in enumerate(LL):
            BB = SkewBrace(GG, lam)
            for D in subs:
                if len(D) in (1, GG.n) or not is_ideal(BB, D) or not trivial_brace(BB, D): continue
                r = analyse(Extension(BB, D)); tot += 1
                if r['Z'] == r['Z_via_310'] and r['eq310']: agree += 1
                else: bad.append((bi, sorted(D), r['Z'], r['Z_via_310']))
        print(f"  {GG.name}: {agree}/{tot} instances agree", bad[:5])

# ---------------------------------------------------------------------------
# Structural description suggested in the review:
#   ker Phi  ~=  C / (Z^1_+ + Z^1_o),   C = { delta : (0, d_o delta) in Z^2_Sb }.
def kerphi_formula(E):
    zero = {k: E.B.G.zero for k in E.beta}
    ths = list(E.thetas())
    C = [d for d in ths if E.is_cocycle(zero, E.ry_tbar_formula(d))]
    Zp = [t for t in ths if all(v == E.B.G.zero for v in E.lemma2_formula(t).values())]
    Zo = [t for t in ths if all(v == E.B.G.zero for v in E.ry_tbar_formula(t).values())]
    S = {tuple(E.dp(a[h], b[h]) for h in E.H) for a in Zp for b in Zo}
    assert S <= {tuple(c) for c in C}
    return len(C) // len(S)

if __name__ == "__main__":
    print("\nker Phi via C/(Z1_+ + Z1_o) vs direct, all trivial-brace-ideal instances:")
    for GG in (z2z4(), d4()):
        LL = skew_braces(GG); subs = subgroups(GG.n, GG.add, GG.zero); ok = tot = 0
        for lam in LL:
            BB = SkewBrace(GG, lam)
            for D in subs:
                if len(D) in (1, GG.n) or not is_ideal(BB, D) or not trivial_brace(BB, D): continue
                EE = Extension(BB, D); tot += 1
                ok += (kerphi_formula(EE) == analyse(EE)['kerPhi'])
        print(f"  {GG.name}: formula agrees in {ok}/{tot} instances")

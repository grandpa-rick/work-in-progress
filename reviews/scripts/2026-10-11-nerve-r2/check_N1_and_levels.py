#!/usr/bin/env python3
"""Round-2 checks: N1 witness and the Remark 8 level count (MacBeth 6089268).

Built from printed definitions (p.3-4):
  * U(p) = sum_a y^{E_a}, E_a = arrows out of a  -> the multiset of out-degrees.
  * tr1 N(C) = reflexive graph (objects, arrows, source, target, identities).
  * tr2 N(C) adds composition d1 : N2 -> N1.
  * 'nerve level' of a class in H^n := n+1 (p.4, l.180-181).
Claims checked:
  (N1a) [1] and B(Z/2) + pt both have U = y^2 + y.
  (N1b) their reflexive graphs are not isomorphic; H^0 (components) already separates them.
  (N1c) brute force: up to isomorphism there are exactly 3 categories with U = y^2 + y.
  (R8)  H^n(-;Z) and H^n(-;F2) of Z/2 and {1,e}, n<=3 (normalized bar complex, Smith normal form).
Negative controls BREAK a claim:
  * Z/2 vs {1,e}: U equal AND tr1 isomorphic, so 'tr1 distinguishes' must FAIL for this pair
    (the tr1 test is not a constant 'different'); tr2 must then separate them.
"""
import itertools, sys
from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_form
from sympy import ZZ

FAIL = []
def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        FAIL.append(name)

class Cat:
    """objs; arrows: name -> (src, tgt); ident: obj -> arrow; comp: (f,g) -> f.g (f then g)."""
    def __init__(self, objs, arrows, ident, comp):
        self.objs, self.arrows, self.ident, self.comp = objs, arrows, ident, comp
    def composable(self):
        return [(f, g) for f in self.arrows for g in self.arrows if self.arrows[f][1] == self.arrows[g][0]]
    def is_category(self):
        A = self.arrows
        for o in self.objs:
            i = self.ident[o]
            if A[i] != (o, o): return False
        for (f, g) in self.composable():
            if (f, g) not in self.comp: return False
            h = self.comp[(f, g)]
            if A[h] != (A[f][0], A[g][1]): return False
        for f in A:
            s, t = A[f]
            if self.comp[(self.ident[s], f)] != f or self.comp[(f, self.ident[t])] != f: return False
        for (f, g) in self.composable():
            for k in A:
                if A[g][1] == A[k][0]:
                    if self.comp[(self.comp[(f, g)], k)] != self.comp[(f, self.comp[(g, k)])]: return False
        return True

def U(C):
    """Out-degree multiset = the polynomial sum_a y^{|E_a|}, as a sorted tuple of exponents."""
    return tuple(sorted(sum(1 for f in C.arrows if C.arrows[f][0] == o) for o in C.objs))

def U_poly(C):
    from collections import Counter
    c = Counter(U(C))
    return " + ".join((f"{m}" if m > 1 else "") + (f"y^{d}" if d > 1 else ("y" if d == 1 else "1"))
                      for d, m in sorted(c.items(), reverse=True))

def iso(C, D, level):
    """level 1: reflexive-graph iso; level 2: plus composition."""
    if len(C.objs) != len(D.objs) or len(C.arrows) != len(D.arrows): return False
    Ca, Da = list(C.arrows), list(D.arrows)
    for po in itertools.permutations(D.objs):
        fo = dict(zip(C.objs, po))
        if any(fo[o] is None for o in C.objs): continue
        for pa in itertools.permutations(Da):
            fa = dict(zip(Ca, pa))
            if any(D.arrows[fa[f]] != (fo[C.arrows[f][0]], fo[C.arrows[f][1]]) for f in Ca): continue
            if any(fa[C.ident[o]] != D.ident[fo[o]] for o in C.objs): continue
            if level >= 2 and any(fa[C.comp[(f, g)]] != D.comp[(fa[f], fa[g])] for (f, g) in C.composable()):
                continue
            return True
    return False

def H0_components(C):
    parent = {o: o for o in C.objs}
    def find(o):
        while parent[o] != o: o = parent[o]
        return o
    for f, (s, t) in C.arrows.items():
        parent[find(s)] = find(t)
    return len({find(o) for o in C.objs})

def comp_table(objs, arrows, ident, extra):
    comp = {}
    for f, (s, t) in arrows.items():
        comp[(ident[s], f)] = f; comp[(f, ident[t])] = f
    comp.update(extra)
    return comp

# [1] = 0 -> 1
A = {'i0': (0, 0), 'i1': (1, 1), 'f': (0, 1)}
I = {0: 'i0', 1: 'i1'}
interval = Cat([0, 1], A, I, comp_table([0, 1], A, I, {}))
# B(Z/2) + pt : object a with loop t, t.t = id; object b alone
A = {'ia': ('a', 'a'), 't': ('a', 'a'), 'ib': ('b', 'b')}
I = {'a': 'ia', 'b': 'ib'}
BZ2pt = Cat(['a', 'b'], A, I, comp_table(['a', 'b'], A, I, {('t', 't'): 'ia'}))
# B{1,e} + pt : e.e = e
BEpt = Cat(['a', 'b'], A, I, comp_table(['a', 'b'], A, I, {('t', 't'): 't'}))
# one-object Z/2 and {1,e}
A1 = {'1': ('*', '*'), 'a': ('*', '*')}
I1 = {'*': '1'}
Z2 = Cat(['*'], A1, I1, comp_table(['*'], A1, I1, {('a', 'a'): '1'}))
E1 = Cat(['*'], A1, I1, comp_table(['*'], A1, I1, {('a', 'a'): 'a'}))

for nm, C in [('[1]', interval), ('B(Z/2)+pt', BZ2pt), ('B{1,e}+pt', BEpt), ('Z/2', Z2), ('{1,e}', E1)]:
    check(f"{nm} satisfies the category axioms", C.is_category())
    print(f"     U({nm}) = {U_poly(C)},  components = {H0_components(C)}")

check("(N1a) U([1]) = U(B(Z/2)+pt) = y^2 + y", U(interval) == U(BZ2pt) == (1, 2) and U_poly(interval) == "y^2 + y")
check("(N1b) tr1 N([1]) NOT iso to tr1 N(B(Z/2)+pt)", not iso(interval, BZ2pt, 1))
check("(N1b) H^0 separates them: 1 vs 2 components", H0_components(interval) == 1 and H0_components(BZ2pt) == 2)

# NEGATIVE CONTROL that breaks the claim 'U-equal pairs are tr1-separated':
check("NEG CONTROL: Z/2 vs {1,e}: U equal", U(Z2) == U(E1))
check("NEG CONTROL: Z/2 vs {1,e}: tr1 ISOMORPHIC (so tr1 does not separate this pair)", iso(Z2, E1, 1))
check("NEG CONTROL: Z/2 vs {1,e}: tr2 NOT isomorphic (composition separates)", not iso(Z2, E1, 2))
check("control: B(Z/2)+pt vs B{1,e}+pt tr1 iso, tr2 not", iso(BZ2pt, BEpt, 1) and not iso(BZ2pt, BEpt, 2))
check("control: relabelled [1] is tr2-iso to [1] (iso test accepts isomorphic input)",
      iso(interval, Cat(['p', 'q'], {'x': ('p', 'p'), 'y': ('q', 'q'), 'g': ('p', 'q')}, {'p': 'x', 'q': 'y'},
                        comp_table(['p', 'q'], {'x': ('p', 'p'), 'y': ('q', 'q'), 'g': ('p', 'q')}, {'p': 'x', 'q': 'y'}, {})), 2))

# (N1c) brute force all categories with out-degree family {2,1}: objects a (deg 2), b (deg 1).
found = []
objs = ['a', 'b']
I = {'a': 'ia', 'b': 'ib'}
for tgt in ['a', 'b']:  # the one non-identity arrow g out of a
    arrows = {'ia': ('a', 'a'), 'ib': ('b', 'b'), 'g': ('a', tgt)}
    base = comp_table(objs, arrows, I, {})
    free_pairs = [(f, h) for f in arrows for h in arrows
                  if arrows[f][1] == arrows[h][0] and (f, h) not in base]
    for vals in itertools.product(list(arrows), repeat=len(free_pairs)):
        comp = dict(base); comp.update(dict(zip(free_pairs, vals)))
        C = Cat(objs, arrows, I, comp)
        if C.is_category() and U(C) == (1, 2):
            if not any(iso(C, D, 2) for D in found):
                found.append(C)
print(f"     categories with U = y^2 + y up to iso: {len(found)}")
check("(N1c) exactly 3 categories with U = y^2+y; exactly 2 tr1-classes",
      len(found) == 3 and len({min(j for j in range(3) if iso(found[i], found[j], 1)) for i in range(3)}) == 2)

# (R8) cohomology of one-object monoids via normalized bar complex.
def bar_matrix(mul, nonid, n, ring_mod=None):
    """Matrix of d^n : C^n -> C^{n+1}, trivial coefficients, normalized (nondegenerate simplices)."""
    rows = list(itertools.product(nonid, repeat=n + 1))
    cols = list(itertools.product(nonid, repeat=n))
    idx = {c: j for j, c in enumerate(cols)}
    Mx = [[0] * len(cols) for _ in rows]
    for i, s in enumerate(rows):
        faces = []
        faces.append((1, s[1:]))
        for k in range(n):
            prodk = mul(s[k], s[k + 1])
            faces.append(((-1) ** (k + 1), s[:k] + (prodk,) + s[k + 2:]))
        faces.append(((-1) ** (n + 1), s[:-1]))
        for sign, f in faces:
            if all(x in nonid for x in f):  # degenerate faces vanish in the normalized complex
                Mx[i][idx[f]] += sign
    return Matrix(len(rows), len(cols), lambda i, j: Mx[i][j]) if rows and cols else Matrix.zeros(len(rows), len(cols))

def H_Z(mul, nonid, n):
    """H^n(;Z) = ker d^n / im d^{n-1}; returns (free rank, torsion list)."""
    dn = bar_matrix(mul, nonid, n)
    rank_n = dn.rank()
    dim = len(nonid) ** n
    ker = dim - rank_n
    if n == 0:
        return (ker, [])
    dm = bar_matrix(mul, nonid, n - 1)
    snf = smith_normal_form(dm, domain=ZZ)
    diag = [abs(snf[i, i]) for i in range(min(snf.shape)) if snf[i, i] != 0]
    tors = [d for d in diag if d != 1]
    return (ker - len(diag), tors)

def H_F2(mul, nonid, n):
    def rk2(Mx):
        A = [[int(x) % 2 for x in Mx.row(i)] for i in range(Mx.rows)]
        r = 0; cols = Mx.cols
        for c in range(cols):
            p = next((i for i in range(r, len(A)) if A[i][c]), None)
            if p is None: continue
            A[r], A[p] = A[p], A[r]
            for i in range(len(A)):
                if i != r and A[i][c]:
                    A[i] = [(a + b) % 2 for a, b in zip(A[i], A[r])]
            r += 1
        return r
    dim = len(nonid) ** n
    ker = dim - rk2(bar_matrix(mul, nonid, n))
    im = rk2(bar_matrix(mul, nonid, n - 1)) if n > 0 else 0
    return ker - im

mulZ2 = lambda a, b: '1' if (a, b) == ('a', 'a') else (b if a == '1' else a)
mulE = lambda a, b: 'a' if 'a' in (a, b) else '1'
for nm, mul in [('Z/2', mulZ2), ('{1,e}', mulE)]:
    print(f"     {nm}: H^n(;Z) n=0..3:", [H_Z(mul, ['a'], n) for n in range(4)],
          "  dim H^n(;F2):", [H_F2(mul, ['a'], n) for n in range(4)])
check("(R8) H^1(;Z) = 0 for both; H^2(Z/2;Z) = Z/2; H^2({1,e};Z) = 0",
      H_Z(mulZ2, ['a'], 1) == (0, []) and H_Z(mulE, ['a'], 1) == (0, [])
      and H_Z(mulZ2, ['a'], 2) == (0, [2]) and H_Z(mulE, ['a'], 2) == (0, []))
check("(R8) over F2 the first separating degree is 1", H_F2(mulZ2, ['a'], 1) == 1 and H_F2(mulE, ['a'], 1) == 0)
# By the note's own definition (nerve level of an H^n class = n+1), the separating classes have
# nerve level 2 over F2 but 3 over Z, so 'nerve level is coefficient-free' (Remark 8) is false as
# defined; only the level of the raw datum a.a (on N2) is coefficient-free.
first_Z = next(n for n in range(4) if H_Z(mulZ2, ['a'], n) != H_Z(mulE, ['a'], n))
first_F2 = next(n for n in range(4) if H_F2(mulZ2, ['a'], n) != H_F2(mulE, ['a'], n))
print(f"     first separating degree: Z -> {first_Z} (nerve level {first_Z+1}), F2 -> {first_F2} (nerve level {first_F2+1})")
check("(R8/N3) nerve level (as defined p.4) of the first separating class depends on coefficients",
      first_Z + 1 != first_F2 + 1)

print()
print("ALL PASS" if not FAIL else f"{len(FAIL)} FAILURES: {FAIL}")
sys.exit(1 if FAIL else 0)

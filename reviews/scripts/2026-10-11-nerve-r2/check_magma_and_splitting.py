#!/usr/bin/env python3
"""Round-2 checks on MacBeth's 'Composition is a nerve-level-two datum' (content commit 6089268).

Independent of MacBeth's scripts. Everything here is built from the printed definitions:
  * Theorem 12 (p.9): the unital magma on {1,x,y} with unit 1 and
        x.x = y, x.y = x, y.x = y, y.y = y.
  * Remark 13 / p.5: 'd^2 omega = 0 is associativity'; '[omega] = 0 <=> the associative weld splits'.
    We model a weld as the twisted composition (Eilenberg-MacLane / Baues-Wirsching linear extension)
        (x,u)*(y,v) = (xy, u + x.v + omega(x,y))
    over a base monoid S with S-module A, and test:
      (T1) omega is a 2-cocycle  <=>  the twisted composition is associative      (both directions)
      (T2) for cocycles: [omega] = 0  <=>  the extension splits (a multiplicative section exists)
      (T3) over Theorem 12's magma the bar 'differential' has d o d != 0, so there is no H^2 at all:
           the cocycle condition is NOT 'the associativity equation of Theorem 12'.
      (T4) the magma is not a twisted composition over any base with nontrivial coefficients.
      (T5) Eilenberg-MacLane 'split' = semidirect; an exact factorization with neither factor normal
           (S4 = S3 . C4) is not an extension of either factor by the other.
Negative controls BREAK claims (they must flip a verdict), they do not relabel.
"""
import itertools, sys

FAIL = []
def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        FAIL.append(name)

# ---------------------------------------------------------------- Theorem 12 magma
E = ['1', 'x', 'y']
M = {('x', 'x'): 'y', ('x', 'y'): 'x', ('y', 'x'): 'y', ('y', 'y'): 'y'}
for a in E:
    M[('1', a)] = a; M[(a, '1')] = a

def is_total(T, X):
    return all((a, b) in T and T[(a, b)] in X for a in X for b in X)
def unit_of(T, X):
    return [e for e in X if all(T[(e, a)] == a and T[(a, e)] == a for a in X)]
def assoc_failures(T, X):
    return [(a, b, c) for a in X for b in X for c in X if T[(T[(a, b)], c)] != T[(a, T[(b, c)])]]

check("T12: composition defined on all composable pairs (total table)", is_total(M, E))
check("T12: unital with unit 1", unit_of(M, E) == ['1'])
fails = assoc_failures(M, E)
print("     associativity failures:", fails)
check("T12: (x,x,x) fails: (xx)x = y, x(xx) = x", M[(M[('x','x')],'x')] == 'y' and M[('x',M[('x','x')])] == 'x')
check("T12: non-associative", len(fails) > 0)

# Separation 'composition exists (N2) vs associative (N3)': the magma and an associative unital
# table on the same 3-element set agree on everything that is pair-level-structural (totality, unit),
# and differ only by the N3 predicate. Enumerate all unital tables on {1,x,y} with unit 1.
unital_tables = []
for vals in itertools.product(E, repeat=4):
    T = {('x','x'): vals[0], ('x','y'): vals[1], ('y','x'): vals[2], ('y','y'): vals[3]}
    for a in E:
        T[('1', a)] = a; T[(a, '1')] = a
    unital_tables.append(T)
n_assoc = sum(1 for T in unital_tables if not assoc_failures(T, E))
print(f"     unital tables on {{1,x,y}} with unit 1: {len(unital_tables)}, associative: {n_assoc}")
check("separation is strict: both associative and non-associative unital tables exist",
      0 < n_assoc < len(unital_tables))

# NEGATIVE CONTROL (breaks 'non-associative'): one-entry mutation x.y = x -> x.y = y gives a monoid?
# Search the one-entry mutations of M and require at least one to be associative; the checker must
# then report NO failures for it (otherwise the checker is broken).
mut_ok = []
for k in [('x','x'), ('x','y'), ('y','x'), ('y','y')]:
    for v in E:
        if v == M[k]:
            continue
        T = dict(M); T[k] = v
        if not assoc_failures(T, E):
            mut_ok.append((k, v))
print("     one-entry mutations of the T12 table that ARE associative:", mut_ok)
check("NEG CONTROL: a one-entry mutation is associative (checker can say 'associative')", len(mut_ok) > 0)

# ---------------------------------------------------------------- generic twisted composition
def twisted(S, mul, A_elems, add, act, omega):
    """Return (elements, product) of the twisted composition on S x A."""
    X = [(s, a) for s in S for a in A_elems]
    def prod(p, q):
        (x, u), (y, v) = p, q
        return (mul(x, y), add(add(u, act(x, v)), omega[(x, y)]))
    return X, prod

def is_assoc(X, prod):
    return all(prod(prod(p, q), r) == prod(p, prod(q, r)) for p in X for q in X for r in X)

def d2(S, mul, add, neg, act, omega):
    """(d omega)(x,y,z) = x.omega(y,z) - omega(xy,z) + omega(x,yz) - omega(x,y)."""
    return {(x, y, z): add(add(add(act(x, omega[(y, z)]), neg(omega[(mul(x, y), z)])), omega[(x, mul(y, z))]), neg(omega[(x, y)]))
            for x in S for y in S for z in S}

def is_cocycle(S, mul, A0, add, neg, act, omega):
    return all(v == A0 for v in d2(S, mul, add, neg, act, omega).values())

def coboundary(S, mul, add, neg, act, beta):
    """(d beta)(x,y) = x.beta(y) - beta(xy) + beta(x)."""
    return {(x, y): add(add(act(x, beta[y]), neg(beta[mul(x, y)])), beta[x]) for x in S for y in S}

def normalized_cochains(S, e, A_elems, A0):
    keys = [(x, y) for x in S for y in S if x != e and y != e]
    for vals in itertools.product(A_elems, repeat=len(keys)):
        om = {(x, y): A0 for x in S for y in S}
        om.update(dict(zip(keys, vals)))
        yield om

def is_coboundary(S, e, mul, A_elems, A0, add, neg, act, omega):
    for vals in itertools.product(A_elems, repeat=len(S)):
        beta = dict(zip(S, vals))
        if coboundary(S, mul, add, neg, act, beta) == omega:
            return True
    return False

def splits(S, e, mul, A_elems, A0, add, act, omega):
    """Exists a section s(x) = (x, b(x)) with s(e)=(e,0) and s(xy) = s(x)*s(y)."""
    X, prod = twisted(S, mul, A_elems, add, act, omega)
    others = [s for s in S if s != e]
    for vals in itertools.product(A_elems, repeat=len(others)):
        b = dict(zip(others, vals)); b[e] = A0
        if all(prod((x, b[x]), (y, b[y])) == (mul(x, y), b[mul(x, y)]) for x in S for y in S):
            return True
    return False

def run_case(label, S, e, mul, A_elems, A0, add, neg, act):
    n = n_coc = n_cob = 0
    t1 = t2 = True
    for om in normalized_cochains(S, e, A_elems, A0):
        n += 1
        X, prod = twisted(S, mul, A_elems, add, act, om)
        coc = is_cocycle(S, mul, A0, add, neg, act, om)
        asc = is_assoc(X, prod)
        if coc != asc:
            t1 = False
        if coc:
            n_coc += 1
            cob = is_coboundary(S, e, mul, A_elems, A0, add, neg, act, om)
            spl = splits(S, e, mul, A_elems, A0, add, act, om)
            n_cob += cob
            if cob != spl:
                t2 = False
    print(f"     {label}: {n} normalized 2-cochains, {n_coc} cocycles, {n_cob} coboundaries")
    check(f"(T1) {label}: cocycle <=> twisted composition associative (all cochains)", t1)
    check(f"(T2) {label}: for cocycles, [omega]=0 <=> extension splits (both directions)", t2)
    return n, n_coc, n_cob

Zn = lambda m: (list(range(m)), 0, (lambda a, b: (a + b) % m), (lambda a: (-a) % m))

# Case 1: S = Z/2, A = Z/2 trivial action.  H^2 = Z/2 (Z/4 is the non-split weld).
S2, e2, m2 = [0, 1], 0, (lambda a, b: (a + b) % 2)
A, A0, add, neg = Zn(2)
run_case("S=Z/2, A=Z/2 trivial", S2, e2, m2, A, A0, add, neg, lambda x, v: v)
om = {(0,0):0,(0,1):0,(1,0):0,(1,1):1}
X, prod = twisted(S2, m2, A, add, lambda x, v: v, om)
def order(p, X, prod, unit):
    k, q = 1, p
    while q != unit:
        q = prod(q, p); k += 1
    return k
check("     omega(1,1)=1 gives an element of order 4 (Z/4): associative, non-split, [omega]!=0",
      is_assoc(X, prod) and order((1, 0), X, prod, (0, 0)) == 4)

# Case 2: S = Z/3, A = Z/2 trivial: H^2 = 0, and NON-cocycles exist -> non-associative welds.
S3, e3, m3 = [0, 1, 2], 0, (lambda a, b: (a + b) % 3)
n, nc, nb = run_case("S=Z/3, A=Z/2 trivial", S3, e3, m3, A, A0, add, neg, lambda x, v: v)
check("     Z/3 case has non-cocycles (so 'no class' really occurs) and H^2=0", nc < n and nc == nb)

# Case 3: S = Z/2 x Z/2 (Klein), A = Z/2 trivial: H^2 = (Z/2)^3.
K4 = [(0,0),(0,1),(1,0),(1,1)]
mK = lambda a, b: ((a[0]+b[0]) % 2, (a[1]+b[1]) % 2)
n, nc, nb = run_case("S=Klein, A=Z/2 trivial", K4, (0,0), mK, A, A0, add, neg, lambda x, v: v)
check("     Klein: #cocycles/#coboundaries = 8 = |H^2(K4;F2)|", nc // nb == 8)

# Case 4: nontrivial action. S = Z/2 acting on A = Z/4 by inversion: H^2 = Z/2 (D8 split, Q8 non-split).
A4, A40, add4, neg4 = Zn(4)
inv = lambda x, v: v if x == 0 else (-v) % 4
n, nc, nb = run_case("S=Z/2, A=Z/4 sign action", S2, e2, m2, A4, A40, add4, neg4, inv)
check("     sign-action case: H^2 = Z/2", nc // nb == 2)
omQ = {(0,0):0,(0,1):0,(1,0):0,(1,1):2}
X, prod = twisted(S2, m2, A4, add4, inv, omQ)
inv_count = sum(1 for p in X if p != (0,0) and prod(p, p) == (0,0))
check("     omega(1,1)=2 gives Q8 (associative; exactly one involution; non-split)",
      is_assoc(X, prod) and inv_count == 1 and not splits(S2, e2, m2, A4, A40, add4, inv, omQ))

# NEGATIVE CONTROL (breaks (T1) if the base is not associative): over T12's magma the twisted
# composition with omega = 0 and A = Z/2 is non-associative even though omega = 0 trivially
# satisfies 'd omega = 0'.  So 'cocycle <=> associative' needs the base to be a category.
MAG = E
mulM = lambda a, b: M[(a, b)]
om0 = {(a, b): 0 for a in MAG for b in MAG}
X, prod = twisted(MAG, mulM, A, add, lambda x, v: v, om0)
coc0 = is_cocycle(MAG, mulM, A0, add, neg, lambda x, v: v, om0)
check("NEG CONTROL (T1 breaks over a non-associative base): omega=0 'cocycle' yet weld non-associative",
      coc0 and not is_assoc(X, prod))

# (T3) over the magma, d o d != 0 on 1-cochains: no cochain complex, no H^2, no class to speak of.
def dd_nonzero(S, mul, A_elems, add, neg):
    for vals in itertools.product(A_elems, repeat=len(S)):
        beta = dict(zip(S, vals))
        db = coboundary(S, mul, add, neg, lambda x, v: v, beta)
        ddb = d2(S, mul, add, neg, lambda x, v: v, db)
        if any(v != 0 for v in ddb.values()):
            return beta, {k: v for k, v in ddb.items() if v != 0}
    return None
A3, A30, add3, neg3 = Zn(3)
res = dd_nonzero(MAG, mulM, A3, add3, neg3)
print("     magma, Z/3 coeffs: beta with d(d beta) != 0:", res)
check("(T3) bar 'complex' over T12 magma has d o d != 0 (no H^2 exists over it)", res is not None)
check("(T3) control: over Z/3 group, d o d = 0 for every 1-cochain", dd_nonzero(S3, m3, A3, add3, neg3) is None)

# (T4) the magma is not a twisted composition with nontrivial coefficients: |magma| = 3 is prime, so
# base has 1 element and A has 3; every such twisted composition (any omega, any A of order 3) is a
# group isomorphic to Z/3, hence associative and not isomorphic to the magma.
def iso_magmas(X1, T1, X2, T2):
    for perm in itertools.permutations(X2):
        f = dict(zip(X1, perm))
        if all(f[T1(a, b)] == T2(f[a], f[b]) for a in X1 for b in X1):
            return True
    return False
all_assoc = True; any_iso = False
for c in A3:
    X, prod = twisted(['*'], lambda a, b: '*', A3, add3, lambda x, v: v, {('*','*'): c})
    all_assoc &= is_assoc(X, prod)
    any_iso |= iso_magmas(X, prod, MAG, mulM)
check("(T4) every order-3 twisted composition over a 1-element base is associative", all_assoc)
check("(T4) none is isomorphic to the T12 magma (its failure is not a cocycle failure)", not any_iso)

# (T5) Eilenberg-MacLane 'split' = semidirect.  S4 = S3 * C4 is an exact factorization (Zappa-Szep)
# with neither factor normal, hence not an extension of either factor by the other.
def compose(p, q):  # apply p then q (diagrammatic)
    return tuple(q[p[i]] for i in range(4))
S4 = list(itertools.permutations(range(4)))
S3sub = [p for p in S4 if p[3] == 3]
c = (1, 2, 3, 0)
C4 = [(0,1,2,3)]
while True:
    nxt = compose(C4[-1], c)
    if nxt == (0,1,2,3): break
    C4.append(nxt)
prods = {compose(a, b) for a in S3sub for b in C4}
def inverse(p):
    q = [0]*4
    for i, pi in enumerate(p): q[pi] = i
    return tuple(q)
def normal(H):
    Hs = set(H)
    return all(compose(compose(inverse(g), h), g) in Hs for g in S4 for h in H)
check("(T5) S4 = S3 . C4 exactly (|S3||C4| = 24, every element once)",
      len(prods) == 24 and len(set(S3sub) & set(C4)) == 1)
check("(T5) neither S3 nor C4 is normal in S4 (a ZS product that is no EM extension)",
      not normal(S3sub) and not normal(C4))

print()
print("ALL PASS" if not FAIL else f"{len(FAIL)} FAILURES: {FAIL}")
sys.exit(1 if FAIL else 0)

# Referee checks for MacBeth h2cluster (ddbe7c0) and change-of-base (21aef97).
import itertools, numpy as np

def rank2(M):
    M = (np.array(M, dtype=np.int64) % 2).copy()
    if M.size == 0: return 0
    r = 0; rows, cols = M.shape
    for c in range(cols):
        piv = [i for i in range(r, rows) if M[i, c]]
        if not piv: continue
        M[[r, piv[0]]] = M[[piv[0], r]]
        for i in range(rows):
            if i != r and M[i, c]: M[i] ^= M[r]
        r += 1
        if r == rows: break
    return r

def h2_F2(elems, mul):
    # H^2(G;F2), trivial action, inhomogeneous bar complex
    n = len(elems); idx = {g: i for i, g in enumerate(elems)}
    C1 = list(elems); C2 = list(itertools.product(elems, repeat=2)); C3 = list(itertools.product(elems, repeat=3))
    i2 = {p: k for k, p in enumerate(C2)}
    d1 = np.zeros((len(C2), len(C1)), dtype=np.int64)
    for k, (g, h) in enumerate(C2):
        for x in (g, h, mul(g, h)): d1[k, idx[x]] ^= 1
    d2 = np.zeros((len(C3), len(C2)), dtype=np.int64)
    for k, (g, h, l) in enumerate(C3):
        for p in ((h, l), (mul(g, h), l), (g, mul(h, l)), (g, h)): d2[k, i2[p]] ^= 1
    return len(C2) - rank2(d2) - rank2(d1)

cyc = lambda n: (list(range(n)), lambda a, b: (a + b) % n)
V4 = ([(a, b) for a in range(2) for b in range(2)], lambda x, y: ((x[0]+y[0]) % 2, (x[1]+y[1]) % 2))
# Q8 as pairs (sign, unit) with units 1,i,j,k
tab = {('1','1'):(1,'1'),('1','i'):(1,'i'),('1','j'):(1,'j'),('1','k'):(1,'k'),
       ('i','1'):(1,'i'),('i','i'):(-1,'1'),('i','j'):(1,'k'),('i','k'):(-1,'j'),
       ('j','1'):(1,'j'),('j','i'):(-1,'k'),('j','j'):(-1,'1'),('j','k'):(1,'i'),
       ('k','1'):(1,'k'),('k','i'):(1,'j'),('k','j'):(-1,'i'),('k','k'):(-1,'1')}
Q8 = ([(s, u) for s in (1, -1) for u in '1ijk'], lambda x, y: (x[0]*y[0]*tab[(x[1], y[1])][0], tab[(x[1], y[1])][1]))
for name, G in [('Z/2', cyc(2)), ('Z/4', cyc(4)), ('V4', V4), ('Z/6', cyc(6)), ('Q8', Q8)]:
    print(f"dim H^2({name};F2) = {h2_F2(*G)}")
# Q8 has a unique involution, so no (Z/2)^2 subgroup: label '(Z/2)^2.Z/2' is wrong
print("Q8 involutions:", [g for g in Q8[0] if Q8[1](g, g) == (1, '1') and g != (1, '1')])

# Additive-side check: loop set L (0 in L) as transversal of subgroup B in abelian A.
def sub(A_add, L): return all(A_add(a, b) in L for a in L for b in L)
print("Z/4 two-element sets containing 0:", {tuple(sorted(L)): sub(lambda a,b:(a+b)%4, set(L)) for L in [(0,1),(0,2),(0,3)]})
# Z/6, B={0,3}: extension split (Z/6 = Z/2 x Z/3) but transversal {0,1,2} is not a subgroup
print("Z/6: {0,1,2} subgroup?", sub(lambda a,b:(a+b)%6, {0,1,2}), " {0,2,4} subgroup?", sub(lambda a,b:(a+b)%6, {0,2,4}))
# normalized cocycle of section s:Z/3->Z/6, s(k)=k, valued in B={0,3}
f = {(a, b): (a + b - ((a + b) % 3)) % 6 for a in range(3) for b in range(3)}
print("Z/6 section {0,1,2} cocycle:", f)

# change-of-base: 'support' monad L(X)=1 if X nonempty, L(0)=0.  Kl(L) is the poset 0<1.
# theta_T : T x L(1) -> L(T) surjective for all T (so extension (C) full), injective only for |T|<=1.
for T in range(4):
    print(f"support monad |T|={T}: |T x L1|={T*1}, |L T|={1 if T else 0}")

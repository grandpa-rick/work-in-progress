"""Part D: mirror-(LA) / two-sided congruence as the unifying sufficient condition; iso classes at order 3; S4 non-normal cases."""
import itertools, importlib.util, sys, io, contextlib
spec = importlib.util.spec_from_file_location("c1", "2026-10-22-macbeth-closure-criterion-check.py")
# re-use functions without re-running heavy part: exec only defs we need
src = open("2026-10-22-macbeth-closure-criterion-check.py").read().split("print(\"== Part A")[0]
ns = {}; exec(src, ns)
monoids, normalizers, analyse = ns["monoids"], ns["normalizers"], ns["analyse"]
print("== Part D1: order-3 witnesses: is ker K a two-sided congruence (mirror (LA): K(c K(x)) = K(c x))? ==")
W = []
for T in monoids(3):
    for K in normalizers(T):
        if analyse(T, K) == (False, True):
            mirror = all(K[T[c][K[x]]] == K[T[c][x]] for c in range(3) for x in range(3))
            W.append((tuple(map(tuple, T)), K)); print(" table", T, "K", K, "mirror(LA) holds:", mirror)
# iso classes (relabel 1<->2)
def relabel(T, K, p):
    n = len(T); q = {p[i]: i for i in range(n)}
    T2 = tuple(tuple(p[T[q[i]][q[j]]] for j in range(n)) for i in range(n)); K2 = tuple(p[K[q[i]]] for i in range(n))
    return T2, K2
classes = set()
for T, K in W:
    classes.add(min([(T, K), relabel(T, K, (0, 2, 1))]))
tabs = set(min([T, relabel(T, (0,1,2), (0,2,1))[0]]) for T, K in W)
print(f" {len(W)} labelled witnesses on {len(set(T for T,K in W))} labelled tables; {len(classes)} up to iso, on {len(tabs)} iso classes of monoid")
# mirror (LA) over ALL order<=3 normalizers: does mirror => assoc always (it must: proof)
bad = 0; cnt = 0
for n in (1,2,3):
    for T in monoids(n):
        for K in normalizers(T):
            if all(K[T[c][K[x]]] == K[T[c][x]] for c in range(n) for x in range(n)):
                cnt += 1; bad += (not analyse(T, K)[1])
print(f" mirror(LA) normalizers order<=3: {cnt}, non-assoc among them: {bad} (expect 0)")

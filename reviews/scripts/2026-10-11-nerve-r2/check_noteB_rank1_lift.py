#!/usr/bin/env python3
"""Round-2 check on MacBeth's uniqueness note (6089268): what pins the rank-1 vertical lift?

Theorem 18's closing and Remark 23 say the lift on M = y is 'eps -> eps (x) eps' (W's own lift),
'pinned by eps -> eps1 eps2 and c l = l'.  A general structure has delta(eps) = k eps(x)eps, k in N.
We check, over N, on the rank-1 fibre:
  * flip-compatibility c.l = l holds for EVERY k (c swaps tensor factors; eps(x)eps is fixed);
  * coassociativity holds for EVERY k (both sides k^2 eps^{(x)3});
  * universality v(a;u,w) = (a, u, k w) on N^3 is a bijection iff k = 1.
So only Theorem 13 step (1) (universality) pins k = 1; 'as one exhibits directly' / 'pinned by c l = l'
does not.  NEGATIVE CONTROL that breaks: k = 2 must FAIL universality (and k = 0 too).
Also: footnote 1 (new) -- no bijection K -> K x K for 2 <= |K| <= 4; Cantor pairing beta1 o beta1 != beta1.
"""
import itertools, sys
FAIL = []
def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond: FAIL.append(name)

# Elements of D^{(x)r} over N as dicts {word: coeff}, words over {'1','e'} (e = eps); D = W = N[eps]/eps^2.
def delta(k):
    # structure map delta : D -> D (x) D on basis: 1 -> 1(x)1 (step (0)/S2), eps -> k eps(x)eps
    return {'1': {('1', '1'): 1}, 'e': {('e', 'e'): k}}
def apply_at(vec, pos, dl):
    out = {}
    for w, c in vec.items():
        for img, c2 in dl[w[pos]].items():
            nw = w[:pos] + img + w[pos + 1:]
            out[nw] = out.get(nw, 0) + c * c2
    return {w: c for w, c in out.items() if c}
def flip_ok(k):      # c = id (x) gamma with gamma = swap of the two factors; check c(l(eps)) = l(eps)
    l = apply_at({('e',): 1}, 0, delta(k))
    return {(w[1], w[0]): c for w, c in l.items()} == l
def coassoc_ok(k):   # (delta (x) id) delta (eps) = (id (x) delta) delta (eps)
    l = apply_at({('e',): 1}, 0, delta(k))
    return apply_at(l, 0, delta(k)) == apply_at(l, 1, delta(k))
def universal_ok(k, box=6):
    # v : N^3 -> N^3, (a,u,w) -> (a,u,k*w). Bijective on N^3 iff surjective onto every target in a box
    # and injective; check on the box [0,box)^3 of targets.
    image = {(a, u, k * w) for a in range(box) for u in range(box) for w in range(box * 2)}
    surj = all((a, u, w) in image for a in range(box) for u in range(box) for w in range(box))
    inj = len({(a, u, k * w) for a in range(box) for u in range(box) for w in range(box)}) == box ** 3
    return surj and inj

for k in range(4):
    print(f"     k={k}: flip {flip_ok(k)}, coassoc {coassoc_ok(k)}, universality {universal_ok(k)}")
check("flip c.l = l holds for all k in 0..3 (does not pin k)", all(flip_ok(k) for k in range(4)))
check("coassociativity holds for all k in 0..3 (does not pin k)", all(coassoc_ok(k) for k in range(4)))
check("universality holds for k = 1", universal_ok(1))
check("NEG CONTROL: universality FAILS for k = 2 and k = 0", not universal_ok(2) and not universal_ok(0))

for n in range(2, 5):
    K = list(range(n)); KK = list(itertools.product(K, K))
    check(f"footnote 1: no bijection K -> KxK for |K|={n} (|K| < |K|^2)", len(K) != len(KK))

def cantor_inv(z):
    w = int(((8 * z + 1) ** 0.5 - 1) // 2)
    t = w * (w + 1) // 2
    y = z - t
    return (w - y, y)
b1 = lambda k: cantor_inv(k)[0]
bad = [k for k in range(50) if b1(b1(k)) != b1(k)]
check("Cantor pairing: beta1 o beta1 != beta1 (first failure at small k)", len(bad) > 0)
print("     first k with beta1(beta1(k)) != beta1(k):", bad[:3])
print()
print("ALL PASS" if not FAIL else f"{len(FAIL)} FAILURES: {FAIL}")
sys.exit(1 if FAIL else 0)

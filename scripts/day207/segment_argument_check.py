"""Day 207: test Rick's Newton-segment argument (re-proving Rado's direction) and
Clio's Prop 6.2 (a6c83ed build) by brute force.

Segment argument: for sigma <| lamhat (partitions of d, <= ell parts), run Clio's Lemma 4.1
from nu = lamhat towards target sigma; each step tau = nu - e_a + e_b equals
(1-s) nu + s (a b)nu with s = 1/(nu_a - nu_b). Carry an explicit convex combination of
permutations of lamhat (exact Fractions) and check it lands on sigma, then on every
permutation alpha of sigma.
"""
from fractions import Fraction as Fr
from itertools import permutations, combinations, product
from collections import defaultdict

def partitions(d, maxpart=None, maxlen=None):
    if maxpart is None: maxpart = d
    if d == 0: yield (); return
    if maxlen == 0: return
    for p in range(min(d, maxpart), 0, -1):
        for rest in partitions(d - p, p, None if maxlen is None else maxlen - 1):
            yield (p,) + rest

def pad(p, ell): return tuple(p) + (0,) * (ell - len(p))
def psums(p):
    s, out = 0, []
    for x in p: s += x; out.append(s)
    return out
def dominated(s, n, ell):  # s <| n (weak), both padded to ell
    return all(a <= b for a, b in zip(psums(pad(s, ell)), psums(pad(n, ell))))

def lemma41(nu, sigma):
    """Clio's Lemma 4.1 exactly as written (1-indexed in text, 0-indexed here)."""
    L = len(nu)
    i = min(c for c in range(L) if nu[c] > sigma[c])
    j = min(c for c in range(i + 1, L) if nu[c] < sigma[c])
    a = max(c for c in range(i, j) if nu[c] == nu[i])
    b = min(c for c in range(i + 1, j + 1) if nu[c] == nu[j])
    return i, j, a, b

stats = defaultdict(int)
fails = []
ELL, DMAX = 5, 10
for d in range(1, DMAX + 1):
    for ell in range(1, ELL + 1):
        parts = [pad(p, ell) for p in partitions(d, maxlen=ell)]
        for lam in parts:
            for sig in parts:
                if not dominated(sig, lam, ell): continue
                # convex combination: dict permutation-tuple -> weight
                comb = {lam: Fr(1)}
                nu, steps = lam, 0
                while nu != sig:
                    i, j, a, b = lemma41(nu, sig)
                    gap = nu[a] - nu[b]
                    if gap < 2: fails.append(('gap', lam, sig, nu)); break      # leak (a)
                    tau = list(nu); tau[a] -= 1; tau[b] += 1; tau = tuple(tau)
                    if list(tau) != sorted(tau, reverse=True): fails.append(('notpart', lam, sig, nu)); break
                    if not (dominated(tau, nu, ell) and tau != nu and dominated(sig, tau, ell)):
                        fails.append(('order', lam, sig, nu)); break
                    s = Fr(1, gap)
                    new = defaultdict(Fr)
                    for v, w in comb.items():
                        vv = list(v); vv[a], vv[b] = vv[b], vv[a]
                        new[v] += (1 - s) * w; new[tuple(vv)] += s * w
                    comb = dict(new)
                    point = tuple(sum(w * v[k] for v, w in comb.items()) for k in range(ell))
                    if point != tau: fails.append(('segment', lam, sig, nu)); break
                    nu, steps = tau, steps + 1
                    stats['steps_total'] += 1
                    if steps > 10**4: fails.append(('noterm', lam, sig)); break   # leak (b)
                else:
                    # leak (c): landed exactly on sigma; all support points are perms of lam
                    assert all(sorted(v, reverse=True) == list(lam) for v in comb)
                    assert abs(sum(comb.values()) - 1) == 0 and all(w >= 0 for w in comb.values())
                    # every permutation alpha of sigma: permute the combination
                    for perm in set(permutations(range(ell))):
                        alpha = tuple(sig[perm[k]] for k in range(ell))
                        pt = tuple(sum(w * v[perm[k]] for v, w in comb.items()) for k in range(ell))
                        if pt != alpha: fails.append(('perm', lam, sig, alpha))
                    stats['pairs'] += 1
                    stats['maxsteps'] = max(stats['maxsteps'], steps)
print('segment argument:', dict(stats), 'failures:', len(fails), fails[:3])

# Sorted form vs subset form (first sentence of Prop 6.2 proof) and repaired exchange.
stats2 = defaultdict(int); f2 = []
for d in range(1, 8):
    for ell in range(1, 5):
        for lam in [pad(p, ell) for p in partitions(d, maxlen=ell)]:
            Lam = [0] + psums(lam)
            subsets = [frozenset(S) for r in range(ell + 1) for S in combinations(range(ell), r)]
            box = [a for a in product(range(d + 1), repeat=ell) if sum(a) == d]
            Jsort = {a for a in box if dominated(tuple(sorted(a, reverse=True)), lam, ell)}
            Jsub = {a for a in box if all(sum(a[k] for k in S) <= Lam[len(S)] for S in subsets)}
            if Jsort != Jsub: f2.append(('ident', lam))
            J = Jsub
            for al in J:
                tight = [S for S in subsets if sum(al[k] for k in S) == Lam[len(S)]]
                for be in J:
                    D = {k for k in range(ell) if al[k] < be[k]}
                    for i in range(ell):
                        if al[i] <= be[i]: continue
                        stats2['triples'] += 1
                        # hypothetical S* of the contradiction: must never exist
                        if any(D <= S and i not in S for S in tight): f2.append(('Sstar', al, be, i))
                        ok = [j for j in D if tuple(al[k] - (k == i) + (k == j) for k in range(ell)) in J]
                        if not ok: f2.append(('exch', al, be, i))
print('Prop 6.2:', dict(stats2), 'failures:', len(f2), f2[:3])

"""Day 217: test MacBeth's corrected basepoint criterion (10-02 PDF, Prop. 1):
  basepoint ([beta]=0 realised)  <=>  exists a section s with delta_.(d) in Hom(H,D) for all d in D.
H = Z/2, (G,+) abelian, D any lambda-invariant index-2 subgroup.  Section s(1)=u not in D, b=2u.
delta_1(d) = lambda_d(u) - u ; delta_0 = 0 ; hom  <=>  2 delta_1(d) = 0 for all d.
Per-brace basepoint: [beta]=0 in H^2_+ = D/2D  <=>  b in 2D.
Class-level: pool all (brace, D) with the same triplet (lambda^D, nu, sigmahat) up to
Aut(D,+) relabelling and section change; basepoint(class) = exists member with b in 2D."""
import sys, itertools, os
COARSE = bool(os.environ.get("COARSE"))  # drop sigmahat from the triplet key (coarser classes)
from collections import defaultdict
from braces import AbGroup, regular_subgroups, brace_ops

def index2_subgroups(A):
    subs = set()
    for gens in itertools.combinations(range(A.n), min(len(A.mods), 2)):
        pass
    # brute: subgroups = closures of subsets of size<=2 generators (enough for rank<=2 groups here)
    def close(S):
        S = set(S) | {A.zero}
        while True:
            new = {A.add[x][y] for x in S for y in S} - S
            if not new: return frozenset(S)
            S |= new
    for x in range(A.n):
        for y in range(x, A.n):
            S = close([x, y])
            if len(S) * 2 == A.n: subs.add(S)
    if len(A.mods) >= 3:
        for x, y, z in itertools.combinations(range(A.n), 3):
            S = close([x, y, z])
            if len(S) * 2 == A.n: subs.add(S)
    return subs

def run(mods):
    A = AbGroup(mods); R = regular_subgroups(A)
    out = []
    for bi, lam in enumerate(R):
        circ, inv = brace_ops(A, lam)
        for D in index2_subgroups(A):
            if not all(lam[a][d] in D for a in range(A.n) for d in D): continue
            Dl = sorted(D)
            twoD = {A.add[x][x] for x in D}
            L = {tuple(lam[d][x] for x in Dl) for d in Dl}
            fixL = {x for x in D if all(lam[d][x] == x for d in D)}
            secs = []
            for u in range(A.n):
                if u in D: continue
                b = A.add[u][u]
                good = all(A.add[A.add[lam[d][u]][A.neg[u]]][A.add[lam[d][u]][A.neg[u]]] == A.zero for d in D)
                secs.append((u, b, good))
            based = secs[0][1] in twoD          # b mod 2D is section-independent
            assert all((s[1] in twoD) == based for s in secs)
            crit = any(s[2] for s in secs)
            # lock sanity: good(u) <=> b in Fix(L)
            assert all(s[2] == (s[1] in fixL) for s in secs)
            out.append(dict(brace=bi, D=D, Lsize=len(L), based=based, crit=crit,
                            ngood=sum(s[2] for s in secs), nsec=len(secs),
                            fix_in_2D=fixL <= twoD, inv=triplet_key(A, lam, circ, inv, D)))
    return A, R, out

def inv_factors(A, D):
    """invariant factors of the abelian 2-group D (orders of cyclic factors)"""
    from collections import Counter
    cnt = Counter(A.order(x) for x in D)
    # number of elements of order dividing 2^k determines the type
    n_le = lambda k: sum(c for o, c in cnt.items() if (2**k) % o == 0)
    facs = []
    k = 1
    import math
    r_prev = None
    ranks = []
    while n_le(k-1) < len(D):
        ranks.append(round(math.log2(n_le(k) // n_le(k-1))))
        k += 1
    # ranks[i] = number of cyclic factors of order >= 2^(i+1)
    for i, r in enumerate(ranks):
        nxt = ranks[i+1] if i+1 < len(ranks) else 0
        facs += [2**(i+1)] * (r - nxt)
    return tuple(sorted(facs))

def triplet_key(A, lam, circ, inv, D):
    """canonical form of (lambda^D, nu_1, sigmahat_1) transported to the model D0 = prod Z/f,
    minimised over all additive isos D0 -> D and all sections u (section change)."""
    facs = inv_factors(A, D)
    D0 = AbGroup(facs)
    basis = []
    for j in range(len(facs)):
        e = [0]*len(facs); e[j] = 1; basis.append(D0.idx[tuple(e)])
    Dl = sorted(D)
    sections = [u for u in range(A.n) if u not in D]
    best = None
    cands = [[x for x in Dl if f % A.order(x) == 0] for f in facs]
    for imgs in itertools.product(*cands):
        psi = []
        for e in D0.elts:
            y = A.zero
            for c, im in zip(e, imgs):
                y = A.add[y][A.mul(c, im)]
            psi.append(y)
        if len(set(psi)) != len(D): continue
        pinv = {v: i for i, v in enumerate(psi)}
        lamD = tuple(tuple(pinv[lam[psi[d]][psi[x]]] for x in range(D0.n)) for d in range(D0.n))
        for u in sections:
            nu = tuple(pinv[lam[u][psi[x]]] for x in range(D0.n))
            sh = tuple(pinv[circ[circ[inv[u]][psi[x]]][u]] for x in range(D0.n))
            key = (lamD, nu, sh) if not COARSE else (lamD, nu)
            if best is None or key < best: best = key
    return (facs, best)

if __name__ == "__main__":
    cases = [tuple(map(int, a.split('x'))) for a in sys.argv[1:]]
    for mods in cases:
        A, R, out = run(mods)
        print(f"=== (G,+)=Z{mods}: {len(R)} braces, {len(out)} (brace,D) pairs, H=Z/2")
        by = defaultdict(lambda: defaultdict(int))
        for r in out:
            Dt = tuple(sorted(A.order(x) for x in r['D']))
            by[(Dt, r['Lsize'])][(r['based'], r['crit'])] += 1
        for k, d in sorted(by.items()):
            tot = sum(d.values()); agree = d[(True, True)] + d[(False, False)]
            print(f"  D-orders-profile={k[0][-1]}-max,|D|={len(k[0])} |L|={k[1]}: pairs {tot}; "
                  f"(based,crit) counts {dict(d)}; agree {agree}/{tot}")
        import pickle
        pickle.dump([(mods, {k: v for k, v in r.items() if k != 'D'}) for r in out],
                    open(f"out_{'x'.join(map(str,mods))}.pkl", "wb"))

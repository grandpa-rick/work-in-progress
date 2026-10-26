"""(c) Corrected item-X sentence (CC Sec 4.2, TR Rem 9.1, 6df0043). Rick, Wake 258.
Biconditional B(G,H): for every right transversal T (id in T): (T,#) is a group <=> T is a subgroup.
Claim: B fails for (G,H) <=> some transversal is associative but is not a complement of H.
Complements are computed independently: subgroups C of G with |C| = [G:H] and C cap H = 1.
Also: subgroup T => associative (always); non-split normal H => every T associative, none a complement;
and the OLD sentence ('false in the presence of a nontrivial normal subgroup') is tested. All 24 groups of order <= 12."""
import time
from common import *

c = cyclic
GROUPS = [("C1", c(1)), ("C2", c(2)), ("C3", c(3)), ("C4", c(4)), ("C2^2", direct(c(2), c(2))), ("C5", c(5)),
          ("C6", c(6)), ("S3", sym(3)), ("C7", c(7)), ("C8", c(8)), ("C4xC2", direct(c(4), c(2))),
          ("C2^3", direct(direct(c(2), c(2)), c(2))), ("D4", dihedral(4)), ("Q8", dicyclic(2)), ("C9", c(9)),
          ("C3^2", direct(c(3), c(3))), ("C10", c(10)), ("D5", dihedral(5)), ("C11", c(11)), ("C12", c(12)),
          ("C6xC2", direct(c(6), c(2))), ("A4", alt4()), ("D6", dihedral(6)), ("Dic3", dicyclic(3))]

def is_group_table(G):
    n = len(G); e = ident(G)
    return is_assoc(G, n) and all(any(G[x][y] == e for y in range(n)) for x in range(n))

t0 = time.time()
tot = dict(groups=0, GH_pairs=0, transversals=0, B_fails=0, claim_mismatch=0, subgrp_not_assoc=0, subgrp_vs_compl_mismatch=0,
           nonsplit_normal=0, nonsplit_normal_violations=0, nontriv_normal_proper=0, old_sentence_counterex=0,
           B_fails_H_nonnormal=0, B_fails_H_normal=0)
old_cex = []; fails = []
assert len({len(G) for _, G in GROUPS} | set()) and all(is_group_table(G) for _, G in GROUPS)
for name, G in GROUPS:
    n = len(G); e = ident(G); inv = inverses(G, e); subs = subgroups(G); tot['groups'] += 1
    for H in subs:
        idx = n // len(H)
        compl = {C for C in subs if len(C) == idx and len(C & H) == 1}
        nor = is_normal(G, H, inv)
        B = True; claim_rhs = False; allassoc = True; anycompl = False
        for T in transversals(G, H, e):
            tot['transversals'] += 1; K = rep_map(G, H, T); Ts = frozenset(T)
            sh = lambda x, y: K[G[x][y]]
            assoc = all(sh(sh(f, g), h) == sh(f, sh(g, h)) for f in T for g in T for h in T)
            subgrp = all(G[a][b] in Ts for a in T for b in T)
            iscompl = Ts in compl
            tot['subgrp_vs_compl_mismatch'] += subgrp != iscompl
            tot['subgrp_not_assoc'] += subgrp and not assoc
            if assoc != subgrp: B = False
            if assoc and not iscompl: claim_rhs = True
            allassoc &= assoc; anycompl |= iscompl
        tot['GH_pairs'] += 1
        if not B:
            tot['B_fails'] += 1; tot['B_fails_H_normal' if nor else 'B_fails_H_nonnormal'] += 1; fails.append((name, len(H), nor))
        tot['claim_mismatch'] += (not B) != claim_rhs
        if nor and not compl:          # non-split normal H
            tot['nonsplit_normal'] += 1; tot['nonsplit_normal_violations'] += not (allassoc and not anycompl)
        if nor and 1 < len(H) < n:
            tot['nontriv_normal_proper'] += 1
            if B: tot['old_sentence_counterex'] += 1; old_cex.append((name, len(H)))
print(tot)
print("old-sentence counterexamples (G,|H|) with H proper nontrivial normal but B holds:", sorted(set(old_cex)))
print("B fails at (G,|H|,normal):", sorted(set(fails)))
# D4 sanity for TR's '(12 transversals)'
D4 = dihedral(4); e = ident(D4); inv = inverses(D4, e)
Z = [H for H in subgroups(D4) if len(H) == 2 and is_normal(D4, H, inv)][0]
print("D4, H = centre: transversals with id:", sum(1 for _ in transversals(D4, Z, e)),
      "| all order-2 H, id-transversals total:", sum(sum(1 for _ in transversals(D4, H, e)) for H in subgroups(D4) if len(H) == 2))
print(f"[{time.time()-t0:.1f}s]")
# strictly-between count (assoc, T not a subgroup, H not normal) per group, |G| <= 12
for name, G in GROUPS:
    n = len(G); e = ident(G); inv = inverses(G, e); cnt = 0; cores = set()
    for H in subgroups(G):
        if is_normal(G, H, inv): continue
        for T in transversals(G, H, e):
            K = rep_map(G, H, T); Ts = frozenset(T); sh = lambda x, y: K[G[x][y]]
            if all(sh(sh(f, g), h) == sh(f, sh(g, h)) for f in T for g in T for h in T) and not all(G[a][b] in Ts for a in T for b in T):
                cnt += 1; cores.add((len(H), len(core(G, H, inv))))
    if cnt: print(f"strictly-between in {name}: {cnt} transversals; (|H|,|core|) = {sorted(cores)}")

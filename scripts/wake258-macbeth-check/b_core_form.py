"""(b) mirror-LA Rem 4.3 core form and Thm 4.2 / Prop 5.1 over groups. Rick, Wake 258.
For every subgroup H and every right transversal T (id in T), K = right-coset rep map:
 assoc(T,#)  vs  [Thm 4.2] f phi(g,h) f^-1 in H all f,g,h  vs  phi(TxT) in core_G(H)  vs  image of T in G/core is a subgroup
 vs Prop 5.1(3) with theta_K computed from the printed formula; also theta_K == core-coset relation; (mLA) <=> H normal."""
import time
from common import *

def run(name, G):
    t0 = time.time(); n = len(G); e = ident(G); inv = inverses(G, e)
    st = dict(pairs=0, pairs_proper_nontriv=0, pairs_excl_H1=0, assoc=0, m_thm42=0, m_core=0, m_img=0, m_p51=0,
              m_theta_core=0, m_mLA_normal=0, LA_fail=0, between=0, ctrl_H_in_place_of_core=0)
    for H in subgroups(G):
        C = core(G, H, inv); nor = is_normal(G, H, inv)
        for T in transversals(G, H, e):
            K = rep_map(G, H, T); Ts = set(T)
            st['pairs'] += 1; st['pairs_proper_nontriv'] += 1 < len(H) < n; st['pairs_excl_H1'] += len(H) > 1
            sharp = lambda x, y: K[G[x][y]]
            assoc = all(sharp(sharp(f, g), h) == sharp(f, sharp(g, h)) for f in T for g in T for h in T)
            phi = {(g, h): G[G[g][h]][inv[sharp(g, h)]] for g in T for h in T}
            thm42 = all(G[G[f][p]][inv[f]] in H for p in phi.values() for f in T)
            corec = all(p in C for p in phi.values())
            ctrl = all(p in H for p in phi.values())          # control: H instead of core (always true -> must mismatch)
            img = {frozenset(G[c][t] for c in C) for t in T}   # cosets C t
            prods = {frozenset(G[c][G[a][b]] for c in C) for a in T for b in T}
            img_sub = prods <= img                             # finite: closed under product => subgroup of G/C
            st['LA_fail'] += not is_normalizer(G, n, e, K)
            th = theta_explicit(G, n, K)
            st['m_theta_core'] += th != [[G[x][inv[y]] in C for y in range(n)] for x in range(n)]
            p51 = p51_three(G, n, K, th)
            st['assoc'] += assoc; st['m_thm42'] += assoc != thm42; st['m_core'] += assoc != corec
            st['m_img'] += assoc != img_sub; st['m_p51'] += assoc != p51
            st['m_mLA_normal'] += mLA_full(G, n, K) != nor
            subgrp = all(G[a][b] in Ts for a in T for b in T)
            st['between'] += assoc and not subgrp and not nor
            st['ctrl_H_in_place_of_core'] += assoc != ctrl
    print(f"{name} (|G|={n}): {st}  [{time.time()-t0:.1f}s]", flush=True)

S3 = sym(3); run("S3", S3)
run("S3xC2", direct(S3, cyclic(2)))
run("S4", sym(4))

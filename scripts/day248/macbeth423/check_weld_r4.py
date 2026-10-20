#!/usr/bin/env python3
"""Rick round-4 check of MacBeth UID 423 (content commit 40f5c95).
New claims in the diff a452538..40f5c95:
  (P1) closed forms of lambda_A, rho_A (+ back maps) for pi_!A, all S,T;
  (P1') lambda_B(s,const_t) = (s, a -> (a,const_t)); rho_B fixes the whole position;
  (T7) e agrees with r_s0 on the fibre over s0 (positions (s0,g)), not elsewhere.
Reuses the round-2 engine unchanged. Negative controls: a wrong back map must fail."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'day245', 'macbeth412'))
from check_weld_worker import *

def r4(nS, nT, out):
    S, T = list(range(nS)), [chr(97+k) for k in range(nT)]
    StS, StT = St(S, 'StS'), St(T, 'StT'); dS, dT = delta_St(S), delta_St(T)
    B = comp(StS, StT); PB = B.positions()
    const = {t: fz({a: t for a in S}) for t in T}
    gof = lambda p: th(p[1])
    lamB = assoc_LtoR(whisk_right(dS, StS, StT), StS, StS, StT)
    rhoB = assoc_RtoL(whisk_left(StS, dT, StT, comp(StT, StT)), StS, StT, StT)
    SxT = [(s, t) for s in S for t in T]; A = St(SxT, 'A'); dA = delta_St(SxT)
    piS = lambda st: (st[0], {a: (a, st[1]) for a in S})
    piT = lambda st: (st[1], {t: (st[0], t) for t in T})
    lamA = compose(whisk_right(piS, A, A), dA); rhoA = compose(whisk_left(A, piT, A, StT), dA)
    LA = bicomodule_laws(A, StS, StT, lamA, rhoA, SxT, eps_St, dS, eps_St, dT); assert all(LA.values())
    # (P1) printed closed forms
    lamA_cf = lambda st: ((st[0], fz({a: (a, st[1]) for a in S})),
                          {(a, (a2, t2)): (a2, t2) for a in S for (a2, t2) in SxT})
    rhoA_cf = lambda st: ((st, fz({(a, t2): t2 for (a, t2) in SxT})),
                          {((a, t2), t3): (a, t3) for (a, t2) in SxT for t3 in T})
    ok_l = all(lamA(x) == lamA_cf(x) for x in SxT); ok_r = all(rhoA(x) == rhoA_cf(x) for x in SxT)
    # negative controls: wrong back maps must FAIL (when there is room to fail)
    bad_l = lambda st: (lamA_cf(st)[0], {(a, (a2, t2)): (a, t2) for a in S for (a2, t2) in SxT})
    bad_r = lambda st: (rhoA_cf(st)[0], {((a, t2), t3): (a, t2) for (a, t2) in SxT for t3 in T})
    nc_l = any(lamA(x) != bad_l(x) for x in SxT); nc_r = any(rhoA(x) != bad_r(x) for x in SxT)
    # (P1') B side at g = const_t
    ok_lB = all(lamB((s, const[t]))[0] == (s, fz({a: (a, const[t]) for a in S})) for s in S for t in T)
    ok_rB = all(rhoB((s, const[t]))[0][0] == (s, const[t]) for s in S for t in T)
    sig = lambda st: ((st[0], const[st[1]]), {(a, t): (a, t) for a in S for t in T})
    okl = same(compose(whisk_left(StS, sig, A, B), lamA), compose(lamB, sig), SxT)[0]
    okr = same(compose(whisk_right(sig, A, StT), rhoA), compose(rhoB, sig), SxT)[0]
    # (T7) e vs r_s0 on fibre over s0
    e = lambda p: ((p[0], const[gof(p)[p[0]]]), {d: d for d in B.dirs(p)})
    t7 = []
    for s0 in S:
        r = lambda p, s0=s0: ((p[0], const[gof(p)[s0]]), {d: d for d in B.dirs(p)})
        on = all(e(p) == r(p) for p in PB if p[0] == s0)
        off = sum(1 for p in PB if p[0] != s0 and e(p) != r(p))
        t7.append((on, off))
    T7on = all(o for o, _ in t7); T7off = [f for _, f in t7]
    out(f'({nS},{nT}): P1 lamA cf {ok_l} rhoA cf {ok_r} | neg-ctrl wrong-back fails: lam {nc_l} rho {nc_r} | '
        f'lamB/rhoB at const {ok_lB}/{ok_rB} | sigma squares {okl}/{okr} | T7 e=r_s0 on fibre s0 {T7on}, '
        f'disagree off-fibre {T7off[0]}')
    assert ok_l and ok_r and ok_lB and ok_rB and okl and okr and T7on
    if nS >= 2 and nT >= 2: assert all(f > 0 for f in T7off)
    if nS >= 1 and len(SxT) >= 2: assert nc_l or nS == 1
    return nc_l, nc_r

if __name__ == '__main__':
    lines = []
    def out(x): print(x); lines.append(x)
    shapes = [(n, m) for n in (1, 2, 3) for m in (1, 2, 3)] + [(1, 4), (4, 1), (3, 4), (4, 3), (2, 5)]
    ncs = [r4(*sh, out) for sh in shapes]
    assert any(a for a, _ in ncs) and any(b for _, b in ncs)
    out('ALL ROUND-4 ASSERTIONS PASSED')

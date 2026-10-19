#!/usr/bin/env python3
"""Rick round-3 check of MacBeth UID 418 (content commit a452538) at GENERIC shapes.
Reuses the round-2 engine (day245/macbeth412/check_weld_worker.py) unchanged."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'day245', 'macbeth412'))
from check_weld_worker import *

def r3(nS, nT, out):
    S, T = list(range(nS)), [chr(97+k) for k in range(nT)]
    StS, StT = St(S, 'StS'), St(T, 'StT'); dS, dT = delta_St(S), delta_St(T)
    B = comp(StS, StT); PB = B.positions()
    const = {t: fz({a: t for a in S}) for t in T}
    gof = lambda p: th(p[1])
    lamB = assoc_LtoR(whisk_right(dS, StS, StT), StS, StS, StT)
    rhoB = assoc_RtoL(whisk_left(StS, dT, StT, comp(StT, StT)), StS, StT, StT)
    out(f'\n=== ({nS},{nT}) |B(1)|={len(PB)}')
    L = bicomodule_laws(B, StS, StT, lamB, rhoB, PB, eps_St, dS, eps_St, dT); assert all(L.values())
    Bp = sorted({(s, const[t]) for s in S for t in T})
    # Prop 1 / Remark 3: pi_! A = B' coaction-by-coaction, all positions, both squares
    SxT = [(s, t) for s in S for t in T]; A = St(SxT, 'A'); dA = delta_St(SxT)
    piS = lambda st: (st[0], {a: (a, st[1]) for a in S})
    piT = lambda st: (st[1], {t: (st[0], t) for t in T})
    lamA = compose(whisk_right(piS, A, A), dA); rhoA = compose(whisk_left(A, piT, A, StT), dA)
    LA = bicomodule_laws(A, StS, StT, lamA, rhoA, SxT, eps_St, dS, eps_St, dT); assert all(LA.values())
    sig = lambda st: ((st[0], const[st[1]]), {(a, t): (a, t) for a in S for t in T})
    okl = same(compose(whisk_left(StS, sig, A, B), lamA), compose(lamB, sig), SxT)[0]
    okr = same(compose(whisk_right(sig, A, StT), rhoA), compose(rhoB, sig), SxT)[0]
    out(f'  Prop1/Rmk3: pi_!A laws {all(LA.values())}; sigma lam-square {okl}, rho-square {okr}'); assert okl and okr
    # Thm 5 / clause-2 readings
    for V in (2, 3):
        s0 = S[0]
        FB = Fbar(B, lamB, rhoB, PB, s0, T, V); FA = Fbar(A, lamA, rhoA, SxT, s0, T, V)
        assert FB == nT**nS * V**nS and FA == nT * V**nS
    wc, ww = (FA == V**nS), (FA == FB)
    out(f'  F_A={FA}, F_B={FB}, content={V**nS} (V=3): worker=content {wc} [pred |T|=1: {nT==1}]; '
        f'worker=weld {ww} [pred: {nS==1 or nT==1}]; carriers A=B {len(SxT)==len(PB)}')
    assert wc == (nT == 1) and ww == (nS == 1 or nT == 1) and (len(SxT) == len(PB)) == ww
    # Prop 2: e fails lam-square exactly at the non-constant g
    e = lambda p: ((p[0], const[gof(p)[p[0]]]), {d: d for d in B.dirs(p)})
    le, re_ = compose(whisk_left(StS, e, B, B), lamB), compose(lamB, e)
    fails = {p for p in PB if le(p) != re_(p)}
    NCp = {p for p in PB if len(set(gof(p).values())) > 1}
    rok = same(compose(whisk_right(e, B, StT), rhoB), compose(rhoB, e), PB)[0]
    out(f'  Prop2: e rho-square {rok}; lam-square fails on exactly the non-constant positions: {fails == NCp} ({len(fails)})')
    assert rok and fails == NCp
    # Thm 7 at EVERY base point s0 (not just s0=0)
    for s0 in S:
        r = lambda p, s0=s0: ((p[0], const[gof(p)[s0]]), {d: d for d in B.dirs(p)})
        ok = same(compose(whisk_left(StS, r, B, B), lamB), compose(lamB, r), PB)[0] and \
             same(compose(whisk_right(r, B, StT), rhoB), compose(rhoB, r), PB)[0]
        idem = all(r(r(p)[0])[0] == r(p)[0] for p in PB)
        img = {r(p)[0] for p in PB} == set(Bp); sec = all(r(p)[0] == p for p in Bp)
        fib = {q: sum(1 for p in PB if r(p)[0] == q) for q in Bp}
        assert ok and idem and img and sec and set(fib.values()) == {nT**(nS-1)}
    out(f'  Thm7: r_s0 bicomodule retraction, idempotent, image B\', fibre {nT**(nS-1)}: all s0 in S')
    # Rmk 4 last sentence: "at |S|=1 the worker is the content, on the nose"
    if nS == 1:
        out(f'  Rmk4 claim "|S|=1 => worker is the content": {FA == V**nS}  (V=3: F_A={FA} vs V^S={V**nS})')

if __name__ == '__main__':
    lines = []
    def out(x): print(x); lines.append(x)
    for sh in [(1, 2), (1, 4), (4, 1), (3, 4), (4, 3), (2, 5)]:
        r3(*sh, out)
    out('\nALL ROUND-3 ASSERTIONS PASSED')

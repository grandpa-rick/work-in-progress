#!/usr/bin/env python3
"""Simpler witness: n=2 variables, alpha=(1,1). Both strips are single boxes in h=2k>=2 sites, so nu(P)=2
for every P; Theorem 11 then needs ord_{t=1} c_(1,1) >= 2.  Compute c_(1,1)(t) and HKKO t=1 signed count."""
import sympy as sp
from warnaar_calpha import calpha, ord1, t
from hkko_paths import paths, cminus
for k in range(1,5):
    for l in range(0,4):
        c=calpha(k,l,[1,1]); s=sum(cminus(xf,k,2*l+2) for xf,_ in paths(k,l,[1,1]))
        print(f"k={k} l={l}: c_(1,1) = {sp.factor(c)}  ord={ord1(c)}  HKKO signed count={s}  c(1)={c.subs(t,1)}")

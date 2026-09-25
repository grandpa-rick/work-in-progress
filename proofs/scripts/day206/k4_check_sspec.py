"""Check q-specialised runs (q=1/s0, t symbolic) against the closed form."""
import json, glob, sympy as sp
from k4_analyze import qint, q, t, u
from k4_closed_form_check import Phat3  # noqa (module prints its own check; harmless)
for fn in sorted(glob.glob('/home/agent/projects/proofs/scripts/day206/k4_r*_m*_s*.json')):
    d=json.load(open(fn)); tau=sp.sympify(d['tau']); s0=d['s0']
    r=int(fn.split('_r')[1].split('_')[0])
    pred=(-(q**4-1)*qint(r+4)*Phat3.subs(u,t**r)/qint(4)/q**10).subs(q,sp.Rational(1,s0))
    print(fn.split('/')[-1], ': tau - closed form at q=1/%d ='%s0, sp.cancel(tau-pred), '[held-out]' if r not in (5,6,7,8) else '[fit r]')

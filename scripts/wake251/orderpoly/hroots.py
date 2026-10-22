# Root data of the TRUE H_kappa(j) on LG(m,2m) via the verified LS model, m<=7.
import sys, itertools
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
from lsmodel import model
from roots import roots_of, strict_in_staircase
from orderpoly import shifted, boxless, longest_chain
import sympy as sp
j = sp.Symbol('j')
for m in range(1, 8):
    viol = {'HP*(Re j>=-m)': 0, 'Re j<=-1/2?': 0, 'OC': 0}; tot = 0; lines = []
    gmin = 1e9; gmax = -1e9
    for lam in strict_in_staircase(m):
        if not lam: continue
        d = sum(lam); vals = [model(tuple(lam), x) for x in range(d + 2)]
        H = sp.expand(sp.interpolate(list(zip(range(d + 1), vals[:d + 1])), j))
        assert H.subs(j, d + 1) == vals[d + 1]  # out-of-fit check
        rs = roots_of(H, j); re = [z.real for z in rs]; tot += 1
        r = longest_chain(shifted(lam), boxless)
        mn, mx = min(re), max(re); gmin = min(gmin, mn); gmax = max(gmax, mx)
        oc = [z for z in rs if abs(z.imag) < 1e-9 and abs(2 * z.real - round(2 * z.real)) < 1e-9 and round(2 * z.real) % 2 and z.real <= -m - 0.5 + 1e-9]
        viol['HP*(Re j>=-m)'] += mn < -m - 1e-9; viol['OC'] += bool(oc)
        lines.append(f'  {tuple(lam)} d={d} r={r} lam1={lam[0]} len={len(lam)} Re j in [{mn:.4f},{mx:.4f}]  minRe+lam1+? -> minRe={mn:.3f} vs -(lam1)={-lam[0]}, -(len+?)')
    print(f'm={m}: #kappa={tot}  global Re j in [{gmin:.4f},{gmax:.4f}]  HP* violations {viol["HP*(Re j>=-m)"]}  OC violations {viol["OC"]}')
    if m <= 5 or True:
        for L in lines: print(L)
    sys.stdout.flush()

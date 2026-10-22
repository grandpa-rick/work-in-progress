import sys, itertools
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake251/oc_roots')
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
import sympy as sp
import oc_roots as O
from orderpoly import *
j = O.j
def Hs_for(m):
    subs = [k for r in range(m + 1) for k in itertools.combinations(range(1, m + 1), r)]
    if m <= 4: return O.src_dem(m, subs)
    return O.src_file('/home/agent/projects/work-in-progress/scripts/wake251/oc_roots/hp0big2_n5_rerun.txt')(m, subs)
def conj(lam):
    return [sum(1 for x in lam if x > i) for i in range(lam[0])] if lam else []
def compl(kap, m):  # complement in [m]
    return tuple(x for x in range(1, m + 1) if x not in kap)
variants = {
  'shifted(kappa) t=j+1': lambda k, m: (shifted(sorted(k, reverse=True)), j + 1),
  'young(kappa) t=j+1': lambda k, m: (young(sorted(k, reverse=True)), j + 1),
  'young(conj kappa) t=j+1': lambda k, m: (young(conj(sorted(k, reverse=True))), j + 1),
  'shifted(kappa) t=2j+1': lambda k, m: (shifted(sorted(k, reverse=True)), 2 * j + 1),
  'young(kappa) t=2j+1': lambda k, m: (young(sorted(k, reverse=True)), 2 * j + 1),
}
def main():
  for m in (3, 4, 5):
    Hs = Hs_for(m)
    print(f'===== m={m}: {len(Hs)} kappas')
    for name, V in variants.items():
        ok = 0; bad = []
        for k, H in Hs.items():
            P, sub = V(k, m)
            Om = omega_poly(P, boxless).subs(t, sub)
            if sp.expand(Om - H) == 0: ok += 1
            else: bad.append(k)
        print(f'  {name}: match {ok}/{len(Hs)}; first fails {bad[:6]}')
    # raw data: H(1) and leading coefficient*dim! (degree)
    for k, H in sorted(Hs.items(), key=lambda x: (sum(x[0]), x[0])):
        d = sum(k); lc = sp.Poly(H, j).LC() * sp.factorial(d) if d else 1
        print('   ', k, 'dim', d, 'H(0..3)=', [H.subs(j, x) for x in range(4)], 'deg*=', lc, 'factor', sp.factor(H))

if __name__ == '__main__': main()

import sys, itertools
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake251/reciprocity')
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
from recip import setup, F, Fint
from orderpoly import shifted, boxless, longest_chain
import sympy as sp
x = sp.Symbol('x')
def polys(kap):
    S = setup(kap); d = sum(kap)
    vals = [F(S, N) for N in range(2*d + 4)]
    E = sp.expand(sp.interpolate([(N, vals[N]) for N in range(0, len(vals), 2)], x))
    O = sp.expand(sp.interpolate([(N, vals[N]) for N in range(1, len(vals), 2)], x))
    return S, E, O, (E+O)/2, (E-O)/2
if __name__ == '__main__':
    m = int(sys.argv[1])
    for r in range(1, m+1):
        for k in itertools.combinations(range(1, m+1), r):
            S, E, O, p, q = polys(k)
            print(k, 'deg p,q', sp.degree(p,x), sp.degree(q,x), ' q=', sp.factor(q), ' q(-M) M=1..9:', [q.subs(x,-M) for M in range(1,10)])

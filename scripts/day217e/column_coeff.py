"""[P_{1^n}(x;q,T)] e_mu  ==  [n;mu]_T * prod (q;T)_{mu_i} / (q;T)_n ?"""
import sys,itertools,sympy as sp
sys.argv=['x','4']; exec(open('scripts/day217e/estar_edges.py').read().split('ok=0;bad=0')[0])
def poch(a,n): return sp.prod([1-a*T**i for i in range(n)])
def tfac(n): return sp.prod([(1-T**i)/(1-T) for i in range(1,n+1)])
for n in range(1,5):
    P,xs=macP(n); Ps=plist(n)
    # solve e_mu = sum c_nu P_nu (triangular): use full linear solve in m-basis
    Mat=sp.Matrix([[P[nu].get(m,0) for nu in Ps] for m in Ps])
    for mu in Ps:
        em=m_of(sp.prod([sum(sp.prod(c) for c in itertools.combinations(xs,r)) for r in mu]),xs)
        c=Mat.LUsolve(sp.Matrix([em.get(m,0) for m in Ps]))
        got=sp.cancel(c[Ps.index(tuple([1]*n))])
        pred=tfac(n)/sp.prod([tfac(r) for r in mu])*sp.prod([poch(q,r) for r in mu])/poch(q,n)
        print(mu, sp.cancel(got-pred)==0)

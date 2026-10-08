# direct subset formula (Prop 2.1 as printed), N=3, e*_{111} and e*_{21}; read off c_{lam,mu}
import sympy as sp, itertools
s,t=sp.symbols('s t'); N=3; x=sp.symbols('x1:4')
def E(k,F):
    tot=0
    for A in itertools.combinations(range(N),k):
        Ac=[j for j in range(N) if j not in A]
        cA=sp.prod([(x[i]-t*x[j])/(x[i]-x[j]) for i in A for j in Ac])
        sub={x[i]:s*x[i] for i in A}
        tot+=cA*sp.prod([x[i] for i in A])*F.subs(sub,simultaneous=True)
    return sp.factor(sp.cancel(tot))
e=[1,sum(x),x[0]*x[1]+x[0]*x[2]+x[1]*x[2],x[0]*x[1]*x[2]]
F=E(1,E(1,E(1,sp.Integer(1))))
a,b,c=sp.symbols('a b c')
ans=sp.solve(sp.Poly(sp.expand(F-(a*e[1]**3+b*e[2]*e[1]+c*e[3])),*x).coeffs(),[a,b,c])
for k,v in ans.items(): print(k,sp.factor(v), '| (s-1)-series:',sp.series(v.subs(s,1+sp.Symbol('h')),sp.Symbol('h'),0,3))

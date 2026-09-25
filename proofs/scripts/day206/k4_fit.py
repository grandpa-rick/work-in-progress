"""Fit Q_r = q^10 tau_r^(4) [4]_t / ((q^4-1)[r+4]_t) as sum_{j<=J} p_j(q,t) u^j, u=t^r, p_j poly (q-deg<=6, t-deg<=D).
Fit on FIT_R, check on all other r (held out)."""
import sympy as sp, sys
from k4_analyze import load, qint, q, t, u
data=load(); rs=sorted({r for r,m in data})
def Q(r):
    tau=data[(r,min(m for rr,m in data if rr==r))]
    return sp.expand(sp.cancel(q**10*tau*qint(4)/((q**4-1)*qint(r+4))))
J=int(sys.argv[1]) if len(sys.argv)>1 else 3
D=int(sys.argv[2]) if len(sys.argv)>2 else 8
FIT=[int(x) for x in sys.argv[3].split(',')] if len(sys.argv)>3 else [5,6,7,8]
cs={}; P=0
for j in range(J+1):
    pj=0
    for a in range(0,7):
        for b in range(0,D+1):
            c=sp.Symbol(f'c_{j}_{a}_{b}'); cs[c]=1; pj+=c*q**a*t**b
    P+=pj*u**j
eqs=[]
for r in FIT:
    diff=sp.Poly(sp.expand(P.subs(u,t**r)-Q(r)),q,t)
    eqs+=diff.coeffs()
sol=sp.solve(eqs,list(cs),dict=True)
if not sol: print('NO SOLUTION', J, D, FIT); sys.exit()
sol=sol[0]; free=[c for c in cs if c not in sol]
Pf=sp.expand(P.subs(sol).subs({c:0 for c in free}))
print(f'J={J} D={D} fit on {FIT}; free params left: {len(free)}')
Pc=sp.collect(Pf,u,evaluate=False)
for j in sorted(Pc, key=lambda e: sp.degree(e,u)):
    print(f'  [{j}]:', sp.factor(Pc[j]))
for r in rs:
    tag='FIT' if r in FIT else 'held-out'
    print(f'r={r}: Q_r - P(t^r) =', sp.expand(Pf.subs(u,t**r)-Q(r)), f'[{tag}]')
import pickle; pickle.dump(Pf,open('P3.pkl','wb'))

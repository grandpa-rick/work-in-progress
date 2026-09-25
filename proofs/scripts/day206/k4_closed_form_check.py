"""Day 206: check q^10 tau_r^(4) = -(q^4-1)[r+4]_t Phat_3(t^r,t)/[4]_t against every computed (r,m)."""
import sympy as sp
from k4_analyze import load, qint, q, t, u
N=lambda j: sp.prod([(t**i*u-1) for i in range(1,j+1)])
Phat3=(q**6*N(3) + q**5*qint(4)*N(2)
       + q**3*(t**2+1)*(q*t**3+q*t**2-q+t+1)*N(1)
       + qint(4)*(q*(t-1)+1)*(q**2*(t**2-1)+1))
FIT={5,6,7,8}
for (r,m),tau in sorted(load().items()):
    d=sp.cancel(q**10*tau + (q**4-1)*qint(r+4)*Phat3.subs(u,t**r)/qint(4))
    print(f'r={r:2d} m={m:2d}: diff = {d}   [{"FIT" if r in FIT else "held-out"}]')
print('Phat3 cyclotomic check: Phat3(i^r, i) for r=0..3 :', [sp.simplify(Phat3.subs({u:sp.I**r,t:sp.I})) for r in range(4)])
print('Phat3(u=(-1)^r, t=-1):', [sp.simplify(Phat3.subs({u:(-1)**r,t:-1})) for r in range(2)])

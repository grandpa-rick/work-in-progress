from psi import *
import random
z=sp.Symbol('z')
def S0(H,M):
    e=sp.Integer(1); xs=[w(m) for m in H]
    for x in xs: e*=(1-x**-2)
    for i in range(len(xs)):
        for k in range(i+1,len(xs)): e*=(1-1/(xs[i]*xs[k]))
    for t in M:
        if t in H: continue
        for x in xs: e*=(x-1/w(t))/(x-w(t))
    return e
def Phi(r,f,M): return sum(S0(H,M)*f([w(m) for m in H]) for H in itertools.combinations(M,r))
def G(S,zz): return sp.Mul(*[(zz-1/w(t))/(zz-w(t)) for t in S])
def res0inf(expr):
    # Res_0 + Res_inf = -(sum of residues at other poles) ; compute directly
    r0=sp.residue(expr,z,0); rinf=-sp.residue(expr.subs(z,1/z)/z**2,z,0)
    return r0+rinf
M=(1,2,3,4); vals={w(k):sp.Rational(random.randint(2,30),random.randint(2,30)) for k in range(1,8)}
for r,f in [(2,lambda x:x[0]*x[1]+x[0]**2+x[1]**2),(2,lambda x:(x[0]*x[1])**-1*(x[0]+x[1])),(3,lambda x:x[0]*x[1]*x[2]*(x[0]+x[1]+x[2]))]:
    lhs=Phi(r,f,M).subs(vals)
    rhs=0
    for Hp in itertools.combinations(M,r-1):
        S=[t for t in M if t not in Hp]
        integrand=(z**-1*G(S,z)*f([w(m) for m in Hp]+[z])).subs(vals)
        rhs+=-S0(Hp,M).subs(vals)*res0inf(sp.together(integrand))
    print(r, sp.simplify(lhs-rhs))

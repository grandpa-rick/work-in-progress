import sympy as sp, pickle
q,t,u=sp.symbols('q t u')
P3=pickle.load(open('P3.pkl','rb'))
def qi(n,x=t): return sum(x**i for i in range(n))
def prodk(k): return sp.prod([(t**i*u-1) for i in range(1,k)])
P1=q*(t*u-1)+qi(2)
P2=q**3*(t*u-1)*(t**2*u-1)+qi(3)*(q**2*(t*u-1)+q*(t-1)+1)
print('P3 irreducible? factor:', sp.factor(P3) == P3 or sp.factor_list(P3))
for sgn in (1,-1):
    R=sp.expand(P3 - sgn*q**6*prodk(4))
    Rq,rem=sp.div(sp.Poly(R,u,t,q),sp.Poly(qi(4),u,t,q))
    print(f'sign {sgn}: P3 - sgn q^6 prod = [4]_t * R + rem, rem =', rem.as_expr())
    if rem.is_zero:
        R=Rq.as_expr()
        print('  R =', sp.collect(sp.expand(R),u))
        # compare with lower levels: try R = a*q^?*(tu-1)(t^2u-1) + ...
        for sg2 in (1,-1):
            R2=sp.expand(R - sg2*q**5*prodk(3))
            print(f'  R - {sg2} q^5 (tu-1)(t^2u-1) =', sp.collect(R2,u), ' factor:', sp.factor(R2))
# P at u = t^{-k} (i.e. r=-k): 
for k,P in ((2,P1),(3,P2),(4,P3)):
    print(f'k={k}: P(u=t^-{k}) =', sp.factor(P.subs(u,t**-k)), '| P(u=1) =', sp.factor(P.subs(u,1)), '| P at t=1:', sp.factor(P.subs(t,1)))

"""Day 209: machine-check each displayed intermediate of the (TC) proof writeup."""
import sympy as sp
s,t,x,y,z,w,A,B,X,W=sp.symbols('s t x y z w A B X W')
ok=True
def chk(name,e):
    global ok; r=sp.simplify(sp.together(e))==0; ok&=r; print(name,r)
g,d=A*z,B*w; c=s**2*z*w/(g*d)
kg=(g-s*z)*(g-s*w)/(g*(g-d)); kd=(d-s*z)*(d-s*w)/(d*(d-g))
F=(y-s*z)*(y-s*w)/(y*(y-g)*(y-d))
chk('(P1) partial fractions of F', F-(c/y+kg/(y-g)+kd/(y-d)))
U=1-kg*x/(g-c*x)-kd*x/(d-c*x)
chk('(P2) U = xF(cx)', U-x*F.subs(y,c*x))
Xs, Ws = s*x/(A*B*z), s*x/(A*B*w)
chk('(P3) U=(1-X)(1-W)/((1-sX/A)(1-sW/B))', U-((1-Xs)*(1-Ws)/((1-s*Xs/A)*(1-s*Ws/B))))
# shifted multiplier
chk('(P4) x/(t^{i-1}z - c t x) = (tBX/s)/(1-st^2X/A)', x/(A/t*z-c*t*x)-(t*B*Xs/s)/(1-s*t**2*Xs/A))
chk('(P5) x/(g-cx) = (BX/s)/(1-sX/A)', x/(g-c*x)-(B*Xs/s)/(1-s*Xs/A))
# K-ratio times kappa'
rho=sp.Symbol('rho')
rKi=B*(A/t*rho-B)*(A/t*rho-1/t)/((A/t*rho-s)*(A/t*rho-B/t))
kgp=(A/t-s)*(A/t*rho-s)/((A/t)*(A/t*rho-B))
chk('(P6) rKi*kappa_gamma(i-1,j) = B(A-st)(A rho-1)/(A(A rho-B))', rKi*kgp-B*(A-s*t)*(A*rho-1)/(A*(A*rho-B)))
# final polynomial identity
u,v=s*t*X,s*t*W
chk('(P7) poly identity', (A*W-B*X)*((1-u)*(1-v)-(A-u)*(B-v))-((1-A)*(B-v)*(A*W-X)+(1-B)*(A-u)*(W-B*X)))
# K ratio with actual products for small i,j
def Kp(i,j): return sp.prod([(t**p*rho-s)/(t**p*rho-t**j) for p in range(i)])*sp.prod([(s*rho-t**r)/(t**i*rho-t**r) for r in range(j)])
for i in range(1,4):
    for j in range(0,4):
        chk(f'(P8) K ratio i={i} j={j}', Kp(i-1,j)/Kp(i,j)-rKi.subs({A:t**i,B:t**j}))
print('ALL OK' if ok else 'FAIL')

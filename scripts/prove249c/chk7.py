from hp import *
n=4; a=[Fr(k*k+1,7) for k in range(1,n+1)]
import random
w={i:Fr(random.randint(2,50),random.randint(2,50)) for i in range(1,n+1)}
print('X exact at c=7, random point:', float(Cfun([],[],[1,3,4],7,n,w)))
c_=sp.Symbol('c')
pts=[(c,Plim([],[],[1,3,4],c,n,a)) for c in range(-2,20,2)]
P=sp.expand(sp.interpolate(pts,c_)); print(sp.factor(P)); print('P(7)=',P.subs(c_,7))
print('check odd 9:', Plim([],[],[1,3,4],9,n,a), P.subs(c_,9))
print([complex(r) for r in sp.Poly(P,c_).nroots()])

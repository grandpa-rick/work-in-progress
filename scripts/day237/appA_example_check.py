# Day 237 (c): print check of the worked example after Thm A.4: m=3, lambda=(2,0,1), Y°_1 = t^2 pi° T_2^{-1} T_1^{-1}
import sympy as sp
s,t=sp.symbols('s t'); x=sp.symbols('x1:4'); m=3
def sw(f,i): return f.subs({x[i]:x[i+1],x[i+1]:x[i]},simultaneous=True)
def T(f,i): g=sw(f,i); return sp.expand(t*g+sp.cancel((t-1)*x[i+1]*(g-f)/(x[i]-x[i+1])))
def Ti(f,i): return sp.expand(T(f,i)/t-(1-1/t)*f)
def pic(f): return sp.expand(f.subs({x[0]:x[1],x[1]:x[2],x[2]:s*x[0]},simultaneous=True))
f=x[0]**2*x[2]
out=sp.expand(t**2*pic(Ti(Ti(f,0),1)))   # T_1^{-1} first (rightmost), then T_2^{-1}, then pi°
claim=s**2*x[0]**2*x[2]-s*(t-1)*x[0]*x[2]**2-t*(t-1)*x[1]**2*x[2]-s*(t-1)*x[0]*x[1]*x[2]
print('example printed (gr + extra term) ==', sp.simplify(out-claim)==0)
print(sp.factor_terms(out))

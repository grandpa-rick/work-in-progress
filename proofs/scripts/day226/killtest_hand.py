# Day 226 cold recheck: sum hand-derived pieces of the (3,3,3)->(7,2) kill test. sympy = calculator only.
from sympy import *
t=symbols('t')
def br(n): return sum(t**i for i in range(n))
br_={}
curly = t**-2/((1-t)*(1-t**2)) + t**-4*(1+t**3)**2 - t**-4*(1+t**3+t**6)**2*(1+t**2+t**4)/(1-t**3)
Phi = -(1-t**7)*(1-t**2)*curly
U = -Phi - (1+t**3+t**6)**3
Lam = 2*br(7)*(1+t**3+t**6)
D = (1+t**3)*(br(7)*(1-t**2+t**4) + br(5)*(1+t**3+t**6))
tot = expand(cancel(U+Lam+D))
print('U =', expand(cancel(U)))
print('total =', tot)
claim = 2*t**13+3*t**12+3*t**11+6*t**10+6*t**9+6*t**8+9*t**7+6*t**6+3*t**5+9*t**4+6*t**3+3*t+4
print('match:', expand(tot-claim)==0)

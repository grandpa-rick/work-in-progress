# Second proof route of the (2,3) tie from Rick's own proved lemmas: label lowering Psi_{(2),0}=Psi_{(1),0}-w1 Psi_{(1),-1}
# (rho_2(x)=1-w1/x) and the (1)-identity Psi_{(1),0}=1-e^-2 (P1, r=1). Symbolic check of each step.
import sympy as sp
from tie23_defs import Psi1
W=sp.symbols('w1:7'); w={i+1:W[i] for i in range(6)}
for S in [(2,3),(1,2,3),(2,3,4),(1,2,4,5)]:
    e=sp.Mul(*[w[m] for m in S])
    ll=sp.cancel(Psi1(2,0,S,w)-(Psi1(1,0,S,w)-w[1]*Psi1(1,-1,S,w)))
    one=sp.cancel(Psi1(1,0,S,w)-(1-e**-2))
    print(S,'label lowering diff',ll,' (1)-identity diff',one)

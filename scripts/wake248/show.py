import sys
from eng import *
from sympy import factor,sympify
N=int(sys.argv[1]); C=setup(N)
for spec in sys.argv[2:]:
    a,b,l,m,k,c=eval(spec)
    r=A(Sys(N,a,b),a,b,l,m,k,c)
    print(spec,'\n  N =',factor(sympify(str(r.n))),'\n  D =',factor(sympify(str(r.d))))

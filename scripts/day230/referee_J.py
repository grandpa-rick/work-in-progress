import sympy as sp, itertools
t=sp.symbols('t')
def K(lam):
    l=len(lam); E=list(itertools.combinations(range(l),2)); tot=0
    for k in range(len(E)+1):
        for H in itertools.combinations(E,k):
            par=list(range(l))
            def f(x):
                while par[x]!=x: x=par[x]
                return x
            for i,j in H: par[f(i)]=f(j)
            if len({f(i) for i in range(l)})==1: tot+=sp.prod([t**(lam[i]*lam[j])-1 for i,j in H])
    return tot
def Lead(lam):
    n=sum(lam); return sp.cancel((1-t**n)/sp.prod([1-t**x for x in lam])*K(lam))
J=sp.factor(sp.cancel(Lead((1,1,1,1))*Lead((2,2))/Lead((2,1,1))**2)); print('J=',J)
print('paper:',sp.simplify(J-(t**2+1)*(t**3+3*t**2+6*t+6)/((t+1)*(t**2+t+2)**2)))
print('L(111)=',sp.factor(Lead((1,1,1))))

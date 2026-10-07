import sympy as sp, itertools
t=sp.symbols('t')
def K(lam):
    l=len(lam); E=list(itertools.combinations(range(l),2)); tot=0
    for k in range(l-1,len(E)+1):
        for H in itertools.combinations(E,k):
            # connectivity
            p=list(range(l))
            def f(x):
                while p[x]!=x: x=p[x]
                return x
            for i,j in H: p[f(i)]=f(j)
            if len({f(i) for i in range(l)})==1:
                tot+=sp.prod([t**(lam[i]*lam[j])-1 for i,j in H])
    return tot
ok=0;bad=[]
for lam in [(1,1),(2,1),(1,1,1),(3,1,1),(2,2,1),(3,2,1,1),(2,2,2,2),(4,3,1),(1,1,1,1,1)]:
    n=sum(lam);l=len(lam)
    L=sp.cancel((1-t**n)/sp.prod([1-t**x for x in lam])*K(lam))
    c=L.subs(t,1); b=L.subs(t,0)
    good = c==(-1)**(l-1)*n**(l-1) and b==(-1)**(l-1)*sp.factorial(l-1)
    sg = all(sp.sign(L.subs(t,v))==(-1)**(l-1) for v in [sp.Rational(1,3),sp.Rational(1,2),sp.Rational(9,10),2,3])
    print(lam,c,b,good,sg)

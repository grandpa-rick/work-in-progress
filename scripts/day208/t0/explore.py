import sympy as sp
s,t=sp.symbols('s t')
def poch(a,n): return sp.prod([1-a*t**i for i in range(n)])
def alpha(j):
    if j<0: return 0
    return sp.prod([s-t**i for i in range(1,j+1)])/poch(t,j)
def c(n,j):
    if j<0 or n-j<0: return 0
    return poch(s,n-j)/poch(t,n-j)*(alpha(j)-s*t**(n-j)*alpha(j-1))
def F(n,w): return sum(c(n,j)*w**j for j in range(n+1))
for k in range(1,6):
  for r in range(0,7):
    coll={}
    for b in range(k+1):
        key=tuple(sorted((b,r+k-b)))
        coll[key]=coll.get(key,0)+s**b*F(k-b,t**(r-b))
    out={}
    for key,v in coll.items():
        v=sp.factor(sp.cancel(v))
        lim=sp.limit(v,t,0)
        out[key]=sp.factor(lim)
    print(k,r,out, sp.simplify(sum(out.values())))

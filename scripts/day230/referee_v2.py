# Referee hand-check: implement Thm 6.6 (closed l=3,kappa=1 leads) LITERALLY from the FPSAC text, compare with engine log.
import sympy as sp, re, itertools
t,w=sp.symbols('t w')
def br(m,q=None):
    q=t if q is None else q
    return sum(q**i for i in range(m))
def Lam(a,r,q): return (-1)**(r+q)*br(a+r+q)/br(a)*br(a,t**r)*br(a,t**q)
def L(a,b): return (1-t**(a+b))*(t**(a*b)-1)/((1-t**a)*(1-t**b))
def M(k,r):  # dict partition(tuple, zeros dropped)->coef
    if k<r: k,r=r,k
    d={}
    def add(p,c):
        p=tuple(sorted([x for x in p if x>0],reverse=True)); d[p]=d.get(p,0)+c
    if r>0: add((k,r),r)
    for j in range(1,r+1): add((k+j,r-j),L(k-r+j,j))
    return d
def Da_of(a,poly):
    out={}
    for p,c in poly.items():
        for i in range(len(p)):
            rest=p[:i]+p[i+1:]
            for q,c2 in M(a,p[i]).items():
                key=tuple(sorted(rest+q,reverse=True)); out[key]=out.get(key,0)+c*c2
    return out
def Phi(a,b,c,x,y):
    n=a+b+c
    g=lambda v: sum(z**b for z in v)*sum(z**c for z in v)
    tot=0
    for A in range(1,a):
        B=a-A
        GA=sp.Poly(sp.expand(g([t**i for i in range(A)]+[w*t**i for i in range(B)])),w)
        co=lambda k: GA.coeff_monomial(w**k) if k>=0 else 0
        term=co(y-B)/((1-t**A)*(1-t**B))
        term+=sum((t**(-A*m)-t**(B*m))*co(y-B-m) for m in range(1,y-B+1))/(1-t**a)
        tot+=t**(-A*B)*term
    tot-=g([t**i for i in range(a)])/(1-t**a)*sum(t**(-j*y) for j in range(a))
    return (-1)**a*(1-t**x)*(1-t**y)*tot
def lead(a,b,c,x,y):
    n=x+y; m=2 if x==y else 1
    U=((-1)**n*Phi(a,b,c,x,y)-Lam(a,b,c))/m
    res=(-1)**(b+c)*U
    res+=sum((-1)**(r+c)*Lam(a,r,c) for r in range(1,b) if sorted([b-r,a+r+c])==sorted([x,y]))
    res+=sum((-1)**(b+q)*Lam(a,b,q) for q in range(1,c) if sorted([c-q,a+b+q])==sorted([x,y]))
    res+=Da_of(a,M(b,c)).get(tuple(sorted([x,y],reverse=True)),0)
    return sp.factor(sp.cancel(res))
ok=bad=0
for f in ['newcases_n10.log','table3_n8.log']:
  for line in open('/home/agent/projects/proofs/scripts/day224/'+f):
    mm=re.match(r'(\([\d, ]*\)) -> (\([\d, ]*\)) : (.*?)  \|',line)
    if not mm: continue
    lam,mu,p=eval(mm.group(1)),eval(mm.group(2)),sp.sympify(mm.group(3))
    if len(mu)!=2 or len(lam)!=3: continue
    for (a,b,c) in set(itertools.permutations(lam)):
        v=lead(a,b,c,*mu); g=sp.expand(v-p)==0; ok+=g; bad+=not g
        if not g: print('MISMATCH',(a,b,c),mu,v,p)
print('ok',ok,'bad',bad)
print('(3,3,3)->(7,2):',sp.expand(lead(3,3,3,7,2)))
print('(4,4,2)->(7,3):',lead(4,4,2,7,3),'|',lead(2,4,4,7,3))

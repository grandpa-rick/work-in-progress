"""Day 223: B(e_k,e_r) = d/ds|_{s=1} of 207b Pieri  sum_b s^b F_{k-b}(t^{r-b}) e_b e_{r+k-b}; compare with
wake223 computed support formula r e_k e_r + sum_{j=1}^r L(k-r+j,j) e_{k+j} e_{r-j}. Symbolic. Grade: computed."""
from sympy import symbols, prod, diff, cancel, simplify
s,t,w=symbols('s t w')
def poch(a,n): return prod([(1-a*t**i) for i in range(n)]) if n>0 else 1
def alpha(j): return 0 if j<0 else prod([(s-t**i) for i in range(1,j+1)])/poch(t,j)
def c(n,j): return 0 if j<0 or j>n else poch(s,n-j)/poch(t,n-j)*(alpha(j)-s*t**(n-j)*alpha(j-1))
def F(n,x): return sum(c(n,j)*x**j for j in range(n+1))
L=lambda a,b:(1-t**(a+b))*(t**(a*b)-1)/((1-t**a)*(1-t**b))
ok=bad=0
for n in range(1,8):  # closed form of F'_n
    d=cancel(diff(F(n,w),s).subs(s,1)); cf=-(1-w**n)*(1-t**n*w)/((1-t**n)*(1-w))
    g=cancel(d-cf)==0; ok+=g; bad+=not g
for k in range(1,8):
  for r in range(0,k+1):
    B={}
    for b in range(0,k+1):
        key=tuple(sorted([x for x in (b,r+k-b) if x>0],reverse=True))
        B[key]=B.get(key,0)+diff(s**b*F(k-b,t**(r-b)),s).subs(s,1)
    B={kk:cancel(v) for kk,v in B.items() if cancel(v)!=0}
    pred={}
    if r>0: pred[tuple(sorted([k,r],reverse=True))]=r
    for j in range(1,r+1):
        key=tuple(sorted([x for x in (k+j,r-j) if x>0],reverse=True)); pred[key]=pred.get(key,0)+L(k-r+j,j)
    pred={kk:cancel(v) for kk,v in pred.items() if cancel(v)!=0}
    g=set(B)==set(pred) and all(cancel(B[x]-pred[x])==0 for x in B); ok+=g; bad+=not g
    print(k,r,'OK' if g else ('FAIL',B,pred),flush=True)
print('TOTAL OK',ok,'FAIL',bad)

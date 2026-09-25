# psi(p_m) in QH*(Gr_{k,n}) = Z[q][e_1..e_k]/<h_{n-k+1},...,h_{n-1}, h_n + (-1)^k q>  (BCF / Siebert-Tian presentation)
import sympy as sp
def check(n,k):
    q=sp.Symbol('q'); e=sp.symbols('e1:%d'%(k+1)); E=[1]+list(e)+[0]*(3*n)
    h=[sp.Integer(1)]
    for a in range(1,3*n): h.append(sp.expand(sum((-1)**(i-1)*E[i]*h[a-i] for i in range(1,a+1))))
    p=[None]
    for m in range(1,2*n+2): p.append(sp.expand(sum((-1)**(i-1)*E[i]*p[m-i] for i in range(1,m)) + (-1)**(m-1)*m*E[m]))
    I=[h[j] for j in range(n-k+1,n)]+[h[n]+(-1)**k*q]
    G=sp.groebner(I,*e,q,order='grevlex')
    red=lambda f: G.reduce(sp.expand(f))[1]
    pn=red(p[n]); hn=red(h[n])
    per=all(sp.expand(red(p[n+m]-(-1)**(k-1)*q*p[m]))==0 for m in range(1,n))
    return pn, hn, sp.expand(pn-(-1)**(k-1)*k*q)==0, per
for n in range(3,7):
    for k in range(1,n): print(n,k,check(n,k))

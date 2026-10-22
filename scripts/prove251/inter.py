import sys, itertools, sympy as sp
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/prove249d')
from dem import Di
n=int(sys.argv[1]); a=int(sys.argv[2]); Y=sp.symbols('y1:%d'%(n+1))
f={}
for S in itertools.product((0,1),repeat=n):
    k=tuple(-a-s for s in S); f[k]=f.get(k,0)+1
for i in range(n,0,-1):
    f=Di(f,i,n)
    P=sp.expand(sum(c*sp.Mul(*[Y[j]**(-k[j]) for j in range(n)]) for k,c in f.items()))
    pre=sp.Mul(*[Y[j]**a*(1+Y[j]) for j in range(i-1)])
    G=sp.factor(sp.cancel(P/pre))
    print('after D_%d:'%i, G)

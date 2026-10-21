from eng import *
import sys
n=int(sys.argv[1]); cs=eval(sys.argv[2]); cases=eval(sys.argv[3])
def at(eps):
    return {i: 1+Fr(i*i+1,7)*eps for i in range(1,n+1)}
full=list(range(1,n+1))
for lam,mu,kap in cases:
  row=[]
  for c in cs:
    vals=[]
    for eps in [Fr(1,10**4),Fr(1,10**6)]:
        w=at(eps)
        vals.append(Cfun(list(lam),list(mu),list(kap),c,n,w)/Cfun([],[],full,2*n,n,w))
    row.append((c,round(float(vals[0]),3),round(float(vals[1]),4)))
  print(lam,mu,kap,row)

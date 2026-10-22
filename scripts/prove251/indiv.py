import sys, itertools
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/prove249d')
from dem import Di, word
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
from orderpoly import shifted, boxless, longest_chain
def w0A(n):  # reduced word of longest element of S_n (generators 1..n-1), as list applied left-to-right in product
    return [j for i in range(1,n) for j in range(i,0,-1)]
def chimax(kap,n,f):
    ws=word(set(kap),n)+w0A(n)   # product s_{ws[0]}...s_{ws[-1]}; apply rightmost first
    for i in reversed(ws): f=Di(f,i,n)
    return sum(f.values())
n=int(sys.argv[1])
for r in range(1,n+1):
  for kap in itertools.combinations(range(1,n+1),r):
    if max(kap)!=n: continue
    rest=sorted(kap,reverse=True)[1:]; rp=longest_chain(shifted(rest),boxless) if rest else 0
    out=[]
    for a in range(0,(rp+1)//2+2):
        vals={S:chimax(kap,n,{tuple(-a-s for s in S):1}) for S in itertools.product((0,1),repeat=n)}
        out.append((a,{''.join(map(str,S)):v for S,v in vals.items() if v}))
    print(kap,"r'",rp,out)

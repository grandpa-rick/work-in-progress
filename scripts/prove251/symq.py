import sys, itertools
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/prove249d')
from dem import Di, word
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
from orderpoly import shifted, boxless, longest_chain
def chiw(kap,n,f):
    for i in reversed(word(set(kap),n)): f=Di(f,i,n)
    return sum(f.values())
def hn_times(n,m,shift):  # e^{shift} * h_n(e^{eps}) as dict ; shift tuple
    f={}
    for c in itertools.combinations_with_replacement(range(m),n):
        w=list(shift)
        for i in c: w[i]+=1
        f[tuple(w)]=f.get(tuple(w),0)+1
    return f
m=int(sys.argv[1])
for r in range(1,m+1):
  for kap in itertools.combinations(range(1,m+1),r):
    lam=sorted(kap,reverse=True); rr=longest_chain(shifted(lam),boxless); W=(rr+1)//2
    row=[]
    for b in range(1,W+1):
        row.append([chiw(kap,m,hn_times(n,m,tuple([-b]*m))) for n in range(0,5)])
    print(kap,'window',W,' chi(Sym^n Q (-b)) b=1..W, n=0..4:',row)

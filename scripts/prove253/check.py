# PROVE 253: independent Demazure check of (Z_int) on LG and (Z_OG) on OG; full-diagram chain drop over BR p=1 supports.
import sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'../prove249d'); sys.path.insert(0,'../wake253'); sys.path.insert(0,'../wake251/orderpoly')
from dem import Di, word
from orderpoly import shifted, boxless, longest_chain
N=int(sys.argv[1])
def sp(n): return [tuple(sorted(c,reverse=True)) for L in range(n+1) for c in itertools.combinations(range(1,n+1),L)]
def DiB(f,i,n):  # type B, doubled coords (spin weight = all 1's)
    from collections import defaultdict
    g=defaultdict(int)
    for mu,a in f.items():
        if i<n: k=(mu[i-1]-mu[i])//2; al=[0]*n; al[i-1]=2; al[i]=-2
        else: k=mu[n-1]; al=[0]*n; al[n-1]=2
        if k>=0:
            for t in range(k+1): g[tuple(mu[j]-t*al[j] for j in range(n))]+=a
        elif k<=-2:
            for t in range(1,-k): g[tuple(mu[j]+t*al[j] for j in range(n))]-=a
    return {k:v for k,v in g.items() if v}
def dem(k,n,j,typ):
    w=word(set(k),n); f={tuple([j]*n):1}
    for i in reversed(w): f=(Di if typ=='C' else DiB)(f,i,n)
    return sum(f.values()), w
def interp(vals,x):  # vals at 0..d
    d=len(vals)-1; s=Fr(0)
    for i,v in enumerate(vals):
        t=Fr(v)
        for l in range(d+1):
            if l!=i: t*=Fr(x-l,i-l)
        s+=t
    return s
def r(k): return longest_chain(shifted(list(k)),boxless) if k else 0
def strict_count(k,t):  # # strictly order-preserving maps shifted diagram -> [t]
    P=sorted(shifted(list(k))); idx={p:i for i,p in enumerate(P)}
    # brute DP over boxes in row-major order: value must exceed left and upper neighbours
    def rec(i,assign):
        if i==len(P): return 1
        (a,b)=P[i]; lo=1
        for q in ((a,b-1),(a-1,b)):
            if q in assign: lo=max(lo,assign[q]+1)
        s=0
        for v in range(lo,t+1):
            assign[(a,b)]=v; s+=rec(i+1,assign)
        assign.pop((a,b),None); return s
    return rec(0,{})
if __name__=="__main__":
 okL=okO=okE=okA=tot=0; badL=[];badO=[]
 for k in sp(N):
     d=sum(k); tot+=1
     for typ in 'CB':
         vals=[];w=None
         for j in range(d+1):
             v,w=dem(k,N,j,typ); vals.append(v)
         if k: okA += (w[0]==N)   # leftmost letter of word list = last applied? record both ends
         Z=[b for b in range(1,2*d+4) if interp(vals,-b)==0]
         if typ=='C':
             pred=list(range(1,(r(k)+1)//2+1))
             if Z==pred: okL+=1
             else: badL.append((k,Z,pred))
         else:
             pred=list(range(1,r(k)+1))
             rec_ok=all(interp(vals,-b)==(-1)**d*strict_count(k,b-1) for b in range(1,r(k)+3))
             if Z==pred and rec_ok: okO+=1
             else: badO.append((k,Z,pred,rec_ok))
 print(f'm={N}: LG integer zero set == [1..ceil(r/2)] (Demazure, LS-free): {okL}/{tot}', badL[:5])
 print(f'm={N}: OG integer zero set == [1..r] and H(-b)=(-1)^|k| strict(P,b-1): {okO}/{tot}', badO[:5])

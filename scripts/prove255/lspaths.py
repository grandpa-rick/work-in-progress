# PROVE 255: LS paths of shape j*omega_n for Sp_{2n}, initial direction <= kappa, from ROOT DATA
# (no box translation). Compare (i) full character sum e^{pi(1)} with Demazure character (prove249d/dem.py),
# (ii) count with F_kappa(2j) (wake255 count), (iii) bond census: bond 1 <=> cover adds diagonal box.
import sys, itertools
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0,'../prove249d'); sys.path.insert(0,'../wake255')
from dem import Di, word
from models_from_reciprocity import count

def roots(n):  # positive roots as (root vector, coroot vector)
    R=[]
    for i in range(n):
        for k in range(i+1,n):
            v=[0]*n; v[i]=1; v[k]=-1; R.append((tuple(v),tuple(v)))
            v=[0]*n; v[i]=1; v[k]=1; R.append((tuple(v),tuple(v)))
        v=[0]*n; v[i]=2; c=[0]*n; c[i]=1; R.append((tuple(v),tuple(c)))
    return R
def pair(a,b): return sum(x*y for x,y in zip(a,b))
def dim(mu,R): return sum(1 for (a,c) in R if pair(mu,c)<0)
def refl(mu,a,c): p=pair(mu,c); return tuple(m-p*x for m,x in zip(mu,a))

def setup(n):
    R=roots(n); V=list(itertools.product([1,-1],repeat=n))
    D={v:dim(v,R) for v in V}
    covers=defaultdict(list)  # covers[v] = list of (u, bond) with u > v cover, u = s_beta v
    for v in V:
        for (a,c) in R:
            u=refl(v,a,c)
            if D[u]==D[v]+1: covers[v].append((u,abs(pair(v,c))))
    # Bruhat order via closure of covers
    up={v:{v} for v in V}
    for v in sorted(V,key=lambda x:-D[x]):
        for (u,b) in covers[v]: up[v]|=up[u]
    return R,V,D,covers,up

def chain_ok(hi,lo,a,j,covers,D):
    # exists saturated chain lo=k_r < ... < k_0=hi with a*j*bond in Z at each step
    if D[hi]<=D[lo]: return False
    cur={lo}
    for _ in range(D[hi]-D[lo]):
        nxt=set()
        for v in cur:
            for (u,b) in covers[v]:
                if (a*j*b).denominator==1: nxt.add(u)
        cur=nxt
    return hi in cur

def lspaths(kap_vec,n,j,S):
    R,V,D,covers,up=S
    below=[v for v in V if kap_vec in up[v]]  # v <= kappa
    times=sorted({Fr(p,q) for q in range(1,2*j+1) for p in range(1,q) })  # candidate denominators
    times=[t for t in times if (2*j*t).denominator==1]  # bonds<=2 force 2j*t integer (asserted by census)
    out=[]
    def rec(path,ts):
        out.append((tuple(path),tuple(ts)))
        last=path[-1]; t0=ts[-1] if ts else Fr(0)
        for t in times:
            if t<=t0: continue
            for v in V:
                if v!=last and last in up[v] and chain_ok(last,v,t,j,covers,D):
                    rec(path+[v],ts+[t])
    for v in below: rec([v],[])
    return out

def weight(path,ts,j):
    T=[Fr(0)]+list(ts)+[Fr(1)]; w=[Fr(0)]*len(path[0])
    for i,v in enumerate(path):
        for k in range(len(v)): w[k]+= (T[i+1]-T[i])*j*v[k]
    assert all(x.denominator==1 for x in w); return tuple(int(x) for x in w)

def demchar(kap,n,m):
    ws=word(kap,n); f={tuple([m]*n):1}
    for i in reversed(ws): f=Di(f,i,n)
    return f

if __name__=='__main__':
    N=int(sys.argv[1]); J=int(sys.argv[2])
    S=setup(N); R,V,D,covers,up=S
    # bond census against boxes: v <-> parts {N+1-p : v_p=-1}
    bad=0
    for v in V:
        for (u,b) in covers[v]:
            pv={N-i for i in range(N) if v[i]==-1}; pu={N-i for i in range(N) if u[i]==-1}
            diag = len(pu)==len(pv)+1
            if (b==1)!=diag or b not in (1,2): bad+=1
    print('bond census: covers',sum(len(c) for c in covers.values()),'violations',bad)
    tot=0; fails=0
    for L in range(N+1):
        for c in itertools.combinations(range(1,N+1),L):
            kap=tuple(sorted(c,reverse=True))
            kv=tuple(-1 if (N-i) in c else 1 for i in range(N))
            for j in range(1,J+1):
                P=lspaths(kv,N,j,S)
                ch=defaultdict(int)
                for (p,t) in P: ch[weight(p,t,j)]+=1
                dc={k:v for k,v in demchar(set(c),N,j).items() if v}
                F=count(kap,range(0,2*j+1),0,False)
                ok = dict(ch)==dc and len(P)==F
                tot+=1; fails+= (not ok)
                if not ok: print('FAIL',kap,j,len(P),sum(dc.values()),F, dict(ch)==dc)
    print(f'n={N}, j<= {J}: {tot} (kappa,j) cases, failures {fails}')

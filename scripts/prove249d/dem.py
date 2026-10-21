# Demazure characters on LG(n,2n): chi(X_kappa, O(m)) for Schubert varieties, type C_n
import itertools, sympy as sp
from collections import defaultdict
def Di(f,i,n):
    # f: dict weight(tuple)->coeff ; C_n simple roots: a_i=e_i-e_{i+1} (i<n), a_n=2e_n
    g=defaultdict(int)
    for mu,a in f.items():
        mu=list(mu)
        if i<n:
            k=mu[i-1]-mu[i]; al=[0]*n; al[i-1]=1; al[i]=-1
        else:
            k=mu[n-1]; al=[0]*n; al[n-1]=2
        if k>=0:
            for t in range(k+1): g[tuple(mu[j]-t*al[j] for j in range(n))]+=a
        elif k<=-2:
            for t in range(1,-k): g[tuple(mu[j]+t*al[j] for j in range(n))]-=a
    return {k:v for k,v in g.items() if v}
def word(kap,n):
    # reduced word of minimal w with w.omega_n having -1 at positions i=n+1-k, k in kap
    # build by acting: start (1,..,1); to make position p negative: move it to n, flip, move back
    # use greedy: path in W^P via sequence of simple reflections increasing length
    target=tuple(-1 if (n+1-p) in kap else 1 for p in range(1,n+1))
    # BFS over sign vectors with length function
    from collections import deque
    start=tuple([1]*n); prev={start:None}; dq=deque([start])
    while dq:
        v=dq.popleft()
        if v==target: break
        for i in range(1,n+1):
            w=list(v)
            if i<n:
                if not(w[i-1]==1 and w[i]==-1) and not(w[i-1]==1 and w[i]==1 and False): pass
                w[i-1],w[i]=w[i],w[i-1]
            else: w[n-1]=-w[n-1]
            w=tuple(w)
            if w not in prev and sum(1 for x in w if x<0)>=0:
                prev[w]=(v,i); dq.append(w)
    # reconstruct
    ws=[]; v=target
    while prev[v] is not None: v,i=prev[v]; ws.append(i)
    return ws  # word applied right-to-left? ws[0] is last applied
def chi(kap,n,m):
    ws=word(kap,n)
    f={tuple([m]*n):1}  # m*omega_n, omega_n=(1,..,1)
    # X_w character = D_{i1}...D_{ik} e^lambda for w=s_{i1}..s_{ik}; reduced word from BFS: target = s_{ws[0]} ... s_{ws[-1]} start
    for i in reversed(ws): f=Di(f,i,n)
    return sum(f.values()), len(ws)
if __name__=="__main__":
    import sys
    n=int(sys.argv[1]); x=sp.Symbol('m')
    for r in range(n+1):
        for kap in itertools.combinations(range(1,n+1),r):
            vals=[chi(set(kap),n,m)[0] for m in range(0,sum(kap)+3)]
            L=chi(set(kap),n,0)[1]
            P=sp.factor(sp.interpolate(list(zip(range(0,sum(kap)+3),vals)),x))
            print(list(kap),'len',L,'chi',P)

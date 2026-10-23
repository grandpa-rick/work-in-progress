# Q552 light census: u->oo single-child chain for F_lambda (twist lam_j-2j), order bookkeeping from Cor I1:
# head j order th_j-(r-j), tail order -r. Slope -1 consecutive blocks normalised by P-infinity first (formal head/head ties).
# Genuine (head/tail) tie at a node iff tail order == max head order. Compare with predicted discriminant.
import itertools
def norm(kap,th):
    th=list(th); r=len(kap); i=0
    while i<r:
        j=i
        while j+1<r and kap[j+1]==kap[j]+1: j+=1
        if j>i and all(th[k]==th[i]-(k-i) for k in range(i,j+1)):
            base=th[i]-(j-i)
            for k in range(i,j+1): th[k]=base
        i=j+1
    return tuple(th)
def chain(lam,maxsteps=40):
    kap=tuple(lam); th=tuple(lam[j]-2*(j+1) for j in range(len(lam)))
    ties=[]; step=0
    while kap and step<maxsteps:
        step+=1; th=norm(kap,th); r=len(kap)
        heads=[th[j]-(r-1-j) for j in range(r)]
        hm=max(heads); tail=-r
        if hm==tail: ties.append((step,kap,th)); return ties  # stop at first genuine tie (single-child fails)
        if hm>tail:
            js=[j for j in range(r) if heads[j]==hm]
            j=js[-1]   # head/head ties formal after norm; take top
            kap=kap[:j]+kap[j+1:]; th=th[:j]+tuple(t+1 for t in th[j+1:])
        else:
            th=tuple(t+1 for t in th)
    return ties
def pred(lam):
    # top block p..r (1-indexed p), length b
    r=len(lam); p=r
    while p>1 and lam[p-2]==lam[p-1]-1: p-=1
    b=r-p+1; d=lam[p-1]-p
    return d==0 or 1<=d<=b-1
agree=dis=0; bad=[]
for n in range(1,10):
  for L in range(1,6):
    for lam in itertools.combinations(range(1,n+1),L):
      if lam[-1]!=n: continue
      t=bool(chain(lam)); pr=pred(lam)
      if t==pr: agree+=1
      else: dis+=1; bad.append((lam,t,pr,chain(lam)))
print('agree',agree,'disagree',dis)
for x in bad[:25]: print(x)
for lam in [(2,3),(3,4),(2,3,4),(3,4,5),(4,5,6),(1,2,3),(1,3,4),(1,4,5),(2,4,5),(1,2,4,5),(2,5,6,7)]:
    print(lam,'tie:',chain(lam),'pred',pred(lam))

print('=== refined predicate: some maximal block B=(p..q) (1-indexed positions) has lam_p - p < |B|  (i.e. lam_p <= q) ===')
def blocks(lam):
    out=[];i=0;r=len(lam)
    while i<r:
        j=i
        while j+1<r and lam[j+1]==lam[j]+1: j+=1
        out.append((i+1,j+1)); i=j+1
    return out
def pred2(lam): return any(lam[p-1]<=q for p,q in blocks(lam))
def pred2_nobase(lam): return any(lam[p-1]<=q for p,q in blocks(lam) if not (p==1 and lam[0]==1 and q==1))
ag=ds=0; bad=[]; tot=0; Gtie=[]
for n in range(1,12):
  for L in range(1,7):
    for lam in itertools.combinations(range(1,n+1),L):
      if lam[-1]!=n: continue
      tot+=1; t=bool(chain(lam)); pr=pred2(lam)
      if t==pr: ag+=1
      else: ds+=1; bad.append((lam,t,pr,chain(lam)))
      if all(q==p for p,q in blocks(lam)) and t: Gtie.append(lam)
print('lam subset [11], len<=6:',tot,'agree',ag,'disagree',ds); print(bad[:10])
print('G lambdas with a tie:',len(Gtie),'all have lam_1=1:',all(l[0]==1 for l in Gtie), Gtie[:6])

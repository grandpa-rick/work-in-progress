"""Independent S_n check of claimed no-twin shapes. Usage: python3 main16.py n l1,l2,...  [nrandQ]"""
import sys, random, json
from vi import *
n=int(sys.argv[1]); lam=tuple(int(a) for a in sys.argv[2].split(',')); nQ=int(sys.argv[3]) if len(sys.argv)>3 else 0
assert sum(lam)==n
bydes={}
cnt=0
for rw in syt_rowwords(lam):
    cnt+=1; bydes.setdefault(rw_desmask(rw),[]).append(rw)
print('shape',lam,'f=',cnt,'distinct Des:',len(bydes))
def m2s(m): return frozenset(i for i in range(1,n) if m>>i&1)
def s2m(S): return sum(1<<i for i in S)
pairs=[]; viol=[]
for Am in bydes:
    A=m2s(Am); I=interior(A,n); Im=s2m(I)
    Zs=[]
    sub=Im
    while sub:
        if sub in bydes: Zs+=bydes[sub]
        sub=(sub-1)&Im
    if Zs:
        assert I!=A
        viol.append(sorted(A))
        sup=[d for d in bydes if d&Am==Am and d!=Am]
        for X in bydes[Am]:
            for Z in Zs: pairs.append((A,I,X,Z,len(sup)))
print('(lam,A) violations:',len(viol),' candidate pairs:',len(pairs))
Q0=[list(range(sum(lam[:i])+1,sum(lam[:i+1])+1)) for i in range(len(lam))]
def allQ(k):
    out=[]
    allrw=[rw for v in bydes.values() for rw in v]
    random.seed(282)
    return [rw_to_tab(rw,lam) for rw in random.sample(allrw,k)]
Qs=[Q0]+allQ(nQ)
res=[]
for A,I,Xr,Zr,nsup in pairs:
    Px=rw_to_tab(Xr,lam); Pz=rw_to_tab(Zr,lam)
    print('\nA=',sorted(A),'intA=',sorted(I),'#Des>A in shape:',nsup)
    print(' Px=',Px); print(' Pz=',Pz,'Des Pz=',sorted(tab_des(Pz)))
    for qi,Q in enumerate(Qs):
        x=inv_rsk(Px,Q); z=inv_rsk(Pz,Q)
        assert rsk(x)==(Px,Q) and rsk(z)==(Pz,Q) and Ldes(x)==A
        lx,lz=length(x),length(z); par=(lx-lz)%2==0
        line=f'  Q#{qi}: l(x)={lx} l(z)={lz} parity {"OK" if par else "FAIL"}'
        if qi==0: line+=f' x={x}'
        # star images in cell with L not subset A
        kills=[]; info=[]
        for s in range(1,n):
            y=lmul(s,x)
            if rsk(y)[1]!=Q: continue
            Ly=Ldes(y)
            if Ly<=A: continue
            assert abs(length(y)-lx)==1 and (bruhat_le(y,x) or bruhat_le(x,y))  # mu(y,x)=1
            le1=bruhat_le(y,z); le2=bruhat_le(z,y)
            ly=length(y)
            if not le1 and not le2: verdict=('incomparable',0,[])
            else:
                m,path=star_search(y,z)
                verdict=('star' if m is not None else 'undecided',m,path)
            info.append((s,ly,le1,le2,verdict[0],verdict[1],len(verdict[2]) if verdict[2] else 0))
            if verdict[1]==0: kills.append(s)
        line+=f' | (i)-elements s,l(sx),sx<=z,z<=sx,how,mu,pathlen: {info} | kills by s: {kills}'
        print(line)
        res.append(dict(A=sorted(A),Px=Px,Pz=Pz,Q=Q,lx=lx,lz=lz,parity=par,info=info,kills=kills))
json.dump(res,open('res_%d_%s.json'%(n,''.join(map(str,lam))),'w'))

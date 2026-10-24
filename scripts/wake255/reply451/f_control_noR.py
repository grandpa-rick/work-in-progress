# Clio uid451 closed form (BR codim indexing, applied to theta_BR = R(theta)) vs Rick Lemma R (dimension indexing).
# Clio: |c| = 2^{f}, f = #boxes alone in row & column & off-diagonal; 0 if some box is pinned to 1 (right of a row pair,
# bottom of a column pair, or diagonal) AND to 1' (left of a row pair or top of a column pair). Non-rim -> 0 (BR Thm 5.5).
# Rick: V = 0 if 3-chain or diagonal (c,c) with (c-1,c) in theta; else 2^{#off-diag boxes incomparable to all others}.
import sys, itertools
m=int(sys.argv[1])
def D(k): return {(r,c) for r in range(len(k)) for c in range(r,r+k[r])}
def strict(m): 
    for n in range(m+1):
        for S in itertools.combinations(range(m,0,-1),n): yield S
le=lambda x,y: x[0]<=y[0] and x[1]<=y[1]
def rick(th):
    th=list(th)
    for x,y,z in itertools.permutations(th,3):
        if x!=y!=z and le(x,y) and le(y,z) and x!=z: return 0
    S=set(th)
    if any(b[0]==b[1] and (b[0]-1,b[1]) in S for b in th): return 0
    i=sum(1 for b in th if b[0]!=b[1] and not any(c!=b and (le(b,c) or le(c,b)) for c in th))
    return 2**i
def clio(phi):
    phi=set(phi)
    if any(a[0]<b[0] and a[1]<b[1] for a in phi for b in phi): return 0   # not a rim
    f=0
    for b in phi:
        row=[x for x in phi if x[0]==b[0]]; col=[x for x in phi if x[1]==b[1]]
        if len(row)>2 or len(col)>2: return 0
        one = (len(row)==2 and b[1]==max(x[1] for x in row)) or (len(col)==2 and b[0]==max(x[0] for x in col)) or b[0]==b[1]
        prm = (len(row)==2 and b[1]==min(x[1] for x in row)) or (len(col)==2 and b[0]==min(x[0] for x in col))
        if one and prm: return 0
        if len(row)==1 and len(col)==1 and b[0]!=b[1]: f+=1
    return 2**f
R=lambda b:(m-1-b[1],m-1-b[0])
ok=bad=nz=0
for k in strict(m):
    Dk=D(k)
    for rho in strict(m):
        if len(rho)>len(k) or any(rho[i]>k[i] for i in range(len(rho))) or rho==k: continue
        th=Dk-D(rho); a=rick(th); b=clio(th)
        if a==b: ok+=1; nz+=(a>0)
        else: bad+=1; print('MISMATCH',k,rho,a,b)
print('m',m,'pairs rho<kappa:',ok,'agree',bad,'disagree;',nz,'nonzero')

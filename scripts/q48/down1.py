import sys
src=open('../q47/ls12_cor34_check.py').read().split('allok = True')[0]
sys.argv=['x',sys.argv[1]]
exec(src)
def inv(w):
    p,s=w; ip=[0]*m; isg=[0]*m
    for i in range(m): ip[p[i]]=i; isg[p[i]]=s[i]
    return (tuple(ip),tuple(isg))
cpart=lambda k: tuple(-x for x in v_of(k))
zz={k:minrep[cpart(k)] for k in strict}
for r in strict:
  for l in strict:
    if not set(D(r))<=set(D(l)): continue
    d=down(zz[r],cpart(l))
    u=mulk(inv(zz[l]),d)
    print('rho',r,'lam',l,'u=',u[0], 'lenexcess',L[d]-L[zz[l]], 'skew',sorted(D(l)-D(r)))

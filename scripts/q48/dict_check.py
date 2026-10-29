# check z_kappa == w_mu, mu = 180-rotated complement of kappa in the staircase, w_mu = box word (colour of diag d = s_{m-d}), w_{mu+b}=s_c w_mu
import sys
src=open('../q47/ls12_cor34_check.py').read().split('allok = True')[0]
sys.argv=['x',sys.argv[1]]; exec(src)
def s(i):  # 1-based simple reflection as signed perm key
    p=list(range(m)); sg=[1]*m
    if i<m: p[i-1],p[i]=p[i],p[i-1]
    else: sg[m-1]=-1
    return (tuple(p),tuple(sg))
stair={(a,b) for a in range(m) for b in range(a,m)}
e=(tuple(range(m)),tuple([1]*m)); bad=0
for k in strict:
    mu={(m-1-b,m-1-a) for (a,b) in stair-D(k)}
    w=e
    for x in sorted(mu, key=lambda x:(x[0]+x[1],x)):  # linear extension
        w=mulk(s(m-(x[1]-x[0])),w)
    z=minrep[tuple(-x for x in v_of(k))]
    if w!=z or L[w]!=len(mu): bad+=1
print('m',m,'dictionary failures',bad)

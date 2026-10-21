import scert8 as S
S.LIMS=(1,); S.IOTA=True
bad=[x for x in S.allA(5) if not S.start('A',5,*x)]
x=bad[0]; print('triple',x)
orig=S.certify; LOG=[]
def cert(Sx,n,a,b,cr,depth=0):
    r=orig(Sx,n,a,b,cr,depth)
    Sc=S.clean(Sx,n,a,b)
    LOG.append((depth,n,cr,r,Sc))
    return r
S.certify=cert
S.start('A',5,*x)
for depth,n,cr,r,Sc in LOG[::-1][:12]:
    print('depth',depth,'n',n,cr,'cert',r)
    for k,v in Sc.items(): print('     ',v,k[0],'lam',k[1],k[2],'mu',k[3],k[4],'dual',S.comp(n,k[5]),k[6],'d',k[7],'e',k[8][1:])
